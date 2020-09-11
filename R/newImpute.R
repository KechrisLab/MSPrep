#' Function for imputing missing values in data.
#'
#' Replaces missing values with non-zero estimates calculated using a
#' selected method.
#'
#' @param data Data set as either a data frame or `SummarizedExperiement`.
#' @param imputeMethod String specifying imputation method. 
#' Options are "halfmin" (half the minimum value), "bpca" (Bayesian PCA), 
#' and "knn" (k-nearest neighbors).
#' @param kKnn Number of clusters for 'knn' method.
#' @param nPcs Number of  principle components used for re-estimation for 
#' 'bpca' method.
#' @param compoundsAsNeighbors For KNN imputation. If TRUE, compounds will be 
#' used as neighbors rather than samples. Note that using compounds as 
#' neighbors is significantly slower than using samples.
#' @param compVars Vector of the columns which identify compounds. If a 
#' `SummarizedExperiment` is used for `data`, row variables will be used.
#' @param sampleVars Vector of the ordered sample variables found in each sample 
#' column.
#' @param colExtraText Any extra text to ignore at the beginning of the sample 
#' columns names. Unused for `SummarizedExperiments`.
#' @param separator Character or text separating each sample variable in sample
#' columns. Unused for `SummarizedExperiment`.
#' @param missingValue Specifies the abundance value which indicates missing 
#' data. May be a numeric or `NA`.
#' @param returnToSE Logical value indicating whether to return as 
#' `SummarizedExperiment`
#' @param returnToDF Logical value indicating whether to return as data frame.
#' 
#' @return A data frame or `SummarizedExperiment` with missing data imputed.
#' Default return type is set to match the data input but may be altered with 
#' the `returnToSE` or `returnToDF` arguments.
#' 
#' @references 
#'   Oba, S.et al.(2003) A Bayesian missing value estimation for gene
#'   expression profile data. Bioinformatics, 19, 2088-2096
#'
#'   Stacklies, W.et al.(2007) pcaMethods A bioconductor package providing
#'   PCA methods for incomplete data. Bioinformatics, 23, 1164-1167.
#'   
#' @examples
#' # Load, tidy, summarize, and filter example dataset
#' data(msquant)
#' 
#' summarizedDF <- msSummarize(msquant,
#'                             compVars = c("mz", "rt"),
#'                             sampleVars = c("spike", "batch", "replicate", 
#'                             "subject_id"),
#'                             cvMax = 0.50,
#'                             minPropPresent = 1/3,
#'                             colExtraText = "Neutral_Operator_Dif_Pos_",
#'                             separator = "_",
#'                             missingValue = 1)
#'                             
#' filteredDF <- msFilter(summarizedDF,
#'                        filterPercent = 0.8,
#'                        compVars = c("mz", "rt"),
#'                        sampleVars = c("spike", "batch", "subject_id"),
#'                        separator = "_")
#' 
#'                            
#' # Impute dataset using 3 possible options
#' hmImputedDF <- msImpute(filteredDF, imputeMethod = "halfmin",
#'                         compVars = c("mz", "rt"),
#'                         sampleVars = c("spike", "batch", "subject_id"),
#'                         separator = "_",
#'                         missingValue = 0)
#' 
#' bpcaImputedDF <- msImpute(filteredDF, imputeMethod = "bpca",
#'                           compVars = c("mz", "rt"),
#'                           sampleVars = c("spike", "batch", "subject_id"),
#'                           separator = "_",
#'                           missingValue = 0)
#' 
#' knnImputedDF <- msImpute(filteredDF, imputeMethod = "knn",
#'                          compVars = c("mz", "rt"),
#'                          sampleVars = c("spike", "batch", "subject_id"),
#'                          separator = "_",
#'                         missingValue = 0)                                
#'
#' @export
msImpute <- function(data,
                     imputeMethod = c("halfmin", "bpca", "knn"),
                     kKnn = 5, 
                     nPcs = 3, 
                     compoundsAsNeighbors = FALSE,
                     compVars = c("mz", "rt"),
                     sampleVars = c("subject_id"),
                     colExtraText = NULL,
                     separator = NULL,
                     missingValue = NA,
                     returnToSE = FALSE,
                     returnToDF = FALSE) {
    
    imputeMethod <- match.arg(imputeMethod)
    
    if (is(data, "SummarizedExperiment")) {
        return <- .seImpute(data, imputeMethod, kKnn, nPcs, 
                            compoundsAsNeighbors, missingValue)
        if (returnToDF) {
            return <- .seToDF(return)
        }
    } else if (is(data, "data.frame")) {
        return <- .dfImpute(data, imputeMethod, kKnn, nPcs, 
                            compoundsAsNeighbors, compVars, sampleVars, 
                            missingValue)
        if (returnToSE) {
            return <- .dfToSE(return, compVars, sampleVars, separator,
                              colExtraText)
        }
    } else {
        stop("'data' must be a data frame or SummarizedExperiment")
    }
    
    return(return)
}

.seImpute <- function(SE, imputeMethod, kKnn, nPcs, compoundsAsNeighbors,
                      missingValue) {
    
    ## Get assay from SE
    assayData <- as_tibble(assay(SE))
    
    ## Select abundance columns to impute, replace missingValue w/ NA
    assayData <- mutate_all(assayData, .replaceMissing, missingValue,
                            setMissing = NA)
    
    ## Impute data
    imputedData <- switch(imputeMethod,
                          "halfmin" = .imputeHalfmin(assayData),
                          "bpca" = .imputeBpca(assayData, nPcs),
                          "knn" = .imputeKnn(assayData, kKnn, 
                                             compoundsAsNeighbors),
                          stop("Invalid impute method - you should never see 
                               this warning."))
    
    ## Return imputed data to SE
    assay(SE, 1, withDimnames = FALSE) <- as.matrix(imputedData)
    
    return(SE)
    
}

#' @importFrom dplyr mutate_all
#' @importFrom dplyr summarise_all
#' @importFrom dplyr bind_cols
.dfImpute <- function(data, imputeMethod, kKnn, nPcs, compoundsAsNeighbors,
                      compVars, sampleVars, missingValue) {
    
    ## Select abundance columns to impute, replace missingValue w/ NA
    abundanceColumns <- select(data, -compVars) %>%
        mutate_all(.replaceMissing, missingValue, setMissing = NA)
    
    missingCount <- sum(abundanceColumns %>% 
                            summarise_all(function(x) sum(is.na(x))))
    
    ## Select compound columns to reattach later
    compColumns <- select(data, compVars)
    
    ## Impute data
    imputedData <- switch(imputeMethod,
                          "halfmin" = .imputeHalfmin(abundanceColumns),
                          "bpca" = .imputeBpca(data = abundanceColumns, 
                                               nPcs = nPcs),
                          "knn" = .imputeKnn(abundanceColumns, kKnn, 
                                             compoundsAsNeighbors),
                          stop("Invalid impute method - you should never see 
                               this warning."))
    
    ## Recombine data and return
    data <- bind_cols(compColumns, imputedData)
}

#' @importFrom VIM kNN
.imputeKnn <- function(data, kKnn = 5, compoundsAsNeighbors) {
    
    if (compoundsAsNeighbors == TRUE) {
        data <- VIM::kNN(as.data.frame(t(data)), k = kKnn, imp_var = FALSE)
        data <- t(data)
    }
    else {
        data <- VIM::kNN(as.data.frame(data), k = kKnn, imp_var = FALSE)
        
    }

    data <- .halfminIfAnyNegative(data, "knn")
    
}



#' @importFrom pcaMethods pca
#' @importFrom pcaMethods completeObs
.imputeBpca <- function(data, nPcs = 3) {
    
    data <- t(data)
    data <- pca(data, nPcs = nPcs, method = "bpca")
    data <- as.data.frame(t(completeObs(data))) # extract imputed dataset
    data <- .halfminIfAnyNegative(data, "BPCA")
    
}

.imputeHalfmin <- function(data) {
    
    halfmin <- function(x) {
        ifelse(is.na(x), min(x, na.rm = TRUE) / 2, x)
    }
    
    data <- as.data.frame(apply(apply(data, 1, halfmin), 1, identity))
    
}

.halfminIfAnyNegative <-  function(data, method) {
    
    numNeg <- sum(data %>% summarise_all(function(x) sum(x < 0)))
    
    if (any(numNeg)) {
        message("Found ", numNeg, " negative imputed values using ", method, 
                ", reverting to half-min imputation for these values")
        data <- as.data.frame(data) %>% mutate_all(.setMissingNegativeVals)
        data <- .imputeHalfmin(data)
    }
    
    return(data)
    
}

.setMissingNegativeVals <- function(var) ifelse(var < 0, NA, var)
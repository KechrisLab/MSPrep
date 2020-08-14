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
    } else if (is(data, "data.frame")) {
        return <- .dfImpute(data, imputeMethod, kKnn, nPcs, 
                            compoundsAsNeighbors, compVars, sampleVars, 
                            colExtraText, separator, missingValue, returnToSE,
                            returnToDF)
    } else {
        stop("'data' must be a data frame or SummarizedExperiment")
    }
}

.seImpute <- function(SE, imputeMethod, kKnn, nPcs, compoundsAsNeighbors,
                      missingValue) {
    
    ## Get assay from SE
    assayData <- as_tibble(assay(SE))
    
    ## Select abundance columns to impute, replace missingValue w/ NA
    assayData <- mutate_all(assayData, .replaceMissing, missingValue)
    
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
                      compVars, sampleVars, colExtraText, separator, missingValue,
                      returnToSE, returnToDF) {
    
    ## Select abundance columns to impute, replace missingValue w/ NA
    abundanceColumns <- select(data, -compVars) %>%
        mutate_all(.replaceMissing, missingValue)
    
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
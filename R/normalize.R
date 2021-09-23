#' Function for performing normalization and batch corrections on imputed data.
#' 
#' Perform normalization and batch corrections on specified imputed dataset.
#' Routines included are quantile, RUV (remove unwanted variation), SVA 
#' (surrogate variable analysis), median, CRMN (cross-contribution 
#' compensating multiple standard normalization), and ComBat to remove batch 
#' effects in raw, quantile, and median normalized data. Generates data 
#' driven controls if none exist.
#' 
#' @param data Data set as either a data frame or `SummarizedExperiement`.
#' @param normalizeMethod  Name of normalization method.
#' "ComBat" (only ComBat batch correction), "quantile" (only quantile 
#' normalization), "quantile + ComBat" (quantile with ComBat batch correction),
#' "median" (only median normalization), "median + ComBat" (median with ComBat
#' batch correction), "CRMN" (cross-contribution compensating multiple 
#' standard normalization), "RUV" (remove unwanted variation), "SVA" (surrogate 
#' variable analysis)
#' @param nControl Number of controls to estimate/utilize (for CRMN and RUV).
#' @param controls Vector of control identifiers.  Leave blank for data driven
#' controls. Vector of column numbers from metafin dataset of that control (for
#' CRMN and RUV).
#' @param nComp Number of factors to use in CRMN algorithm. 
#' @param kRUV Number of factors to use in RUV algorithm.
#' @param batch Name of the sample variable identifying batch.
#' @param covariatesOfInterest Sample variables used as covariates in
#' normalization algorithms (required for ComBat, CRMN, and SVA).
#' @param transform  Select transformation to apply to data prior to 
#' normalization. Options are "log10", "log2", "ln" and "none".
#' @param compVars Vector of the columns which identify compounds. If a 
#' `SummarizedExperiment` is used for `data`, row variables will be used.
#' @param sampleVars Vector of the ordered sample variables found in each sample
#' column.
#' @param colExtraText Any extra text to ignore at the beginning of the sample 
#' columns names. Unused for `SummarizedExperiments`.
#' @param separator Character or text separating each sample variable in sample
#' columns. Unused for `SummarizedExperiment`.
#' @param returnToSE Logical value indicating whether to return as 
#' `SummarizedExperiment`
#' @param returnToDF Logical value indicating whether to return as data frame.
#' 
#' @return  A data frame or `SummarizedExperiment` with transformed and
#' normalized data. Default return type is set to match the data input but may 
#' be altered with the `returnToSE` or `returnToDF` arguments.
#' 
#' @references 
#' Bolstad, B.M.et al.(2003) A comparison of normalization methods for high
#' density oligonucleotide array data based on variance and bias.
#' Bioinformatics, 19, 185-193
#' 
#' DeLivera, A.M.et al.(2012) Normalizing and Integrating Metabolomic Data.
#' Anal. Chem, 84, 10768-10776.
#' 
#' Gagnon-Bartsh, J.A.et al.(2012) Using control genes to correct for unwanted
#' variation in microarray data. Biostatistics, 13, 539-552.
#' 
#' Johnson, W.E.et al.(2007) Adjusting batch effects in microarray expression
#' data using Empirical Bayes methods. Biostatistics, 8, 118-127.
#' 
#' Leek, J.T.et al.(2007) Capturing Heterogeneity in Gene Expression Studies by
#' Surrogate Variable Analysis. PLoS Genetics, 3(9), e161
#' 
#' Wang, W.et al.(2003) Quantification of Proteins and Metabolites by Mass
#' Spectrometry without Isotopic Labeling or Spiked Standards. Anal. Chem., 75,
#' 4818-4826.
#' 
#' @examples
#' # Load, tidy, summarize, filter, and impute example dataset
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
#' hmImputedDF <- msImpute(filteredDF, imputeMethod = "halfmin",
#'                         compVars = c("mz", "rt"),
#'                         sampleVars = c("spike", "batch", "subject_id"),
#'                         separator = "_",
#'                         missingValue = 0)
#' 
#' # Normalize data set
#' medianNormalizedDF <- msNormalize(hmImputedDF, normalizeMethod = "median",
#'                                   compVars = c("mz", "rt"),
#'                                   sampleVars = c("spike", "batch", 
#'                                   "subject_id"),
#'                                   separator = "_")
#'
#' @export
msNormalize <- function(data, 
                        normalizeMethod = c("median", "ComBat", "quantile",
                                            "quantile + ComBat", 
                                            "median + ComBat", "CRMN", "RUV",
                                            "SVA"),
                        nControl = 10, controls  = NULL, nComp = 2, kRUV = 3,
                        batch = "batch", covariatesOfInterest = NULL,
                        transform = c("log10", "log2", "ln", "none"),
                        compVars = c("mz", "rt"),
                        sampleVars = c("subject_id"),
                        colExtraText = NULL,
                        separator = NULL,
                        returnToSE = FALSE,
                        returnToDF = FALSE) {
    
    normalizeMethod <- match.arg(normalizeMethod)
    transform <- match.arg(transform)
    
    .normalizeParamValidation(data, normalizeMethod, nControl, controls, nComp,
                              kRUV, batch, covariatesOfInterest, transform, 
                              compVars, sampleVars, colExtraText, separator, 
                              returnToSE, returnToDF)
    
    if (is(data, "SummarizedExperiment")) {
        return <- .seNormalize(data, normalizeMethod, nControl, controls, nComp,
                               kRUV, batch, transform, covariatesOfInterest)
        if (returnToDF) {
            return <- .seToDF(return)
        }
    } else if (is(data, "data.frame")) {
        return <- .dfNormalize(data, normalizeMethod, nControl, controls, nComp,
                               kRUV, batch, transform, compVars, 
                               sampleVars, colExtraText, separator, returnToSE, 
                               covariatesOfInterest)
        if (returnToSE) {
            return <- .dfToSE(return, compVars, sampleVars, separator,
                              colExtraText)
        }
    } else {
        stop("'data' must be a data frame or SummarizedExperiment")
    }
    
    return(return)
}

.dfNormalize <- function(data, normalizeMethod, nControl, controls, nComp, kRUV,
                         batch, transform, compVars, sampleVars, 
                         colExtraText, separator, returnToSE, 
                         covariatesOfInterest) {
    
    ## Get abundance and compound columns, transform data
    abundanceColumns <- select(data, -compVars) %>%
        .dfLogTransform(transform)
    compColumns <- select(data, compVars)
    
    ## Check for NA values
    if (any(apply(abundanceColumns, 2, function(x) any(is.na(x))))) {
        "NA values present in data"
    }
    
    normalizedData <-
        switch(normalizeMethod,
               "ComBat" = .dfNormalizeCombat(abundanceColumns, batch, 
                                             sampleVars, colExtraText, 
                                             separator, covariatesOfInterest),
               "quantile" = .normalizeQuantile(abundanceColumns),
               "quantile + ComBat" = 
                   .dfNormalizeQuantileCombat(abundanceColumns, batch, 
                                              sampleVars, colExtraText, 
                                              separator, covariatesOfInterest),
               "median" = .normalizeMedian(abundanceColumns),
               "median + ComBat" = 
                   .dfNormalizeMedianCombat(abundanceColumns, batch, sampleVars,
                                            colExtraText, separator,
                                            covariatesOfInterest),
               "CRMN" = .normalizeCRMN(abundanceColumns, compColumns, compVars, 
                                       sampleVars, colExtraText, separator, 
                                       covariatesOfInterest, nComp, nControl, 
                                       controls),
               "RUV" = .normalizeRUV(abundanceColumns, compColumns, compVars, 
                                     sampleVars, kRUV, nControl, controls, 
                                     colExtraText, separator),
               "SVA" = .normalizeSVA(abundanceColumns, sampleVars, 
                                     colExtraText, separator, 
                                     covariatesOfInterest),
               stop("Invalid normalize method - provide argument from list in",
                    "function definition and help file"))
    
    ## Recombine data
    if(normalizeMethod != "CRMN") {
        normalizedData <- bind_cols(compColumns, normalizedData)
    }
    
    return(normalizedData)
}

.normalizeSVA <- function(data, sampleVars, colExtraText, separator, 
                          covariatesOfInterest) {
    
    svaFactors <- .svaFactors(data, sampleVars, colExtraText, separator, 
                              covariatesOfInterest)
    
    normalizedData <- .genAdj(data, svaFactors)
}

#' @importFrom stats model.matrix
#' @importFrom sva sva
#' @importFrom tidyr separate
.svaFactors <- function(data, sampleVars, colExtraText, separator, 
                        covariatesOfInterest) {
    
    colNames <- data.frame("colNames" = colnames(data))
    colData <- separate(colNames, col = "colNames", into = sampleVars, 
                        sep = separator, remove = FALSE)
    
    sampleInfo <- select(colData, covariatesOfInterest)
    sampleInfo <- model.matrix(~ ., data = sampleInfo)
    
    svaFactors <- sva(as.matrix(data), sampleInfo, method = "irw")
    svaFactors <- svaFactors$sv
}

#' @importFrom stats lm
.genAdj <- function(data, factors) {
    matData <- as.matrix(data)
    normalizedData <- vapply(seq_len(nrow(matData)),
                            function(j) {
                                fit <- lm(matData[j, ] ~
                                              as.matrix(factors, ncol = 1))
                                fit$fitted.values
                                },
                            numeric(ncol(data)))
    
    normalizedData <- t(normalizedData)
    
    rownames(normalizedData) <- rownames(data)
    as.data.frame(normalizedData)
}

#' @importFrom sva ComBat
#' @importFrom dplyr mutate_if
.dfNormalizeCombat <- function(data, batch, sampleVars, colExtraText, separator,
                               covariatesOfInterest) {
    
    if (is.null(batch) & "batch" %in% sampleVars) {
        batch <- "batch"
    } 
    
    if (is.null(batch)) {
        stop("'batch' must be included for ComBat normalization methods")
    }
    
    if (length(batch) > 1) stop("only one batch variable allowed for ComBat") 
    
    colNames <- data.frame("colNames" = colnames(data))
    colData <- separate(colNames, col = "colNames", into = sampleVars, 
                     sep = separator, remove = FALSE)

    # Create model matrix of ColData
    sampleInfo <- select(colData, covariatesOfInterest)
    sampleInfo <- model.matrix(~ ., data = sampleInfo)
    
    batchInfo <- select(colData, batch) %>%
        mutate_if(is.character, as.factor)
    batchVec <- batchInfo[, 1]
    
    normalizedData <- sva::ComBat(data, mod = sampleInfo, batch = batchVec)
    
    as.data.frame(normalizedData)
}

#' @importFrom preprocessCore normalize.quantiles
.normalizeQuantile <- function(data) {
    colNames <- colnames(data)
    normalizedData <- as.data.frame(normalize.quantiles(as.matrix(data)))
    colnames(normalizedData) <- colNames
    return(normalizedData)
}

.dfNormalizeQuantileCombat <- function(data, batch, sampleVars, colExtraText, 
                                       separator, covariatesOfInterest) {
    quantNormalized <- .normalizeQuantile(data)
    quantCombNormalize <- .dfNormalizeCombat(quantNormalized, batch, sampleVars,
                                             colExtraText, separator,
                                             covariatesOfInterest)
    return(quantCombNormalize)
}

.normalizeMedian <- function(data) {
    
    ## Subtracts median of each COMPOUND
    ## normalizedData <- sweep(data, 1, apply(data, 1, median), "-")
    
    ## Subtracts median of each SAMPLE
    normalizedData <- as.data.frame(vapply(data, function(x) x - median(x),
                                    FUN.VALUE = numeric(nrow(data))))
    
}

.dfNormalizeMedianCombat <- function(data, batch, sampleVars, colExtraText, 
                                   separator, covariatesOfInterest) {
    medNormalized <- .normalizeMedian(data)
    medCombNormalized <- .dfNormalizeCombat(data, batch, sampleVars,
                                           colExtraText, separator,
                                           covariatesOfInterest)
    return(medCombNormalized)
}

#' @importFrom crmn normalize
.normalizeCRMN <- function(abCols, compCols, compVars, sampleVars, colExtraText,
                           separator, covariatesOfInterest, nComp, nControl, 
                           controls) {
    
    crmnInputs <- .createCRMNInputs(abCols, compCols, compVars, sampleVars, 
                                    colExtraText, separator, 
                                    covariatesOfInterest, nComp, nControl, 
                                    controls)
    
    crmnNormalized  <- normalize(as.matrix(crmnInputs$data),
                                 method = "crmn",
                                 factors = crmnInputs$modelMatrix,
                                 standards = crmnInputs$ISVec,
                                 lg = FALSE)
    
    compCols <- compCols[!crmnInputs$ISVec, ]
    
    cbind(compCols, as.data.frame(crmnNormalized))
    
}

.createCRMNInputs <- function(abCols, compCols, compVars, sampleVars, 
                              colExtraText, separator, covariatesOfInterest, 
                              nComp, nControl, controls) {
    
    svaFactors <- .svaFactors(abCols, sampleVars, colExtraText, separator, 
                              covariatesOfInterest)
    
    nCompounds <- nrow(abCols)
    
    colNames <- data.frame("colNames" = colnames(abCols))
    colData <- separate(colNames, col = "colNames", into = sampleVars, 
                        sep = separator, remove = FALSE)
    sampleInfo <- select(colData, covariatesOfInterest)
    sampleInfo <- model.matrix(~ -1 + ., data = sampleInfo)
    
    # Transform data to tidy
    fullData <- bind_cols(compCols, abCols)
    data <- .msTidy(fullData, compVars, sampleVars, colExtraText, separator)
    
    ISVec <- .isIS(data, compVars, sampleVars, svaFactors, nCompounds, 
                   nControl, controls)
    
    list("data" = abCols, "modelMatrix" = sampleInfo, "ISVec" = ISVec)
    
}

## Requires tidy data
.isIS <- function(data, compVars, sampleVars, svaFactors, nComp, nControl, 
                  controls) {
    
    controlSummary <- .controlSummary(data, compVars, sampleVars)
    ctl <- .ctlCompounds(controlSummary, nControl, controls)
    ctlo <- ctl[order(ctl)]
    j <- 1
    ISvec <- rep(FALSE, nComp)
    
    for (i in seq_len(nComp)) {
        if (j <= 10) {
            if (ctlo[j] == i) {
                ISvec[i] <- TRUE
                j       <- j + 1
            } else {
                ISvec[i] <- FALSE
            }
        }
    }
    
    return(ISvec)
    
}

#' @importFrom dplyr rowwise
## Requires tidy data
.controlSummary <- function(data, compVars, sampleVars) {
    
    compVars <- syms(compVars)
    
    counteq0 <- function(x) sum(x == 0)
    
    rtn <- group_by(data, `!!!`(compVars)) %>% 
        summarise(mean = mean(.data$abundance), sd = sd(.data$abundance),
                     counteq0 = counteq0(.data$abundance))
    rtn <- ungroup(rtn)
    rtn <- mutate(rtn, cv = sd / mean, number = seq_len(n()))
    rtn <- select(rtn, `!!!`(compVars), .data$number, .data$mean, .data$sd,
                  .data$cv, .data$counteq0)
    
}

#' @importFrom dplyr arrange
.ctlCompounds <- function(data, nControl, controls) {
    
    if (length(controls) > 0) { 
        dat <- controls 
    } else {
        dat <- arrange(data, .data$counteq0, .data$cv)
        dat <- dat[seq_len(nControl), "number"]
        dat <- dat[["number"]]
    }
    
    return(dat)
}

.normalizeRUV <- function(abCols, compCols, compVars, sampleVars, kRUV, 
                          nControl, controls, colExtraText, separator) {
    
    if(kRUV > nControl){
        stop ("kRUV must be less than or equal to nControl")
    }
    
    ruvFactors <- .ruvFactors(abCols, compCols, compVars, sampleVars, kRUV,
                              nControl, controls, colExtraText, separator)
    
    adjusted <- .genAdj(abCols, ruvFactors)
    
    cat("\n")
    
    return(adjusted)
}

.ruvFactors <- function(abCols, compCols, compVars, sampleVars, kRUV, nControl, 
                        controls, colExtraText, separator) {
    
    # ctl -- n compounds w/ lowest CV and no missing
    fullData <- bind_cols(compCols, abCols)
    data <- .msTidy(fullData, compVars, sampleVars, colExtraText, separator)
    controlSummary <- .controlSummary(data, compVars, sampleVars)
    ctl <- .ctlCompounds(controlSummary, nControl, controls)
    
    Y <- as.matrix(abCols) # Only data Used in calculation
    
    # Calculate RUV 
    #   Uses: Y, ctl, sva_factors
    Z   <- matrix(rep(1, ncol(Y)))
    RZY <- Y - Y %*% Z %*% solve(t(Z) %*% Z) %*% t(Z)
    W   <- svd(RZY[ctl, ])$v
    W   <- W[, seq_len(kRUV)]
    
    # Format output
    rtn <- as.matrix(W, ncol = kRUV)
    colnames(rtn) <- paste0("f", seq_len(ncol(rtn)))
    
    return(rtn)
}

.dfLogTransform <- function(data, transform) {
    data <- switch(transform,
                   "log10" = log10(data),
                   "log2" = log2(data),
                   "ln" = log(data),
                   "none" = data)
}

.normalizeParamValidation <- function(data, normalizeMethod, nControl, controls,
                                      nComp, kRUV, batch, covariatesOfInterest,
                                      transform, compVars, sampleVars, 
                                      colExtraText, separator, returnToSE, 
                                      returnToDF) {
    
    c("median", "ComBat", "quantile", "quantile + ComBat", "median + ComBat", 
      "CRMN", "RUV", "SVA")
    
    if (returnToSE && returnToDF) {
        stop("Only one of returnToSE and returnToDF may be TRUE")
    }
    
    if(is.null(covariatesOfInterest) && normalizeMethod %in% 
       c("quantile + ComBat", "median + ComBat", "CRMN", "SVA")) {
        stop("covariatesOfInterest must be included for ComBat, CRMN, and SVA")
    }
    
    if (is(data, "data.frame")) {
        
        .dfParamValidation(data, compVars, sampleVars, colExtraText, separator)
        
    } else if (is(data, "SummarizedExperiment")) {
        
        if (length(assays(data)) != 1) {
            stop("Current version of MSPrep only supports one assay")
        }
        
    }
}
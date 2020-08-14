msNormalize <- function(data, 
                        normalizeMethod = c("median", "ComBat", "quantile",
                                            "quantile + ComBat", 
                                            "median + ComBat", "CRMN", "RUV",
                                            "SVA"),
                        nControl = 10, controls  = NULL, nComp = 2, kRUV = 3,
                        batch = NULL, covariatesOfInterest = NULL,
                        transform = c("log10", "log2", "none"),
                        compVars = c("mz", "rt"),
                        sampleVars = c("subject_id"),
                        colExtraText = NULL,
                        separator = NULL,
                        missingValue = NA,
                        returnToSE = FALSE,
                        returnToDF = FALSE) {
    
    normalizeMethod <- match.arg(normalizeMethod)
    transform <- match.arg(transform)
    

    
    if (is(data, "SummarizedExperiment")) {
        return <- .seNormalize( )
    } else if (is(data, "data.frame")) {
        return <- .dfNormalize(data, normalizeMethod, nControl, controls, nComp, 
                               kRUV, batch, transform, compVars, 
                               sampleVars, colExtraText, separator, 
                               missingValue, returnToSE, covariatesOfInterest)
    } else {
        stop("'data' must be a data frame or SummarizedExperiment")
    }
}

.dfNormalize <- function(data, normalizeMethod, nControl, controls, nComp, kRUV, 
                         batch, transform, compVars, sampleVars, 
                         colExtraText, separator, missingValue, returnToSE, 
                         covariatesOfInterest) {
    
    ## Get abundance and compound columns, transform data
    abundanceColumns <- select(data, -compVars) %>%
        .dfLogTransform(transform)
    compColumns <- select(data, compVars)
    
    ## Check for NA values
    if (any(apply(abundanceColumns, 2, function(x) any(is.na(x))))) {
        "NA values present in data"
    }
    
    ## Normalize with appropriate method
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

.seNormalize <- function(data, normalizeMethod, nControl, controls, nComp, kRUV,
                         transform, returnToDF) {
    return(0)
}

.normalizeSVA <- function(data, sampleVars, colExtraText, separator, 
                          covariatesOfInterest) {
    # if (is.null(spike) & "spike" %in% sampleVars) {
    #     spike <- "spike"
    # } 
    # 
    # if (is.null(spike)) {
    #     stop("'spike' must be included for SVA normalization")
    # }
    
    svaFactors <- .svaFactors(data, sampleVars, colExtraText, separator, 
                              covariatesOfInterest)
    
    normalizedData <- .genAdj(data, svaFactors)
}

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

.genAdj <- function(data, factors) {
    matData <- as.matrix(data)
    normalizedData <- sapply(seq_len(nrow(matData)),
                            function(j) {
                                fit <- lm(matData[j, ] ~
                                              as.matrix(factors, ncol = 1))
                                fit$fitted.values
                                })
    # for(j in 1:nrow(matData)) {
    #     fit <- lm(matData[j, ] ~ as.matrix(factors, ncol = 1))
    #     matData[j, ] <- fit$fitted.values
    # }
    normalizedData <- t(normalizedData)
    
    rownames(normalizedData) <- rownames(data)
    as.data.frame(normalizedData)
}

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
    
    as.data.frame(crmnNormalized)
    
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
    data <- msTidy(fullData, compVars, sampleVars, colExtraText, separator)
    
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
    
    for (i in 1:nComp) {
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
    rtn <- mutate(rtn, cv = sd / mean, number = 1:n())
    rtn <- select(rtn, `!!!`(compVars), .data$number, .data$mean, .data$sd,
                  .data$cv, .data$counteq0)
    
}

#' @importFrom dplyr arrange
.ctlCompounds <- function(data, nControl, controls) {
    
    if (length(controls) > 0) { 
        dat <- controls 
    } else {
        dat <- arrange(data, .data$counteq0, .data$cv)
        dat <- dat[1:nControl, "number"]
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
    data <- msTidy(fullData, compVars, sampleVars, colExtraText, separator)
    controlSummary <- .controlSummary(data, compVars, sampleVars)
    ctl <- .ctlCompounds(controlSummary, nControl, controls)
    
    Y <- as.matrix(abCols) # Only data Used in calculation
    
    # Calculate RUV 
    #   Uses: Y, ctl, sva_factors
    Z   <- matrix(rep(1, ncol(Y)))
    RZY <- Y - Y %*% Z %*% solve(t(Z) %*% Z) %*% t(Z)
    W   <- svd(RZY[ctl, ])$v
    W   <- W[, 1:kRUV]
    
    # Format output
    rtn <- as.matrix(W, ncol = kRUV)
    colnames(rtn) <- paste0("f", 1:ncol(rtn))
    
    return(rtn)
}

.dfLogTransform <- function(data, transform) {
    data <- switch(transform,
                   "log10" = log10(data),
                   "log2" = log2(data),
                   "none" = data)
}

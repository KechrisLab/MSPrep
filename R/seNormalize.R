.seNormalize <- function(SE, normalizeMethod, nControl, controls, nComp, kRUV, 
                         batch, transform, covariatesOfInterest) {
    
    ## Get abundance columns, transform data
    abundanceColumns <- as_tibble(assay(SE)) %>%
        .dfLogTransform(transform)
    
    assay(SE, 1, withDimnames = FALSE) <- as.matrix(abundanceColumns)
    
    ## Check for NA values
    if (any(apply(abundanceColumns, 2, function(x) any(is.na(x))))) {
        "NA values present in data"
    }
    
    ## Normalize with appropriate method
    normalizedSE <-
        switch(normalizeMethod,
               "ComBat" = .seNormalizeCombat(SE, batch, covariatesOfInterest),
               "quantile" = .seNormalizeQuantile(SE),
               "quantile + ComBat" = 
                   .seNormalizeQuantileCombat(SE, batch, covariatesOfInterest),
               "median" = .seNormalizeMedian(SE),
               "median + ComBat" = 
                   .seNormalizeMedianCombat(SE, batch, covariatesOfInterest),
               "CRMN" = .seNormalizeCRMN(SE, covariatesOfInterest, nComp, 
                                         nControl, controls),
               "RUV" = .seNormalizeRUV(SE, kRUV, nControl, controls),
               "SVA" = .seNormalizeSVA(SE, covariatesOfInterest),
               stop("Invalid normalize method - provide argument from list in",
                    "function definition and help file"))
    
    return(normalizedSE)
}

.seNormalizeCombat <- function(SE, batch, covariatesOfInterest) {
    
    if (is.null(batch) & "batch" %in% colnames(colData(SE))) {
        batch <- "batch"
    } 
    
    if (is.null(batch)) {
        stop("'batch' must be included for ComBat normalization methods")
    }
    
    if (length(batch) > 1) stop("only one batch variable allowed for ComBat") 
    
    colData <- as.data.frame(colData(SE))
    sampleInfo <- select(colData, covariatesOfInterest)
    sampleInfo <- model.matrix(~ ., data = sampleInfo)
    
    batchInfo <- select(colData, batch) %>%
        mutate_if(is.character, as.factor)
    batchVec <- batchInfo[, 1]
    
    abundanceData <- assay(SE)
    normalizedData <- sva::ComBat(abundanceData, mod = sampleInfo, 
                                  batch = batchVec)
    
    assay(SE, 1, withDimnames = FALSE) <- as.matrix(normalizedData)
    
    return(SE)
    
}

#' @importFrom preprocessCore normalize.quantiles
.seNormalizeQuantile <- function(SE) {
    abundanceData <- assay(SE)
    normalizedData <- normalize.quantiles(abundanceData)
    
    assay(SE, 1, withDimnames = FALSE) <- as.matrix(normalizedData)
    return(SE)
}

.seNormalizeQuantileCombat <- function(SE, batch, covariatesOfInterest) {
    
    quantNormalized <- .seNormalizeQuantile(SE)
    
    quantCombNormalize <- .seNormalizeCombat(quantNormalized, batch,
                                             covariatesOfInterest)
    
    return(quantCombNormalize)
}

.seNormalizeMedian <- function(SE) {
    
    abundanceData <- as.data.frame(assay(SE))
    
    normalizedData <- vapply(abundanceData, function(x) x - median(x),
                             FUN.VALUE = numeric(nrow(abundanceData)))
    
    assay(SE, 1, withDimnames = FALSE) <- as.matrix(normalizedData)
    return(SE)
}

.seNormalizeMedianCombat <- function(SE, batch, covariatesOfInterest) {
    medianNormalized <- .seNormalizeMedian(SE)
    
    medCombatNormalized <- .seNormalizeCombat(SE, batch,covariatesOfInterest)
    
    return(medCombatNormalized)
}

#' @importFrom crmn normalize
.seNormalizeCRMN <- function(SE, covariatesOfInterest, nComp, nControl,
                             controls) {
    
    rowData <- rowData(SE)
    colData <- colData(SE)
    metaData <- metadata(SE)
    
    
    crmnInputs <- .seCreateCRMNInputs(SE, covariatesOfInterest, nComp, nControl,
                                      controls)
    
    rowData <- rowData[!crmnInputs$ISVec, ]
    
    crmnNormalized  <- normalize(as.matrix(crmnInputs$data),
                                 method = "crmn",
                                 factors = crmnInputs$modelMatrix,
                                 standards = crmnInputs$ISVec,
                                 lg = FALSE)
    
    rtn <- SummarizedExperiment(assays = 
                                    list(abundance = 
                                             as.matrix((crmnNormalized))),
                                colData = colData,
                                rowData = rowData,
                                metadata = metaData)
    
    
}

.seCreateCRMNInputs <- function(SE, covariatesOfInterest, nComp, nControl,
                                controls) {
    
    svaFactors <- .seSVAFactors(SE, covariatesOfInterest)
    
    nCompounds <- nrow(assay(SE))
    
    colData <- as_tibble(colData(SE))
    sampleInfo <- select(colData, covariatesOfInterest)
    sampleInfo <- model.matrix(~ -1 + ., data = sampleInfo)
    
    # Transform data to tidy
    data <- .msTidy(SE)
    
    ISVec <- .isIS(data, compVars=colnames(rowData(SE)), 
                   sampleVars=colnames(colData(SE)), svaFactors, nCompounds, 
                   nControl, controls)
    
    list("data" = assay(SE), "modelMatrix" = sampleInfo, "ISVec" = ISVec)
    
}

.seNormalizeRUV <- function(SE, kRUV, nControl, controls) {
    
    if(kRUV > nControl){
        stop ("kRUV must be less than or equal to nControl")
    }
    
    ruvFactors <- .seRUVFactors(SE, kRUV, nControl, controls)
    
    normalizedData <- .genAdj(assay(SE), ruvFactors)
    
    cat("\n")
    
    assay(SE, 1, withDimnames = FALSE) <- as.matrix(normalizedData)
    
    return(SE)
}

.seRUVFactors <- function(SE, kRUV, nControl, controls) {
    
    # ctl -- n compounds w/ lowest CV and no missing
    data <- .msTidy(SE)
    compVars <- colnames(rowData(SE))
    sampleVars <- colnames(colData(SE))
    controlSummary <- .controlSummary(data, compVars, sampleVars)
    ctl <- .ctlCompounds(controlSummary, nControl, controls)
    
    Y <- as.matrix(assay(SE)) # Only data Used in calculation
    
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


.seNormalizeSVA <- function(SE, covariatesOfInterest) {
    
    svaFactors <- .seSVAFactors(SE, covariatesOfInterest)
    abundanceData <- assay(SE)
    
    normalizedData <- .genAdj(abundanceData, svaFactors)
    
    assay(SE, 1, withDimnames = FALSE) <- as.matrix(normalizedData)
    
    return(SE)
}

.seSVAFactors <- function(SE, covariatesOfInterest) {
    
    colData <- as_tibble(colData(SE))
    sampleInfo <- select(colData, covariatesOfInterest)
    sampleInfo <- model.matrix(~ ., data = sampleInfo)
    
    svaFactors <- sva(as.matrix(assay(SE)), sampleInfo, method = "irw")
    svaFactors <- svaFactors$sv
    
}
#' Summarize, filter, impute, transform and normalize metabolomics dataset
#' 
#' Wrapper function for the entire MSPrep pre-analytics pipeline. Calls 
#' msSummarize(), msFilter, msImpute(), and msNormalize().
#' 
#' @param data Data set as either a data frame or `SummarizedExperiement`.
#' @param cvMax Decimal value from 0 to 1 representing the acceptable level of 
#' coefficient of variation between replicates.
#' @param minPropPresent  Decimal value from 0 to 1 representing the 
#' minimum proportion present to summarize with median or mean. Below this the 
#' compound will be set to 0.
#' @param filterPercent Decimal value indicating filtration threshold. 
#' Compounds which are present in fewer samples than the specified proportion 
#' will be removed. 
#' @param imputeMethod String specifying imputation method. 
#' Options are "halfmin" (half the minimum value), "bpca" (Bayesian PCA), 
#' and "knn" (k-nearest neighbors), or "none" to skip imputation.
#' @param kKnn Number of clusters for 'knn' method.
#' @param nPcs Number of  principle components used for re-estimation for 
#' 'bpca' method.
#' @param maxIterRf Maximum number of iterations to be performed given the
#' stopping criterion is not met beforehand for 'rf' method.
#' @param nTreeRf Number of trees to grow in each forest for 'rf' method.
#' @param compoundsAsNeighbors For KNN imputation. If TRUE, compounds will be 
#' used as neighbors rather than samples. Note that using compounds as 
#' neighbors is significantly slower than using samples.
#' @param normalizeMethod  Name of normalization method.
#' "ComBat" (only ComBat batch correction), "quantile" (only quantile 
#' normalization), "quantile + ComBat" (quantile with ComBat batch correction),
#' "median" (only median normalization), "median + ComBat" (median with ComBat
#' batch correction), "CRMN" (cross-contribution compensating multiple 
#' standard normalization), "RUV" (remove unwanted variation), "SVA" (surrogate 
#' variable analysis), or "none" to skip normalization.
#' @param nControl Number of controls to estimate/utilize (for CRMN and RUV).
#' @param controls Vector of control identifiers.  Leave blank for data driven
#' controls. Vector of column numbers from metafin dataset of that control (for
#' CRMN and RUV).
#' @param nComp Number of factors to use in CRMN algorithm. 
#' @param kRUV Number of factors to use in RUV algorithm.
#' @param covariatesOfInterest Sample variables used as covariates in
#' normalization algorithms (required for ComBat, CRMN, and SVA).
#' @param batch Name of the sample variable identifying batch.
#' @param transform  Select transformation to apply to data prior to 
#' normalization. Options are "log10", "log2", and "none".
#' @param replicate Name of sample variable specifying replicate. Must match an
#' element in `sampleVars` or a column in the column data of a 
#' `SummarizedExperiment`.
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
#' @param returnSummaryDetails Logical value specifying whether to return
#' details of replicate summarization.
#' @param returnToSE Logical value specifying whether to return as 
#' `SummarizedExperiment`
#' @param returnToDF Logical value specifying whether to return as data frame.
#' 
#' @examples 
#' # Load example data
#' data(msquant)
#' 
#' # Call function to tidy, summarize, filter, impute, and normalize data
#' peparedData <- msPrepare(msquant, cvMax = 0.50, minPropPresent = 1/3,
#'                          filterPercent = 0.8, imputeMethod = "halfmin",
#'                          normalizeMethod = "quantile",
#'                          compVars = c("mz", "rt"),
#'                          sampleVars = c("spike", "batch", "replicate", 
#'                                         "subject_id"),
#'                          colExtraText = "Neutral_Operator_Dif_Pos_",
#'                          separator = "_", missingValue = 1, 
#'                          returnToSE = FALSE)
#' 
#' @return A data frame or `SummarizedExperiment` with summarized technical 
#' replicates (if present), filtered compounds, missing values imputed, and 
#' transformed and normalized abundances. Default return type is set to match 
#' the data input but may be altered with the `returnToSE` or `returnToDF` 
#' arguments.
#' 
#' @export
msPrepare <- function(data, cvMax = 0.50, minPropPresent = 1/3,
                      filterPercent = 0.8, imputeMethod = c("halfmin", "bpca", 
                                                          "knn", "rf", "none"),
                      kKnn = 5, nPcs = 3, maxIterRf = 10, nTreeRf = 100,
                      compoundsAsNeighbors = FALSE,
                      normalizeMethod = c( "median", "ComBat", "quantile", 
                                          "quantile + ComBat",
                                          "median + ComBat", "CRMN", "RUV",
                                          "SVA", "none"),
                      nControl = 10, controls  = NULL, nComp = 2, kRUV = 3,
                      covariatesOfInterest = NULL, batch = NULL,
                      transform = c("log10", "log2", "none"),
                      replicate = "replicate", compVars = c("mz", "rt"),
                      sampleVars = c("subject_id"), colExtraText = NULL,
                      separator = NULL, missingValue = NA,
                      returnSummaryDetails = FALSE, returnToSE = FALSE,
                      returnToDF = FALSE) {
    
    imputeMethod <- match.arg(imputeMethod)
    normalizeMethod <- match.arg(normalizeMethod)
    transform <- match.arg(transform)
    
    
    if (replicate %in% sampleVars) {
        cat("Summarizing\n")
        data <- msSummarize(data, cvMax, minPropPresent, replicate, compVars,
                            sampleVars, colExtraText, separator, missingValue,
                            returnSummaryDetails)
        
        sampleVars <- sampleVars[sampleVars != replicate]
    }
    
    cat("Filtering\n")  
    data <- msFilter(data, filterPercent, compVars, sampleVars, colExtraText,
                     separator, missingValue = 0)
    
    if (imputeMethod != "none") {
        cat("Imputing\n")
        data <- msImpute(data, imputeMethod, kKnn, nPcs, maxIterRf, nTreeRf, 
                         compoundsAsNeighbors, compVars, sampleVars, 
                         colExtraText, separator, missingValue = 0, returnToSE, 
                         returnToDF)
    }
    
    if (normalizeMethod != "none") {
        cat("Normalizing\n")
        data <- msNormalize(data, normalizeMethod, nControl, controls, nComp,
                            kRUV, batch, covariatesOfInterest, transform,
                            compVars, sampleVars, colExtraText, separator,
                            returnToSE, returnToDF)
    }
    
    return(data)
}
#' Function for summarizing technical replicates.
#'
#' Reads data and summarizes technical replicates as the
#' mean of observations for compounds found in 2 or 3 replicates and with
#' coefficient of variation below specified level, or median for those found in
#' 3 replicates but with excessive coefficient of variation (CV). Compounds 
#' found in only 1 replicate are assigned as missing.
#' 
#' @param data Data set as either a data frame or `SummarizedExperiement`.
#' @param cvMax Decimal value from 0 to 1 representing the acceptable level of 
#' coefficient of variation between replicates.
#' @param minPropPresent  Decimal value from 0 to 1 representing the 
#' minimum proportion present to summarize with median or mean. Below this the 
#' compound will be set to 0.
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
#' @return A data frame or `SummarizedExperiment` containing abundance data 
#' with summarized technical replicates. Default return type is set to match the
#' data nput but may be altered with the `returnToSE` or `returnToDF` 
#' arguments. If `returnSummaryDetails` is selected, a list will be returned 
#' containing the summarized data and a separate tidy data frame with 
#' summarization details included for each compound/sample pair.
#' 
#' @examples
#' # Read in data file
#' data(msquant)
#' 
#' # Summarize technical replicates
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
#' @importFrom methods is
#' 
#' 
#' @export
msSummarize <- function(data,
                        cvMax = 0.50,
                        minPropPresent = 1/3,
                        replicate = "replicate",
                        compVars = c("mz", "rt"),
                        sampleVars = c("subject_id"),
                        colExtraText = NULL,
                        separator = NULL,
                        missingValue = NA,
                        returnSummaryDetails = FALSE,
                        returnToSE = FALSE,
                        returnToDF = FALSE) {
    
    if (is(data, "SummarizedExperiment")) {
        return <- .seSummarize(data, cvMax, minPropPresent, replicate,
                               missingValue, returnSummaryDetails, returnToDF)
    } else if (is(data, "data.frame")) {
        return <- .dfSummarize(data, cvMax, minPropPresent, replicate,
                               compVars, sampleVars, colExtraText, separator,
                               missingValue, returnSummaryDetails, returnToSE)
    } else {
        stop("'data' must be a data frame or SummarizedExperiment")
    }
    
}

#' @importFrom tidyr as_tibble
#' @importFrom dplyr bind_cols
#' @importFrom S4Vectors metadata
.seSummarize <- function(SE, cvMax, minPropPresent, replicate, missingValue, 
                         returnSummaryDetails, returnToDF) {
    
    ## Check for presence of 'replicate' in column data
    if(replicate %notin% colnames(colData(SE))) {
        stop("'replicate' must correspond to a column name in 'colData'")
    }
    
    ## Store existing metadata, get compVars and sampleVars (w/o replicate)
    metaData <- metadata(SE)
    compVars <- colnames(rowData(SE))
    sampleVars <- colnames(colData(SE)) 
    sampleVars <- sampleVars[sampleVars != replicate]
    
    ## Tidy Data
    tidyData <- .msTidy(data = SE, missingValue = missingValue, setMissing = NA)
    
    ## Summarize Replicates
    summarizedData <- .tidySummarize(tidyData = tidyData, compVars = compVars,
                                     sampleVars = sampleVars, 
                                     replicate = replicate,
                                     cvMax = cvMax, 
                                     minPropPresent = minPropPresent,
                                     returnSummaryDetails = 
                                         returnSummaryDetails)
    
    ## Return data according to user preference
    rtn <- .returnSummarized(data = summarizedData, compVars = compVars,
                             sampleVars = sampleVars, metaData = metaData,
                             returnSummaryDetail = returnSummaryDetails,
                             toSE = !returnToDF)
}

.dfSummarize <- function(data, cvMax, minPropPresent, replicate, compVars, 
                         sampleVars, colExtraText, separator, missingValue, 
                         returnSummaryDetails, returnToSE) {
    
    tidyData <- .msTidy(data, compVars, sampleVars, colExtraText, separator,
                        missingValue, setMissing = NA) 
    
    sampleVars <- sampleVars[sampleVars != replicate]
    
    summarizedData <- .tidySummarize(tidyData, compVars, sampleVars, replicate, 
                                     cvMax, minPropPresent, 
                                     returnSummaryDetails)
    
    rtn <- .returnSummarized(data = summarizedData, compVars = compVars,
                             sampleVars = sampleVars, metaData = NULL,
                             returnSummaryDetail = returnSummaryDetails,
                             toSE = returnToSE)
}



#' @importFrom dplyr mutate_at
#' @importFrom dplyr group_by 
#' @importFrom dplyr summarise 
#' @importFrom dplyr mutate 
#' @importFrom dplyr ungroup 
#' @importFrom dplyr case_when
#' @importFrom dplyr distinct
#' @importFrom dplyr filter
#' @importFrom dplyr vars
#' @importFrom tibble as_tibble
#' @importFrom tidyr replace_na
#' @importFrom rlang .data
#' @importFrom rlang UQ
#' @importFrom rlang sym
#' @importFrom rlang syms
#' @importFrom rlang as_string
#' @importFrom rlang !!!
#' @importFrom stats median
#' @importFrom stats sd
.tidySummarize <- function(tidyData, compVars, sampleVars, replicate, cvMax,
                           minPropPresent, returnSummaryDetails) {
    replicateCount <- length(unique(tidyData[[replicate]]))
    
    sampleVars <- syms(sampleVars)
    compVars <- syms(compVars)
    
    if (length(sampleVars) == 0) {
        tidyData <- group_by(tidyData, `!!!`(compVars))
    } else {
        tidyData <- group_by(tidyData, `!!!`(sampleVars), `!!!`(compVars))
    }
    
    ## Determine appropriate summary measure, apply summary
    tidyData <- summarise(tidyData,
                          nPresent = sum(!is.na(.data$abundance)),
                          propPresent = UQ(sym("nPresent")) / replicateCount,
                          meanAbundance = ifelse(.data$nPresent == 0, NA,
                                                 mean(.data$abundance,
                                                      na.rm = TRUE)),
                          sdAbundance = ifelse(.data$nPresent == 0, NA,
                                               sd(.data$abundance,
                                                  na.rm = TRUE)),
                          medianAbundance = ifelse(.data$nPresent == 0, NA,
                                                   median(.data$abundance,
                                                          na.rm = TRUE))) %>%
        mutate(cvAbundance = ifelse(is.na(.data$meanAbundance), NA, 
                                    .data$sdAbundance / 
                                        .data$meanAbundance)) %>%
        ungroup() %>%
        mutate(summaryMeasure = .selectSummaryMeasure(.data$nPresent,
                                                      .data$cvAbundance,
                                                      replicateCount,
                                                      minPropPresent,
                                                      cvMax)) %>%
        mutate(abundance = case_when(.data$summaryMeasure == "median" ~
                                                 .data$medianAbundance,
                                             .data$summaryMeasure == "mean" ~
                                                 .data$meanAbundance,
                                             TRUE ~ 0)) %>%
        mutate_at(vars("summaryMeasure", `!!!`(sampleVars)), factor)
    
    summaryData <- tidyData
    tidyData <- select(tidyData, c(!!!sampleVars, !!!compVars), .data$abundance)

    ifelse(returnSummaryDetails, 
           return(list("data" = tidyData, "summaryDetails" = summaryData)), 
           return(tidyData))
}

#' @importFrom dplyr case_when
.selectSummaryMeasure <- function(nPresent,
                                  cvAbundance,
                                  nReplicates,
                                  minProportionPresent,
                                  cvMax) {
    
    case_when((nPresent / nReplicates) <= minProportionPresent ~ 
                  "none: proportion present <= minProportionPresent",
              cvAbundance > cvMax & (nReplicates == 3 & nPresent == 2) ~ 
                  "none: cv > cvMax & 2 present",
              cvAbundance > cvMax & (nPresent == nReplicates) ~ "median",
              TRUE ~ "mean")
}

# Internal function to return data in right format and with requested metadata
.returnSummarized <- function(data, compVars, sampleVars, metaData, 
                              returnSummaryDetail, toSE) {
    if (returnSummaryDetail & toSE) {
        returnSE <- .tidyReturn(tidyData = data$data, compVars = compVars,
                                metaData = metaData, sampleVars = sampleVars,
                                toSE = toSE)
        
        metadata(returnSE)$summaryDetails <- data$summaryDetails
        return(returnSE)
    } else if (returnSummaryDetail) {
        summaryData <- .tidyReturn(tidyData = data$data, compVars = compVars, 
                                   sampleVars = sampleVars, toSE = toSE)
        summaryDetails <- data$summaryDetails
        return(list("data" = summaryData, "summaryDetails" = summaryDetails))
    } else {
        summaryData <- .tidyReturn(tidyData = data, compVars = compVars,
                                   sampleVars = sampleVars, toSE = toSE)
        return(summaryData)
    }
}

`%notin%` <- Negate(`%in%`)

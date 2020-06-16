msSummarize <- function(data,
                        cvMax = 0.50,
                        minPropPresent = 1/3,
                        replicate = "replicate",
                        compVars = c("mz", "rt"),
                        techVars = c("subject_id"),
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
                               compVars, techVars, colExtraText, separator,
                               missingValue, returnSummaryDetails, returnToSE)
    } else {
        stop("'data' must be a data frame or SummarizedExperiment")
    }
    
}

#' Internal function to summarize replicates in a SummarizedExperiment
#' @importFrom tidyr as_tibble
#' @importFrom dplyr bind_cols
#' @importFrom S4Vectors metadata
.seSummarize <- function(SE, cvMax, minPropPresent, replicate, missingValue, 
                         returnSummaryDetails, returnToDF) {
    
    ## Check for presence of 'replicate' in column data
    if(replicate %notin% colnames(colData(SE))) {
        stop("'replicate' must correspond to a column name in 'colData'")
    }
    
    ## Store existing metadata, get compVars and techVars (w/o replicate)
    metaData <- metadata(SE)
    compVars <- colnames(rowData(SE))
    techVars <- colnames(colData(SE)) 
    techVars <- techVars[techVars != replicate]
    
    ## Tidy Data
    tidyData <- msTidy(data = SE, missingValue = missingValue)
    
    ## Summarize Replicates
    summarizedData <- .tidySummarize(tidyData = tidyData, compVars = compVars,
                                     techVars = techVars, replicate = replicate,
                                     cvMax = cvMax, 
                                     minPropPresent = minPropPresent,
                                     returnSummaryDetails = 
                                         returnSummaryDetails)
    
    ## Return data according to user preference
    rtn <- .returnSummarized(data = summarizedData, compVars = compVars,
                             techVars = techVars, metaData = metaData,
                             returnSummaryDetail = returnSummaryDetails,
                             toSE = !returnToDF)
}

.dfSummarize <- function(data, cvMax, minPropPresent, replicate, compVars, 
                         techVars, colExtraText, separator, missingValue, 
                         returnSummaryDetails, returnToSE) {
    
    tidyData <- msTidy(data, compVars, techVars, colExtraText, separator,
                       missingValue) 
    
    techVars <- techVars[techVars != replicate]
    
    summarizedData <- .tidySummarize(tidyData, compVars, techVars, replicate, 
                                     cvMax, minPropPresent, 
                                     returnSummaryDetails)
    
    rtn <- .returnSummarized(data = summarizedData, compVars = compVars,
                             techVars = techVars, metaData = NULL,
                             returnSummaryDetail = returnSummaryDetails,
                             toSE = returnToSE)
}

.tidySummarize <- function(tidyData, compVars, techVars, replicate, cvMax,
                           minPropPresent, returnSummaryDetails) {
    replicateCount <- length(unique(tidyData[[replicate]]))
    
    techVars <- syms(techVars)
    compVars <- syms(compVars)
    
    if (length(techVars) == 0) {
        tidyData <- group_by(tidyData, `!!!`(compVars))
    } else {
        tidyData <- group_by(tidyData, `!!!`(techVars), `!!!`(compVars))
    }
    
    ## Determine approprate summary measure, apply summary
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
        mutate_at(vars(subject_id, "summaryMeasure", `!!!`(techVars)), factor)
    
    summaryData <- tidyData
    tidyData <- select(tidyData, c(!!!techVars, !!!compVars, abundance))

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

#' Internal function to return data in right format and with requested metadata
.returnSummarized <- function(data, compVars, techVars, metaData, 
                              returnSummaryDetail, toSE) {
    if (returnSummaryDetail & toSE) {
        returnSE <- .msReturn(tidyData = data$data, compVars = compVars, 
                              metaData = metaData, techVars = techVars, 
                              toSE = toSE)
        
        metadata(returnSE)$summaryDetails <- data$summaryDetails
        return(returnSE)
    } else if (returnSummaryDetail) {
        summaryData <- .msReturn(tidyData = data$data, compVars = compVars, 
                                 techVars = techVars, toSE = toSE)
        summaryDetails <- data$summaryDetails
        return(list("data" = summaryData, "summaryDetails" = summaryDetails))
    } else {
        summaryData <- .msReturn(tidyData = data, compVars = compVars,  
                                 techVars = techVars, toSE = toSE)
        return(summaryData)
    }
}
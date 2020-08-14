msFilter <- function(data,
                     filterPercent = 0.8,
                     compVars = c("mz", "rt"),
                     sampleVars = c("subject_id"),
                     colExtraText = NULL,
                     separator = NULL,
                     missingValue = NA,
                     returnToSE = FALSE,
                     returnToDF = FALSE) {
    
    if (is(data, "SummarizedExperiment")) {
        rtn <- .seFilter(data, filterPercent, missingValue, returnToDF)
    } else if (is(data, "data.frame")) {
        rtn <- .dfFilter(data, filterPercent, compVars, 
                         sampleVars, colExtraText, separator, missingValue, 
                         returnToSE)
    } else {
        stop("'data' must be a data frame or SummarizedExperiment")
    }
    
}

.seFilter <- function(SE, filterPercent, missingValue, returnToDF) {
    
    ## Store existing metadata, get compVars and sampleVars 
    metaData <- metadata(SE)
    compVars <- colnames(rowData(SE))
    sampleVars <- colnames(colData(SE))
    
    ## Tidy Data
    tidyData <- msTidy(data = SE, missingValue = missingValue)
    
    ## Filter Data
    filteredData <- .tidyFilter(tidyData, compVars, filterPercent)
    
    ## Return Data
    rtn <- .msReturn(filteredData, compVars, sampleVars, metaData, 
                     toSE = !returnToDF)
}

.dfFilter <- function(data, filterPercent, compVars, sampleVars, colExtraText,
                      separator, missingValue, returnToSE) {
    ## Tidy Data
    tidyData <- msTidy(data = data, compVars = compVars, 
                       sampleVars = sampleVars, colExtraText = colExtraText, 
                       separator = separator, missingValue = missingValue)
    
    ## Filter Data
    filteredData <- .tidyFilter(tidyData, compVars, filterPercent)
    
    ## Return Data
    rtn <- .msReturn(filteredData, compVars, sampleVars, toSE = returnToSE)
}

.tidyFilter <- function(tidyData, compVars, filterPercent) {
    
    compVarsSyms <- syms(compVars)
    
    filterStatus <- group_by(tidyData, `!!!`(compVarsSyms)) %>%
        summarise(percPresent = sum(.data$abundance != 0) / n()) %>%
        mutate(keep = .data$percPresent >= filterPercent) %>% 
        ungroup
    
    filteredData <- full_join(tidyData, select(filterStatus, 
                                               `!!!`(compVarsSyms), .data$keep),
                              by = compVars) %>% 
        filter(.data$keep) %>% 
        select(-.data$keep)
}
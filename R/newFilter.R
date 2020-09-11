#' Function for filtering abundance data set.
#'
#' Filters compounds to those found in specified proportion of subjects.
#' 
#' @param data Data set as either a data frame or `SummarizedExperiement`.
#' @param filterPercent Decimal value indicating filtration threshold. 
#' Compounds which are present in fewer samples than the specified proportion 
#' will be removed. 
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
#' @return A data frame or `SummarizedExperiment` with filtered abundance data.
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
#'
#' # Load example data set, summarize replicates
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
#' # Filter the dataset using a 80% filter rate
#' filteredDF <- msFilter(summarizedDF,
#'                        filterPercent = 0.8,
#'                        compVars = c("mz", "rt"),
#'                        sampleVars = c("spike", "batch", "subject_id"),
#'                        separator = "_")
#'
#'
#' @export
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

#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
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
    rtn <- .tidyReturn(filteredData, compVars, sampleVars, metaData,
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
    rtn <- .tidyReturn(filteredData, compVars, sampleVars, toSE = returnToSE)
}

#' @importFrom dplyr mutate 
#' @importFrom dplyr group_by 
#' @importFrom dplyr summarise 
#' @importFrom dplyr ungroup 
#' @importFrom dplyr filter
#' @importFrom dplyr full_join
#' @importFrom dplyr n
#' @importFrom rlang .data
#' @importFrom rlang syms
#' @importFrom rlang !!!
#' @importFrom magrittr %>%
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
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr separate
#' @importFrom dplyr select
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment SummarizedExperiment
## Returns tidy data to SE or DF
.tidyReturn <- function(tidyData, compVars, sampleVars, metaData = NULL, toSE) {
    
    ## Return data to wide format data frame
    rtn <- tidyData %>%
        pivot_wider(id_cols = compVars,
                    names_from = sampleVars,
                    values_from = "abundance")
    
    ## If selected, convert data to SummarizedExperiment
    if(toSE) {
        ## Get row data
        seRowData <- select(rtn, compVars)

        ## Get assay data
        seAssay <- select(rtn, -compVars) %>%
            as.matrix()

        ## Get column data
        seColumnData <- tibble("samples" = colnames(seAssay)) %>%
            separate(col = "samples", into = sampleVars, sep = "_")

        ## Build SummarizedExperiment
        rtn <- SummarizedExperiment(assays = list(abundance = seAssay),
                                    colData = seColumnData,
                                    rowData = seRowData,
                                    metadata = metaData)
    }
    
    return(rtn)
}

#' @importFrom tidyr pivot_wider
#' @importFrom tidyr separate
#' @importFrom dplyr select
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
## Converts DF to SE
.dfToSE <- function(DF, compVars, sampleVars, separator, colExtraText = NULL, 
                    metaData = NULL) {
    
    ## Get column data
    sampleCols <- select(DF, -compVars)
    seColumnData <- tibble("samples" = colnames(sampleCols))
    
    ## Remove colExtraText text if present
    if(!is.null(colExtraText)){
        seColumnData <- mutate(seColumnData, "samples" = 
                                   str_replace_all(.data$samples,
                                                   colExtraText, ""))
    }
        
    seColumnData <- separate(seColumnData, col = "samples", into = sampleVars, 
                             sep = separator)
    
    ## Get row data
    seRowData <- select(DF, compVars)
    
    ## Get abundance data
    seAssay <-  select(DF, -compVars) %>%
        as.matrix()
    
    ## Build SummarizedExperiment
    rtn <- SummarizedExperiment(assays = list(abundance = seAssay),
                                colData = seColumnData,
                                rowData = seRowData,
                                metadata = metaData)
}

## Converts SE to DF
## Note: Currently metadata is lost in conversion
.seToDF <- function(SE, colExtraText = NULL, seperator = "_") {
    
    ## Get row data, compound variables, and technical variables
    rowData <- as_tibble(rowData(SE))
    colData <- as_tibble(colData(SE))
    abundanceData <- as.data.frame(assay(SE))
    
    if (!is.null(colExtraText)) {
        colData <- tibble::add_column(colData, 
                                      "colExtraText" = 
                                          rep(colExtraText, nrow(colData)),
                                      .before = 1)
    }
    
    columnNames <- apply(colData, 1, paste, collapse=seperator)
    colnames(abundanceData) <- columnNames
     
    cbind(rowData, abundanceData)
    
}
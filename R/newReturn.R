.msReturn <- function(tidyData, compVars, sampleVars, metaData = NULL, toSE) {
    
    ## Return data to wide format data frame
    rtn <- tidyData %>%
        pivot_wider(id_cols = compVars,
                    names_from = sampleVars,
                    values_from = "abundance")
    
    ## If selected, convert data to SummarizedExperiment
    if(toSE) {
        # Get row data
        seRowData <- select(rtn, compVars)

        # Get assay data
        seAssay <- select(rtn, -compVars) %>%
            as.matrix()

        # Get column data
        seColumnData <- tibble("samples" = colnames(seAssay)) %>%
            separate(col = "samples", into = sampleVars, sep = "_")

        # Build SummarizedExperiment
        rtn <- SummarizedExperiment(assays = list(abundance = seAssay),
                                    colData = seColumnData,
                                    rowData = seRowData,
                                    metadata = metaData)
    }
    
    return(rtn)
}
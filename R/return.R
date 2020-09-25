# 
# ms_return <- function(msprepObject, returnToSE = FALSE){
# 
#     # Validate inputs
#     stopifnot(class(msprepObject) == "msprep")
#     
#     # Return data to wide format data frame
#     msprepObject$data <- msprepObject$data %>%
#         pivot_wider(id_cols = met_vars(msprepObject),
#                     names_from = col_order(msprepObject),
#                     values_from = "abundance_summary")
# 
#     # If selected, convert data to SummarizedExperiment
#     if(returnToSE) {
#         msprepObject$data <- .msprepToSE(msprepObject)
#     }
#   
#     return(msprepObject)
# }
# 
# 
# 
# .msprepToSE <- function(msprepObject) {
#     
#     # Get row data
#     seRowData <- msprepObject$data %>%
#         select(met_vars(msprepObject))
#     
#     # Get assay data
#     seAssay <- msprepObject$data %>%
#         select(-met_vars(msprepObject)) %>%
#         as.matrix()
#     
#     # Get column data
#     seColumnData <- tibble("samples" = colnames(seAssay)) %>%
#         separate(col = "samples", into = col_order(msprepObject), sep = "_")
#     
#     # Build SummarizedExperiment
#     rtn <- SummarizedExperiment(assays = list(abundance = seAssay), 
#                                 colData = seColumnData, 
#                                 rowData = seRowData)
# }
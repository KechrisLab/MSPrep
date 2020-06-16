#' Function to return data prepared by MSPrep pipeline to compound by matrix
#' data frame or to a SummarizedExperiment
#' 
#' @param msprep_obj An object of class `msprep`
#' 
#' @importFrom tidyr pivot_wider
#' 
#' @examples 
#' # Load, tidy, summarize, filter, impute, and normalize example dataset
#' data(msquant)
#' 
#' tidied_data <- ms_tidy(msquant, mz = "mz", rt = "rt",
#'                        col_extra_txt = "Neutral_Operator_Dif_Pos_",
#'                        separator = "_", 
#'                        col_names = c("spike", "batch", "replicate", 
#'                                      "subject_id"))
#' 
#' summarized_data <- ms_summarize(tidied_data, 
#'                                 mz = "mz", 
#'                                 rt = "rt", 
#'                                 replicate = "replicate", 
#'                                 batch = "batch", 
#'                                 groupingvars = "spike", 
#'                                 subject_id = "subject_id", 
#'                                 cvmax = 0.50, 
#'                                 min_proportion_present = 1/3, 
#'                                 missing_val = 1)
#' 
#' filtered_data <- ms_filter(summarized_data, 
#'                            filter_percent =  0.80)
#' 
#' imputed_data <- ms_impute(filtered_data, 
#'                           imputeMethod = "halfmin")
#' 
#' normalized_data <- ms_normalize(imputed_data,
#'                                 normalizeMethod = "median")
#' 
#' # Return dataset to original form
#' returned_data <- ms_return(normalized_data)
#' 
#' @return An msprep object with data formatted with metabolites as rows and
#' and samples as columns, or with data as a SummarizedExperiment
#' 
#' @export
ms_return <- function(msprepObject, returnToSE = FALSE){

    # Validate inputs
    stopifnot(class(msprepObject) == "msprep")
    
    # Return data to wide format data frame
    msprepObject$data <- msprepObject$data %>%
        pivot_wider(id_cols = met_vars(msprepObject),
                    names_from = col_order(msprepObject),
                    values_from = "abundance_summary")

    # If selected, convert data to SummarizedExperiment
    if(returnToSE) {
        msprepObject$data <- .msprepToSE(msprepObject)
    }
  
    return(msprepObject)
}


#' Function to convert MSPrep object to SE
#' 
#' @importFrom dplyr select
#' @importFrom tidyr separate
#' @importFrom tibble tibble
#' @import SummarizedExperiment
.msprepToSE <- function(msprepObject) {
    
    # Get row data
    seRowData <- msprepObject$data %>%
        select(met_vars(msprepObject))
    
    # Get assay data
    seAssay <- msprepObject$data %>%
        select(-met_vars(msprepObject)) %>%
        as.matrix()
    
    # Get column data
    seColumnData <- tibble("samples" = colnames(seAssay)) %>%
        separate(col = "samples", into = col_order(msprepObject), sep = "_")
    
    # Build SummarizedExperiment
    rtn <- SummarizedExperiment(assays = list(abundance = seAssay), 
                                colData = seColumnData, 
                                rowData = seRowData)
}
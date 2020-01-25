#' Function to return data prepared by MSPrep pipeline to compound by sample format
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
#'                        col_names = c("spike", "batch", "replicate", "subject_id"))
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
#' @export

ms_return <- function(msprep_obj){
  
  # Validate inputs
  stopifnot(class(msprep_obj) == "msprep")
  #stopifnot(stage(msprep_obj) %in% c("normalized"))
  
  # Get column order from object attribute
  col_order <- col_order(msprep_obj)
  
  # Get met_id columns from object attribute
  met_vars <- met_vars(msprep_obj)
  
  # Return data to wide format
  msprep_obj$data <- msprep_obj$data %>% pivot_wider(id_cols = met_vars, 
                                                        names_from = col_order,
                                                        values_from = "abundance_summary")
  
  return(msprep_obj)
}
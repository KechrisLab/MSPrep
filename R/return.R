#' Function to return data prepared by MSPrep pipeline to compound by sample format
#' 
#' @param msprep_obj An object of class `msprep`
#' 
#' @importFrom tidyr pivot_wider
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
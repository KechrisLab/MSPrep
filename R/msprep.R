# #' 
# #' The primary class object created by the MSPrep pipeline functions
# #'
# #'
# #' @param data A tidy dataframe of quantification data.
# #' @param replicate_info Placeholder
# #' @param medians Placeholder
# #' @param replicate_count Placeholder
# #' @param cvmax Placeholder
# #' @param min_proportion_present Placeholder
# #' @param stage Attribute holding the stage of the msprep object in the
# #' pipeline.  Valid values can be viewed with the `stages()` function.
# #' @return An `msprep` object
# as.msprep <- function(data, 
#                       replicate_info = NULL, 
#                       medians = NULL,
#                       replicate_count = NULL, 
#                       cvmax = NULL,
#                       min_proportion_present = NULL,
#                       stage = NULL) {
# 
#   # Create return object & return
#   structure(list("data" = summary_data,
#                  replicate_info  = replicate_info,
#                  medians         = medians),
#             replicate_count = replicate_count,
#             cvmax           = cvmax,
#             min_proportion_present = min_proportion_present,
#             class = "msprep",
#             stage = "prepared")
# 
# }
# 
# 

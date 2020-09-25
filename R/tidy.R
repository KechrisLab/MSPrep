# 
# ms_tidy <- function(quantification_data,
#                     met_id = NULL,
#                     mz = NULL, 
#                     rt = NULL,
#                     col_extra_txt = NULL,
#                     col_names = c("subject_id"),
#                     separator = NULL) {
#   
#   # Check that at either met_id or both mz and rt are included
#     if(is.null(met_id) & (is.null(mz) | is.null(rt))){
#         stop("Must include met_id or both mz and rt")
#     }
#     
#     # Store whatever metabolite id args are present
#     met_vars <- syms(c(met_id, mz, rt))
#     
#     # Gather data to long format (adds id/varnames as column), ensure mz and rt
#     # are numeric if present
#     rtn <- as_tibble(quantification_data) %>%
#         gather(key = "id_col", value = "abundance", -c(!!!met_vars))
#     
#     
#     # Ensure mz and rt are numeric if present
#     if (!is.null(mz) & !is.null(rt)) {
#         rtn <- mutate_at(rtn, vars(mz, rt), as.numeric)
#     }
#     
#     # Remove col_extra_txt text if present
#     if (!is.null(col_extra_txt)) {
#         rtn <- mutate(rtn, id_col = str_replace_all(.data$id_col,
#                                                     col_extra_txt, ""))
#     }
#   
#     # If only one column name, rename id_col appropriately. Otherwise split
#     # 'id_col' into new variable columns designated by col_names 
#     if (length(col_names) == 1) {
#         colnames(rtn)[colnames(rtn) == "id_col"] <- col_names[1]
#     }
#     else if (!is.null(separator)) {
#         rtn <- separate(rtn, .data$id_col, sep = separator, into = col_names)
#     }
#     else {
#         stop("Must include 'separator' if multiple 'col_names'")
#     }
# 
#     return(rtn)
# }
# 

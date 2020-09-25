# ms_filter <- function(msprep_obj, filter_percent = 0.8) {
#   
#   stopifnot(class(msprep_obj) == "msprep")
#   stopifnot(stage(msprep_obj) == "summarized")
#   
#   met_vars <- met_vars(msprep_obj)
#   met_syms <- syms(met_vars)
#   
#   filter_status <- group_by(msprep_obj$data, `!!!`(met_syms)) %>% 
#     summarise(perc_present = sum(.data$abundance_summary != 0)/n()) %>% 
#     mutate(keep = .data$perc_present >= filter_percent) %>% 
#     ungroup
#   
#   filtereddata <- full_join(msprep_obj$data, select(filter_status, 
#                                                     `!!!`(met_syms), 
#                                                     .data$keep), 
#                                                     by = met_vars) %>% 
#   filter(.data$keep) %>% 
#     select(-.data$keep)
#   
#   msprep_obj$data <- filtereddata
#   attr(msprep_obj, "filter_status") <- filter_status
#   attr(msprep_obj, "filter_percent") <- filter_percent
#   stage(msprep_obj) <- "filtered"
#   
#   return(msprep_obj)
#   
# }
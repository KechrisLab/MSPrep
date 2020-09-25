# ms_prepare <- function(quantification_data,
#                        met_id = NULL,
#                        mz = NULL,
#                        rt = NULL,
#                        col_extra_txt = NULL,
#                        col_names = c("subject_id"),
#                        separator = NULL,
#                        abundance = "abundance",
#                        subject_id = "subject_id",
#                        replicate = NULL,
#                        batch = NULL,
#                        groupingvars = NULL,
#                        cvmax = 0.50,
#                        missing_val = 1,
#                        min_proportion_present = 1/3,
#                        filter_percent = .8,
#                        imputeMethod = c("halfmin",
#                                         "bpca", 
#                                         "knn", 
#                                         "none"),
#                        k_knn = 5, 
#                        n_pcs = 3, 
#                        compoundsAsNeighbors = FALSE,
#                        normalizeMethod = c("ComBat",
#                                            "quantile",
#                                            "quantile + ComBat",
#                                            "median",
#                                            "median + ComBat",
#                                            "CRMN",
#                                            "RUV",
#                                            "SVA",
#                                            "none"),
#                        n_control = 10,
#                        controls  = NULL,
#                        n_comp    = 2,
#                        k_ruv     = 3,
#                        transform = c("log10",
#                                      "log2",
#                                      "none")) {
#   
#   cat("Tidying\n")
#   data <- ms_tidy(quantification_data,
#                   met_id = met_id,
#                   mz = mz,
#                   rt = rt,
#                   col_extra_txt = col_extra_txt,
#                   col_names = col_names,
#                   separator = separator)
#   
#   cat("Summarizing\n")
#   data <- ms_summarize(data,
#                        abundance     = "abundance",
#                        met_id        = met_id,
#                        mz            = mz,
#                        rt            = rt,
#                        subject_id    = "subject_id",
#                        replicate     = replicate,
#                        batch         = batch,
#                        groupingvars = groupingvars,
#                        cvmax         = cvmax,
#                        missing_val   = missing_val,
#                        min_proportion_present = min_proportion_present)
#   
#   cat("Filtering\n")  
#   data <- ms_filter(data,
#                     filter_percent = filter_percent)
#   
#   if(imputeMethod != "none") {
#     cat("Imputing\n")
#     data <- ms_impute(data,
#                       imputeMethod = imputeMethod,
#                       k_knn = k_knn, 
#                       n_pcs = n_pcs, 
#                       compoundsAsNeighbors = compoundsAsNeighbors)
#   }
#   
#   if(normalizeMethod != "none"){
#     cat("Normalizing\n")
#     data <- ms_normalize(data,
#                          normalizeMethod = normalizeMethod,
#                          n_control = n_control,
#                          controls  = controls,
#                          n_comp    = n_comp,
#                          k_ruv     = k_ruv,
#                          transform = transform)
#   }
#   
#   cat("Returning")
#   data <- ms_return(data)
#           
#   return(data)
# }
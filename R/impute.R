# ms_impute <- function(msprep_obj,
#                       imputeMethod = c("halfmin", "bpca", "knn"),
#                       k_knn = 5, 
#                       n_pcs = 3, 
#                       compoundsAsNeighbors = FALSE) {
# 
#   # Validate inputs
#   stopifnot(class(msprep_obj) == "msprep")
#   stopifnot(stage(msprep_obj) %in% c("filtered", "normalized"))
#   imputeMethod <- match.arg(imputeMethod) # requires 1 argument from vec in function arg
# 
#   # Prep data - replace 0's with NA's -- for minval and bpca() (all methods?)
#   data <- mutate_at(msprep_obj$data, vars("abundance_summary"), 
#                     replace_missing, 0)
#   grp <- grouping_vars(msprep_obj)
#   batch <- batch_var(msprep_obj)
#   met_vars <- met_vars(msprep_obj)
#   
#   # Store number of missing values
#   missingCount <- sum(is.na(data$abundance_summary))
# 
#   # Impute data
#   data <-
#     switch(imputeMethod,
#            "halfmin" = impute_halfmin(data, grp, batch, met_vars),
#            "bpca"    = impute_bpca(data, grp, batch, met_vars, n_pcs),
#            "knn"     = impute_knn(data, grp, batch, met_vars, k_knn, compoundsAsNeighbors),
#            stop("Invalid impute method - you should never see this warning."))
# 
#   # Prep output object
#   msprep_obj$data  <- data
#   attr(msprep_obj, "impute_method") <- imputeMethod
#   attr(msprep_obj, "missing_count") <- missingCount
#   stage(msprep_obj) <- "imputed"
# 
#   # ...and:
#   return(msprep_obj)
# 
# }
# 
# # Old return
# #list(minval = as.data.frame(minval), withzero = as.data.frame(metafin), 
# #     bpca = as.data.frame(bpca), count = as.data.frame(count))
# 
# 
# 
# 
# impute_halfmin <- function(data, groupingvars, batch, met_vars) {
# 
#   sym_mz <- sym("mz")
#   sym_rt <- sym("rt")
#   sym_met_vars <- syms(met_vars)
# 
#   # NOTE: other imputation methods operate across all batchs/groups for a given
#   # mz_rt
#   #   grp  <- syms(c(groupingvars, batch))
#   # data <- group_by(data, `!!`(sym_mz), `!!!`(sym_met_vars))
#   data <- group_by(data, `!!!`(sym_met_vars))
# 
#   halfmin <- function(x) {
#     ifelse(is.na(x), min(x, na.rm = TRUE)/2, x)
#   }
# 
#   data <- mutate_at(data, vars("abundance_summary"), halfmin)
#   data <- ungroup(data)
# 
#   return(data)
# 
# }
# 
# 
# 
# 
# impute_bpca <- function(data, groupingvars, batch, met_vars, n_pcs = 3) {
# 
#   # 1. Bayesian pca imputation
#   data <- data_to_wide_matrix(data, groupingvars, batch, met_vars) 
#   data <- pca(data, nPcs = n_pcs, method = "bpca")
#   data <- completeObs(data) # extract imputed dataset
#   data <- wide_matrix_to_data(data, groupingvars, batch, met_vars)
#   data <- halfmin_if_any_negative(data, groupingvars, batch, met_vars)
# 
#   return(data)
# 
# }
# 
# 
# 
# impute_knn <- function(data, groupingvars, batch, met_vars, k_knn = 5, compoundsAsNeighbors) {
#   
#   data <- data_to_wide_matrix(data, groupingvars, batch, met_vars) 
#   rwnm <- rownames(data)
#   cnm <- colnames(data)
#   
#   # Perform kNN imputation using compounds or samples as neighbors
#   if (compoundsAsNeighbors == TRUE) {
#     data <- VIM::kNN(as.data.frame(data), k = k_knn, imp_var = FALSE)
#   }
#   else {
#     ### transpose
#     data <- VIM::kNN(as.data.frame(t(data)), k = k_knn, imp_var = FALSE)
#     
#     ## transpose again (to 'untranspose' basically)
#     data <- t(data)
#   }
#   
#   ### set column names same as prior to knn function
#   colnames(data) <- cnm
#   
#   rownames(data) <- rwnm
#   data <- wide_matrix_to_data(data, groupingvars, batch, met_vars)
#   data <- halfmin_if_any_negative(data, groupingvars, batch, met_vars)
# 
#   return(data)
# 
# }
# 
# 
# 
# halfmin_if_any_negative <-  function(data, groupingvars, batch, met_vars) {
# 
#   num_neg <- sum(data$abundance_summary < 0)
#   if (any(num_neg)) {
#     message("Found ", num_neg, " negative imputed values using KNN, reverting to half-min imputation for these values")
#     data$abundance_summary <- setmissing_negative_vals(data$abundance_summary)
#     data <- impute_halfmin(data, groupingvars, batch, met_vars)
#   }
# 
#   return(data)
# 
# }
# 
# truncate_negative_vals   <- function(var) ifelse(var < 0, 0, var)
# setmissing_negative_vals <- function(var) ifelse(var < 0, NA, var)
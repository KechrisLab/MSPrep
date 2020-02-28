#' Function for imputing missing values in data.
#'
#' Replaces missing values with non-zero estimates calculated using a
#' selected method.
#'
#' @param msprep_obj Filtered MSPrep object.
#' @param imputeMethod String specifying imputation method choice. 
#' Options are "halfmin" (half the minimum value), "bpca" (Bayesian PCA), 
#' and "knn" (k-nearest neighbors)
#' @param k_knn Number of clusters for 'knn' method.
#' @param n_pcs Number of  principle components used for re-estimation for 
#' 'bpca' method.
#' @param compoundsAsNeighbors For KNN imputation. If TRUE, compounds will be 
#' used as neighbors rather than samples. Note that using compounds as 
#' neighbors is significantly slower than using samples.
#' @return An `msprep` object with missing data imputed.
#' @references 
#'   Oba, S.et al.(2003) A Bayesian missing value estimation for gene
#'   expression profile data. Bioinformatics, 19, 2088-2096
#'
#'   Stacklies, W.et al.(2007) pcaMethods A bioconductor package providing
#'   PCA methods for incomplete data. Bioinformatics, 23, 1164-1167.
#' @examples
#' # Load, tidy, summarize, and filter example dataset
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
#' # Impute dataset using 3 possible options
#' imputed_data_hm <- ms_impute(filtered_data, 
#'                              imputeMethod = "halfmin")
#' imputed_data_knn <- ms_impute(filtered_data, 
#'                               imputeMethod = "knn",
#'                               k_knn = 5)
#' imputed_data_bpca <- ms_impute(filtered_data,
#'                                imputeMethod = "bpca",
#'                                n_pcs = 3)
#'
#' @importFrom dplyr case_when
#' @importFrom dplyr mutate_at
#' @export
ms_impute <- function(msprep_obj,
                      imputeMethod = c("halfmin", "bpca", "knn"),
                      k_knn = 5, 
                      n_pcs = 3, 
                      compoundsAsNeighbors = FALSE) {

  # Validate inputs
  stopifnot(class(msprep_obj) == "msprep")
  stopifnot(stage(msprep_obj) %in% c("filtered", "normalized"))
  imputeMethod <- match.arg(imputeMethod) # requires 1 argument from vec in function arg

  # Prep data - replace 0's with NA's -- for minval and bpca() (all methods?)
  data <- mutate_at(msprep_obj$data, vars("abundance_summary"), 
                    replace_missing, 0)
  grp <- grouping_vars(msprep_obj)
  batch <- batch_var(msprep_obj)
  met_vars <- met_vars(msprep_obj)
  
  # Store number of missing values
  missingCount <- sum(is.na(data$abundance_summary))

  # Impute data
  data <-
    switch(imputeMethod,
           "halfmin" = impute_halfmin(data, grp, batch, met_vars),
           "bpca"    = impute_bpca(data, grp, batch, met_vars, n_pcs),
           "knn"     = impute_knn(data, grp, batch, met_vars, k_knn, compoundsAsNeighbors),
           stop("Invalid impute method - you should never see this warning."))

  # Prep output object
  msprep_obj$data  <- data
  attr(msprep_obj, "impute_method") <- imputeMethod
  attr(msprep_obj, "missing_count") <- missingCount
  stage(msprep_obj) <- "imputed"

  # ...and:
  return(msprep_obj)

}

# Old return
#list(minval = as.data.frame(minval), withzero = as.data.frame(metafin), 
#     bpca = as.data.frame(bpca), count = as.data.frame(count))



#' @importFrom dplyr group_by
#' @importFrom dplyr mutate_at
#' @importFrom dplyr ungroup
#' @importFrom dplyr vars
#' @importFrom dplyr funs
#' @importFrom rlang sym
#' @importFrom rlang syms
#' @importFrom rlang !!
#' @importFrom rlang !!!
# Imputation using half of the minimum value
impute_halfmin <- function(data, groupingvars, batch, met_vars) {

  sym_mz <- sym("mz")
  sym_rt <- sym("rt")
  sym_met_vars <- syms(met_vars)

  # NOTE: other imputation methods operate across all batchs/groups for a given
  # mz_rt
  #   grp  <- syms(c(groupingvars, batch))
  # data <- group_by(data, `!!`(sym_mz), `!!!`(sym_met_vars))
  data <- group_by(data, `!!!`(sym_met_vars))

  halfmin <- function(x) {
    ifelse(is.na(x), min(x, na.rm = TRUE)/2, x)
  }

  data <- mutate_at(data, vars("abundance_summary"), halfmin)
  data <- ungroup(data)

  return(data)

}



#' @importFrom pcaMethods pca
#' @importFrom pcaMethods completeObs
impute_bpca <- function(data, groupingvars, batch, met_vars, n_pcs = 3) {

  # 1. Bayesian pca imputation
  data <- data_to_wide_matrix(data, groupingvars, batch, met_vars) 
  data <- pca(data, nPcs = n_pcs, method = "bpca")
  data <- completeObs(data) # extract imputed dataset
  data <- wide_matrix_to_data(data, groupingvars, batch, met_vars)
  data <- halfmin_if_any_negative(data, groupingvars, batch, met_vars)

  return(data)

}


#' @importFrom VIM kNN
impute_knn <- function(data, groupingvars, batch, met_vars, k_knn = 5, compoundsAsNeighbors) {
  
  data <- data_to_wide_matrix(data, groupingvars, batch, met_vars) 
  rwnm <- rownames(data)
  cnm <- colnames(data)
  
  # Perform kNN imputation using compounds or samples as neighbors
  if (compoundsAsNeighbors == TRUE) {
    data <- VIM::kNN(as.data.frame(data), k = k_knn, imp_var = FALSE)
  }
  else {
    ### transpose
    data <- VIM::kNN(as.data.frame(t(data)), k = k_knn, imp_var = FALSE)
    
    ## transpose again (to 'untranspose' basically)
    data <- t(data)
  }
  
  ### set column names same as prior to knn function
  colnames(data) <- cnm
  
  rownames(data) <- rwnm
  data <- wide_matrix_to_data(data, groupingvars, batch, met_vars)
  data <- halfmin_if_any_negative(data, groupingvars, batch, met_vars)

  return(data)

}



halfmin_if_any_negative <-  function(data, groupingvars, batch, met_vars) {

  num_neg <- sum(data$abundance_summary < 0)
  if (any(num_neg)) {
    message("Found ", num_neg, " negative imputed values using KNN, reverting to half-min imputation for these values")
    data$abundance_summary <- setmissing_negative_vals(data$abundance_summary)
    data <- impute_halfmin(data, groupingvars, batch, met_vars)
  }

  return(data)

}

truncate_negative_vals   <- function(var) ifelse(var < 0, 0, var)
setmissing_negative_vals <- function(var) ifelse(var < 0, NA, var)
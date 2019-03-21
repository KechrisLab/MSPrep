#' Function for imputing on filtered data.
#'
#' Performs data imputation for a given mz_rt value (not separately for batches
#' or grouping variables).
#'
#' @param msprep_obj Filtered MSPrep object.
#' @param method Name of imputation method to use. 
#' Options are:
#' - halfmin (half the minimum value)
#' - bpca (Bayesian PCA)
#' - knn (k-nearest neighbors)
#' @param k Number of clusters for 'knn' method.
#' @return An msprep object with missing data imputed.
#' @details minval Filtered dataset with missing values replaced by 1/2 minimum
#' observed value for that compound.
#' @details bpca Filtered dataset with missing values imputed by a Bayesian PCA
#' from PCAMethods package.
#' @references 
#'   Oba, S.et al.(2003) A Bayesian missing value estimation for gene
#'   expression profile data. Bioinformatics, 19, 2088-2096
#'
#'   Stacklies, W.et al.(2007) pcaMethods A bioconductor package providing
#'   PCA methods for incomplete data. Bioinformatics, 23, 1164-1167.
#' @examples
#' library(magrittr)
#'
#' # Load object generated from readdata() function
#' data(msquant)
#'
#' filtered_data <- msquant %>% ms_tidy %>% 
#'   ms_prepare(replicate = "replicate", 
#'              batch = "batch", 
#'              groupingvars = "spike") %>% 
#'   ms_filter(0.80)
#' imputed_data  <- ms_impute(filtered_data, "halfmin")
#'
#' @importFrom dplyr case_when
#' @importFrom dplyr mutate_at
#' @export
ms_impute <- function(msprep_obj,
                      method = c("halfmin", "bpca", "knn"),
                      k = 5) {

  # Validate inputs
  stopifnot(class(msprep_obj) == "msprep")
  stopifnot(stage(msprep_obj) %in% c("filtered", "normalized"))
  method <- match.arg(method) # requires 1 argument from vec in function arg

  # Prep data - replace 0's with NA's -- for minval and bpca() (all methods?)
  data <- mutate_at(msprep_obj$data, vars("abundance_summary"), 
                    replace_missing, 0)
  grp <- grouping_vars(msprep_obj)
  batch <- batch_var(msprep_obj)

  # Impute data
  data <-
    switch(method,
           "halfmin" = impute_halfmin(data, grp, batch),
           "bpca"    = impute_bpca(data, grp, batch),
           "knn"     = impute_knn(data, grp, batch, k),
           stop("Invalid impute method - you should never see this warning."))

  # Prep output object
  msprep_obj$data  <- data
  attr(msprep_obj, "impute_method") <- method
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
impute_halfmin <- function(data, groupingvars, batch) {

  sym_mz <- sym("mz")
  sym_rt <- sym("rt")

  # NOTE: other imputation methods operate across all batchs/groups for a given
  # mz_rt
  #   grp  <- syms(c(groupingvars, batch))
  data <- group_by(data, `!!`(sym_mz), `!!`(sym_rt))

  halfmin <- function(x) {
    ifelse(is.na(x), min(x, na.rm = TRUE)/2, x)
  }

  data <- mutate_at(data, vars("abundance_summary"), funs(halfmin))
  data <- ungroup(data)

  return(data)

}



#' @importFrom pcaMethods pca
#' @importFrom pcaMethods completeObs
impute_bpca <- function(data, groupingvars, batch) {

  # 1. Bayesian pca imputation
  data <- data_to_wide_matrix(data, groupingvars, batch) 
  data <- pca(data, nPcs = 3, method = "bpca")
  data <- completeObs(data) # extract imputed dataset
  data <- wide_matrix_to_data(data, groupingvars, batch)
  data <- halfmin_if_any_negative(data, groupingvars, batch)

  return(data)

}


#' @importFrom VIM kNN
impute_knn <- function(data, groupingvars, batch, k = 5) {

  data <- data_to_wide_matrix(data, groupingvars, batch) 
  rwnm <- rownames(data)
  cnm<- colnames(data)
  
  ### transpose
  data <- VIM:::kNN(as.data.frame(t(data)), k = k, imp_var = FALSE)
  
  ## transpose again (to 'untranspose' bsaically)
  data<- t(data)
  
  ### set column names same as prior to knn function
  colnames(data) <- cnm
  
  rownames(data) <- rwnm
  data <- wide_matrix_to_data(data, groupingvars, batch)
  data <- halfmin_if_any_negative(data)

  return(data)

}



halfmin_if_any_negative <-  function(data, groupingvars, batch) {

  num_neg <- sum(data$abundance_summary < 0)
  if (any(num_neg)) {
    message("Found ", num_neg, " negative imputed values using KNN, reverting to half-min imputation for these values")
    data$abundance_summary <- setmissing_negative_vals(data$abundance_summary)
    data <- impute_halfmin(data, groupingvars, batch)
  }

  return(data)

}


truncate_negative_vals   <- function(var) ifelse(var < 0, 0, var)
setmissing_negative_vals <- function(var) ifelse(var < 0, NA, var)

# OTHER METHODS
# 4. zero imputation -- 0 becomes 0.0001 (for normalization methods)
# 5. median imputation
# 6. mean imputation








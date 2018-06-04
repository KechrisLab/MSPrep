
#' Imputes dataset
#'
#' Performs data imputation.
#'
#' @param msprep_obj Placeholder
#' @param method Placeholder
#' @return Placeholder
#' @details minval Filtered dataset with missing values replaced by 1/2 minimum
#' observed value for that compound.
#' @details bpca Filtered dataset with missing values imputed by a Bayesian PCA
#' from PCAMethods package.
#' @details withzero Filtered dataset with no imputation.
#' @details count List of all compounds and the percent present for each
#' compound.
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
#' filtered_data <- msquant %>% ms_tidy %>% ms_prepare %>% ms_filter(0.80)
#' imputed_data <- ms_impute(filtered_data, "halfmin")
#'
#' @importFrom dplyr case_when
#' @export
ms_impute <- function(msprep_obj,
                      method = c("halfmin", "bpca", "knn")) {

  # Validate inputs
  stopifnot(class(msprep_obj) == "msprep")
  stopifnot(stage(msprep_obj) %in% c("filtered", "normalized"))
  match.arg(method) # requires 1 argument from vec in function arg

  # Prep data - replace 0's with NA's -- for minval and bpca() (all methods?)
  data <- mutate_at(msprep_obj$data, vars("abundance_summary"), 
                    replace_missing, 0)

  # Impute data
  data <-
    switch(method,
           "halfmin" = impute_halfmin(data),
           "bpca"    = impute_bpca(data),
           "knn"     = impute_knn(data),
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
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr vars
#' @importFrom dplyr funs
#' @importFrom rlang sym
#' @importFrom rlang UQ
# Imputation using half of the minimum value
impute_halfmin <- function(data) {

  sym_mz <- sym("mz")
  sym_rt <- sym("rt")

  data <- group_by(data, UQ(sym_mz), UQ(sym_rt))

  halfmin <- function(x) {
    ifelse(is.na(x), min(x, na.rm = TRUE)/2, x)
  }

  data <- mutate(data, vars("sym_abundance"), funs(halfmin))
  data <- ungroup(data)

  return(data)

}



#' @importFrom pcaMethods pca
#' @importFrom pcaMethods completeObs
impute_bpca <- function(data) {

  # 1. Bayesian pca imputation
  data <- data_to_wide_matrix(data) 
  data <- pca(data, nPcs = 3, method = "bpca")
  data <- completeObs(data) # extract imputed dataset
  data <- wide_matrix_to_data(data)
  data <- halfmin_if_any_negative(data)

  return(data)

}



#' @importFrom VIM kNN
impute_knn <- function(data, k = 5) {

  data <- data_to_wide_matrix(data) 
  rwnm <- rownames(data)
  data <- kNN(as.data.frame(data), k = k, imp_var = FALSE)
  rownames(data) <- rwnm
  data <- wide_matrix_to_data(data)
  data <- halfmin_if_any_negative(data)

  return(data)

}



halfmin_if_any_negative <-  function(data) {

  num_neg <- sum(data$abundance_summary < 0)
  if (any(num_neg)) {
    message("Found ", num_neg, " negative imputed values using KNN, reverting to half-min imputation for these values")
    data$abundance_summary <- setmissing_negative_vals(data$abundance_summary)
    data <- impute_halfmin(data)
  }

  return(data)

}


truncate_negative_vals   <- function(var) ifelse(var < 0, 0, var)
setmissing_negative_vals <- function(var) ifelse(var < 0, NA, var)

# OTHER METHODS
# 4. zero imputation -- 0 becomes 0.0001 (for normalization methods)
# 5. median imputation
# 6. mean imputation








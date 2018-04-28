
#' Imputes dataset
#'
#' Performs data imputation.
#'
#' @param msprepped Placeholder
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
ms_impute <- function(msprepped,
                      method = c("halfmin", "bpca", "bpca + halfmin", "knn")) {

  # Validate inputs
  stopifnot("msprepped" %in% class(msprepped))
  match.arg(method) # requires 1 argument from vec in function arg
  #check_ms_flow() # for checking valid path through pipeline

  # Prep data - replace 0's with NA's -- for minval and bpca() (all methods?)
  data <- mutate_at(msprepped$summary_data, vars(abundance_summary), 
                    replace_missing, 0)

  # Impute data
  data <-
    switch(method,
           "halfmin"        = impute_halfmin(data),
           "bpca"           = impute_bpca_zero(data),
           "bpca + halfmin" = impute_bpca_halfmin(data),
           "knn"            = impute_knn(data),
           stop("Invalid impute method - you should never see this warning."))

  # Prep output object
  msprepped$summary_data  <- data
  msprepped$impute_method <- method
  class(msprepped) <- c("msimputed", class(msprepped))
  #msprepped <- mark_completed(msprepped, "imputation") # mark pipeline step as complete

  # ...and:
  return(msprepped)

}

# Old return
#list(minval = as.data.frame(minval), withzero = as.data.frame(metafin), 
#     bpca = as.data.frame(bpca), count = as.data.frame(count))



#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
# Imputation using half of the minimum value 
# should this be within subject or across all datapoints?  currently across
# all
impute_halfmin <- function(data) {

  data <- group_by(data, spike, mz, rt)
  data <- mutate(data, 
                 abundance_summary = ifelse(is.na(abundance_summary),
                                            min(abundance_summary, na.rm = TRUE)/2,
                                            abundance_summary))
  data <- ungroup(data)

  return(data)

}



#' @importFrom pcaMethods pca
#' @importFrom pcaMethods completeObs
impute_bpca <- function(data) {

  # 1. Bayesian pca imputation
  browser()
  data <- data_to_wide_matrix(data) 
  data <- pca(data, nPcs = 3, method = "bpca")
  data <- completeObs(data) # extract imputed dataset
  data <- wide_matrix_to_data(data)

  return(data)

}


impute_bpca_halfmin <- function(data) {

  data <- impute_bpca(data)
  data <- impute_halfmin(data)

  return(data)

}


impute_bpca_zero <- function(data) {

  data <- impute_bpca(data)
  data$abundance_summary <- truncate_negative_vals(data$abundance_summary)

  return(data)

}


#' @importFrom VIM kNN
impute_knn <- function(data, k = 5) {

  data <- data_to_wide_matrix(data) 
  rwnm <- rownames(data)
  data <- kNN(as.data.frame(data), k = k, imp_var = FALSE)
  rownames(data) <- rwnm

  data <- wide_matrix_to_data(data)
  data$abundance_summary <- truncate_negative_vals(data$abundance_summary)

  return(data)

}

truncate_negative_vals <- function(var) ifelse(var < 0, 0, var)

# OTHER METHODS
# 4. zero imputation -- 0 becomes 0.0001 (why?)
# 5. median imputation
# 6. mean imputation








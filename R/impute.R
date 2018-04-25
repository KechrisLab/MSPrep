# #' Imputes dataset
# #'
# #' and
# #' performs data imputation.
# #'
# #' @param msprepped Placeholder
# #' @param impute_function Placeholder
# #' @return Placeholder
# #' @details minval Filtered dataset with missing values replaced by 1/2 minimum
# #' observed value for that compound.
# #' @details bpca Filtered dataset with missing values imputed by a Bayesian PCA
# #' from PCAMethods package.
# #' @details withzero Filtered dataset with no imputation.
# #' @details count List of all compounds and the percent present for each
# #' compound.
# #' @references 
# #'   Oba, S.et al.(2003) A Bayesian missing value estimation for gene
# #'   expression profile data. Bioinformatics, 19, 2088-2096
# #'
# #'   Stacklies, W.et al.(2007) pcaMethods A bioconductor package providing
# #'   PCA methods for incomplete data. Bioinformatics, 23, 1164-1167.
# #' @examples
# #' library(magrittr)
# #'
# #' # Load object generated from readdata() function
# #' data(msquant)
# #'
# #' filtered_data <- msquant %>% ms_tidy %>% ms_prepare %>% ms_filter(0.80)
# #' imputed_data <- ms_impute(filtered_data, 0.80)
# #'
# #' @export
# ms_impute <- function(msprepped, impute_function = ~ minval(.x)) {
# 
#   stopifnot("msprepped" %in% class(msprepped))
# 
#   # replace 0's with NA's for minval and pca() (all methods?)
#   data <- msprepped$summary_data %>% 
#     mutate_at(vars(abundance_summary), replace_missing, 0) 
# 
#   data %>% impute_function
# 
# }



ms_impute <- function(msprepped, impute_function = ~ minval(.x)) {

  stopifnot("msprepped" %in% class(msprepped))

  # replace 0's with NA's for minval and pca() (all methods?)
  data <- msprepped$summary_data %>% 
    mutate_at(vars(abundance_summary), replace_missing, 0)
  # pca imputation
  metabpca <- data_to_wide_matrix(data) %>% pcaMethods::pca(., nPcs = 3, method = "bpca")
  # k-nearest-neighbors imputation
  knnimpute <- data_to_wide_matrix(data) %>% as.data.frame %>% VIM::kNN(., k=5)
  #[1:ncol(data), 1:nrow(data)]

  data %>% impute_function

}

# Filterft to only include compounds that are found in specified percentage of
# subjects and perform imputation of missing data

# Manuscript description: 
# Filtering: The resulting summarized dataset contains all compounds with one
# observation per subject (or sample). The next processing step filters the data
# to only compounds found in a user-specified percentage of subjects.

data_to_wide_matrix <- function(data) {

  data %>% 
    unite(col = sample, subject_id, spike) %>% 
    unite(col = compound, mz, rt) %>% 
    spread(key = compound, value = abundance_summary)  %>%
    as.data.frame %>%
    column_to_rownames(var = "sample") %>%
    as.matrix

}



#  colnames(toss) <- "toss"
#  tests <- cbind(toss, t(metaf))
#  metafin <- t(subset(tests, toss != 0)[, -1])
#  metaimp <- matrix(NA, nrow = nrow(metafin), ncol = ncol(metafin))
#  for (i in 1:nrow(metafin)) {
#    for (j in 1:ncol(metafin)) {
#      # sean jacobson added this so that missing data would be ignored
#      # if(is.na(metafin[i, j])) next
#      if (metafin[i, j] == 0) {
#        metaimp[i, j] <- NA
#      }
#      else {
#        metaimp[i, j] <- metafin[i, j]
#      }
#    }
#  }
#
#}
#
#  present <- matrix(1, nrow = nrow(metaf), ncol = ncol(metaf))
#  colnames(present) <- colnames(metaf)
#  rownames(present) <- rownames(metaf)
#  for (j in 1:eval(ncol(metaf))) {
#    for (i in 1:dim(metaf)[1]) {
#      # sean jacobson added this for missing data
#      # if(is.na(metaf[i, j])){present[i,j] <- NA; next}
#      if (metaf[i, j] == 0) {
#        present[i, j] <- 0
#      }
#    }
#  }
#
#  minval <- metafin
#  for (i in 1:nrow(metaf)) {
#    for (j in 1:ncol(metafin)) {
#      # sean jacobson added this for missing data
#      # if(is.na(metafin[i, j])){minval[i,j] <- NA; next}
#      if (metafin[i, j] == 0) {
#        minval[i, j] <- summary(metaimp[, j])[1]/2
#      }
#      else {
#        minval[i, j] <- metafin[i, j]
#      }
#    }
#  }
#
#  metabpca <- pcaMethods::pca(metaimp, nPcs = 3, method = "bpca")
#  bpca <- completeObs(metabpca)
#  colnames(bpca) <- colnames(metafin)
#  rownames(bpca) <- rownames(metafin)
#  coldrop <- 1
#  for (i in 1:nrow(bpca)) {
#    for (j in 1:ncol(bpca)) {
#      if (bpca[i, j] < 0) {
#        bpca[i, j] <- minval[i, j]
#      }
#    }
#  }
#  list(minval = as.data.frame(minval), withzero = as.data.frame(metafin), 
#    bpca = as.data.frame(bpca), count = as.data.frame(count))
#}
#
#

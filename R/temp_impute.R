
impute_ms <- function(msprepped, impute_function = ~ minval(.x)) {

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

data_to_wide_matrix <- function(data) {

  data %>% 
    unite(col = sample, subject_id, spike) %>% 
    unite(col = compound, mz, rt) %>% 
    spread(key = compound, value = abundance_summary)  %>%
    as.data.frame %>%
    column_to_rownames(var = "sample") %>%
    as.matrix

}


# TODO: 
# Add a dataframe object to first step (ms_prepare) that has a column for
# pipline step (prepare, filter, impute, etc.) and a logical column for
# conducted yes/no.  Then create a function that defines acceptable paths
# through the pipeline -- give error when invalid path.


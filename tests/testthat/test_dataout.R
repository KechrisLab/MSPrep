
context("ms_impute()")

set.seed(9999)
data(msquant_subject1)
tidy_data       <- ms_tidy(msquant_subject1, mz = "mz", rt = "rt",
                           col_extra_txt = "Neutral_Operator_Dif_Pos_", 
                           separator = "_", 
                           col_names = c("spike", "batch", "replicate", "subject_id"))
prepped_data    <- tidy_data %>% 
  ms_prepare(mz = "mz", rt = "rt", replicate = "replicate", batch = "batch", groupingvars = "spike") %>%
  ms_filter(0.8)

# impute_methods <- c("halfmin", "knn", "bpcs")
# imputed_lst <- impute_methods %>%
#     mclapply(., function(x) prepped_data %>% ms_impute(x),
#              mc.cores = 3)
# halfmin_imputed <- imputed_lst[[1]]
# knn_imputed     <- imputed_lst[[2]]
# bpca_imputed    <- imputed_lst[[3]]

halfmin_imputed <- prepped_data %>% ms_impute("halfmin")
knn_imputed     <- prepped_data %>% ms_impute("knn")
bpca_imputed    <- prepped_data %>% ms_impute("bpca")

test_that("Check dataset columns", {

    halfmin_colnames <- sort(colnames(halfmin_imputed$data))
    knn_colnames     <- sort(colnames(knn_imputed$data))
    bpca_colnames    <- sort(colnames(bpca_imputed$data))
    halfmin_extravars <- c(grouping_vars(halfmin_imputed),
                           batch_var(halfmin_imputed))
    knn_extravars     <- c(grouping_vars(knn_imputed), batch_var(knn_imputed))
    bpca_extravars    <- c(grouping_vars(bpca_imputed), batch_var(bpca_imputed))
    
    expect_equal(halfmin_colnames, valid_cols(halfmin_extravars))
    expect_equal(knn_colnames, valid_cols(knn_extravars))
    expect_equal(bpca_colnames, valid_cols(bpca_extravars))

})

test_that("Check imputations", {

    missvals <- which(prepped_data$data$abundance_summary == 0)
    halfmin_miss <- halfmin_imputed$data[missvals, ]
    knn_miss     <- knn_imputed$data[missvals, ]
    bpca_miss    <- bpca_imputed$data[missvals, ]
    expect_true(all(halfmin_miss$abundance_summary != 0))
    expect_true(all(knn_miss$abundance_summary != 0))
    expect_true(all(bpca_miss$abundance_summary != 0))
    expect_true(sum(halfmin_miss$abundance_summary == bpca_miss$abundance_summary) == 2)
    halfmin_imputed$data == bpca_imputed$data

})


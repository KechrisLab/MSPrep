
context("ms_impute()")

set.seed(9999)
data(msquant)
tidy_data       <- ms_tidy(msquant, mz = "mz", rt = "rt")
prepped_data    <- tidy_data %>% ms_prepare %>% ms_filter(0.8)
halfmin_imputed <- prepped_data %>% ms_impute("halfmin")
knn_imputed     <- prepped_data %>% ms_impute("knn")
bpca_imputed    <- prepped_data %>% ms_impute("bpca")

test_that("Check dataset columns", {

   halfmin_colnames <- sort(colnames(halfmin_imputed$data))
   knn_colnames     <- sort(colnames(knn_imputed$data))
   bpca_colnames    <- sort(colnames(bpca_imputed$data))
   expect_equal(halfmin_colnames, valid_cols(grouping_vars(halfmin_imputed)))
   expect_equal(knn_colnames, valid_cols(grouping_vars(knn_imputed)))
   expect_equal(bpca_colnames, valid_cols(grouping_vars(bpca_imputed)))

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

})


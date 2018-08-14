context("utility functions")

set.seed(9999)
data(msquant_subject1)
tidy_data       <- ms_tidy(msquant_subject1, mz = "mz", rt = "rt")
prepped_data    <- tidy_data %>% 
  ms_prepare(replicate = "replicate", batch = "batch", grouping_vars = "spike")

test_that("Check dataset transformations give back original dataset.", {

    grps  <- grouping_vars(prepped_data)
    batch <- batch_var(prepped_data)
    original_data <- prepped_data$data
    xfrmd_data <-
      prepped_data$data %>%
      data_to_wide_matrix(., grps, batch) %>%
      wide_matrix_to_data(., grps, batch)

    original_data == xfrmd_data

})

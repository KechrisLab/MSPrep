context("utility functions")

set.seed(9999)
data(msquant_subject1)
tidy_data       <- ms_tidy(msquant_subject1, mz = "mz", rt = "rt",
                           col_extra_txt = "Neutral_Operator_Dif_Pos_", 
                           separator = "_", 
                           col_names = c("spike", "batch", "replicate", "subject_id"))
prepped_data    <- tidy_data %>% 
  ms_prepare(mz = "mz", rt = "rt", replicate = "replicate", batch = "batch", groupingvars = "spike")

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

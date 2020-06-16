context("check that new equals old")

data("msquant")

tidied_data <- ms_tidy(msquant, mz = "mz", rt = "rt",
                       col_extra_txt = "Neutral_Operator_Dif_Pos_",
                       separator = "_", 
                       col_names = c("spike", "batch", "replicate", 
                                     "subject_id"))

summarized_data <- ms_summarize(tidied_data, mz = "mz", rt = "rt", 
                                replicate = "replicate", 
                                batch = "batch", groupingvars = "spike", 
                                subject_id = "subject_id", 
                                cvmax = 0.50, min_proportion_present = 1/3, 
                                missing_val = 1)

filtered_data <- ms_filter(summarized_data, filter_percent = .8)

summarized_data <- summarized_data %>% ms_return()
filtered_data <- filtered_data %>% ms_return()

summarizedDF <- msSummarize(msquant,
                            compVars = c("mz", "rt"),
                            techVars = c("spike", "batch", "replicate", 
                                         "subject_id"),
                            cvMax = 0.50,
                            minPropPresent = 1/3,
                            returnSummaryDetails = FALSE,
                            colExtraText = "Neutral_Operator_Dif_Pos_",
                            separator = "_",
                            missingValue = 1)

filteredDF <- msFilter(summarizedDF,
                       filterPercent = 0.8,
                       compVars = c("mz", "rt"),
                       techVars = c("spike", "batch", "subject_id"),
                       separator = "_",
                       returnToSE = FALSE)

test_that("msSummarize() == ms_summarize()", {
    expect_true(all(summarized_data$data[["1x_O1_01"]] == 
                        summarizedDF[["1x_O1_01"]],
                    summarized_data$data[["2x_O1_01"]] == 
                        summarizedDF[["2x_O1_01"]],
                    summarized_data$data[["4x_O1_01"]] == 
                        summarizedDF[["4x_O1_01"]],
                    summarized_data$data[["1x_O2_01"]] == 
                        summarizedDF[["1x_O2_01"]],
                    summarized_data$data[["2x_O2_01"]] == 
                        summarizedDF[["2x_O2_01"]],
                    summarized_data$data[["4x_O2_01"]] == 
                        summarizedDF[["4x_O2_01"]],
                    summarized_data$data[["1x_O3_01"]] == 
                        summarizedDF[["1x_O3_01"]],
                    summarized_data$data[["2x_O3_01"]] == 
                        summarizedDF[["2x_O3_01"]],
                    summarized_data$data[["4x_O3_01"]] == 
                        summarizedDF[["4x_O3_01"]],
                    summarized_data$data[["1x_O1_02"]] == 
                        summarizedDF[["1x_O1_02"]],
                    summarized_data$data[["2x_O1_02"]] == 
                        summarizedDF[["2x_O1_02"]],
                    summarized_data$data[["4x_O1_02"]] == 
                        summarizedDF[["4x_O1_02"]],
                    summarized_data$data[["1x_O2_02"]] == 
                        summarizedDF[["1x_O2_02"]],
                    summarized_data$data[["2x_O2_02"]] == 
                        summarizedDF[["2x_O2_02"]],
                    summarized_data$data[["4x_O2_02"]] == 
                        summarizedDF[["4x_O2_02"]],
                    summarized_data$data[["1x_O3_02"]] == 
                        summarizedDF[["1x_O3_02"]],
                    summarized_data$data[["2x_O3_02"]] == 
                        summarizedDF[["2x_O3_02"]],
                    summarized_data$data[["4x_O3_02"]] == 
                        summarizedDF[["4x_O3_02"]]))
})

test_that("msFilter() == ms_filter()", {
    expect_true(all(filtered_data$data[["1x_O1_01"]] == 
                        filteredDF[["1x_O1_01"]],
                    filtered_data$data[["2x_O1_01"]] == 
                        filteredDF[["2x_O1_01"]],
                    filtered_data$data[["4x_O1_01"]] == 
                        filteredDF[["4x_O1_01"]],
                    filtered_data$data[["1x_O2_01"]] == 
                        filteredDF[["1x_O2_01"]],
                    filtered_data$data[["2x_O2_01"]] == 
                        filteredDF[["2x_O2_01"]],
                    filtered_data$data[["4x_O2_01"]] == 
                        filteredDF[["4x_O2_01"]],
                    filtered_data$data[["1x_O3_01"]] == 
                        filteredDF[["1x_O3_01"]],
                    filtered_data$data[["2x_O3_01"]] == 
                        filteredDF[["2x_O3_01"]],
                    filtered_data$data[["4x_O3_01"]] == 
                        filteredDF[["4x_O3_01"]],
                    filtered_data$data[["1x_O1_02"]] == 
                        filteredDF[["1x_O1_02"]],
                    filtered_data$data[["2x_O1_02"]] == 
                        filteredDF[["2x_O1_02"]],
                    filtered_data$data[["4x_O1_02"]] == 
                        filteredDF[["4x_O1_02"]],
                    filtered_data$data[["1x_O2_02"]] == 
                        filteredDF[["1x_O2_02"]],
                    filtered_data$data[["2x_O2_02"]] == 
                        filteredDF[["2x_O2_02"]],
                    filtered_data$data[["4x_O2_02"]] == 
                        filteredDF[["4x_O2_02"]],
                    filtered_data$data[["1x_O3_02"]] == 
                        filteredDF[["1x_O3_02"]],
                    filtered_data$data[["2x_O3_02"]] == 
                        filteredDF[["2x_O3_02"]],
                    filtered_data$data[["4x_O3_02"]] == 
                        filteredDF[["4x_O3_02"]]))
})
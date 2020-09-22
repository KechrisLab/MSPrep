context("msSummarize() related functions")

# ### Read in data files
# pathclinical <- system.file("data-raw", "Clinical.csv", 
#                             package = "MSPrep")
# pathquant <- system.file("extdata", "Quantification.csv", package = "MSPrep")
# pathlink <- system.file("data-raw", "SubjectLinks.csv", package = "MSPrep")
# path_olddata <- system.file("extdata", "old_object.Rda", package = "MSPrep")
# data(msquant_subject1)
# quant        <- msquant_subject1
# load(path_olddata)
# 
# Generate tidy dataset from wide quant data
# tidy_data    <- ms_tidy(quant, mz = "mz", rt = "rt",
#                         col_extra_txt = "Neutral_Operator_Dif_Pos_",
#                         separator = "_",
#                         col_names = c("spike", "batch", "replicate", 
#                                       "subject_id"))
# prepped_data <- tidy_data %>% ms_summarize(mz = "mz",
#                                          rt = "rt",
#                                          replicate = "replicate",
#                                          batch = "batch",
#                                          groupingvars = "spike")

data(msquant)

preppedData <- msSummarize(msquant,
                           returnSummaryDetails = TRUE,
                           compVars = c("mz", "rt"),
                           sampleVars = c("spike", "batch", "replicate",
                                         "subject_id"),
                           colExtraText = "Neutral_Operator_Dif_Pos_",
                           separator = "_",
                           missingValue = 1)
summaryDetails <- preppedData$summaryDetails
preppedData <- preppedData$data



test_that(".replaceMissing()", {

    # numeric to NA
    testVec <- c(4, 6, 1, 9, 5, 8, 3, 7, 10, 1)
    replacedVec <- MSPrep:::.replaceMissing(testVec, missingValue = 1,
                                            setMissing = NA)
    expect_equal(sum(is.na(replacedVec)), 2)
    expect_equal(replacedVec[c(3, 10)], c(NA_real_, NA_real_))
    
    # numeric to numeric
    replacedVec <- MSPrep:::.replaceMissing(testVec, missingValue = 1,
                                            setMissing = 0)
    expect_equal(sum(replacedVec == 0), 2)
    expect_equal(replacedVec[c(3, 10)], c(0, 0))

    # character to NA
    testVec     <- as.character(testVec)
    replacedVec <- MSPrep:::.replaceMissing(testVec, 
                                            missingVal = as.character(1),
                                            setMissing = NA)
    expect_equal(sum(is.na(replacedVec)), 2)
    expect_equal(replacedVec[c(3, 10)], c(NA_character_, NA_character_))
    
})

test_that(".selectSummaryMeasure", {

  nReplicates <- 3
  cvMax <- 0.50
  minPropPresent <- 1/3
  
  dat <- expand.grid(nPresent = c(0, 1, 2, 3), cvAbundance = c(0.1, 0.5, 1.0))
  
  results <- MSPrep:::.selectSummaryMeasure(dat$nPresent, dat$cvAbundance, 
                                            nReplicates, minPropPresent, cvMax)
  
  expect_equal(results, c("none: proportion present <= minProportionPresent",
                          "none: proportion present <= minProportionPresent",
                          "mean", "mean",
                          "none: proportion present <= minProportionPresent",
                          "none: proportion present <= minProportionPresent",
                          "mean", "mean",
                          "none: proportion present <= minProportionPresent",
                          "none: proportion present <= minProportionPresent",
                          "none: cv > cvMax & 2 present",
                          "median"))

})

test_that("Summary measure count remains constant in test data", {

    expect_equal(46, sum(summaryDetails$summaryMeasure == "median"))
    expect_equal(23872, sum(summaryDetails$summaryMeasure == "mean"))
    expect_equal(23638, 
                 sum(summaryDetails$summaryMeasure == 
                         "none: proportion present <= minProportionPresent"))
    expect_equal(36, 
                 sum(summaryDetails$summaryMeasure == 
                         "none: cv > cvMax & 2 present"))

})

# test_that("New version of summarized dataset matches old version", {
# 
#   # Add old readdata() function that creates old_readdata_result 
#   # in inst/extdata/
#   # folder
# 
#   # read in old dataset
#   # load("R/sysdata.rda")
# 
#   sum_data <-
#     prepped_data$data %>%
#     dplyr::select(batch, spike, mz, rt, abundance_summary) %>%
#     tidyr::unite(id, spike, batch, sep = "_") %>%
#     tidyr::unite(metabolite, mz, rt, sep = "_")
# 
#   #new_sum_data <- sum_data %>%
#   #  tidyr::spread(key = id, value = abundance_summary)
# 
#   old_sum <- old_readdata_result$sum_data1
#   old_sum <- old_sum %>% as.data.frame
#   old_sum <- old_sum %>% tibble::rownames_to_column(var = "id")
#   old_sum <- old_sum %>% tidyr::gather(key = metabolite, 
#                                        value = abundance_summary, -id)
#   old_sum <- old_sum %>% tibble::as_tibble(.)
#   old_sum <- old_sum %>% dplyr::select(id, metabolite, abundance_summary)
#   old_sum <- old_sum %>% tidyr::separate(metabolite, into = c("mz", "rt"), 
#                                          sep = "_")
#   old_sum <- old_sum %>% dplyr::mutate(mz = as.numeric(mz))
# #   old_sum <- mutate(rt = as.numeric(rt))
#   old_sum <- old_sum %>% dplyr::arrange(mz)
#   old_sum <- old_sum %>% tidyr::unite(metabolite, mz, rt, sep = "_")
# 
#   new <- sum_data %>% tidyr::separate(metabolite, into = c("mz", "rt"), 
#                                       sep = "_") %>%
#     dplyr::arrange(id, mz, rt)
#   old <- old_sum %>% tidyr::separate(metabolite, into = c("mz", "rt"), 
#                                      sep = "_") %>%
#     dplyr::arrange(id, mz, rt)
#   diffrows <- !(new == old)[, 4]
#   newdiffs <- new[diffrows, ] %>% dplyr::rename(new_summary = 
#                                                     abundance_summary)
#   olddiffs <- old[diffrows, ] %>% dplyr::rename(old_summary = 
#                                                     abundance_summary)
# 
#   expect_true(nrow(dplyr::anti_join(old_sum, sum_data)) == 0)
# 
#   #comparing <-
#   #  .data %>% mutate(mz = as.character(mz), rt = as.character(rt)) %>%
#   #    unite(id, spike, subject_id) %>%
#   #    right_join(., newdiffs) %>%
#   #    left_join(., olddiffs)
# 
#   expect_true(identical(new, old))
# 
# })




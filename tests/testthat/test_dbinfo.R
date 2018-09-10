

context("ms_prepare() related functions")


### Read in data files
# pathclinical <- system.file("data-raw", "Clinical.csv", package = "MSPrep")
pathquant    <- system.file("extdata", "Quantification.csv", package = "MSPrep")
# pathlink     <- system.file("data-raw", "SubjectLinks.csv", package = "MSPrep")
path_olddata <- system.file("extdata", "old_object.Rda", package = "MSPrep")
data(msquant_subject1)
quant        <- msquant_subject1 
load(path_olddata)

# Generate tidy dataset from wide quant data
tidy_data    <- ms_tidy(quant, mz = "mz", rt = "rt")
prepped_data <- tidy_data %>% ms_prepare(replicate = "replicate", 
                                         batch = "batch",
                                         groupingvars = "spike")



test_that("replace_missing replaces missing val w/ NA ", {

    testvec     <- c(4, 6, 1, 9, 5, 8, 3, 7, 10, 1)
    replacedvec <- MSPrep:::replace_missing(testvec, missing_val = 1)
    expect_equal(sum(is.na(replacedvec)), 2)
    expect_equal(replacedvec[3], NA_real_)

    testvec     <- as.character(testvec)
    replacedvec <- MSPrep:::replace_missing(testvec, missing_val = as.character(1))
    expect_equal(sum(is.na(replacedvec)), 2)
    expect_equal(replacedvec[3], NA_character_)

})



test_that("Check default select_summary_measure works", {

  n_replicates           <- 3
  cv_max                 <- 0.50
  min_proportion_present <- 1/3
  dat <- expand.grid(n_present    = c(0, 1, 2, 3),
                     cv_abundance = c(0.1, 0.5, 1.0))
  results <- MSPrep:::select_summary_measure(dat$n_present,
                                             dat$cv_abundance,
                                             n_replicates,
                                             min_proportion_present,
                                             cv_max)
  expect_equal(results, c("none: proportion present <= min_proportion_present",
                          "none: proportion present <= min_proportion_present",
                          "mean", "mean", 
                          "none: proportion present <= min_proportion_present", 
                          "none: proportion present <= min_proportion_present", 
                          "mean", "mean", 
                          "none: proportion present <= min_proportion_present", 
                          "none: proportion present <= min_proportion_present", 
                          "none: cv > cvmax & 2 present",
                          "median"))

})


test_that("Number summarised by median stays constant in test data", {

   expect_equal(23, nrow(prepped_data$median))

})

test_that("New version of summarized dataset matches old version", {

  # Add old readdata() function that creates old_readdata_result in inst/extdata/
  # folder

  # read in old dataset
  # load("R/sysdata.rda")

  sum_data <-
    prepped_data$data %>%
    dplyr::select(batch, spike, mz, rt, abundance_summary) %>%
    tidyr::unite(id, spike, batch, sep = "_") %>% 
    tidyr::unite(metabolite, mz, rt, sep = "_")

  #new_sum_data <- sum_data %>% 
  #  tidyr::spread(key = id, value = abundance_summary)

  old_sum <- old_readdata_result$sum_data1 
  old_sum <- old_sum %>% as.data.frame 
  old_sum <- old_sum %>% tibble::rownames_to_column(var = "id") 
  old_sum <- old_sum %>% tidyr::gather(key = metabolite, value = abundance_summary, -id) 
  old_sum <- old_sum %>% tibble::as_data_frame(.)
  old_sum <- old_sum %>% dplyr::select(id, metabolite, abundance_summary)
  old_sum <- old_sum %>% tidyr::separate(metabolite, into = c("mz", "rt"), sep = "_")
  old_sum <- old_sum %>% dplyr::mutate(mz = as.numeric(mz))
#   old_sum <- mutate(rt = as.numeric(rt))
  old_sum <- old_sum %>% dplyr::arrange(mz)
  old_sum <- old_sum %>% tidyr::unite(metabolite, mz, rt, sep = "_")

  new <- sum_data %>% tidyr::separate(metabolite, into = c("mz", "rt"), sep = "_") %>%
    dplyr::arrange(id, mz, rt)
  old <- old_sum %>% tidyr::separate(metabolite, into = c("mz", "rt"), sep = "_") %>% 
    dplyr::arrange(id, mz, rt)
  diffrows <- !(new == old)[, 4]
  newdiffs <- new[diffrows, ] %>% dplyr::rename(new_summary = abundance_summary)
  olddiffs <- old[diffrows, ] %>% dplyr::rename(old_summary = abundance_summary)

  expect_true(nrow(dplyr::anti_join(old_sum, sum_data)) == 0)

  #comparing <-
  #  .data %>% mutate(mz = as.character(mz), rt = as.character(rt)) %>% 
  #    unite(id, spike, subject_id) %>% 
  #    right_join(., newdiffs) %>%
  #    left_join(., olddiffs)

  expect_true(identical(new, old))

})




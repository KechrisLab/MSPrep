

context("prepare_msprep() related functions")


### Read in data files
clinical <- read.csv("./data-raw/Clinical.csv")
quant    <- read.csv("./data-raw/Quantification.csv")
link     <- read.csv("./data-raw/SubjectLinks.csv")
clinical_data <- clinical
link_data     <- link
cvmax         <- 0.50
missing       <- 1
linktxt       <- "LCMS_Run_ID"

.data         <- tidy_quant(quant, mz = "mz", rt = "rt")

test_that("replace_missing replaces missing val w/ NA ", {


    testvec     <- c(4, 6, 1, 9, 5, 8, 3, 7, 10, 1)
    replacedvec <- replace_missing(testvec, missing_val = missing)
    expect_equal(sum(is.na(replacedvec)), 2)
    expect_equal(replacedvec[3], NA_real_)

    testvec     <- as.character(testvec)
    replacedvec <- replace_missing(testvec, missing_val = as.character(missing))
    expect_equal(sum(is.na(replacedvec)), 2)
    expect_equal(replacedvec[3], NA_character_)

})

test_that("Check default select_summary_measure works", {

  n_replicates           <- 3
  cv_max                 <- 0.50
  min_proportion_present <- 1/3
  dat <- expand.grid(n_present    = c(0, 1, 2, 3),
                     cv_abundance = c(0.1, 0.5, 1.0))
  results <- select_summary_measure(dat$n_present,
                                    dat$cv_abundance,
                                    n_replicates,
                                    min_proportion_present,
                                    cv_max)
  expect_equal(results, c(NA, NA, "mean", "mean", NA, NA, "mean", "mean", NA,
                          NA, NA, "median"))

})


test_that("New version of summarized dataset matches old version", {

  # Add old readdata() function that creates old_readdata_result in data-raw/
  # folder

  sum_data <-
    quant_summary %>% 
    select(subject_id, spike, mz, rt, abundance_summary, summary_measure) %>%
    unite(id, spike, subject_id, sep = "_") %>% 
    unite(metabolite, mz, rt, sep = "_")

  new_sum_data <- sum_data %>% spread(key = id, value = abundance_summary)

  old_sum_data <- 
    old_readdata_result$sum_data1 %>% 
      as.data.frame %>% 
      rownames_to_column(var = "id") %>%
      gather(key = metabolite, value = abundance_summary, -id) %>%
      as_data_frame  %>%
      select(id, metabolite, abundance_summary) %>%
      separate(metabolite, into = c("mz", "rt"), sep = "_") %>%
      mutate(mz = as.numeric(mz)) %>%
      #mutate(rt = as.numeric(rt)) %>%
      arrange(mz) %>%
      unite(metabolite, mz, rt, sep = "_")

  # 
  new <- sum_data %>% separate(metabolite, into = c("mz", "rt"), sep = "_") %>% arrange(id, mz, rt)
  old <- old_sum_data %>% separate(metabolite, into = c("mz", "rt"), sep = "_") %>% arrange(id, mz, rt)
  diffrows <- !(select(new, -summary_measure) == old)[, 4]
  newdiffs <- new[diffrows, ] %>% rename(new_summary = abundance_summary)
  olddiffs <- old[diffrows, ] %>% rename(old_summary = abundance_summary)

  nrow(anti_join(old_sum_data, sum_data)) == 0

  comparing <-
    .data %>% mutate(mz = as.character(mz), rt = as.character(rt)) %>% 
      unite(id, spike, subject_id) %>% 
      right_join(., newdiffs) %>%
      left_join(., olddiffs)

}



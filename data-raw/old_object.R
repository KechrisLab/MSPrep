################################################################################
# clean old object data for tests
################################################################################

# file.rename("inst/extdata/old_object.Rda", "data-raw/old_object.Rda")

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

load("data-raw/old_object.Rda")
str(old_readdata_result)


# Remove columns w/ ":" in the title/non-numeric rt (retention time) portion of
# column names 
old_readdata_result$sum_data1 <- 
  colnames(old_readdata_result$sum_data1) %>% 
  str_detect(":") %>%  `!` %>% old_readdata_result$sum_data1[, .]

# Other list elements don't need to be changed
old_readdata_result$clinical
old_readdata_result$medians[, 4] %>% str_detect(":") %>% `!`

save(old_readdata_result, file = "inst/extdata/old_object.Rda")



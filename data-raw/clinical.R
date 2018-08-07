################################################################################
# Clean clinical, quantification, and subjectlinks datasets for use in tests,
# examples, and vignettes
#
# NOTE: only quantification needed -- other datasets info in its column names
################################################################################


file.rename("inst/extdata/Clinical.csv", "data-raw/Clinical.csv")
file.rename("inst/extdata/Quantification.csv", "data-raw/Quantification.csv")
file.rename("inst/extdata/SubjectLinks.csv", "data-raw/SubjectLinks.csv")


library(dplyr)
library(tidyr)
library(readr)
library(stringr)

clinical <- read_csv("data-raw/Clinical.csv")
msquant  <- read_csv("data-raw/Quantification.csv")
subject  <- read_csv("data-raw/SubjectLinks.csv")

msquant <- msquant %>% filter(!str_detect(rt, ":"))

devtools::use_data(msquant, overwrite = TRUE)




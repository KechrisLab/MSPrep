################################################################################
# Clean clinical, quantification, and subjectlinks datasets for use in tests,
# examples, and vignettes
#
# Notes:
#   - only quantification needed -- other datasets info in its column names
#   - generating second subject from first, for testing multiple subjects 
################################################################################


# file.rename("inst/extdata/Clinical.csv", "data-raw/Clinical.csv")
# file.rename("inst/extdata/Quantification.csv", "data-raw/Quantification.csv")
# file.rename("inst/extdata/SubjectLinks.csv", "data-raw/SubjectLinks.csv")


library(dplyr)
library(tidyr)
library(readr)
library(stringr)

clinical <- read_csv("data-raw/Clinical.csv")
msquant  <- read_csv("data-raw/Quantification.csv")
subject  <- read_csv("data-raw/SubjectLinks.csv")

# Remove duplicated columns
msquant_subject1 <- msquant %>% filter(!str_detect(rt, ":"))

# Add random noise to second subject and add subject id to col names
msquant_2nd <- msquant_subject1 %>% 
  mutate_at(vars(-rt, -mz), ~ ifelse(.x == 1, 1, .x + rnorm(1, 0, 100))) %>%
  rename_at(vars(-mz, -rt), ~ paste(.x, "02", sep = "_"))

# Add subject ID
msquant_subject1 <- msquant_subject1 %>%
  rename_at(vars(-mz, -rt), ~ paste(.x, "01", sep = "_"))

# Save 1-subject dataset
devtools::use_data(msquant_subject1, overwrite = TRUE)


# msquant_2nd$Neutral_Operator_Dif_Pos_1x_O1_A == msquant_subject1$Neutral_Operator_Dif_Pos_1x_O1_A

# Combine into one 2-subject dataset
msquant <-
  cbind(msquant_subject1, select(msquant_2nd, -rt, -mz))

# Save 2-subject dataset
devtools::use_data(msquant, overwrite = TRUE)



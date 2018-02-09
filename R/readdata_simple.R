#' 
#' Prepare a mass spec quantification data frame for filtering, imputation,
#' and normalization. Also provides summaries of data structure (replicates,
#' subjects, spike, etc.) 
#'
#' Function reads in raw data files and summarizes technical replicates as the
#' mean of observations for compounds found in 2 or 3 replicates and with
#' coefficient of variation below specified level, or median for those found in
#' 3 replicates but excess CV.
#'
#' @param clinical_data Name of the clinical data file.
#' @param quantification_data Name of the data file containing the
#' quantification data.
#' @param link_data Name of Subject link file.
#' @param rt Retention time.
#' @param mz Mass-to-charge ratio.
#' @param cvmax Acceptable level of coefficient of variation between replicates.
#' @param missing Value of missing data in the quantification data file.
#' @param linktxt Column name for Run ID field in Subject Link dataset
#' @return sum_data1 Matrix of summarized replicates, one obs per subject per
#' compound
#' @return clinical Clinical dataset
#' @return medians List of compounds that had excess CV and utilized the median
#' @note Must have unique compounds in raw data file.  Raw quantification data
#' should have two columns titled mz and rt that are combined to make the
#' column header. 
#' @examples
#'
#'   # LCMS_Run_ID = operator/replicate (A-C), subject (01-03), concentration  (1x,2x,4x)
#'   # SubjectID   = subject (01-03), concentration  (1x,2x,4x)
#'
#'   ### Read in data files
#'   clinical <- read.csv("./data-raw/Clinical.csv")
#'   quant    <- read.csv("./data-raw/Quantification.csv")
#'   link     <- read.csv("./data-raw/SubjectLinks.csv")
#'
#'   tidyquant  <- tidy_quant(quant, mz = "mz", rt = "rt")
#'   msprep_obj <- msprep(quant, 
#' 
#'   ### Set variables for program
#'   cvmax   <- 0.5
#'   missing <- 1
#'   linktxt <- "LCMS_Run_ID"
#'
#'   # Convert dataset to tidy-ish format
#'   #  need replicate, subject_id, mz, and rt as variables
#'
#'   test <- msprep(clinical_data       = clinical,
#'                  quantification_data = quant,
#'                  link_data           = link,
#'                  cvmax               = 0.50,
#'                  missing             = 1,
#'                  linktxt             = linktxt)
#'   #save(test, file = paste0("./data/test.Rdata"))
#'
#' str(quant)
#' str(test$sum_data)
#'
#'
#' save.image("readdata_workspace.RData")
#'
#'
#' ########################
#' # Start temporary code #
#'
#' # source 3 fns at bottom of script,
#' clinical      <- read.csv("./data-raw/Clinical.csv")
#' quant         <- read.csv("./data-raw/Quantification.csv")
#' link          <- read.csv("./data-raw/SubjectLinks.csv")
#' .data         <- tidy_quant(quant, mz = "mz", rt = "rt")
#' clinical_data <- clinical
#' link_data     <- link
#' cvmax         <- 0.50
#' missing_val   <- 1
#' linktxt       <- "LCMS_Run_ID"
#' # End temporary code #
#' ########################
#'
#'
#'
#'
#'
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @export
summarize_msprep <- function(.data,
                             clinical_data,
                             link_data,
                             subject_id,
                             replicate = "replicate",
                             mz,
                             rt,
                             cvmax   = 0.50,
                             missing_val = 1,
                             linktxt,
                             min_proportion_present = 1/3) {

  # Check args
  stopifnot(is.data.frame(.data))

  # Replace miss val with NAs 
  .data <- .data %>% mutate(abundance = replace_missing(abundance, missing_val))

  # Get replicate count for each mz/rt/spike/subject combo
  # NOTE: verify that replicate is a scalar for a given dataset
  replicate_count <- length(unique(.data[[replicate]]))

  # Check if all
  stopifnot(nrow(.data) %% replicate_count == 0)

  # Calculate initial summary measures
  quant_summary <-
    .data %>%
    group_by(subject_id, spike, mz, rt) %>%
    arrange(mz, rt, subject_id, spike, replicate) %>%
    summarise(n_present        = sum(!is.na(abundance)),
              prop_present     = n_present / replicate_count,
              mean_abundance   = mean(abundance, na.rm = T),
              sd_abundance     = sd(abundance, na.rm = T),
              median_abundance = median(abundance, na.rm = T)) %>%
    mutate(cv_abundance = sd_abundance / mean_abundance) %>%
    ungroup

  # Identify and select summary measure
  quant_summary <-
    quant_summary %>% 
    mutate(summary_measure   = select_summary_measure(n_present,
                                                      cv_abundance,
                                                      replicate_count,
                                                      min_proportion_present,
                                                      cvmax),
           abundance_summary = 
             case_when(
                       summary_measure == "median" ~ median_abundance,
                       summary_measure == "mean"   ~ mean_abundance,
                       is.na(summary_measure)      ~ 0
                       )
           )

    compounds <- quant_summary %>% select(mz, rt) %>% distinct %>% nrow
    # actually count of subject,spike pairs (not subjects only)
    subjects  <- quant_summary %>% select(subject_id, spike) %>% distinct %>% nrow


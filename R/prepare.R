#' 
#' Prepare a mass spec quantification data frame for filtering, imputation,
#' and normalization. Also provides summaries of data structure (replicates,
#' subjects, group_by, etc.)
#'
#' Function reads in raw data files and summarizes technical replicates as the
#' mean of observations for compounds found in 2 or 3 replicates and with
#' coefficient of variation below specified level, or median for those found in
#' 3 replicates but excess CV.
#'
#' @param data A tidy dataframe of quantification data.
#' @param subject_id Name (string) of the subject ID column.
#' @param replicate Name (string) of the replicate column. Set to NULL if no
#' replicates. TODO: test NULL.
#' @param abundance Name (string) of the abundance column.
#' @param group_by Name (string) of the group_by column.
#' @param mz Mass-to-charge ratio variable name.
#' @param rt Retention time variable name.
#' @param cvmax Acceptable level of coefficient of variation between replicates.
#' @param missing_val Value of missing data in the quantification data file.
#' @param min_proportion_present  Minimum proportion present to summarize with
#' median or mean. Below this the compound will be set to 0.
#' @return An `msprep` object with `stage(rtn) == "prepared"` containing
#' summarised quantification data, a dataset of compounds summarised by medians,
#' and other related summaries.
#' @examples
#'
#' # Read in data file
#' data(msquant)
#' 
#' # Convert dataset to tidy format
#' tidy_data    <- ms_tidy(msquant, mz = "mz", rt = "rt")
#' prepped_data <- ms_prepare(tidy_data, replicate = "replicate")
#' 
#' # Or, using tidyverse/magrittr pipes 
#' library(magrittr)
#' prepped_data <- msquant %>% ms_tidy %>% ms_prepare
#'
#' str(prepped_data)
#' str(prepped_data$summary_data)
#'
#' @importFrom dplyr select
#' @importFrom dplyr select_at
#' @importFrom dplyr mutate 
#' @importFrom dplyr mutate_at
#' @importFrom dplyr group_by 
#' @importFrom dplyr arrange 
#' @importFrom dplyr summarise 
#' @importFrom dplyr ungroup 
#' @importFrom dplyr case_when
#' @importFrom dplyr distinct
#' @importFrom dplyr filter
#' @importFrom dplyr vars
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @importFrom rlang UQ
#' @importFrom rlang sym
#' @importFrom rlang as_string
#' @importFrom stats median
#' @importFrom stats sd
#' @export
ms_prepare <- function(data,
                       subject_id  = "subject_id",
                       replicate   = "replicate",
                       abundance   = "abundance",
                       group_by    = "spike", 
                       mz          = "mz",
                       rt          = "rt",
                       cvmax       = 0.50,
                       missing_val = 1,
                       min_proportion_present = 1/3) {

  # Check args
  stopifnot(is.data.frame(data))
  my_args  <- mget(names(formals()), sys.frame(sys.nframe()))

  # rlang magic 
  group_by = sym(group_by)

  # Convert to tibble data frame
  data <- as_tibble(data)

  # Replace provided variable names with standardized ones
  data <- standardize_dataset(data, subject_id, replicate, abundance, group_by, mz, rt)

  # Replace miss val with NAs 
  data <- mutate_at(data, vars("abundance"), replace_missing, missing_val)

  # Get replicate count for each mz/rt/group_by/subject combo
  replicate_count <- length(unique(data[["replicate"]]))

  # Roughly check if all compounds are present in each replicate
  stopifnot(nrow(data) %% replicate_count == 0)

  # Calculate initial summary measures 
  #   Note;(matrix algebra would be faster --
  #     consider later)
  quant_summary <- group_by(data, subject_id, UQ(group_by), mz, rt)
  quant_summary <- arrange(quant_summary, mz, rt, subject_id, UQ(group_by), replicate)
#   quant_summary <- summarise(quant_summary,
#                              "n_present"        = sum(!is.na(abundance)),
#                              "prop_present"     = n_present / replicate_count,
#                              "mean_abundance"   = mean(abundance, na.rm = T),
#                              "sd_abundance"     = sd(abundance, na.rm = T),
#                              "median_abundance" = median(abundance, na.rm = T))
#   quant_summary <- mutate(quant_summary, cv_abundance = sd_abundance / mean_abundance)
#   quant_summary <- ungroup(quant_summary)
# 
#   # Identify and select summary measure -- TODO: decompose this to get rid of 'no
#   #                                              visible binding error'
#   quant_summary <-
#     mutate(quant_summary,
#            summary_measure = 
#              select_summary_measure(.data$n_present, .data$cv_abundance, replicate_count,
#                                     min_proportion_present, cvmax),
#            abundance_summary = 
#              case_when(summary_measure == "median" ~ .data$median_abundance,
#                        summary_measure == "mean"   ~ .data$mean_abundance,
#                        TRUE                        ~ 0))
#   quant_summary <- mutate_at(quant_summary,
#                              vars(subject_id, "summary_measure", UQ(group_by)),
#                              factor)
# 
#   # Extract summarized dataset
#   summary_data  <- select_at(quant_summary, 
#                              vars(subject_id, UQ(group_by), "mz", "rt", "abundance_summary"))
# 
#   # Additional info extracted in summarizing replicates
#   replicate_info <- select_at(quant_summary, 
#                               vars(subject_id, UQ(group_by), "mz", "rt", "n_present",
#                                    "cv_abundance", "summary_measure"))
# 
#   # Summaries that used medians
#   medians        <- filter(quant_summary, summary_measure == "median")
#   medians        <- select_at(medians, vars(subject_id, UQ(group_by), "mz", "rt",
#                                             "abundance_summary"))
# 
#   # Total number of compounds identified
#   n_compounds    <- nrow(distinct(select(quant_summary, "mz", "rt")))
# 
#   # Create return object & return
#   structure(list("data" = summary_data,
#                  replicate_info  = replicate_info,
#                  medians         = medians),
#             replicate_count = replicate_count,
#             cvmax           = cvmax,
#             min_proportion_present = min_proportion_present,
#             group_by = as_string(group_by),
#             class = "msprep",
#             stage = "prepared")

}

print.msprep <- function(x) {

  stage <- stage(x)

  cat("msprep dataset\n")
  cat(paste("    Stage:", stage, "\n"))
  cat("    Replicate count: ", attr(x, "replicate_count"), "\n")
  cat("    Patient count: ", length(unique(x$data$subject_id)), "\n")
  cat("    Count of patient-compounds summarized by median: ", nrow(x$medians), "\n")
  cat("    Prepare summary: \n")
  cat("        User-defined parameters: \n")
  cat("          cvmax = ", attr(x, "cvmax"), "\n")
  cat("          min_proportion_present = ", round(attr(x, "min_proportion_present"), digits=3), "\n")
  if (stage %in% msprep_stages()[2:4]) {
    cat("    Filter summary:\n")
    cat("      User-defined parameters: \n")
    cat("        filter percent = ", attr(x, "filter_percent"), "\n")
    cat("      Resulting stats: \n")
    cat("        Count of compounds present in >= ", 
        round(attr(x, "filter_percent")*100, digits = 3), "% of patients = ", 
        sum(attr(x, "filter_status")$keep), "\n")
  }

  cat("    Dataset:\n")
  print(x$data, n = 6)

}



#
# Internal ms_prepare functions
#

#' @importFrom dplyr case_when
select_summary_measure <- function(n_present,
                                   cv_abundance,
                                   n_replicates,
                                   min_proportion_present,
                                   cv_max) {

  case_when((n_present / n_replicates) <= min_proportion_present
              ~ "none: proportion present <= min_proportion_present",
            cv_abundance > cv_max & (n_replicates == 3 & n_present == 2) 
              ~ "none: cv > cvmax & 2 present",
            cv_abundance > cv_max & (n_present == n_replicates) 
              ~ "median",
            TRUE ~ "mean")

}

#' @importFrom dplyr rename
#' @importFrom rlang sym
#' @importFrom rlang UQ
standardize_dataset <- function(data, subject_id, replicate, abundance, group_by,
                                mz, rt) {

  # Rename required variables
  subject_id = sym(subject_id)
  abundance  = sym(abundance)
  mz         = sym(mz)
  rt         = sym(rt)

  data <- data %>% 
    rename("subject_id" = UQ(subject_id),
           "abundance"  = UQ(abundance),
           "mz"         = UQ(mz),
           "rt"         = UQ(rt))

  # Rename optional variables if present
  if (!is.null(replicate)) {
    replicate <- sym(replicate)
    data      <- data %>% rename("replicate" = UQ(replicate))
  } else {
    data$replicate <- "None"
  }

  #if (!is.null(group_by)) {
  #  group_by <- sym(group_by)
  #  data     <- data %>% rename("group_by" = UQ(group_by)) 
  #} else { 
  #  data$group_by <- "Not provided" 
  #}

  return(data)

}




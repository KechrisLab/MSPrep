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
#' @param data A tidy dataframe of quantification data.
#' @param rt Retention time variable name.
#' @param mz Mass-to-charge ratio variable name.
#' @param cvmax Acceptable level of coefficient of variation between replicates.
#' @param missing Value of missing data in the quantification data file.
#' @param linktxt Column name for Run ID field in Subject Link dataset
#' @return An \code{msprepped} object containing summarised quantification data,
#' a dataset of compounds summarised by medians, and other related summaries.
#' @examples
#'
#' # Read in data file
#' quant <- read.csv("./data-raw/Quantification.csv")
#' 
#' # Convert dataset to tidy format
#' tidy_data    <- tidy_ms(quant, mz = "mz", rt = "rt")
#' prepped_data <- prepare_ms(tidy_data)
#' 
#' # Or, using tidyverse/magrittr pipes 
#' prepped_data <- quant %>% tidy_ms %>% prepare_ms
#'
#' str(prepped_data)
#' str(prepped_data$summary_data)
#'
#' @importFrom dplyr select
#' @importFrom dplyr mutate 
#' @importFrom dplyr mutate_at
#' @importFrom dplyr group_by 
#' @importFrom dplyr arrange 
#' @importFrom dplyr summarise 
#' @importFrom dplyr ungroup 
#' @importFrom dplyr case_when
#' @importFrom dplyr distinct
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
prepare_ms <- function(data,
                       subject_id  = "subject_id",
                       replicate   = NULL,
                       abundance   = "abundance",
                       spike       = "spike",
                       mz          = "mz",
                       rt          = "rt",
                       cvmax       = 0.50,
                       missing_val = 1,
                       min_proportion_present = 1/3) {

  # Check args
  stopifnot(is.data.frame(data))
  my_args  <- mget(names(formals()), sys.frame(sys.nframe()))

  # Replace provided variable names with standardized ones
  data <- standardize_dataset(data, subject_id, replicate, abundance, spike, mz, rt)

  # Replace miss val with NAs 
  data <- data %>% mutate_at(vars(abundance), replace_missing, missing_val)

  # Get replicate count for each mz/rt/spike/subject combo
  replicate_count <- length(unique(data[["replicate"]]))

  # Roughly check if all compounds are present in each replicate
  stopifnot(nrow(data) %% replicate_count == 0)

  # Calculate initial summary measures 
  #   Note;(matrix algebra would be faster --
  #     consider later)

  quant_summary <-
    data %>% 
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
    mutate(summary_measure = 
             select_summary_measure(n_present, cv_abundance, replicate_count,
                                    min_proportion_present, cvmax),
           abundance_summary = 
             case_when(summary_measure == "median" ~ median_abundance,
                       summary_measure == "mean"   ~ mean_abundance,
                       TRUE                        ~ 0)) %>%
    mutate_at(c("subject_id", "summary_measure", "spike"), factor)

  # Extract summarized dataset
  summary_data   <-
    quant_summary %>% select(subject_id, spike, mz, rt, abundance_summary)

  # Additional info extracted in summarizing replicates
  replicate_info <-
    quant_summary %>% select(subject_id, spike, mz, rt, n_present, cv_abundance, summary_measure)

  # Summaries that used medians
  medians        <-
    quant_summary %>% filter(summary_measure == "median") %>%
    select(subject_id, spike, mz, rt, abundance_summary)

  # Total number of compounds identified
  n_compounds    <-
    quant_summary %>% select(mz, rt) %>% distinct %>% nrow

  # Count of subject,spike pairs
  n_subject_spike_pairs  <-
    quant_summary %>% select(subject_id, spike) %>% distinct %>% nrow

  # Subject and spike id combos
  subjects_summary <- 
    quant_summary %>% select(subject_id, spike) %>% distinct

  # Create return object & return
  structure(list(summary_data    = summary_data,
                 replicate_info  = replicate_info,
                 clinical        = subjects_summary,
                 medians         = medians,
                 replicate_count = replicate_count,
                 cvmax           = cvmax,
                 min_proportion_present = min_proportion_present),
            class = "msprepped")

}

#' @rdname prepare_ms
print.msprepped <- function(x) {
  cat("msprepped object\n")
  cat("    Replicate count: ", x$replicate_count, "\n")
  cat("    Patient count: ", length(unique(x$clinical$subject_id)), "\n")
  cat("    Count of spike levels: ", length(unique(x$clinical$spike)), "\n")
  cat("    Count patient-spike combinations: ", nrow(x$clinical), "\n")
  cat("    Count of patient-spike compounds summarized by median: ", nrow(x$medians), "\n")
  cat("    User-defined parameters \n")
  cat("        cvmax = ", x$cvmax, "\n")
  cat("        min_proportion_present = ", round(x$min_proportion_present, digits=3), "\n")
  cat("    Summarized dataset:\n")
  print(x$summary_data, n = 6)
}




#' Function for converting wide mass spec quantification data into a tidy data
#' frame
#'
#' Function reads in wide dataset of mass spec quantification data and converts
#' it to a tidy dataset.  This function assumes that the dataset is in a wide
#' format, with a column representing the retention time (rt), another
#' representing the mass-to-charge ratio (mz), and the remaining columns
#' containing MS quantification data.  
#' 
#' It also assumes that the column names of quantification data start with some
#' consistent, informational but unnecessary text (id_extra_txt), and contain
#' the spike, subject ID, and replicate ID in a consistent position, all
#' separated by a consistent separator.
#'
#' See \code{data(quantification)} for an example.  If your data doesn't fit
#' this format, view the function code for hints on tidying your data.  The core
#' of this function consists of \code{tidyr::gather()}, \code{dplyr::mutate()},
#' and \code{tidyr::separate()}.
#'
#' @param quantification_data Data frame containing the quantification data.
#' @param mz Name of the column containing mass-to-charge ratios.
#' @param rt Name of the column containing retention time.
#' @param id_extra_txt Text to remove when converting quant variable names to
#'   variables.
#' @param separator Character/string separating spike, subject and replicate
#'   ids.
#' @param id_spike_pos Order in which the spike number occurs in the quant
#'   variable names (after \code{id_extra_txt} is removed).
#' @param id_subject_id_pos Order in which the subject id occurs in the quant
#'   variable names.
#' @param id_replicate_id_pos Order in which the replicate id occurs in the quant
#'   variable names.
#'
#' @return A tidy data frame of quant data, with columns mz, rt, spike,
#' subject_id, replicate, and abundance.
#'
#' @examples
#'
#'   quant <- read.csv("./data-raw/Quantification.csv")
#'   quant <- tidy_ms(quant, mz = "mz", rt = "rt")
#'
#' @importFrom tibble as_data_frame
#' @importFrom tibble data_frame
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom tidyr separate
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @export
tidy_ms <- function(quantification_data,
                    mz = "mz", 
                    rt = "rt",
                    id_extra_txt     = "Neutral_Operator_Dif_Pos_",
                    separator        = "_",
                    id_spike_pos     = 1,
                    id_subject_pos   = 2,
                    id_replicate_pos = 3) {

  into <-
    data_frame(name = c("spike", "subject_id", "replicate"),
               order = c(id_spike_pos, id_subject_pos, id_replicate_pos)) %>%
    arrange(order) %>% .[["name"]]

  rtn <- 
    quantification_data %>%
      tibble::as_data_frame(.) %>%
      tidyr::gather(key = id_col, value = abundance, -mz, -rt) %>%
      dplyr::mutate(id_col = stringr::str_replace_all(id_col, 
                                                      id_extra_txt, "")) %>%
      tidyr::separate(id_col, sep = separator, into = into)

  return(rtn)

}




#' @rdname prepare_ms
replace_missing <- function(abundance, missing_val) {
  ifelse(abundance == missing_val, NA, abundance)
}

#' @rdname prepare_ms
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

#' @rdname prepare_ms
#' @importFrom dplyr rename
#' @importFrom rlang sym
#' @importFrom rlang UQ
standardize_dataset <- function(data, subject_id, replicate, abundance, spike,
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
    replicate  = sym(replicate)
    data <- data %>%
      rename("replicate" = UQ(replicate))
  } else {
    data$replicate <- "None"
  }

  if (!is.null(spike)) {
    spike  = sym(spike)
    data <- data %>%
      rename("spike" = UQ(spike))
  } else {
    data$spike <- "Not provided"
  }

  return(data)

}




#' 
#' Prepare a mass spec quantification data frame for filtering, imputation,
#' and normalization. 
#'
#' Also provides summaries of data structure (replicates,
#' subjects, spike, etc.). 
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
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @export









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
#'   quant     <- read.csv("./data-raw/Quantification.csv")
#'   tidyquant <- tidy_quant(quant, mz = "mz", rt = "rt")
#'
#' @importFrom tibble as_data_frame
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom tidyr separate
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @export
tidy_quant <- function(quantification_data, mz, rt,
                       id_extra_txt = "Neutral_Operator_Dif_Pos_",
                       separator    = "_",
                       id_spike_pos        = 1,
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


#' @rdname prepare
replace_missing <- function(abundance, missing_val) {
  ifelse(abundance == missing, NA, abundance)
}

#' @rdname prepare
mean_or_median <- function(cv_abundance, cv_max, presence_count, replicate_count) {
  ifelse(cv_abundance > cvmax & presence_count > ceiling(replicate_count/2), "median", "mean")
}


select_summary_measure <- function(min_proportion_present,
                                   cv_abundance,
                                   cv_max,
                                   n_replicates,
                                   n_present) {

  if ((n_present / n_replicates) <= min_proportion_present) {
    return( NA )
  } else if (cv_abundance > cv_max & (n_replicates == 3 & n_present == 2)) {
    return( NA )
  } else if (cv_abundance > cv_max & (n_present == n_replicates)) {
    return( "median" )
  } else {
    return( "mean" )
  }

}

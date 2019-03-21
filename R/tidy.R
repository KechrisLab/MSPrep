#' Function for converting wide mass spectrometry quantification data into a tidy data
#' frame.
#'
#' Function reads in wide dataset of mass spectrometry quantification data and converts
#' it to a tidy dataset.  This function assumes that the dataset is in a wide
#' format, with a column representing the retention time (rt), another
#' representing the mass-to-charge ratio (mz), and the remaining columns
#' containing MS quantification data.  
#' 
#' It also assumes that the column names of quantification data start with some
#' consistent, informational but unnecessary text (col_extra_txt), and contain
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
#' @param col_extra_txt Text to remove from column names when converting column names to
#'   variables.
#' @param col_names Vector of the ordered ID names to extract from the variable
#' names.
#' @param separator Character/string separating spike, subject, and replicate
#' ids in column names.
#'
#' @return A tidy data frame of quant data, with columns mz, rt,
#' replicate, and abundance.
#'
#' @examples
#'
#'   # load raw abundance/quantification dataset
#'   data(msquant)
#'
#'   # print original dataset (using tibble's nice printing)
#'   tibble::as_data_frame(msquant) 
#'
#'   # tidy the dataset and print it
#'   ms_tidy(msquant, mz = "mz", rt = "rt")
#'
#' @importFrom tibble as_data_frame
#' @importFrom tibble data_frame
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom tidyr separate
#' @importFrom tidyr unite
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @importFrom rlang UQ
#' @export
ms_tidy <- function(quantification_data,
                    mz = "mz", 
                    rt = "rt",
                    col_extra_txt     = "Neutral_Operator_Dif_Pos_",
                    col_names         = c("spike", "batch", "replicate", "subject_id"),
                    separator        = "_") {

  # gather data to long format (adds id/varnames as column), remove col_extra_txt
  # from id column, and convert id column to separate subject,replicate
  # variables
  rtn <-
    as_data_frame(quantification_data) %>%
    gather(key = "id_col", value = "abundance", -mz, -rt) %>%
    mutate_at(vars("mz", "rt"), as.numeric) %>%
    mutate(id_col = str_replace_all(.data$id_col, col_extra_txt, ""))
  # Split and recombine id names 
  rtn <- separate(rtn, .data$id_col, sep = separator, into = col_names)

  return(rtn)

}


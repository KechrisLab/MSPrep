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
#' @param id_subject_pos Order in which the subject id occurs in the quant
#'   variable names.
#' @param id_replicate_pos Order in which the replicate id occurs in the quant
#'   variable names.
#'
#' @return A tidy data frame of quant data, with columns mz, rt, spike,
#' subject_id, replicate, and abundance.
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
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @export
ms_tidy <- function(quantification_data,
                    mz = "mz", 
                    rt = "rt",
                    id_extra_txt     = "Neutral_Operator_Dif_Pos_",
                    separator        = "_",
                    id_spike_pos     = 1,
                    id_subject_pos   = 2,
                    id_replicate_pos = 3) {

  # create vector of order of id parts
  into <-
    data_frame(name = c("spike", "subject_id", "replicate"),
               order = c(id_spike_pos, id_subject_pos, id_replicate_pos)) %>%
    arrange(order) %>% .[["name"]]

  # gather data to long format (adds id/varnames as column), remove id_extra_txt
  # from id column, and convert id column to separate spike,subject,replicate
  # variables
  rtn <-
    tibble::as_data_frame(quantification_data) %>%
    tidyr::gather(key = id_col, value = abundance, -mz, -rt) %>%
    dplyr::mutate(id_col = stringr::str_replace_all(id_col, id_extra_txt, "")) %>%
    tidyr::separate(id_col, sep = separator, into = into)

  return(rtn)

}


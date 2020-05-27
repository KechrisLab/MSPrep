#' Function for converting wide mass spectrometry quantification data into a tidy data
#' frame.
#'
#' Function reads in wide dataset of mass spectrometry quantification data and converts
#' it to a tidy dataset.  This function assumes that the dataset is in a wide
#' format, with a combination of columns representing the retention time (rt), another
#' representing the mass-to-charge ratio (mz), and the remaining columns
#' containing MS quantification data. Optionally, the data may also include a column
#' specifying compound name (met_id) as an addition to or replacement of rt/mz columns. 
#'
#' @param quantification_data Data frame containing the quantification data.
#' @param met_id Name of the column containing compound names.
#' @param mz Name of the column containing mass-to-charge ratios.
#' @param rt Name of the column containing retention time.
#' @param col_extra_txt Text to remove from column names when converting column names to
#'   variables.
#' @param col_names Vector of the ordered ID names to extract from the variable
#' names.
#' @param separator Character or text separating spike, subject, and replicate
#' ids in column names.
#' 
#' @details 
#' Function reads in wide dataset of mass spectrometry quantification data and converts
#' it to a tidy dataset.  This function assumes that the dataset is in a wide
#' format, with a combination of columns representing the retention time (rt), another
#' representing the mass-to-charge ratio (mz), and the remaining columns
#' containing MS quantification data. Optionally, the data may also include a column
#' specifying compound name (met_id) as an addition to or replacement of rt/mz columns.
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
#' @return A tidy data frame of quant data, with columns mz, rt,
#' replicate, and abundance.
#'
#' @examples
#'
#'   # load raw abundance/quantification dataset
#'   data(msquant)
#'   
#'   # print original dataset (using tibble's nice printing)
#'   tibble::as_tibble(msquant) 
#'   
#'   # tidy the dataset and print it
#'   ms_tidy(msquant, mz = "mz", rt = "rt", 
#'           col_extra_txt = "Neutral_Operator_Dif_Pos_", 
#'           separator = "_", 
#'           col_names = c("spike", "batch", "replicate", "subject_id"))
#'
#' @importFrom tibble as_data_frame
#' @importFrom tibble data_frame
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_at
#' @importFrom dplyr arrange
#' @importFrom tidyr separate
#' @importFrom tidyr unite
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @importFrom rlang UQ
#' @importFrom rlang !!!
#' @importFrom rlang syms
#' @export
ms_tidy <- function(quantification_data,
                    met_id = NULL,
                    mz = NULL, 
                    rt = NULL,
                    col_extra_txt = NULL,
                    col_names = c("subject_id"),
                    separator = NULL) {
  
  # Check that at either met_id or both mz and rt are included
    if(is.null(met_id) & (is.null(mz) | is.null(rt))){
        stop("Must include met_id or both mz and rt")
    }
    
    # Store whatever metabolite id args are present
    met_vars <- syms(c(met_id, mz, rt))
    
    # Gather data to long format (adds id/varnames as column), ensure mz and rt
    # are numeric if present
    rtn <- as_tibble(quantification_data) %>%
        gather(key = "id_col", value = "abundance", -c(!!!met_vars))
    
    
    # Ensure mz and rt are numeric if present
    if (!is.null(mz) & !is.null(rt)) {
        rtn <- mutate_at(rtn, vars(mz, rt), as.numeric)
    }
    
    # Remove col_extra_txt text if present
    if (!is.null(col_extra_txt)) {
        rtn <- mutate(rtn, id_col = str_replace_all(.data$id_col,
                                                    col_extra_txt, ""))
    }
  
    # If only one column name, rename id_col appropriately. Otherwise split
    # 'id_col' into new variable columns designated by col_names 
    if (length(col_names) == 1) {
        colnames(rtn)[colnames(rtn) == "id_col"] <- col_names[1]
    }
    else if (!is.null(separator)) {
        rtn <- separate(rtn, .data$id_col, sep = separator, into = col_names)
    }
    else {
        stop("Must include 'separator' if multiple 'col_names'")
    }

    return(rtn)
}


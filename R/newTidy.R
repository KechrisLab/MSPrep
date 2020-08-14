#' Function for converting wide mass spectrometry abundance data into a tidy data
#' frame.
#'
#' Function reads in wide dataset of mass spectrometry abundance data and converts
#' it to a tidy dataset.  This function assumes that the dataset is in a wide
#' format, with a combination of columns representing the retention time (rt), another
#' representing the mass-to-charge ratio (mz), and the remaining columns
#' containing MS abundance data. Optionally, the data may also include a column
#' specifying compound name (compID) as an addition to or replacement of rt/mz columns. 
#'
#' @param data Data frame or SummarizedExperiment containing the abundance data.
#' @param compID Name of the column containing compound names.
#' @param mz Name of the column containing mass-to-charge ratios.
#' @param rt Name of the column containing retention time.
#' @param colExtraText Text to remove from column names when converting column names to
#'   variables.
#' @param colNames Vector of the ordered ID names to extract from the variable
#' names.
#' @param separator Character or text separating spike, subject, and replicate
#' ids in column names.
#' 
#' @details 
#' Function reads in wide dataset of mass spectrometry abundance data and converts
#' it to a tidy dataset.  This function assumes that the dataset is in a wide
#' format, with a combination of columns representing the retention time (rt), another
#' representing the mass-to-charge ratio (mz), and the remaining columns
#' containing MS abundance data. Optionally, the data may also include a column
#' specifying compound name (compID) as an addition to or replacement of rt/mz columns.
#' 
#' It also assumes that the column names of abundance data start with some
#' consistent, informational but unnecessary text (colExtraText), and contain
#' the spike, subject ID, and replicate ID in a consistent position, all
#' separated by a consistent separator.
#'
#' See \code{data(abundance)} for an example.  If your data doesn't fit
#' this format, view the function code for hints on tidying your data.  The core
#' of this function consists of \code{tidyr::gather()}, \code{dplyr::mutate()},
#' and \code{tidyr::separate()}.
#'
#' @return A tidy data frame of quant data, with columns mz, rt,
#' replicate, and abundance.
#'
#' @examples
#'
#'   # load raw abundance/abundance dataset
#'   data(msquant)
#'   
#'   # print original dataset (using tibble's nice printing)
#'   tibble::as_tibble(msquant) 
#'   
#'   # tidy the dataset and print it
#'   ms_tidy(msquant, mz = "mz", rt = "rt", 
#'           colExtraText = "Neutral_Operator_Dif_Pos_", 
#'           separator = "_", 
#'           colNames = c("spike", "batch", "replicate", "subject_id"))
#'
#' @importFrom tibble as_data_frame
#' @importFrom tibble data_frame
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_at
#' @importFrom dplyr arrange
#' @importFrom tidyr separate
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @importFrom rlang UQ
#' @importFrom rlang !!!
#' @importFrom rlang syms
#' @import SummarizedExperiment
msTidy <- function(data,
                   compVars = c("mz", "rt"),
                   sampleVars = c("subject_id"),
                   colExtraText = NULL,
                   separator = NULL,
                   missingValue = NA) {
    
    ## Check data structure
    if (is(data, "SummarizedExperiment")) {
        
        rtn <- .seTidy(data)
        
    } else if (is(data, "data.frame")) {

        rtn <- .dfTidy(data,
                       compVars = compVars,
                       colExtraText = colExtraText,
                       sampleVars = sampleVars,
                       separator = separator)
    } else {
        
        stop("'data' must be a data frame or SummarizedExperiment")
        
    }
    
    ## Replace miss val with NAs (if not already)
    rtn <- mutate_at(rtn, vars("abundance"), .replaceMissing, missingValue)
    
}


#' Reformats SummarizedExperiments into tidy data frames
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr as_tibble
#' @importFrom dplyr bind_cols
#' @importFrom dplyr left_join
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
## Adapted from 
## https://github.com/pkimes/upbm/blob/master/R/tidy-PBMExperiment.R
.seTidy <- function(SE) {
        
        ## Get row data, compound variables, and technical variables
        rowDat <- as_tibble(rowData(SE))
        compVars <- colnames(rowDat)
        sampleVars <- colnames(colData(SE))

        ## Extract abundances, join with row data
        rtn <- SummarizedExperiment::assay(SE) %>%
            as_tibble() %>%
            bind_cols(rowDat)
        
        ## Make long data frame
        rtn <- gather(rtn, column_name, abundance, colnames(SE))
        
        ## Extract column data
        coldat <- as.data.frame(colData(SE), optional = TRUE) %>%
            rownames_to_column("column_name")
        
        ## Join abundance data with column data
        rtn <- left_join(rtn, coldat, by = "column_name")
        
        ## Reorder columns for cleaner appearance
        rtn <- select(rtn, c(column_name, sampleVars, compVars, abundance))
        
        return(rtn)
}

## Function to convert matrix to tidy data frame
.dfTidy <- function(data, compVars, colExtraText, sampleVars,
                    separator) {
    
    # ## Check that at either compID or both mz and rt are included
    # if (is.null(compID) & (is.null(mz) | is.null(rt))){
    #     stop("Must include compID or both mz and rt")
    # }
    
    ## Store whatever metabolite id args are present
    compVars <- syms(compVars)
    
    ## Gather data to long format (adds id/varnames as column), ensure mz and rt
    ## are numeric if present
    rtn <- as_tibble(data) %>%
        gather(key = "column_name", value = "abundance", -c(!!!compVars))
    
    # ## Ensure mz and rt are numeric if present
    # if (!is.null(mz) & !is.null(rt)) {
    #     rtn <- mutate_at(rtn, vars(mz, rt), as.numeric)
    # }
    
    ## Remove colExtraText text if present
    if(!is.null(colExtraText)){
        rtn <- mutate(rtn, id_col = str_replace_all(.data$column_name, 
                                                    colExtraText, ""))
    } else {
        rtn <- mutate(rtn, id_col = .data$column_name)
    }
    
    ## If only one column name, rename id_col appropriately. Otherwise split
    ##  'id_col' into new variable columns designated by colNames 
    if (length(sampleVars) == 1) { 
        colnames(rtn)[colnames(rtn) == "id_col"] <- sampleVars[1]
    } else if (!is.null(separator)) {
        rtn <- separate(rtn, .data$id_col, sep = separator, into = sampleVars)
    } else {
        stop("Must include 'separator' if multiple 'sampleVars'")
    }
    
    ## Reorder columns for cleaner appearance
    rtn <- select(rtn, c(column_name, sampleVars, !!!compVars, abundance))
}

## Replace missing values with NA
.replaceMissing <- function(abundance, missingValue) {
    if (is.na(missingValue)) {
        return(abundance)
    } else {
        ifelse(abundance == missingValue, NA, abundance)
    }
}
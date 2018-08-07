################################################################################
# Internal utility functions
################################################################################


# Three fns for assigning and getting stage attribute (for tracking object
# progress through pipeline)
msprep_stages <- function() {
  c("prepared", "filtered", "imputed", "normalized")
}

stage <- function(x) attr(x, "stage")

`stage<-` <- function(x, value) {
  attr(x,  "stage") <- value
  x
}


# Functions for assigning and getting grouping_vars attribute (e.g. spike)
grouping_vars <- function(x) attr(x, "grouping_vars")

`grouping_vars<-` <- function(x, value) {
  attr(x, "grouping_vars") <- value
  x
}


# Standard arrange for data object used in msprep_obj
ms_arrange <- function(data, ...) arrange(data, .data$subject_id, .data$mz, .data$rt, ...)


valid_cols <- function(grouping_vars)  {
  grouping_vars <- as.character(grouping_vars)
  sort(c("subject_id", "mz", "rt", "abundance_summary", grouping_vars))
}



# Internal function - replace missing val with NA
replace_missing <- function(abundance, missing_val) {
  ifelse(abundance == missing_val, NA, abundance)
}



# Internal function -- spread data
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr unite
#' @importFrom tidyr spread
data_to_wide_matrix <- function(data, grouping_vars) {

  data %>% 
    unite(col = "compound", "mz", "rt") %>% 
    spread(key = "compound", value = "abundance_summary")  %>%
    as.data.frame %>%
    unite(col = "rwnm", "subject_id", grouping_vars) %>%
    column_to_rownames(var = "rwnm") %>%
    as.matrix

}

# Internal function -- undo spread
#' @importFrom rlang sym
#' @importFrom rlang UQ
#' @importFrom tibble rownames_to_column
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
#' @importFrom tidyr separate
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate_at
#' @importFrom dplyr vars
wide_matrix_to_data <- function(wide, grouping_vars) {

  sym_id <- sym("subject_id")
  sym_mz <- sym("mz")
  sym_rt <- sym("rt")

  rtn <- wide %>% as.data.frame
  rtn <- rtn %>% rownames_to_column(var = "rwnm")
  rtn <- rtn %>% separate("rwnm", sep = "_", into = c("subject_id", grouping_vars))
  rtn <- rtn %>% as_tibble
  rtn <- rtn %>% gather(key = "mz_rt", value = "abundance_summary", -"subject_id", -grouping_vars)
  rtn <- rtn %>% separate("mz_rt", sep = "_", into = c("mz", "rt"))
  rtn <- standardize_datatypes(rtn, grouping_vars)
  rtn <- ms_arrange(rtn)

  return(rtn)

}




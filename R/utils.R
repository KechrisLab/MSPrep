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
grouping_vars <- function(x) attr(x, "groupingvars")

`grouping_vars<-` <- function(x, value) {
  attr(x, "groupingvars") <- value
  x
}

# Functions for assigning and getting batch attribute (e.g. operator)
batch_var <- function(x) attr(x, "batch")

`batch_var<-` <- function(x, value) {
  attr(x, "batch_var") <- value
  x
}

# Functions for assigning and getting replicate attribute
replicate_var <- function(x) attr(x, "replicate")

`replicate_var<-` <- function(x, value) {
  attr(x, "replicate_var") <- value
  x
}

# Standard arrange for data object used in msprep_obj
#' @importFrom dplyr arrange
ms_arrange <- function(data, ...) {
  other_vars <- list(...)
  if (("mz" %in% colnames(data)) & ("rt" %in% colnames(data)) & 
      ("met_id" %in% colnames(data))) {
    rtn <- arrange(data, .data$subject_id, .data$mz, .data$rt, .data$met_id, 
                   `!!!`(syms(other_vars)))
  } else if (("mz" %in% colnames(data)) & ("rt" %in% colnames(data))) {
    rtn <- arrange(data, .data$subject_id, .data$mz, .data$rt, 
                   `!!!`(syms(other_vars)))
  } else {
    rtn <- arrange(data, .data$subject_id, .data$met_id, 
                   `!!!`(syms(other_vars)))
  }
  return(rtn)
}


valid_cols <- function(extra_vars)  {
  extra_vars <- as.character(extra_vars)
  sort(c("subject_id", "mz", "rt", "abundance_summary", extra_vars))
}



# Internal function - replace missing val with NA
replace_missing <- function(abundance, missing_val) {
  ifelse(abundance == missing_val, NA, abundance)
}



# Internal function -- spread data
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr unite
#' @importFrom tidyr spread
data_to_wide_matrix <- function(data, groupingvars, batch, asmatrix = TRUE) {

  internal_id <- internal_id_order(groupingvars, batch)

  rtn <- data %>% 
    unite(col = "compound", "mz", "rt") %>% 
    spread(key = "compound", value = "abundance_summary")  %>%
    as.data.frame %>%
    unite(col = "rwnm", internal_id) %>%
    column_to_rownames(var = "rwnm") 

  if (asmatrix) rtn <- as.matrix(rtn)

  return(rtn)

}

internal_id_order <- function(groupingvars = NULL, batch = NULL) {
  c("subject_id", batch, groupingvars)
}

# Internal function -- undo spread
#' @importFrom rlang sym
#' @importFrom rlang UQ
#' @importFrom tibble rownames_to_column
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate_at
#' @importFrom dplyr vars
wide_matrix_to_data <- function(wide, groupingvars, batch) {

  sym_id <- sym("subject_id")
  sym_mz <- sym("mz")
  sym_rt <- sym("rt")

  internal_id <- internal_id_order(groupingvars, batch)

  rtn <- wide %>% as.data.frame
  rtn <- rtn %>% rownames_to_column(var = "rwnm")
  rtn <- rtn %>% separate("rwnm", sep = "_", into = internal_id)
  rtn <- rtn %>% as_tibble
  if(!is.null(groupingvars) & !is.null(batch)){
    rtn <- rtn %>% gather(key = "mz_rt", value = "abundance_summary",
                        -"subject_id", -groupingvars, -batch)
  }
  else
    rtn <- rtn %>% gather(key = "mz_rt", value = "abundance_summary", 
                          -"subject_id")
  rtn <- rtn %>% separate("mz_rt", sep = "_", into = c("mz", "rt"))
  rtn <- standardize_datatypes(rtn, groupingvars = groupingvars, batch = batch)
  if(!is.null(groupingvars) & !is.null(batch)){
    rtn <- ms_arrange(rtn, batch, groupingvars)
  }
  else {
    rtn <- ms_arrange(rtn)
  }

  return(rtn)

}




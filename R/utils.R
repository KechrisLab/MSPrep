
# script of internal utility functions



msprep_stages <- function() {
  c("prepared", "filtered", "imputed", "normalized")
}

stage <- function(x) attr(x, "stage")
`stage<-` <- function(x, value) {
  attr(x,  "stage") <- value
  x
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
# TODO: should rt be a numeric -- several levels with colons in the value

#' @importFrom rlang sym
#' @importFrom rlang UQ
#' @importFrom tibble rownames_to_column
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
#' @importFrom tidyr separate
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate_at
#' @importFrom dplyr vars
wide_matrix_to_data <- function(wide) {

  sym_id <- sym("subject_id")
  sym_mz <- sym("mz")
  sym_rt <- sym("rt")

  wide %>%
    as.data.frame %>%
    rownames_to_column(var = "rwnm") %>%
    separate("rwnm", sep = "_", into = c("subject_id", grouping_vars)) %>%
    as_tibble %>%
    gather(key = "mz_rt", value = "abundance_summary", -"subject_id") %>%
    separate("mz_rt", sep = "_", into = c("mz", "rt")) %>%
    arrange(UQ(sym_id), UQ(sym_mz), UQ(sym_rt)) %>%
    mutate_at(vars("subject_id", "rt"), as.factor) %>%
    mutate_at(vars("mz"), as.double) %>%
    arrange(UQ(sym_id), UQ(sym_mz), UQ(sym_rt))

}




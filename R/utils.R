
# script of internal utility functions



# Internal function - replace missing val with NA
replace_missing <- function(abundance, missing_val) {
  ifelse(abundance == missing_val, NA, abundance)
}

# Internal function -- spread data
data_to_wide_matrix <- function(data) {

  data %>% 
    unite(col = sample, subject_id, spike) %>% 
    unite(col = compound, mz, rt) %>% 
    spread(key = compound, value = abundance_summary)  %>%
    as.data.frame %>%
    column_to_rownames(var = "sample") %>%
    as.matrix

}

# Internal function -- undo spread
wide_matrix_to_data <- function(wide) {

  wide %>%
    as.data.frame %>%
    rownames_to_column(var = "id_col") %>%
    tibble::as_data_frame(.) %>%
    tidyr::gather(key = mz_rt, value = abundance_summary, -id_col) %>%
    tidyr::separate(id_col, sep = "_", into = c("subject_id", "spike")) %>%
    tidyr::separate(mz_rt, sep = "_", into = c("mz", "rt")) %>%
    arrange(subject_id, spike, mz, rt) %>%
    mutate_at(vars("subject_id", "spike", "rt"), as.factor) %>%
    mutate_at(vars("mz"), as.double) %>%
    arrange(subject_id, spike, mz, rt) 

}




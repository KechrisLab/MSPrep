
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

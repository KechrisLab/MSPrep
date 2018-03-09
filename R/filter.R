
#' Filters and imputes dataset
#'
#' Filters compounds to those found in specified percentage of subjects and
#' performs data imputation.
#'
#' @param msprepped Summarized dataset output as sum_data1 from readdata() function
#' @param filter_percent Percent to filter the data
#' @return Placeholder
#' @details minval Filtered dataset with missing values replaced by 1/2 minimum
#' observed value for that compound.
#' @details bpca Filtered dataset with missing values imputed by a Bayesian PCA
#' from PCAMethods package.
#' @details withzero Filtered dataset with no imputation.
#' @details count List of all compounds and the percent present for each
#' compound.
#' @references 
#'   Oba, S.et al.(2003) A Bayesian missing value estimation for gene
#'   expression profile data. Bioinformatics, 19, 2088-2096
#'
#'   Stacklies, W.et al.(2007) pcaMethods A bioconductor package providing
#'   PCA methods for incomplete data. Bioinformatics, 23, 1164-1167.
#' @examples
#'
#' library(magrittr)
#'
#' # Load example dataset, tidy it, prepare it, then filter it
#' data(msquant)
#'
#' prepped_data <- msquant %>% ms_tidy %>% ms_prepare
#' filtered_data <- ms_filter(prepped_data, 0.80)
#'
#' @importFrom dplyr select
#' @importFrom dplyr mutate 
#' @importFrom dplyr mutate_at
#' @importFrom dplyr group_by 
#' @importFrom dplyr summarise 
#' @importFrom dplyr ungroup 
#' @importFrom dplyr filter
#' @importFrom dplyr full_join
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
ms_filter <- function (msprepped, filter_percent = 0.5) {

  stopifnot("msprepped" %in% class(msprepped))
  filter_status <- 
    msprepped$summary_data %>%
      group_by(mz, rt) %>% 
      summarise(perc_present = sum(abundance_summary != 0) / n()) %>%
      mutate(keep = perc_present >= filter_percent) %>%
      ungroup

  filtereddata <- 
    full_join(msprepped$summary_data,
              filter_status %>% select(mz, rt, keep),
              by = c("mz", "rt")) %>% 
    filter(keep) %>% select(-keep)

  msfiltered <- 
    msprepped %>%
      append(list("filter_status" = filter_status), after = 4) %>% 
      append(list("filter_percent" = filter_percent))

  msfiltered$summary_data <- filtereddata
  class(msfiltered) <- c("msfiltered", class(msprepped))

  return(msfiltered)

}

print.msfiltered <- function(x) {
  cat("filtered msprep object\n")
  cat("    Count of compounds present in >= ", round(x$filter_percent*100, digits = 3), "% of patients = ", sum(x$filter_status$keep), "\n")
  cat("    Replicate count: ", x$replicate_count, "\n")
  cat("    Patient count: ", length(unique(x$clinical$subject_id)), "\n")
  cat("    Count of spike levels: ", length(unique(x$clinical$spike)), "\n")
  cat("    Count patient-spike combinations: ", nrow(x$clinical), "\n")
  cat("    Count of patient-spike compounds summarized by median: ", nrow(x$medians), "\n")
  cat("    User-defined parameters \n")
  cat("        cvmax = ", x$cvmax, "\n")
  cat("        min_proportion_present = ", round(x$min_proportion_present, digits=3), "\n")
  cat("        filter percent = ", x$filter_percent, "\n")
  cat("    Dataset:\n")
  print(x$summary_data, n = 6)
}

impute_ms <- function(msprepped, impute_function = ~ minval(.x)) {

  stopifnot("msprepped" %in% class(msprepped))

  # replace 0's with NA's for minval and pca() (all methods?)
  data <- msprepped$summary_data %>% 
    mutate_at(vars(abundance_summary), replace_missing, 0)
  # pca imputation
  metabpca <- data_to_wide_matrix(data) %>% pcaMethods::pca(., nPcs = 3, method = "bpca")
  # k-nearest-neighbors imputation
  knnimpute <- data_to_wide_matrix(data) %>% as.data.frame %>% VIM::kNN(., k=5)
  #[1:ncol(data), 1:nrow(data)]

  data %>% impute_function

}

# Filterft to only include compounds that are found in specified percentage of
# subjects and perform imputation of missing data

# Manuscript description: 
# Filtering: The resulting summarized dataset contains all compounds with one
# observation per subject (or sample). The next processing step filters the data
# to only compounds found in a user-specified percentage of subjects.

data_to_wide_matrix <- function(data) {

  data %>% 
    unite(col = sample, subject_id, spike) %>% 
    unite(col = compound, mz, rt) %>% 
    spread(key = compound, value = abundance_summary)  %>%
    as.data.frame %>%
    column_to_rownames(var = "sample") %>%
    as.matrix

}




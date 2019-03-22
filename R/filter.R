
#' Function for filtering prepared dataset.
#'
#' Filters compounds to those found in specified proportion of subjects.
#'
#' @param msprep_obj Prepared MSPrep object.
#' @param filter_percent Decimal value representing to proportion to filter the data.
#' @return An `msprep` object with `stage(rtn) == "filtered"` containing
#' filtered quantification data.
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
#' prepped_data <- msquant %>% ms_tidy %>%
#'   ms_prepare(replicate = "replicate",
#'              batch = "batch",
#'              groupingvars = "spike")
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
#' @importFrom dplyr n
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
ms_filter <- function (msprep_obj, filter_percent = 0.5) {

  stopifnot(class(msprep_obj) == "msprep")
  stopifnot(stage(msprep_obj) == "prepared")

  filter_status <- 
      group_by(msprep_obj$data, .data$mz, .data$rt) %>%
      summarise(perc_present = sum(.data$abundance_summary != 0) / n()) %>%
      mutate(keep = .data$perc_present >= filter_percent) %>%
      ungroup

  filtereddata <- 
    full_join(msprep_obj$data,
              select(filter_status, .data$mz, .data$rt, .data$keep),
              by = c("mz", "rt")) %>% 
    filter(.data$keep) %>% select(-.data$keep)

  msprep_obj$data <- filtereddata
  attr(msprep_obj, "filter_status")  <- filter_status
  attr(msprep_obj, "filter_percent") <- filter_percent
  stage(msprep_obj)                  <- "filtered"

  return(msprep_obj)

}




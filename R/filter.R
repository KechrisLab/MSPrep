#' Function for filtering prepared dataset.
#'
#' Filters compounds to those found in specified proportion of subjects.
#'
#' @param msprep_obj Summarized MSPrep object.
#' @param filter_percent Decimal value indicating filtration threshold. 
#' Metabolites which are present in fewer samples than the specified proportion 
#' will be removed. 
#' @return An `msprep` object with containing filtered quantification data.
#' @references 
#'   Oba, S.et al.(2003) A Bayesian missing value estimation for gene
#'   expression profile data. Bioinformatics, 19, 2088-2096
#'
#'   Stacklies, W.et al.(2007) pcaMethods A bioconductor package providing
#'   PCA methods for incomplete data. Bioinformatics, 23, 1164-1167.
#' @examples
#'
#' # Load example dataset, tidy it, and summarize it
#' data(msquant)
#' 
#' tidied_data <- ms_tidy(msquant, mz = 'mz', rt = 'rt',
#'                        col_extra_txt = 'Neutral_Operator_Dif_Pos_',
#'                        separator = '_', 
#'                        col_names = c('spike', 'batch', 'replicate', 
#'                        'subject_id'))
#' 
#' summarized_data <- ms_summarize(tidied_data, 
#'                                 mz = 'mz', 
#'                                 rt = 'rt', 
#'                                 replicate = 'replicate', 
#'                                 batch = 'batch', 
#'                                 groupingvars = 'spike', 
#'                                 subject_id = 'subject_id', 
#'                                 cvmax = 0.50, 
#'                                 min_proportion_present = 1/3, 
#'                                 missing_val = 1)
#' 
#' # Filter the dataset using a 80% filter rate
#' filtered_data <- ms_filter(summarized_data, 
#'                           filter_percent =  0.80)
#'                           
#' # Print results
#' print(filtered_data)
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
#' @importFrom rlang syms
#' @importFrom rlang !!!
#' @importFrom magrittr %>%
#' @export
ms_filter <- function(msprep_obj, filter_percent = 0.8) {
  
  stopifnot(class(msprep_obj) == "msprep")
  stopifnot(stage(msprep_obj) == "summarized")
  
  met_vars <- met_vars(msprep_obj)
  met_syms <- syms(met_vars)
  
  filter_status <- group_by(msprep_obj$data, `!!!`(met_syms)) %>% 
    summarise(perc_present = sum(.data$abundance_summary != 0)/n()) %>% 
    mutate(keep = .data$perc_present >= filter_percent) %>% 
    ungroup
  
  filtereddata <- full_join(msprep_obj$data, select(filter_status, 
                                                    `!!!`(met_syms), 
                                                    .data$keep), 
                                                    by = met_vars) %>% 
  filter(.data$keep) %>% 
    select(-.data$keep)
  
  msprep_obj$data <- filtereddata
  attr(msprep_obj, "filter_status") <- filter_status
  attr(msprep_obj, "filter_percent") <- filter_percent
  stage(msprep_obj) <- "filtered"
  
  return(msprep_obj)
  
}
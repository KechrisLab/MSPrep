#' Tidy, summarize, filter, impute, and normalize metabolomics dataset
#' 
#' Calls ms_tidy, ms_summarize, ms_filter, ms_impute, and ms_normalize
#' 
#' @param quantification_data Data frame containing the quantification data.
#' @param met_id Name of the column containing compound names.
#' @param mz Name of the column containing mass-to-charge ratios.
#' @param rt Name of the column containing retention time.
#' @param col_extra_txt Text to remove from column names when converting column names to
#'   variables.
#' @param col_names Vector of the ordered ID names to extract from the existing column
#' names.
#' @param separator Character or text separating the variable names sepecified by
#' \code{col_names}
#' 
#' @param subject_id Name of the subject ID column.
#' @param replicate Name of the replicate column. Set to NULL if no
#' replicates.
#' @param abundance Name of the abundance column.
#' @param groupingvars Variable name or vector of names of the
#' phenotypes or comparison groups. Set to NULL if none are present.
#' @param batch Name of the column representing batches. Set to NULL if no batches present.
#' @param cvmax Decimal value from 0 to 1 representing the acceptable level of coefficient 
#' of variation between replicates.
#' @param missing_val Value of missing data in the quantification data file.
#' @param min_proportion_present  Decimal value from 0 to 1 representing the minimum proportion present 
#' to summarize with median or mean. Below this the compound will be set to 0.
#' 
#' @param filter_percent Decimal value representing to proportion to filter the data.
#' 
#' @param imputeMethod Name of imputation method to use. 
#' Options are:
#' - halfmin (half the minimum value)
#' - bpca (Bayesian PCA)
#' - knn (k-nearest neighbors)
#' @param k_knn Number of clusters for 'knn' method.
#' @param n_pcs Number of  principle components used for re-estimation for 
#' 'bpca' method.
#' @param compoundsAsNeighbors If TRUE, will use compounds as neighbors for KNN imputation
#' rather than samples. Note, using compounds as neighbors is significantly slower than using
#' samples as neighbors.
#' 
#' @param normalizeMethod Name of normalization method.
#' - ComBat (only ComBat batch correction)
#' - quantile (only quantile normalization)
#' - quantile + ComBat (quantile with ComBat batch correction)
#' - median (only median normalization)
#' - median + ComBat (median with ComBat batch correction)
#' - CRMN
#' - RUV
#' - SVA
#' @param n_control Number of controls to estimate/utilize.
#' @param controls Vector of control identifiers.  Leave blank for data driven
#' controls. Vector of column numbers from metafin dataset of that control.
#' @param n_comp Number of factors to use in CRMN algorithm.
#' @param k_ruv Number of factors to use in RUV algorithm.
#' @param transform Select transformation to apply to data prior to normalization.
#' - log10
#' - log2
#' - none
#' 
#' @details Wrapper function for ms_tidy, ms_summarize, ms_filter, ms_impute, and
#' ms_normalize.
#' 
#' @examples 
#' # Load example data
#' data(msquant)
#' 
#' # Call function to tidy, summarize, filter, impute, and normalize data
#' prepared_data <- ms_prepare(msquant,
#'                             mz = "mz",
#'                             rt = "rt",
#'                             col_extra_txt = "Neutral_Operator_Dif_Pos_",
#'                             col_names = c("spike", "batch", "replicate", "subject_id"),
#'                             separator = "_",
#'                             abundance = "abundance",
#'                             subject_id = "subject_id",
#'                             replicate = "replicate",
#'                             batch = "batch",
#'                             groupingvars = "spike",
#'                             cvmax = 0.50,
#'                             missing_val = 1,
#'                             min_proportion_present = 1/3,
#'                             filter_percent = .8,
#'                             imputeMethod = "halfmin",
#'                             normalizeMethod = "median")
#' 
#' # Print summary
#' print(prepared_data)
#' 
#' 
#' @export

ms_prepare <- function(quantification_data,
                       met_id = NULL,
                       mz = NULL,
                       rt = NULL,
                       col_extra_txt = NULL,
                       col_names = c("subject_id"),
                       separator = NULL,
                       abundance = "abundance",
                       subject_id = "subject_id",
                       replicate = NULL,
                       batch = NULL,
                       groupingvars = NULL,
                       cvmax = 0.50,
                       missing_val = 1,
                       min_proportion_present = 1/3,
                       filter_percent = .8,
                       imputeMethod = c("halfmin",
                                        "bpca", 
                                        "knn", 
                                        "none"),
                       k_knn = 5, 
                       n_pcs = 3, 
                       compoundsAsNeighbors = FALSE,
                       normalizeMethod = c("ComBat",
                                           "quantile",
                                           "quantile + ComBat",
                                           "median",
                                           "median + ComBat",
                                           "CRMN",
                                           "RUV",
                                           "SVA",
                                           "none"),
                       n_control = 10,
                       controls  = NULL,
                       n_comp    = 2,
                       k_ruv     = 3,
                       transform = c("log10",
                                     "log2",
                                     "none")) {
  
  cat("Tidying\n")
  data <- ms_tidy(quantification_data,
                  met_id = met_id,
                  mz = mz,
                  rt = rt,
                  col_extra_txt = col_extra_txt,
                  col_names = col_names,
                  separator = separator)
  
  cat("Summarizing\n")
  data <- ms_summarize(data,
                       abundance     = "abundance",
                       met_id        = met_id,
                       mz            = mz,
                       rt            = rt,
                       subject_id    = "subject_id",
                       replicate     = replicate,
                       batch         = batch,
                       groupingvars = groupingvars,
                       cvmax         = cvmax,
                       missing_val   = missing_val,
                       min_proportion_present = min_proportion_present)
  
  cat("Filtering\n")  
  data <- ms_filter(data,
                    filter_percent = filter_percent)
  
  if(imputeMethod != "none") {
    cat("Imputing\n")
    data <- ms_impute(data,
                      imputeMethod = imputeMethod,
                      k_knn = k_knn, 
                      n_pcs = n_pcs, 
                      compoundsAsNeighbors = compoundsAsNeighbors)
  }
  
  if(normalizeMethod != "none"){
    cat("Normalizing\n")
    data <- ms_normalize(data,
                         normalizeMethod = normalizeMethod,
                         n_control = n_control,
                         controls  = controls,
                         n_comp    = n_comp,
                         k_ruv     = k_ruv,
                         transform = transform)
  }
  
  cat("Returning")
  data <- ms_return(data)
          
  return(data)
}
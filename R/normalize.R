# 
# #' Perform normalization and batch corrections
# #' 
# #' Perform normalization and batch corrections on specified imputation dataset.
# #' Routines included are quantile, RUV, SVA, median, and CRMN.  Combat to
# #' remove batch effects in raw, quantile, and median normalized data.
# #' Generates data driven controls if none exist.
# #' 
# #' @param metafin Imputated and filtered dataset to use for normalization
# #' @param clindat Name of clinical datafile
# #' @param link1 Column header in clindat that links to subject IDs.
# #' @param pheno Name of phenotype variable in clindat.
# #' @param batch Name of batch variable in clindat.
# #' @param ncont Number of controls to estimate/utilize.
# #' @param controls Vector of control identifiers.  Leave blank for data driven
# #' controls. Vector of column numbers from metafin dataset of that control.
# #' @param ncomp Number of factors to use in CRMN algorithm.
# #' @return controls List of compounds that were used as controls
# #' @return crmn_adj CRMN adjusted dataset (log2)
# #' @return log_data Log2 data with no adjustment
# #' @return log_data_combat Log2 data with Combat batch adjustment
# #' @return log_quant Log2 Quantile adjusted
# #' @return log_quant_combat Log2 Quantile with Combat data
# #' @return med_adj Median adjusted dataset
# #' @return med_combat Median adjusted with Combat
# #' @return ruv_adj RUV adjusted dataset (log2)
# #' @return ruv_factors Matrix of factors generated in RUV for use as covariates
# #' in later analysis.
# #' @return sva_adj SVA adjusted dataset (log2)
# #' @return sva_factors Matrix of factors generated in SVA for use as covariates
# #' in later analysis.
# #' @references 
# #' Bolstad, B.M.et al.(2003) A comparison of normalization methods for high
# #' density oligonucleotide array data based on variance and bias.
# #' Bioinformatics, 19, 185-193
# #' 
# #' DeLivera, A.M.et al.(2012) Normalizing and Integrating Metabolomic Data.
# #' Anal. Chem, 84, 10768-10776.
# #' 
# #' Gagnon-Bartsh, J.A.et al.(2012) Using control genes to correct for unwanted
# #' variation in microarray data. Biostatistics, 13, 539-552.
# #' 
# #' Johnson, W.E.et al.(2007) Adjusting batch effects in microarray expression
# #' data using Empirical Bayes methods. Biostatistics, 8, 118-127.
# #' 
# #' Leek, J.T.et al.(2007) Capturing Heterogeneity in Gene Expression Studies by
# #' Surrogate Variable Analysis. PLoS Genetics, 3(9), e161
# #' 
# #' Wang, W.et al.(2003) Quantification of Proteins and Metabolites by Mass
# #' Spectrometry without Isotopic Labeling or Spiked Standards. Anal. Chem., 75,
# #' 4818-4826.
# #' 
# #' @examples
# #' library(magrittr)
# #'
# #' # Load object generated from readdata() function
# #' data(msquant)
# #'
# #' imputed_data <- 
# #'   msquant %>% 
# #'     ms_tidy %>% ms_prepare %>% ms_filter(filter_percent = 0.80) %>%
# #'     ms_impute(method = "halfmin")
# #' normalized_data <- ms_normalize(imputed_data)
# #'
# #'
# #' @importFrom dplyr case_when
# #' @export
# ms_normalize <- function(msprepped,
#                          method = c("halfmin", "bpca", "bpca + halfmin", "knn")) {
# 
# }

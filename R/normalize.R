
#' Perform normalization and batch corrections
#' 
#' Perform normalization and batch corrections on specified imputation dataset.
#' Routines included are quantile, RUV, SVA, median, and CRMN.  Combat to
#' remove batch effects in raw, quantile, and median normalized data.
#' Generates data driven controls if none exist.
#' 
#' @param metafin Imputated and filtered dataset to use for normalization
#' @param clindat Name of clinical datafile
#' @param link1 Column header in clindat that links to subject IDs.
#' @param pheno Name of phenotype variable in clindat.
#' @param batch Name of batch variable in clindat.
#' @param ncont Number of controls to estimate/utilize.
#' @param controls Vector of control identifiers.  Leave blank for data driven
#' controls. Vector of column numbers from metafin dataset of that control.
#' @param ncomp Number of factors to use in CRMN algorithm.
#' @return controls List of compounds that were used as controls
#' @return crmn_adj CRMN adjusted dataset (log2)
#' @return log_data Log2 data with no adjustment
#' @return log_data_combat Log2 data with Combat batch adjustment
#' @return log_quant Log2 Quantile adjusted
#' @return log_quant_combat Log2 Quantile with Combat data
#' @return med_adj Median adjusted dataset
#' @return med_combat Median adjusted with Combat
#' @return ruv_adj RUV adjusted dataset (log2)
#' @return ruv_factors Matrix of factors generated in RUV for use as covariates
#' in later analysis.
#' @return sva_adj SVA adjusted dataset (log2)
#' @return sva_factors Matrix of factors generated in SVA for use as covariates
#' in later analysis.
#' @references 
#' Bolstad, B.M.et al.(2003) A comparison of normalization methods for high
#' density oligonucleotide array data based on variance and bias.
#' Bioinformatics, 19, 185-193
#' 
#' DeLivera, A.M.et al.(2012) Normalizing and Integrating Metabolomic Data.
#' Anal. Chem, 84, 10768-10776.
#' 
#' Gagnon-Bartsh, J.A.et al.(2012) Using control genes to correct for unwanted
#' variation in microarray data. Biostatistics, 13, 539-552.
#' 
#' Johnson, W.E.et al.(2007) Adjusting batch effects in microarray expression
#' data using Empirical Bayes methods. Biostatistics, 8, 118-127.
#' 
#' Leek, J.T.et al.(2007) Capturing Heterogeneity in Gene Expression Studies by
#' Surrogate Variable Analysis. PLoS Genetics, 3(9), e161
#' 
#' Wang, W.et al.(2003) Quantification of Proteins and Metabolites by Mass
#' Spectrometry without Isotopic Labeling or Spiked Standards. Anal. Chem., 75,
#' 4818-4826.
#' 
#' @examples
#' library(magrittr)
#'
#' # Load object generated from readdata() function
#' data(msquant)
#'
#' imputed_data <- 
#'   msquant %>% 
#'     ms_tidy %>% ms_prepare %>% ms_filter(filter_percent = 0.80) %>%
#'     ms_impute(method = "halfmin")
#' normalized_data <- ms_normalize(imputed_data)
#'
#'
#' @importFrom dplyr case_when
#' @export
ms_normalize <- function(msprepped,
                         method = c("ComBat",
                                    "quantile + ComBat",
                                    "median + ComBat",
                                    "CRMN",
                                    "RUV",
                                    "SVA"),
                         n_control = NULL,
                         controls  = NULL,
                         n_comp    = NULL) {

  # Validate inputs
  stopifnot(class(msprep_obj) == "msprep")
  # NOTE: check whether all require imputation or if other valid stages
  stopifnot(stage(msprep_obj) %in% c("imputed"))
  match.arg(method) # requires 1 argument from vec in function arg

  # # Prep data - replace 0's with NA's -- for minval and bpca() (all methods?)
  # data <- mutate_at(msprep_obj$data, vars("abundance_summary"), 
  #                   replace_missing, 0)

  # normalize/batch correct data
  data <-
    switch(method,
           "ComBat"            = normalize_combat(data),
           "quantile + ComBat" = normalize_quantile_combat(data),
           "median + ComBat"   = normalize_median_combat(data),
           "CRMN"              = normalize_crmn(data),
           "RUV"               = normalize_ruv(data),
           "SVA"               = normalize_sva(data),
           stop("Invalid normalize method - provide argument from list in",
                "function definition and help file"))

  # Prep output object
  msprep_obj$data  <- data
  attr(msprep_obj, "normalize_method") <- method
  stage(msprep_obj) <- "normalized"

  # ...and:
  return(msprep_obj)


}


#' @importFrom preprocessCore normalize.quantiles
normalize_combat <- function(data) {

  # Combat batch correction
  rtn <- combat(rtn)

  return(rtn)

}

#' @importFrom preprocessCore normalize.quantiles
normalize_quantile_combat <- function(data) {

  # Quantile normalization
  mat   <- data_to_wide_matrix(data)
  rwnm  <- rownames(mat)
  colnm <- colnames(mat)
  rtn   <- normalize.quantiles(mat)
  rownames(rtn) <- rwnm
  colnames(rtn) <- colnm

  # Combat batch correction
  rtn <- combat(rtn)

  return(rtn)

}


normalize_median_combat <- function(data) {

}




normalize_crmn <- function(data) {
}

# wait on this -- need to rewrite genadj and ruv functions
normalize_ruv <- function(data) {

  gen_adj()
}

# wait on this -- need to rewrite genadj function and pick apart related code in
# main  function
normalize_sva <- function(data) {
}




# Utility functions for ms_normalize()
log_base2 <- function(wide_matrix) apply(wide_matrix, 2, log2)

combat <- function(data, phenotype, batch) {


}




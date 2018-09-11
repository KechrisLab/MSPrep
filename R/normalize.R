#' Perform normalization and batch corrections
#' 
#' Perform normalization and batch corrections on specified imputation dataset.
#' Routines included are quantile, RUV, SVA, median, and CRMN.  Combat to
#' remove batch effects in raw, quantile, and median normalized data.
#' Generates data driven controls if none exist.
#' 
#' @param msprep_obj An MSPrep object that has been filtered, normalized, and
#' imputed.
#' @param method One of 6 methods.  3 are normalization and 3 are batch
#' correction + normalization? TODO: updated this.
#' @param n_control Number of controls to estimate/utilize.
#' @param controls Vector of control identifiers.  Leave blank for data driven
#' controls. Vector of column numbers from metafin dataset of that control.
#' @param n_comp Number of factors to use in CRMN algorithm.
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
#'   ms_tidy %>%
#'   ms_prepare %>%
#'   ms_filter(filter_percent = 0.80) %>%
#'   ms_impute(method = "halfmin")
#' normalized_data <- ms_normalize(imputed_data)
#'
#'
#' @importFrom dplyr case_when
#' @export
ms_normalize <- function(msprep_obj,
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
  method <- match.arg(method) # requires 1 argument from vec in function arg

  # # Prep data - replace 0's with NA's -- for minval and bpca() (all methods?)
  # data <- mutate_at(msprep_obj$data, vars("abundance_summary"), 
  #                   replace_missing, 0)

  # normalize/batch correct data
  data <-
    switch(method,
           "ComBat"            = normalize_combat(data),
           "quantile + ComBat" = normalize_quantile_combat(data, grouping_vars(msprep_obj)),
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
normalize_quantile_combat <- function(data, groupingvars) {

  # Quantile normalization
  mat   <- data_to_wide_matrix(data, groupingvars)
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






# # Utility functions for ms_normalize()
log_base2 <- function(wide_matrix) apply(wide_matrix, 2, log2)




# NOTE: returns log base 2 results
combat_new <- function(data, groupingvars, batch) {

  # testing setup
  data         <- dat$data
  groupingvars <- MSPrep:::grouping_vars(dat)
  batch        <- MSPrep:::batch_var(dat)

  # Function dev
  # Convert to wide matrix -- do we want spike in the ID?
  wide <- data_to_wide_matrix(data, groupingvars, batch, asmatrix = FALSE)
  wide <- wide %>% rownames_to_column
  wide <- separate(wide, col = "rowname",
                   into = internal_id_order(groupingvars, batch),
                   remove = FALSE)

  # Generate long matrix  for feeding to ComBat()
  widemat <- data_to_wide_matrix(data, groupingvars, batch, asmatrix = TRUE)
  longmat <- t(widemat)
  loglongmat <- log_base2(longmat)

  # Create model matrix of sample_info variables (soon to be formerly groupingvars)
  sample_info_dat <- select(wide, groupingvars)
  sample_info_dat <- mutate_if(sample_info_dat, is.character, funs(as.factor))
  sample_info_modmat <- model.matrix(~ ., data = sample_info_dat)

  # Create model matrix of batch variable
  if (length(batch) > 1) stop("only one batch variable allowed for ComBat") 
  batch_dat <- select(wide, batch)
  batch_dat <- mutate_if(batch_dat, is.character, funs(as.factor))
  batch_vec <- batch_dat[, 1]

  comadj <- sva::ComBat(loglongmat, mod = sample_info_modmat, batch = batch_vec) 
  comadj <- t(comadj)

  # What form to return this in??
  return(comadj)

}



# combat <- function(data, phenotype, batch) {
# 
#   # data == wide matrix data w/ rownames
#   # batch == "Operator" -- i.e. 
#   # pheno == "Spike"
#   # clindat == Clinical.csv
#   # link1 == "SubjectID"
# 
#   path_olddata <- system.file("extdata", "old_object.Rda", package = "MSPrep")
#   load(path_olddata)
#   clinical     <- read_csv("data-raw/Clinical.csv")
# 
# 
#   phenotype <- "Spike"
#   batch     <- "Operator"
#   link1     <- "SubjectID"
#   olddata   <- old_readdata_result$sum_data %>% as.data.frame
#   clindat   <- clinical
# 
# 
#   #   compounds <- ncol(data)
#   temp      <- olddata
#   temp$mid  <- rownames(temp)
#   # wide matrix format + subject id (or spike_subjectid)
#   test      <- merge(temp, clindat, by.x = "mid", by.y = as.character(link1))
#   d1        <- t(olddata)
# #   colnames(d1) <- olddata[, 1]
#   d2    <- subset(test, select = phenotype) # spike only
#   d2a   <- model.matrix(~as.factor(d2[, 1])) # spike in model matrix form
#   d3    <- subset(test, select = batch) # operator ids
# #   count <- matrix(0, nrow = compounds, ncol = 1)
# #   d1mod <- matrix(0, nrow = compounds, ncol = subjects)
# 
#   d1mod <- log_base2(d1)
# 
#   # Here the element "d2a" is Gold Stage phenotype and "d3" is the batch.
#   comadj <- t(sva::ComBat(d1mod, mod = d2a, batch = d3[,1])) 
#   colnames(comadj) <- colnames(olddata)
#   rownames(comadj) <- rownames(olddata)
# 
#   return(comadj)
# 
# }
# 
# 
# 

# ################################################################################
# # testing setup
# ################################################################################
# # New testing
# dat <- imputed_data
# data         <- dat$data
# groupingvars <- MSPrep:::grouping_vars(dat)
# batch        <- MSPrep:::batch_var(dat)
################################################################################
# Old testing
#   # data == wide matrix data w/ rownames
#   # batch == "Operator" -- i.e. 
#   # pheno == "Spike"
#   # clindat == Clinical.csv
#   # link1 == "SubjectID"
# path_olddata <- system.file("extdata", "old_object.Rda", package = "MSPrep")
# load(path_olddata)
# clinical <- read_csv("data-raw/Clinical.csv")
# pheno    <- "Spike"
# batch    <- "Operator"
# link1    <- "SubjectID"
# olddata  <- as.data.frame(old_readdata_result$sum_data)
# clindat  <- clinical
################################################################################



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
#'   ms_prepare(replicate = "replicate", batch = "batch", groupingvars = "spike") %>%
#'   ms_filter(filter_percent = 0.80) %>%
#'   ms_impute(method = "halfmin")
#' normalized_data_qC  <- ms_normalize(imputed_data, method = "quantile + ComBat")
#' normalized_data_C   <- ms_normalize(imputed_data, method = "ComBat")
#' normalized_data_sva <- ms_normalize(imputed_data, method = "SVA")
#' normalized_data_ruv <- ms_normalize(imputed_data, method = "RUV")
#' normalized_data_crmn <- ms_normalize(imputed_data, method = "CRMN")
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
                         n_control = 10,
                         controls  = NULL,
                         n_comp    = 2) {

  # Validate inputs
  stopifnot(class(msprep_obj) == "msprep")
  # NOTE: check whether all require imputation or if other valid stages
  stopifnot(stage(msprep_obj) %in% c("imputed"))
  method <- match.arg(method) # requires 1 argument from vec in function arg

  # Grab attributes for passing to fns
  groupingvars <- grouping_vars(msprep_obj)
  batch        <- batch_var(msprep_obj)
  data         <- msprep_obj$data

  # normalize/batch correct data
  data <-
    switch(method,
           "ComBat"            = normalize_combat(data, groupingvars, batch),
           "quantile + ComBat" = normalize_quantile_combat(data, groupingvars, batch),
           #"median + ComBat"   = normalize_median_combat(data),
           "CRMN"              = normalize_crmn(data, groupingvars, batch, n_comp, n_control),
           "RUV"               = normalize_ruv(data, groupingvars, batch, n_control),
           "SVA"               = normalize_sva(data, groupingvars, batch),
           stop("Invalid normalize method - provide argument from list in",
                "function definition and help file"))

  # Prep output object
  msprep_obj$data  <- data
  attr(msprep_obj, "normalize_method") <- method
  stage(msprep_obj) <- "normalized"

  # ...and:
  return(msprep_obj)


}

normalize_combat <- function(data, groupingvars, batch) {

  # Combat batch correction
  rtn <- combat(data, groupingvars, batch, log2_transform = TRUE)

  return(rtn)

}


#' @importFrom preprocessCore normalize.quantiles
normalize_quantile_combat <- function(data, groupingvars, batch) {

  # Quantile normalization
  mat   <- data_to_wide_matrix(data, groupingvars, batch)
  rwnm  <- rownames(mat)
  colnm <- colnames(mat)
  rtn   <- t(normalize.quantiles(t(mat)))
  rownames(rtn) <- rwnm
  colnames(rtn) <- colnm

  # This is backtransforming to 'tidy' dataset before combat, which then
  # transforms back to wide matrix and then back to tidy.  This could be skipped
  # if have time to generalize/separate out a 'prep for combat' fuinction or add
  # checks for current format.
  rtn <- wide_matrix_to_data(rtn, groupingvars, batch)

  # Combat batch correction
  rtn <- combat(rtn, groupingvars, batch)

  return(rtn)

}


# TODO: isIS is a mess -- seems to be dependent on a whole bunch of stuff
# including SVA -- is there a simpler way?
# #
# normalize_median_combat <- function(data, groupingvars, batch) {
# 
#   # User inputs: 
#   # - n_control: numbe4r of controls to estimate/utilize?????
#   n_control <- 10
# 
#   compounds       <- ncol(metafin)
#   subjects        <- nrow(metafin)
#   final           <- as.data.frame(metafin)
#   dset    <- final
#   control <- matrix(NA, nrow = compounds, ncol = 6)
# 
#   # NOT NECESSARY???????
#   sva <- as.data.frame(svafactors(log_final, pheno, batch, compounds))
#   for (i in 1:ncol(sva)) {
#     colnames(sva)[i] <- paste("f", i, sep = "")
#   }
#   colnames(sva) <- paste0("f", 1:ncol(sva))
# 
#   # Old versions
# #   for (j in 1:compounds) {
# # 
# #     control[j, 1] <- colnames(dset)[j]                                           # compund mz_rt
# #     control[j, 2] <- j                                                           # compound number
# #     control[j, 3] <- mean(dset[, j], na.rm = TRUE)                               # mean
# #     control[j, 4] <- sd(dset[, j], na.rm = TRUE)                                 # SD
# #     control[j, 5] <- sd(dset[, j], na.rm = TRUE) / mean(dset[, j], na.rm = TRUE) # coef of variation
# #     control[j, 6] <- sum(dset[, j] == 0, na.rm = TRUE)                           # sum == 0
# # 
# #   }
# #   control2 <- control[order(control[, 6], control[, 5]), ] # sort by col 6(sum == 0) then 5(sd/mean - coeffiecient of variation)
# #   ctl      <- as.numeric(control2[, 2][1:max(n_control, ncol(sva) + 5)])
#   # end
# 
#   # New versions
#   control  <- control_summary(data)
#   control2 <- arrange(control, counteq0, cv)
#   ctl      <- control2[1:max(n_control, ncol(sva) + 5), "number"]
#   # end
# 
#   isIS     <- rep(FALSE, compounds)
#   ctlo     <- ctl[order(ctl)]
#   j        <- 1
# 
#   for (i in 1:compounds) {
#     if (j <= 10) {
#       if (ctlo[j] == i) {
#         isIS[i] <- TRUE
#         j       <- j + 1
#       } else {
#         isIS[i] <- FALSE
#       }
#     }
#   }
#   dset     <- final_shift
#   dset$mid <- rownames(dset)
#   test     <- merge(dset, clindat, by.x = "mid", by.y = link1)
#   dset     <- final_shift
#   dset$mid <- rownames(dset)
#   test     <- merge(dset, clindat, by.x = "mid", by.y = link1)
#   X        <- subset(test, select = c(pheno))
#   colnames(X)[1] <- "pheno"
#   G <- model.matrix(~-1 + as.factor(X$pheno)) # Here, "G" comes from the Gold Stage phenotype.
#   Y           <- t(final_shift)
# 
#   # For normalize/median, we need
#   #   - Y     <- transposed (converts to matrix) final_shift <- final <- dataframed <- input dataset in wide_matrix format
#   #   - G     <- cell means model matrix of groupingvars
#   #   - isIS
#   # TODO: Now figure out what these are and simplist way to get them.
#   normed.med   <- crmn::normalize(Y, "median", factors = G, standards = isIS)
#   med_com      <- combat(as.data.frame(t(normed.med)), pheno, batch)
# }
# 


# Not entirely sure what this function is/does
isIS <- function(data, factors, n_control, n_compounds) {
  # compounds -> number of mz_rt
  # ctlo -> ctl[order(ctl)]

  control   <- control_summary(data)
  ctl       <- ctl_compounds(control, factors, n_control)
  ctlo      <- ctl[order(ctl)]
  j         <- 1
  isISvec      <- rep(FALSE, n_compounds)

  for (i in 1:n_compounds) {
    if (j <= 10) {
      cat("sanity check:\n", j, "\n", i, "\n", ctlo[j], "\n")
      if (ctlo[j] == i) {
        isISvec[i] <- TRUE
        j       <- j + 1
      } else {
        isISvec[i] <- FALSE
      }
    }
  }

  return(isISvec)

}



ctl_compounds <- function(control_summary_data, factors, n_control) {

  dat <- arrange(control_summary_data, .data$counteq0, .data$cv)
  dat <- dat[1:max(n_control, ncol(factors) + 5), "number"]
  dat <- dat[["number"]]

  return(dat)

}



normalize_crmn <- function(data, groupingvars, batch, n_comp, n_control) {

  internal_id <- internal_id_order(groupingvars, batch)

  sva         <- svafactors(data, groupingvars, batch)
  widedf      <- data_to_wide_matrix(data, groupingvars, batch, asmatrix = FALSE)

  n_compounds <- ncol(widedf)
  Y    <- t(widedf) # For crmn::normalize()

  wide_plus_rownames <- rownames_to_column(widedf) # for model.matrix
  wide_plus_rownames <- separate(wide_plus_rownames, col = "rowname", into = internal_id, remove = FALSE)
  sample_info_dat    <- select(wide_plus_rownames, groupingvars)
  sample_info_dat    <- mutate_if(sample_info_dat, is.character, funs(as.factor))
  sample_info_modmat <- model.matrix(~ -1 + ., data = sample_info_dat)

  isIS_vec <- isIS(data, sva, n_control, n_compounds)

  # Requires:
  # - Y -> long, not logged, matrix w/ record id as column names and mz_rt as rownames
  # - G -> sample_info_modmat -> cell means model matrix of phenotypes/sample info
  # - isIS -> result of isIS() function
  # - n_comp -> user input to main ms_normalize() function

   normed.crmn  <- crmn::normalize(Y, "crmn", factors = sample_info_modmat, standards = isIS_vec, ncomp = n_comp)
   lnormed.crmn <- log_base2(normed.crmn)
   final_crmn   <- as.data.frame(t(lnormed.crmn))
   final_crmn   <- wide_matrix_to_data(final_crmn, groupingvars, batch)

   return(final_crmn)

}


# done
normalize_ruv <- function(data, groupingvars, batch, n_control) {

  ruv_factors <- ruvfactors(data, groupingvars, batch, n_control)
  adjusted    <- genadj(data, groupingvars, batch, ruv_factors)

  return(adjusted)

}


# done
normalize_sva <- function(data, groupingvars, batch) {

  sva_factors <- svafactors(data, groupingvars, batch)
  adjusted    <- genadj(data, groupingvars, batch, sva_factors)

  return(adjusted)

}


# NOTE: returns log base 2 results

#' @importFrom dplyr mutate_if
#' @importFrom stats model.matrix
combat <- function(data, groupingvars, batch, log2_transform = TRUE) {

  internal_id <- internal_id_order(groupingvars, batch)

  # Function dev
  # Convert to wide matrix -- do we want spike in the ID?
  wide <- data_to_wide_matrix(data, groupingvars, batch, asmatrix = FALSE)
  wide <- rownames_to_column(wide)
  wide <- separate(wide, col = "rowname",
                   into = internal_id,
                   remove = FALSE)

  # Generate long matrix  for feeding to ComBat()
  widemat <- data_to_wide_matrix(data, groupingvars, batch, asmatrix = TRUE)
  longmat <- t(widemat)

  # combat + median is not log2'd, otherwise should be
  if (log2_transform) {
    longmat <- log_base2(longmat)
  }

  # Create model matrix of sample_info variables (soon to be formerly groupingvars)
  sample_info_dat    <- select(wide, groupingvars)
  sample_info_dat    <- mutate_if(sample_info_dat, is.character, funs(as.factor))
  sample_info_modmat <- model.matrix(~ ., data = sample_info_dat)

  # Create model matrix of batch variable
  if (length(batch) > 1) stop("only one batch variable allowed for ComBat") 
  batch_dat <- select(wide, batch)
  batch_dat <- mutate_if(batch_dat, is.character, funs(as.factor))
  batch_vec <- batch_dat[, 1]

  comadj <- sva::ComBat(longmat, mod = sample_info_modmat, batch = batch_vec) 
  comadj <- t(comadj)

  comadj_tidy <- wide_matrix_to_data(comadj, groupingvars, batch)

  # What form to return this in??
  return(comadj_tidy)

}


#' @importFrom dplyr mutate_if
#' @importFrom dplyr select
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr separate
#' @importFrom stats model.matrix
#' @importFrom sva sva
svafactors <- function(data, groupingvars, batch) {

  # This can probably be simplified some more

  # Create long, logged matrix for sva
  wide      <- data_to_wide_matrix(data, groupingvars, batch, asmatrix = FALSE)
  wide_log2 <- log_base2(wide)
  long_log2 <- t(wide_log2)

  # Create model matrix of sample_info variables (soon to be formerly groupingvars)
  internal_id        <- internal_id_order(groupingvars, batch)
  wide               <- rownames_to_column(wide)
  wide               <- separate(wide, col = "rowname", into = internal_id, remove = FALSE)
  sample_info_dat    <- select(wide, groupingvars)
  sample_info_dat    <- mutate_if(sample_info_dat, is.character, funs(as.factor))
  sample_info_modmat <- model.matrix(~ ., data = sample_info_dat)

  svaadj <- sva(long_log2, sample_info_modmat, method = "irw")
  rtn    <- svaadj$sv
  colnames(rtn) <- paste0("f", 1:ncol(rtn))

  return(rtn)

}


ruvfactors <- function(data, groupingvars, batch, n_control) {

  # Prep datasets/matrices
  wide        <- data_to_wide_matrix(data, groupingvars, batch)
  dset        <- log_base2(wide)
  Y           <- t(dset) # Only data Used in calculation

  # Get num factors from sva -- why sva? who knows
  sva_factors <- svafactors(data, groupingvars, batch)

  # ctl -- n compounds w/ lowest CV and no missing
  control  <- control_summary(data)
  ctl      <- ctl_compounds(control, sva_factors, n_control)

  # Calculate RUV 
  #   Uses: Y, ctl, sva_factors
  Z   <- matrix(rep(1, ncol(Y)))
  RZY <- Y - Y %*% Z %*% solve(t(Z) %*% Z) %*% t(Z)
  W   <- svd(RZY[ctl, ])$v
  W   <- W[, 1:ncol(sva_factors)]

  # Format output
  rtn <- as.matrix(W, ncol = ncol(sva))
  colnames(rtn) <- paste0("f", 1:ncol(rtn))

  return(rtn)

}

#' @importFrom stats lm
genadj <- function(data, groupingvars, batch, factors) {
  # dset      == log_final
  # factors   == sva or ruv_raw factors/loadings?
  # compounds == num mz_rt

  wide      <- data_to_wide_matrix(data, groupingvars, batch, asmatrix = FALSE)
  wide_log2 <- log_base2(wide)

  out <- sapply(1:ncol(wide_log2),
    function(j) {
      fit <- lm(wide_log2[, j] ~ as.matrix(factors, ncol = 1))
      fit$fitted.values
    })

  colnames(out) <- colnames(wide_log2)
  rownames(out) <- rownames(wide_log2)

  rtn <- wide_matrix_to_data(out, groupingvars, batch)

  return(rtn)

}


#' Summarise dataset for use in control compounds (ctl_compounds) function
#' 
#' @param data tidy dataset from msprep_obj (msprep_obj$data)
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise_at
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom rlang .data
control_summary <- function(data) {

  counteq0 <- function(x) sum(x == 0)
  rtn <- group_by(data, .data$mz, .data$rt)
  rtn <- summarise_at(rtn, vars(.data$abundance_summary),
                      funs(mean, sd, counteq0)) #, na.rm = TRUE)
  rtn <- ungroup(rtn)
  rtn <- mutate(rtn, cv = sd / mean, number = 1:n())
  rtn <- select(rtn, .data$mz, .data$rt, .data$number, .data$mean, .data$sd,
                .data$cv, .data$counteq0)
  return(rtn)

}


# Utility functions for ms_normalize()
log_base2 <- function(wide_matrix) apply(wide_matrix, 2, log2)




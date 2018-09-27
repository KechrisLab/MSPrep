
#' Perform normalization and batch corrections
#'
#' Perform normalization and batch corrections on specified imputation dataset.
#' Routines included are quantile, RUV, SVA, median, and CRMN.  Combat to
#' remove batch effects in raw, quantile, and median normalized data.
#' Generates data driven controls if none exist.
#'
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3866554/#btt589-B3
#'
#' 1. Normalization - 5 methods
#'    0. Unnormalized ->  ComBat
#'    *1. Quantile    ->  ComBat
#'    *2. median      ->  ComBat
#'    *3. CRMN 
#'    +4. RUV 
#'    +5. SVA 
#' 
#'  * Two choices Normalization & ComBat -- but only for first 3 (or 4 incl.
#'  unnormalized)
#'
#'  * normalized datasets
#'  + supervised factor analysis
#'
#'
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
#'
#'   path_olddata <- system.file("extdata", "old_object.Rda", package = "MSPrep")
#'   load(path_olddata)
#'   source("./old-R-scripts/old_filterft.R")
#'   source("./old-R-scripts/old_normdata.R")
#'   old_filtered_result <- filterft(metaf = old_readdata_result$sum_data1)
#'
#'   metafin  <- old_filtered_result$bpca
#'   clindat  <- old_readdata_result$clinical
#'   #metafin <- test2$bpca
#'   #clindat  <- test$clinical
#'   link1    <- "SubjectID"
#'   pheno    <- "Spike"
#'   batch    <- "Operator"
#'   ncont    <- 10
#'   controls <- c()
#'   ncomp    <- 2
#'   
#'   test3 <- normdata(metafin, clindat, link1, pheno, batch, ncont = 10,
#'                     controls, ncomp)
#'   str(test3, max.level = 1)
#'
#' @export 

# Katerina's function, unedited (formerly normdata.new_kk)
normdata <- function (metafin, # imputed dataset
                      clindat,
                      link1,
                      pheno,
                      batch,
                      ncont = 10,
                      controls = c(),
                      ncomp = 2) {

  compounds       <- ncol(metafin)
  subjects        <- nrow(metafin)

  # qunatile normalization
  metan           <- t(preprocessCore::normalize.quantiles(t(metafin)))
  colnames(metan) <- colnames(metafin)
  rownames(metan) <- rownames(metafin)
  finaln          <- metan

  final           <- as.data.frame(metafin)


  # Convert to log2 
  log_final  <- convert(final) # for combat only, sva
  log_finaln <- convert(finaln) # only for quantile + combat

  # Combat procedure for combat-only and quantile + combat 
  final_rc  <- combat(as.data.frame(final), pheno, batch)
  final_com <- combat(as.data.frame(finaln), pheno, batch)


  sva <- as.data.frame(svaestimate(log_final, pheno, batch, compounds))
  for (i in 1:ncol(sva)) {
    colnames(sva)[i] <- paste("f", i, sep = "")
  }

  dset    <- final
  control <- matrix(NA, nrow = compounds, ncol = 6)

  for (j in 1:compounds) {

    control[j, 1] <- colnames(dset)[j]
    control[j, 2] <- j
    control[j, 3] <- mean(dset[, j], na.rm = TRUE)
    control[j, 4] <- sd(dset[, j], na.rm = TRUE)
    control[j, 5] <-
      sd(dset[, j], na.rm = TRUE) / mean(dset[, j], na.rm = TRUE)
    control[j, 6] <- sum(dset[, j] == 0, na.rm = TRUE)

  }

  control2 <- control[order(control[, 6], control[, 5]), ]
  # Grab top n with no missing (=0) and smallest Coef of Variation
  ctl <- as.numeric(control2[, 2][1:max(ncont, ncol(sva) + 5)])

  if (length(controls) > 0) { ctl <- controls }

  contout <- cbind(ctl, colnames(log_final)[ctl])
  ruv_raw <- as.data.frame(ruv(log_final, ncol(sva), ctl, pheno, compounds))

  for (i in 1:ncol(ruv_raw)) {
    colnames(ruv_raw)[i] <- paste("f", i, sep = "")
  }

  ruv_adj     <- genadj(log_final, ruv_raw, compounds)
  sva_adj     <- genadj(log_final, sva, compounds)
  final_shift <- final
  Y           <- t(final_shift)
  Ya          <- t(final_shift)[ctl, ]
  Ys          <- t(final_shift)[-ctl, ]
  j           <- 1
  isIS        <- rep(FALSE, compounds)
  ctlo        <- ctl[order(ctl)]

  for (i in 1:compounds) {
    if (j <= 10) {
      if (ctlo[j] == i) {
        isIS[i] <- TRUE
        j       <- j + 1
      } else {
        isIS[i] <- FALSE
      }
    }
  }

  dset     <- final_shift
  dset$mid <- rownames(dset)
  test     <- merge(dset, clindat, by.x = "mid", by.y = link1)
  X        <- subset(test, select = c(pheno))
  colnames(X)[1] <- "pheno"
  G <- model.matrix(~-1 + as.factor(X$pheno)) # Here, "G" comes from the Gold Stage phenotype.

  normed.crmn  <- crmn::normalize(Y, "crmn", factors = G, standards = isIS, ncomp = ncomp)
  lnormed.crmn <- convert(normed.crmn)
  normed.med   <- crmn::normalize(Y, "median", factors = G, standards = isIS)
  lnormed.med  <- convert(normed.med)
  final_crmn   <- as.data.frame(t(lnormed.crmn))
  final_med    <- as.data.frame(t(lnormed.med))
  med_com      <- combat(as.data.frame(t(normed.med)), pheno, batch)

  list(#  unnormalized/batch corrected
       log_data         = log_final, # log2'd
       log_data_combat  = final_rc, # log2'd
       # quantile normalized/batch corrected
       log_quant        = log_finaln, # log2'd
       log_quant_combat = final_com, # log2'd
       # median normalized/batch corrected 
       med_adj          = final_med, # log2'd
       med_combat       = med_com, # NOT log2'd
       # SVA factor analysis 
       sva_factors      = sva, # log2'd
       sva_adj          = sva_adj, # log2'd
       # RUV factor analysis
       ruv_factors      = ruv_raw, # log2'd
       ruv_adj          = ruv_adj, # log2'd
       # CRMN normalized
       crmn_adj         = final_crmn, # log2'd
       controls         = contout) # log2'd

}







# what does this do? --- all it does is log2 all of the values............
convert <- function(dset) {

  s_dset <- as.data.frame(matrix(NA, ncol = ncol(dset), 
                                 nrow = nrow(dset)))
  for (j in 1:ncol(dset)) {
    s_dset[, j] <- log2(as.numeric(dset[, j]))
  }
  colnames(s_dset) <- colnames(dset)
  rownames(s_dset) <- rownames(dset)

  return(s_dset)

}

combat <- function(dset, pheno, batch) {

#   compounds <- ncol(dset)
  temp      <- dset
  temp$mid  <- rownames(temp)
  # wide matrix format + subject id (or spike_subjectid)
  test      <- merge(temp, clindat, by.x = "mid", by.y = as.character(link1))
  d1        <- t(dset)
  colnames(d1) <- dset[, 1]
  d2    <- subset(test, select = pheno) # spike only
  d2a   <- model.matrix(~as.factor(d2[, 1])) # spike in model matrix form
  d3    <- subset(test, select = batch) # operator ids
#   count <- matrix(0, nrow = compounds, ncol = 1)
#   d1mod <- matrix(0, nrow = compounds, ncol = subjects)

  d1mod <- log_base2(d1)

  # Here the element "d2a" is Gold Stage phenotype and "d3" is the batch.
  comadj <- t(sva::ComBat(d1mod, mod = d2a, batch = d3[,1])) 
  colnames(comadj) <- colnames(dset)
  rownames(comadj) <- rownames(dset)

  return(comadj)

}

svaestimate <- function(dset, pheno, batch, compounds) {

  d1           <- t(dset[, 1:compounds])
  colnames(d1) <- rownames(dset)
  rownames(d1) <- colnames(dset)
  dset$mid     <- rownames(dset)
  test <- merge(dset, clindat, by.x = "mid", by.y = as.character(link1))
  d2   <- subset(test, select = c(batch))
  d3   <- subset(test, select = c(pheno))
  colnames(d3)[1] <- "pheno"
  d3m    <- model.matrix(~as.factor(pheno), data = d3)
  count  <- matrix(0, nrow = compounds, ncol = 1)
  d1mod  <- d1
  svaadj <- sva::sva(d1mod, d3m, method = "irw")

  return(as.matrix(svaadj$sv, ncol = svaadj$n.sv))

}

ruv <- function(dset, k, ctl, pheno, compounds) {

  d1       <- t(dset[, 1:compounds])
  count    <- matrix(0, nrow = compounds, ncol = 1)
  Y        <- d1
  dset$mid <- rownames(dset)
  test     <- merge(dset, clindat, by.x = "mid", by.y = link1)
  X        <- subset(test, select = pheno)
  Z        <- matrix(rep(1, ncol(Y)))
  RZY      <- Y - Y %*% Z %*% solve(t(Z) %*% Z) %*% t(Z)
  W        <- svd(RZY[ctl, ])$v
  W        <- W[, 1:k]

  return(as.matrix(W, ncol = k))

}

genadj <- function(dset, factors, compounds) {

  out <- sapply(1:compounds,
    function(j) {
      fit <- lm(dset[, j] ~ as.matrix(factors, ncol = 1))
      fit$fitted.values
    })

  colnames(out) <- colnames(dset)
  rownames(out) <- rownames(dset)

  return(out)

}


log_base2 <- function(wide_matrix) apply(wide_matrix, 2, log2)

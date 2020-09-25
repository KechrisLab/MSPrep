# ms_normalize <- function(msprep_obj,
#                          normalizeMethod = c("median",
#                                              "ComBat",
#                                              "quantile",
#                                              "quantile + ComBat",
#                                              "median + ComBat",
#                                              "CRMN",
#                                              "RUV",
#                                              "SVA"),
#                          n_control = 10,
#                          controls  = NULL,
#                          n_comp    = 2,
#                          k_ruv     = 3,
#                          transform = c("log10",
#                                        "log2",
#                                        "none")) {
# 
#   # Validate inputs
#   stopifnot(class(msprep_obj) == "msprep")
#   # NOTE: check whether all require imputation or if other valid stages
#   #         consider checking for 0's and disallow transformation if present
#   stopifnot(stage(msprep_obj) %in% c("imputed"))
#   normalizeMethod <- match.arg(normalizeMethod) # requires 1 argument from vec in function arg
#   transform <- match.arg(transform)
# 
#   # Grab attributes for passing to fns
#   groupingvars <- grouping_vars(msprep_obj)
#   batch        <- batch_var(msprep_obj)
#   met_vars     <- met_vars(msprep_obj)
#   data         <- msprep_obj$data
# 
#   # normalize/batch correct data
#   data <-
#     switch(normalizeMethod,
#            "ComBat"            = normalize_combat(data, groupingvars, batch, met_vars, transform),
#            "quantile"          = normalize_quantile(data, groupingvars, batch, met_vars, transform),
#            "quantile + ComBat" = normalize_quantile_combat(data, groupingvars, batch, met_vars, transform),
#            "median"            = normalize_median(data, groupingvars, batch, met_vars, transform),
#            "median + ComBat"   = normalize_median_combat(data, groupingvars, batch, met_vars, n_control, controls, transform),
#            "CRMN"              = normalize_crmn(data, groupingvars, batch, met_vars, n_comp, n_control, controls, transform),
#            "RUV"               = normalize_ruv(data, groupingvars, batch, met_vars, n_control, controls, k_ruv, transform),
#            "SVA"               = normalize_sva(data, groupingvars, batch, met_vars, transform),
#            stop("Invalid normalize method - provide argument from list in",
#                 "function definition and help file"))
# 
#   # Prep output object
#   msprep_obj$data  <- data
#   attr(msprep_obj, "normalize_method") <- normalizeMethod
#   attr(msprep_obj, "transformation") <- transform
#   stage(msprep_obj) <- "normalized"
# 
#   # ...and:
#   return(msprep_obj)
# 
# 
# }
# 
# normalize_combat <- function(data, groupingvars, batch, met_vars, transform) {
# 
#   # Combat batch correction
#   rtn <- combat(data, groupingvars, batch, met_vars, transform)
# 
#   return(rtn)
# 
# }
# 
# normalize_quantile <- function(data, groupingvars = NULL, batch = NULL, met_vars, transform) {
#   
#   wide   <- data_to_wide_matrix(data, groupingvars, batch, met_vars)
#   
#   if (transform == "log2") {
#     wide_trans <- log_base2(wide)
#   }
#   else if (transform == "log10") {
#     wide_trans <- log_base10(wide)
#   }
#   else if (transform == "none") {
#     wide_trans <- wide
#   }
# 
#   rwnm  <- rownames(wide_trans)
#   colnm <- colnames(wide_trans)
#   rtn   <- t(normalize.quantiles(t(wide_trans)))
#   rownames(rtn) <- rwnm
#   colnames(rtn) <- colnm
#   
#   rtn <- wide_matrix_to_data(rtn, groupingvars, batch, met_vars)
#   
#   return(rtn)
#   
# }
# 
# 
# 
# normalize_quantile_combat <- function(data, groupingvars, batch, met_vars, transform) {
#   
#   # Check that batches are present in data before proceeding
#   if (is.null(batch)) stop("No batches found for ComBat")
# 
#   # Transform data to wide matrix
#   wide   <- data_to_wide_matrix(data, groupingvars, batch, met_vars)
#   
#   if (transform == "log2") {
#     wide_trans <- log_base2(wide)
#   }
#   else if (transform == "log10") {
#     wide_trans <- log_base10(wide)
#   }
#   else if (transform == "none") {
#     wide_trans <- wide
#   }
#   
#   # Quantile normalization
#   rwnm  <- rownames(wide_trans)
#   colnm <- colnames(wide_trans)
#   rtn   <- t(normalize.quantiles(t(wide_trans)))
#   rownames(rtn) <- rwnm
#   colnames(rtn) <- colnm
# 
#   # This is backtransforming to 'tidy' dataset before combat, which then
#   # transforms back to wide matrix and then back to tidy.  This could be skipped
#   # if have time to generalize/separate out a 'prep for combat' fuinction or add
#   # checks for current format.
#   rtn <- wide_matrix_to_data(rtn, groupingvars, batch, met_vars)
#   
#   # Combat batch correction
#   rtn <- combat(rtn, groupingvars, batch, met_vars, transform = "none")
#   
# 
#   return(rtn)
# 
# }
# 
# 
# 
# normalize_crmn <- function(data, groupingvars, batch, met_vars, n_comp, n_control, controls, transform) {
#   
#   # Check to ensure spike-ins/internal standards are present
#   if(is.null(groupingvars)) stop("No internal standards found for CRMN")
# 
#   crmn_inputs <- create_crmn_inputs(data, groupingvars, batch, met_vars, n_control, controls)
#   # Requires:
#   # - Y -> long, not logged, matrix w/ record id as column names and mz_rt as rownames
#   # - G -> sample_info_modmat -> cell means model matrix of phenotypes/sample info
#   # - isIS -> result of isIS() function
#   # - n_comp -> user input to main ms_normalize() function
#   
#   normed.crmn  <- normalize(crmn_inputs$longmat,
#                             method    = "crmn",
#                             factors   = crmn_inputs$model_matrix,
#                             standards = crmn_inputs$isIS_vec,
#                             ncomp     = n_comp,
#                             lg = TRUE)
#   
#   if (transform == "log2") {
#     wide_trans <- log_base2(normed.crmn)
#   }
#   else if (transform == "log10") {
#     wide_trans <- log_base10(normed.crmn)
#   }
#   else if (transform == "none") {
#     wide_trans <- normed.crmn
#   }
# 
#   final_crmn   <- as.data.frame(t(wide_trans))
#   final_crmn   <- wide_matrix_to_data(final_crmn, groupingvars, batch, met_vars)
# 
#   return(final_crmn)
# 
# }
# 
# normalize_ruv <- function(data, groupingvars, batch, met_vars, n_control, controls, k_ruv, transform) {
#   
#   if(k_ruv > n_control){
#     stop ("k_ruv must be less than or equal to n_control")
#   }
#   
#   wide <- data_to_wide_matrix(data, groupingvars, batch, met_vars)
#   if (transform == "log2") {
#     wide_trans <- log_base2(wide)
#   }
#   else if (transform == "log10") {
#     wide_trans <- log_base10(wide)
#   }
#   else if (transform == "none") {
#     wide_trans <- wide
#   }
#   data_trans <- wide_matrix_to_data(wide_trans, groupingvars, batch, met_vars)
# 
#   ruv_factors <- ruvfactors(data_trans, groupingvars, batch, met_vars, n_control, controls, k_ruv)
#   adjusted    <- genadj(data_trans, groupingvars, batch, met_vars, ruv_factors)
# 
#   cat("\n")
#   return(adjusted)
# 
# }
# 
# 
# normalize_sva <- function(data, groupingvars, batch, met_vars, transform) {
#   if(is.null(groupingvars)){
#     stop("Spike-ins missing for SVA normalization.")
#   }
#   
#   wide <- data_to_wide_matrix(data, groupingvars, batch, met_vars)
#   if (transform == "log2") {
#     wide_trans <- log_base2(wide)
#   }
#   else if (transform == "log10") {
#     wide_trans <- log_base10(wide)
#   }
#   else if (transform == "none") {
#     wide_trans <- wide
#   }
#   data_trans <- wide_matrix_to_data(wide_trans, groupingvars, batch, met_vars)
#   
#   sva_factors <- svafactors(data_trans, groupingvars, batch, met_vars)
#   adjusted    <- genadj(data_trans, groupingvars, batch, met_vars, sva_factors)
# 
#   cat("\n")
#   return(adjusted)
# 
# }
# 
# 
# normalize_median_combat <- function(data, groupingvars, batch, met_vars, n_control, controls, transform) {
# 
# #   # For normalize/median, we need
# #   #   - Y     <- transposed (converts to matrix) final_shift <- final <- dataframed <- input dataset in wide_matrix format
# #   #   - G     <- cell means model matrix of groupingvars
# #   #   - isIS
# #   # TODO: Now figure out what these are and simplist way to get them.
# #   normed.med   <- crmn::normalize(Y, "median", factors = G, standards = isIS)
# #   med_com      <- combat(as.data.frame(t(normed.med)), pheno, batch)
# 
#   # from normalize_crmn -- verify if true 
#   # - Y -> long, not logged, matrix w/ record id as column names and mz_rt as rownames
#   # - G -> sample_info_modmat -> cell means model matrix of phenotypes/sample info
#   # - isIS -> result of isIS() function
#   # - n_comp -> user input to main ms_normalize() function
#   
#   #wide <- data_to_wide_matrix(data, groupingvars, batch, met_vars)
#   
#   #if (transform == "log2") {
#     #wide_trans <- log_base2(wide)
#   #}
#   #else if (transform == "log10") {
#     #wide_trans <- log_base10(wide)
#   #}
#   #else if (transform == "none") {
#     #wide_trans <- wide
#   #}
#   
#   #data_trans <- wide_matrix_to_data(wide_trans, groupingvars, batch, met_vars)
#   
#   #crmn_inputs <- create_crmn_inputs(data_trans, groupingvars, batch, met_vars, n_control, controls)
#   
#   # Requires:
#   # - Y -> long, not logged, matrix w/ record id as column names and mz_rt as rownames
#   # - G -> sample_info_modmat -> cell means model matrix of phenotypes/sample info
#   # - isIS -> result of isIS() function
#   # - n_comp -> user input to main ms_normalize() function
# 
#   #normed_med  <- normalize(object    = crmn_inputs$longmat,
#                            #method    = "median",
#                            #factors   = crmn_inputs$model_matrix,
#                            #standards = crmn_inputs$isIS_vec,
#                            #lg = FALSE)
#   #normed_med <- wide_matrix_to_data(t(normed_med), groupingvars, batch, met_vars)
#   
#   wide <- data_to_wide_matrix(data, groupingvars, batch, met_vars)
#   
#   if (transform == "log2") {
#     wide_trans <- log_base2(wide)
#   }
#   else if (transform == "log10") {
#     wide_trans <- log_base10(wide)
#   }
#   else if (transform == "none") {
#     wide_trans <- wide
#   }
#   
#   #if(transform == "log2" | transform == "log10") {
#    # normed_med <- sweep(wide_trans, 2, apply(wide_trans, 2, median), "/")
#   #}
#   #else {
#     normed_med <- sweep(wide_trans, 1, apply(wide_trans, 1, median), "-")
#   #}
#   
#   normed_med <- wide_matrix_to_data(normed_med, groupingvars, batch, met_vars)
#   rtn        <- combat(normed_med, groupingvars, batch, met_vars, transform = "None")
# 
#   return(rtn)
# 
# }
# 
# normalize_median <- function(data, groupingvars, batch, met_vars, transform) {
#   
#   wide <- data_to_wide_matrix(data, groupingvars, batch, met_vars)
#   
#   if (transform == "log2") {
#     wide_trans <- log_base2(wide)
#   }
#   else if (transform == "log10") {
#     wide_trans <- log_base10(wide)
#   }
#   else if (transform == "none") {
#     wide_trans <- wide
#   }
#   
#   #if(transform == "log2" | transform == "log10") {
#    # normed_med <- sweep(wide_trans, 2, apply(wide_trans, 2, median), "/")
#   #}
#   #else {
#     normed_med <- sweep(wide_trans, 1, apply(wide_trans, 1, median), "-")
#   #}
#   
#   normed_med <- wide_matrix_to_data(normed_med, groupingvars, batch, met_vars)
#   
#   return(normed_med)
# }
# 
# 
# create_crmn_inputs <- function(data, groupingvars, batch, met_vars, n_control, controls) {
# 
#   internal_id <- internal_id_order(groupingvars, batch)
# 
#   sva         <- svafactors(data, groupingvars, batch, met_vars)
#   widedf      <- data_to_wide_matrix(data, groupingvars, batch, met_vars, asmatrix = FALSE)
# 
#   n_compounds <- ncol(widedf)
#   long_data   <- t(widedf) # For crmn::normalize()
# 
#   wide_plus_rownames <- rownames_to_column(widedf) # for model.matrix
#   wide_plus_rownames <- separate(wide_plus_rownames, col = "rowname", into = internal_id, remove = FALSE)
#   sample_info_dat    <- select(wide_plus_rownames, groupingvars)
#   sample_info_dat    <- mutate_if(sample_info_dat, is.character, funs(as.factor))
#   sample_info_modmat <- model.matrix(~ -1 + ., data = sample_info_dat)
# 
#   isIS_vec <- isIS(data, sva, n_control, n_compounds, controls, met_vars)
# 
#   return(structure(list("longmat"      = long_data,
#                         "model_matrix" = sample_info_modmat,
#                         "isIS_vec"     = isIS_vec),
#                    class = "crmn_inputs"))
# 
# }
# 
# # Not entirely sure what this function is/does
# isIS <- function(data, factors, n_control, n_compounds, controls, met_vars) {
#   # compounds -> number of mz_rt
#   # ctlo -> ctl[order(ctl)]
# 
#   control_sum   <- control_summary(data, met_vars)
#   ctl       <- ctl_compounds(control_sum, factors, n_control, controls)
#   ctlo      <- ctl[order(ctl)]
#   j         <- 1
#   isISvec      <- rep(FALSE, n_compounds)
# 
#   for (i in 1:n_compounds) {
#     if (j <= 10) {
#       if (ctlo[j] == i) {
#         isISvec[i] <- TRUE
#         j       <- j + 1
#       } else {
#         isISvec[i] <- FALSE
#       }
#     }
#   }
# 
#   return(isISvec)
# 
# }
# 
# 
# 
# ctl_compounds <- function(control_summary_data, factors = NULL, n_control, controls) {
# 
#   if (length(controls) > 0) { 
#     dat <- controls 
#   } else {
#     dat <- arrange(control_summary_data, .data$counteq0, .data$cv)
#     dat <- dat[1:n_control, "number"]
#     dat <- dat[["number"]]
#   }
# 
#   return(dat)
# 
# }
# 
# 
# 
# combat <- function(data, groupingvars, batch, met_vars, transform) {
#   
#   # Check that batches are present in data before proceeding
#   if (is.null(batch)) stop("No batches found for ComBat")
# 
#   internal_id <- internal_id_order(groupingvars, batch)
# 
#   # Function dev
#   # Convert to wide matrix -- do we want spike in the ID?
#   wide <- data_to_wide_matrix(data, groupingvars, batch, met_vars, asmatrix = FALSE)
#   wide <- rownames_to_column(wide)
#   wide <- separate(wide, col = "rowname",
#                    into = internal_id,
#                    remove = FALSE)
# 
#   # Generate long matrix  for feeding to ComBat()
#   widemat <- data_to_wide_matrix(data, groupingvars, batch, met_vars, asmatrix = TRUE)
#   longmat <- t(widemat)
# 
#   # combat + median is not log'd, otherwise should be
#   if (transform == "log2") {
#     longmat <- log_base2(longmat)
#   }
#   else if (transform == "log10") {
#     longmat <- log_base10(longmat)
#   }
# 
#   # Create model matrix of sample_info variables (soon to be formerly groupingvars)
#   sample_info_dat    <- select(wide, groupingvars)
#   sample_info_dat    <- mutate_if(sample_info_dat, is.character, funs(as.factor))
#   sample_info_modmat <- model.matrix(~ ., data = sample_info_dat)
# 
#   # Create model matrix of batch variable
#   if (length(batch) > 1) stop("only one batch variable allowed for ComBat") 
#   batch_dat <- select(wide, batch)
#   batch_dat <- mutate_if(batch_dat, is.character, funs(as.factor))
#   batch_vec <- batch_dat[, 1]
# 
#   comadj <- sva::ComBat(longmat, mod = sample_info_modmat, batch = batch_vec) 
#   comadj <- t(comadj)
# 
#   comadj_tidy <- wide_matrix_to_data(comadj, groupingvars, batch, met_vars)
# 
#   # What form to return this in??
#   return(comadj_tidy)
# 
# }
# 
# 
# svafactors <- function(data, groupingvars, batch, met_vars) {
# 
#   # This can probably be simplified some more
# 
#   # Create long, logged matrix for sva
#   wide      <- data_to_wide_matrix(data, groupingvars, batch, met_vars, asmatrix = FALSE)
#   
#   long <- t(wide)
# 
#   # Create model matrix of sample_info variables (soon to be formerly groupingvars)
#   internal_id        <- internal_id_order(groupingvars, batch)
#   wide               <- rownames_to_column(wide)
#   wide               <- separate(wide, col = "rowname", into = internal_id, remove = FALSE)
#   sample_info_dat    <- select(wide, groupingvars) 
#   sample_info_dat    <- mutate_if(sample_info_dat, is.character, funs(as.factor))
#   sample_info_modmat <- model.matrix(~ ., data = sample_info_dat)
#   
#  
#   svaadj <- sva(long, sample_info_modmat, method = "irw")
#   rtn    <- svaadj$sv
#   colnames(rtn) <- paste0("f", 1:ncol(rtn))
# 
#   return(rtn)
# 
# }
# 
# 
# ruvfactors <- function(data, groupingvars, batch, met_vars, n_control, controls, k_ruv) {
# 
#   # Prep datasets/matrices
#   wide        <- data_to_wide_matrix(data, groupingvars, batch, met_vars)
#   
#   Y           <- t(wide) # Only data Used in calculation
# 
#   # ctl -- n compounds w/ lowest CV and no missing
#   control_sum  <- control_summary(data, met_vars)
#   ctl <- ctl_compounds(control_summary_data = control_sum, n_control = n_control, controls = controls)
# 
#   # Calculate RUV 
#   #   Uses: Y, ctl, sva_factors
#   Z   <- matrix(rep(1, ncol(Y)))
#   RZY <- Y - Y %*% Z %*% solve(t(Z) %*% Z) %*% t(Z)
#   W   <- svd(RZY[ctl, ])$v
#   W   <- W[, 1:k_ruv]
#   
#   # Format output
#   rtn <- as.matrix(W, ncol = k_ruv)
#   colnames(rtn) <- paste0("f", 1:ncol(rtn))
#   
#   return(rtn)
# 
# }
# 
# 
# genadj <- function(data, groupingvars, batch, met_vars, factors) {
#   # dset      == log_final
#   # factors   == sva or ruv_raw factors/loadings?
#   # compounds == num mz_rt
# 
#   wide      <- data_to_wide_matrix(data, groupingvars, batch, met_vars, asmatrix = FALSE)
# 
#   out <- sapply(1:ncol(wide),
#     function(j) {
#       fit <- lm(wide[, j] ~ as.matrix(factors, ncol = 1))
#       fit$fitted.values
#     })
# 
#   colnames(out) <- colnames(wide)
#   rownames(out) <- rownames(wide)
# 
#   rtn <- wide_matrix_to_data(out, groupingvars, batch,  met_vars)
# 
#   return(rtn)
# 
# }
# 
# 
# control_summary <- function(data, met_vars) {
#   
#   # Summarises dataset for use in control compounds (ctl_compounds) function
#   met_syms <- syms(met_vars)
# 
#   counteq0 <- function(x) sum(x == 0)
#   rtn <- group_by(data, `!!!`(met_syms))
#   rtn <- summarise_at(rtn, vars(.data$abundance_summary),
#                       funs(mean, sd, counteq0)) #, na.rm = TRUE)
#   rtn <- ungroup(rtn)
#   rtn <- mutate(rtn, cv = sd / mean, number = 1:n())
#   rtn <- select(rtn, `!!!`(met_syms), .data$number, .data$mean, .data$sd,
#                 .data$cv, .data$counteq0)
#   return(rtn)
# 
# }
# 
# # Utility functions for ms_normalize()
# log_base2 <- function(wide_matrix) apply(wide_matrix, 2, log2)
# log_base10 <- function(wide_matrix) apply(wide_matrix, 2, log10)
# 
# # ################################################################################
# # # testing setup
# # ################################################################################
# # # New testing
# # dat <- imputed_data
# # data         <- dat$data
# # groupingvars <- MSPrep:::grouping_vars(dat)
# # batch        <- MSPrep:::batch_var(dat)
# ################################################################################
# # Old testing
# #   # data == wide matrix data w/ rownames
# #   # batch == "Operator" -- i.e. 
# #   # pheno == "Spike"
# #   # clindat == Clinical.csv
# #   # link1 == "SubjectID"
# # path_olddata <- system.file("extdata", "old_object.Rda", package = "MSPrep")
# # load(path_olddata)
# # clinical <- read_csv("data-raw/Clinical.csv")
# # pheno    <- "Spike"
# # batch    <- "Operator"
# # link1    <- "SubjectID"
# # olddata  <- as.data.frame(old_readdata_result$sum_data)
# # clindat  <- clinical
# ################################################################################
# 
# 
# 

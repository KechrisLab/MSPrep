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
# #'   load("test2.Rdata")
# #'   
# #'   metafin  <- test2$bpca
# #'   clindat  <- test$clinical
# #'   link1    <- "SubjectID"
# #'   pheno    <- "Spike"
# #'   batch    <- "Operator"
# #'   ncont    <- 10
# #'   controls <- c()
# #'   ncomp    <- 2
# #'   
# #'   test3 <- normdata(metafin, clindat, link1, pheno, batch, ncont = 10,
# #'                     controls, ncomp)
# #'
# #' @export 
# 
# # Katerina's function, unedited (formerly normdata.new_kk)
# normdata <- function (metafin, clindat, link1, pheno, batch, ncont = 10, 
#                       controls = c(), ncomp = 2) 
# {
#   compounds <- ncol(metafin)
#   subjects <- nrow(metafin)
#   metan <- t(normalize.quantiles(t(metafin)))
#   colnames(metan) <- colnames(metafin)
#   rownames(metan) <- rownames(metafin)
#   final <- as.data.frame(metafin)
#   finaln <- metan
#   convert <- function(dset) {
#     s_dset <- as.data.frame(matrix(NA, ncol = ncol(dset), 
#                                    nrow = nrow(dset)))
#     for (j in 1:ncol(dset)) {
#       s_dset[, j] <- log2(as.numeric(dset[, j]))
#     }
#     colnames(s_dset) <- colnames(dset)
#     rownames(s_dset) <- rownames(dset)
#     return(s_dset)
#   }
#   log_final <- convert(final)
#   log_finaln <- convert(finaln)
#   combat <- function(dset, pheno, batch) {
#     compounds <- ncol(dset)
#     temp <- dset
#     temp$mid <- rownames(temp)
#     test <- merge(temp, clindat, by.x = "mid", by.y = as.character(link1))
#     d1 <- t(dset)
#     colnames(d1) <- dset[, 1]
#     d2 <- subset(test, select = pheno)
#     d2a <- model.matrix(~as.factor(d2[, 1]))
#     d3 <- subset(test, select = batch)
#     count <- matrix(0, nrow = compounds, ncol = 1)
#     d1mod <- matrix(0, nrow = compounds, ncol = subjects)
#     for (j in 1:compounds) {
#       d1mod[j, ] <- log2(as.numeric(d1[j, ]))
#     }
#     comadj <- t(ComBat(d1mod, mod = d2a, batch = d3[,1])) # Here the element "d2a" is Gold Stage phenotype and "d3" is the batch.
#     colnames(comadj) <- colnames(dset)
#     rownames(comadj) <- rownames(dset)
#     return(comadj)
#   }
#   final_rc <- combat(as.data.frame(final), pheno, batch)
#   final_com <- combat(as.data.frame(finaln), pheno, batch)
#   svaestimate <- function(dset, pheno, batch) {
#     d1 <- t(dset[, 1:compounds])
#     colnames(d1) <- rownames(dset)
#     rownames(d1) <- colnames(dset)
#     dset$mid <- rownames(dset)
#     test <- merge(dset, clindat, by.x = "mid", by.y = as.character(link1))
#     d2 <- subset(test, select = c(batch))
#     d3 <- subset(test, select = c(pheno))
#     colnames(d3)[1] <- "pheno"
#     d3m <- model.matrix(~as.factor(pheno), data = d3)
#     count <- matrix(0, nrow = compounds, ncol = 1)
#     d1mod <- d1
#     svaadj <- sva(d1mod, d3m, method = "irw")
#     return(as.matrix(svaadj$sv, ncol = svaadj$n.sv))
#   }
#   sva <- as.data.frame(svaestimate(log_final, pheno, batch))
#   for (i in 1:ncol(sva)) {
#     colnames(sva)[i] <- paste("f", i, sep = "")
#   }
#   ruv <- function(dset, k, ctl, pheno) {
#     d1 <- t(dset[, 1:compounds])
#     count <- matrix(0, nrow = compounds, ncol = 1)
#     d1mod <- d1
#     Y <- d1mod
#     dset$mid <- rownames(dset)
#     test <- merge(dset, clindat, by.x = "mid", by.y = link1)
#     X <- subset(test, select = pheno)
#     Z = matrix(rep(1, ncol(Y)))
#     RZY = Y - Y %*% Z %*% solve(t(Z) %*% Z) %*% t(Z)
#     W = svd(RZY[ctl, ])$v
#     W = W[, 1:k]
#     return(as.matrix(W, ncol = k))
#   }
#   dset <- final
#   control <- matrix(NA, nrow = compounds, ncol = 6)
#   for (j in 1:compounds) {
#     control[j, 1] <- colnames(dset)[j]
#     control[j, 2] <- j
#     control[j, 3] <- mean(dset[, j], na.rm = TRUE)
#     control[j, 4] <- sd(dset[, j], na.rm = TRUE)
#     control[j, 5] <- sd(dset[, j], na.rm = TRUE)/mean(dset[, 
#                                                       j], na.rm = TRUE)
#     control[j, 6] <- sum(dset[, j] == 0, na.rm = TRUE)
#   }
#   control2 <- control[order(control[, 6], control[, 5]), ]
#   ctl <- as.numeric(control2[, 2][1:max(ncont, ncol(sva) + 
#                                         5)])
#   if (length(controls) > 0) {
#     ctl <- controls
#   }
#   contout <- cbind(ctl, colnames(log_final)[ctl])
#   ruv_raw <- as.data.frame(ruv(log_final, ncol(sva), ctl, pheno))
#   for (i in 1:ncol(ruv_raw)) {
#     colnames(ruv_raw)[i] <- paste("f", i, sep = "")
#   }
#   genadj <- function(dset, factors) {
#     out <- sapply(1:compounds, function(j) lm(dset[, j] ~ 
#                                               as.matrix(factors, ncol = 1))$fitted.values)
#     colnames(out) <- colnames(dset)
#     rownames(out) <- rownames(dset)
#     return(out)
#   }
#   ruv_adj <- genadj(log_final, ruv_raw)
#   sva_adj <- genadj(log_final, sva)
#   final_shift <- final
#   Y <- t(final_shift)
#   Ya <- t(final_shift)[ctl, ]
#   Ys <- t(final_shift)[-ctl, ]
#   j <- 1
#   isIS <- rep(FALSE, compounds)
#   ctlo <- ctl[order(ctl)]
#   for (i in 1:compounds) {
#     if (j <= 10) {
#       if (ctlo[j] == i) {
#         isIS[i] <- TRUE
#         j <- j + 1
#       }
#       else {
#         isIS[i] <- FALSE
#       }
#     }
#   }
#   dset <- final_shift
#   dset$mid <- rownames(dset)
#   test <- merge(dset, clindat, by.x = "mid", by.y = link1)
#   X <- subset(test, select = c(pheno))
#   colnames(X)[1] <- "pheno"
#   G <- model.matrix(~-1 + as.factor(X$pheno)) # Here, "G" comes from the Gold Stage phenotype.
#   normed.crmn <- normalize(Y, "crmn", factors = G, standards = isIS, 
#                            ncomp = ncomp)
#   lnormed.crmn <- convert(normed.crmn)
#   normed.med <- normalize(Y, "median", factors = G, standards = isIS)
#   lnormed.med <- convert(normed.med)
#   final_crmn <- as.data.frame(t(lnormed.crmn))
#   final_med <- as.data.frame(t(lnormed.med))
#   med_com <- combat(as.data.frame(t(normed.med)), pheno, batch)
#   list(log_data = log_final, log_data_combat = final_rc, log_quant = log_finaln, 
#        log_quant_combat = final_com, med_adj = final_med, med_combat = med_com, 
#        sva_factors = sva, sva_adj = sva_adj, ruv_factors = ruv_raw, 
#        ruv_adj = ruv_adj, crmn_adj = final_crmn, controls = contout)
# }
# 
# ## Katerina's function, edited by me for multiple phenotypes and only ComBat
# 
# normdata.new_multpheno <- function (metafin, clindat, link1, pheno, batch, ncont = 10, 
#                                     controls = c(), ncomp = 2, combat_log = F, do_crmn = F) 
# {
#   # if(length(pheno) == 1) stop("if pheno is only 1 element, use normdata.new_kk")
#   if(length(pheno) > 1){
#     remove_subjects <- which(apply(clindat[,pheno], 1, function(x) any(is.na(x))))
#   }else if(length(pheno) == 1){
#     remove_subjects <- which(is.na(clindat[,pheno]))
#   }else remove_subjects <- numeric(0)
#   if(length(remove_subjects) > 1){
#     metafin <- metafin[-remove_subjects, ]
#     clindat <- clindat[-remove_subjects, ]
#   }
#   compounds <- ncol(metafin)
#   subjects <- nrow(metafin)
#   metan <- t(normalize.quantiles(t(metafin)))
#   colnames(metan) <- colnames(metafin)
#   rownames(metan) <- rownames(metafin)
#   final <- as.data.frame(metafin)
#   finaln <- metan
#   convert <- function(dset) {
#     s_dset <- as.data.frame(matrix(NA, ncol = ncol(dset), 
#                                    nrow = nrow(dset)))
#     for (j in 1:ncol(dset)) {
#       s_dset[, j] <- log2(as.numeric(dset[, j]))
#     }
#     colnames(s_dset) <- colnames(dset)
#     rownames(s_dset) <- rownames(dset)
#     return(s_dset)
#   }
#   # log_final <- convert(final)
#   # log_finaln <- convert(finaln)
#   combat <- function(dset, pheno, batch) {
#     compounds <- ncol(dset)
#     temp <- dset
#     temp$mid <- rownames(temp)
#     test <- merge(temp, clindat, by.x = "mid", by.y = as.character(link1))
#     d1 <- t(dset)
#     colnames(d1) <- dset[, 1]
#     d2 <- subset(test, select = pheno)
#     # d2_p <- apply(d2, 1, paste0, collapse = "")
#     pheno_oc <- sapply(pheno, function(x) length(unique(d2[,x])) > 7)
#     if(length(pheno) > 0){
#       d2a <- eval(parse(text = sprintf("model.matrix(~%s,data=d2)", paste0(sapply(pheno, function(x, oc){
#                                                           if(x %in% pheno[oc]){
#                                                             res <- x
#                                                           }else res <- paste0("as.factor(", x, ")")
#                                                           return(res)
# }, pheno_oc), collapse = "+") )))
#     }else d2a <- NULL
#     d3 <- subset(test, select = batch)
#     count <- matrix(0, nrow = compounds, ncol = 1)
#     d1mod <- matrix(0, nrow = compounds, ncol = nrow(d3))#subjects)
# 
#     for (j in 1:compounds) {
#       ifelse(combat_log, d1mod[j, ] <- log2(as.numeric(d1[j, ])), d1mod[j, ] <- as.numeric(d1[j, ]))
#     }
#     comadj <- t(ComBat(d1mod, mod = d2a, batch = d3[,1]))
#     colnames(comadj) <- colnames(dset)
#     rownames(comadj) <- rownames(dset)
#     return(comadj)
#   }
#   final_rc <- combat(as.data.frame(final), pheno, batch)
#   final_com <- combat(as.data.frame(finaln), pheno, batch)
#   # svaestimate <- function(dset, pheno, batch) {
#   #   d1 <- t(dset[, 1:compounds])
#   #   colnames(d1) <- rownames(dset)
#   #   rownames(d1) <- colnames(dset)
#   #   dset$mid <- rownames(dset)
#   #   test <- merge(dset, clindat, by.x = "mid", by.y = as.character(link1))
#   #   d2 <- subset(test, select = c(batch))
#   #   d3 <- subset(test, select = c(pheno))
#   #   colnames(d3)[1] <- "pheno"
#   #   d3m <- model.matrix(~as.factor(pheno), data = d3)
#   #   count <- matrix(0, nrow = compounds, ncol = 1)
#   #   d1mod <- d1
#   #   svaadj <- sva(d1mod, d3m, method = "irw")
#   #   return(as.matrix(svaadj$sv, ncol = svaadj$n.sv))
#   # }
#   # sva <- as.data.frame(svaestimate(log_final, pheno, batch))
#   # for (i in 1:ncol(sva)) {
#   #   colnames(sva)[i] <- paste("f", i, sep = "")
#   # }
#   # ruv <- function(dset, k, ctl, pheno) {
#   #   d1 <- t(dset[, 1:compounds])
#   #   count <- matrix(0, nrow = compounds, ncol = 1)
#   #   d1mod <- d1
#   #   Y <- d1mod
#   #   dset$mid <- rownames(dset)
#   #   test <- merge(dset, clindat, by.x = "mid", by.y = link1)
#   #   X <- subset(test, select = pheno)
#   #   Z = matrix(rep(1, ncol(Y)))
#   #   RZY = Y - Y %*% Z %*% solve(t(Z) %*% Z) %*% t(Z)
#   #   W = svd(RZY[ctl, ])$v
#   #   W = W[, 1:k]
#   #   return(as.matrix(W, ncol = k))
#   # }
#   if(do_crmn){
#     dset <- final
#     control <- matrix(NA, nrow = compounds, ncol = 6)
#     for (j in 1:compounds) {
#       control[j, 1] <- colnames(dset)[j]
#       control[j, 2] <- j
#       control[j, 3] <- mean(dset[, j], na.rm = TRUE)
#       control[j, 4] <- sd(dset[, j], na.rm = TRUE)
#       control[j, 5] <- sd(dset[, j], na.rm = TRUE)/mean(dset[, 
#                                                         j], na.rm = TRUE)
#       control[j, 6] <- sum(dset[, j] == 0, na.rm = TRUE)
#     }
#     control2 <- control[order(control[, 6], control[, 5]), ]
#     ctl <- as.numeric(control2[, 2][1:max(ncont, ncol(sva) + 
#                                           5)])
#     if (length(controls) > 0) {
#       ctl <- controls
#     }
#     # contout <- cbind(ctl, colnames(log_final)[ctl])
#     # ruv_raw <- as.data.frame(ruv(log_final, ncol(sva), ctl, pheno))
#     # for (i in 1:ncol(ruv_raw)) {
#     #   colnames(ruv_raw)[i] <- paste("f", i, sep = "")
#     # }
#     # genadj <- function(dset, factors) {
#     #   out <- sapply(1:compounds, function(j) lm(dset[, j] ~ 
#     #     as.matrix(factors, ncol = 1))$fitted.values)
#     #   colnames(out) <- colnames(dset)
#     #   rownames(out) <- rownames(dset)
#     #   return(out)
#     # }
#     # ruv_adj <- genadj(log_final, ruv_raw)
#     # sva_adj <- genadj(log_final, sva)
#     final_shift <- final
#     Y <- t(final_shift)
#     # Ya <- t(final_shift)[ctl, ]
#     # Ys <- t(final_shift)[-ctl, ]
#     j <- 1
#     isIS <- rep(FALSE, compounds)
#     ctlo <- ctl[order(ctl)]
#     for (i in 1:compounds) {
#       if (j <= 10) {
#         if (ctlo[j] == i) {
#           isIS[i] <- TRUE
#           j <- j + 1
#         }
#         else {
#           isIS[i] <- FALSE
#         }
#       }
#     }
#     dset <- final_shift
#     dset$mid <- rownames(dset)
#     test <- merge(dset, clindat, by.x = "mid", by.y = link1)
#     X <- subset(test, select = c(pheno))
#     pheno_oc <- sapply(pheno, function(x) length(unique(X[,x])) > 7)
#     if(length(pheno) > 0){
#       G <- eval(parse(text = sprintf("model.matrix(~-1 + %s,data=X)", paste0(sapply(pheno, function(x, oc){
#                                                         if(x %in% pheno[oc]){
#                                                           res <- x
#                                                         }else res <- paste0("as.factor(", x, ")")
#                                                         return(res)
# }, pheno_oc), collapse = "+") )))
#     }else G <- model.matrix(~1, data = X)
# 
#     # # colnames(X)[1] <- "pheno"
#     # # d2_p <- apply(d2, 1, paste0, collapse = "")
#     # # d2a <- model.matrix(~as.factor(d2_p))
#     # G <- model.matrix(~-1 + as.factor(apply(X, 1, paste0, collapse = "")))
#     normed.crmn <- normalize(Y, "crmn", factors = G, standards = isIS, 
#                              ncomp = ncomp)
#     lnormed.crmn <- convert(normed.crmn)
#     normed.med <- normalize(Y, "median", factors = G, standards = isIS)
#     lnormed.med <- convert(normed.med)
#     final_crmn <- as.data.frame(t(lnormed.crmn))
#     final_med <- as.data.frame(t(lnormed.med))
#   }else final_crmn <- final_med <- NULL
#   # med_com <- combat(as.data.frame(t(normed.med)), pheno, batch)
#   list(data = final, data_combat = final_rc, quant = finaln, 
#        quant_combat = final_com, crmn_adj = final_crmn, med_adj = final_med)
#   # data_combat is logged if using quant normalization
#   #  med_combat = med_com, 
#   # sva_factors = sva, sva_adj = sva_adj, ruv_factors = ruv_raw, 
#   # ruv_adj = ruv_adj, controls = contout)
# } 
# 
# 
# ## Katerina's function, edited to not use any phenotypes
# 
# normdata.new_nopheno <- function (metafin, clindat, link1, batch, ncont = 10, 
#                                   controls = c(), ncomp = 2) 
# {
#   compounds <- ncol(metafin)
#   subjects <- nrow(metafin)
#   metan <- t(normalize.quantiles(t(metafin)))
#   colnames(metan) <- colnames(metafin)
#   rownames(metan) <- rownames(metafin)
#   final <- as.data.frame(metafin)
#   finaln <- metan
#   convert <- function(dset) {
#     s_dset <- as.data.frame(matrix(NA, ncol = ncol(dset), 
#                                    nrow = nrow(dset)))
#     for (j in 1:ncol(dset)) {
#       s_dset[, j] <- log2(as.numeric(dset[, j]))
#     }
#     colnames(s_dset) <- colnames(dset)
#     rownames(s_dset) <- rownames(dset)
#     return(s_dset)
#   }
#   # log_final <- convert(final)
#   # log_finaln <- convert(finaln)
#   combat <- function(dset, pheno = NULL, mybatch) {
#     compounds <- ncol(dset)
#     temp <- dset
#     temp$mid <- rownames(temp)
#     testqq <- merge(temp, clindat, by.x = "mid", by.y = as.character(link1))
#     d1 <- t(dset)
#     colnames(d1) <- dset[, 1]
#     d2 <- NULL# subset(testqq, select = pheno)
#     d2a <- NULL# model.matrix(~as.factor(d2[, 1]))
#     d3 <- subset(testqq, select = mybatch)
#     # print(dim(d3))
#     count <- matrix(0, nrow = compounds, ncol = 1)
#     d1mod <- matrix(0, nrow = compounds, ncol = subjects)
#     for (j in 1:compounds) {
#       d1mod[j, ] <- as.numeric(d1[j, ]) #log2(as.numeric(d1[j, ]))
#     }
#     comadj <- t(ComBat(d1mod, mod = NULL, batch = d3[,1])) # Here the element "d2a" is Gold Stage phenotype and "d3" is the batch.
#     colnames(comadj) <- colnames(dset)
#     rownames(comadj) <- rownames(dset)
#     return(comadj)
#   }
#   final_rc <- combat(as.data.frame(final), mybatch = batch)
#   final_com <- combat(as.data.frame(finaln), pheno, mybatch = batch)
#   # svaestimate <- function(dset, pheno = NULL, mybatch) {
#   #   d1 <- t(dset[, 1:compounds])
#   #   colnames(d1) <- rownames(dset)
#   #   rownames(d1) <- colnames(dset)
#   #   dset$mid <- rownames(dset)
#   #   test <- merge(dset, clindat, by.x = "mid", by.y = as.character(link1))
#   #   d2 <- subset(test, select = c(batch))
#   #   d3 <- subset(test, select = c(pheno))
#   #   colnames(d3)[1] <- "pheno"
#   #   d3m <- model.matrix(~as.factor(pheno), data = d3)
#   #   count <- matrix(0, nrow = compounds, ncol = 1)
#   #   d1mod <- d1
#   #   svaadj <- sva(d1mod, d3m, method = "irw")
#   #   return(as.matrix(svaadj$sv, ncol = svaadj$n.sv))
#   # }
#   # sva <- as.data.frame(svaestimate(log_final, pheno, batch))
#   # for (i in 1:ncol(sva)) {
#   #   colnames(sva)[i] <- paste("f", i, sep = "")
#   # }
#   # ruv <- function(dset, k, ctl, pheno) {
#   #   d1 <- t(dset[, 1:compounds])
#   #   count <- matrix(0, nrow = compounds, ncol = 1)
#   #   d1mod <- d1
#   #   Y <- d1mod
#   #   dset$mid <- rownames(dset)
#   #   test <- merge(dset, clindat, by.x = "mid", by.y = link1)
#   #   X <- subset(test, select = pheno)
#   #   Z = matrix(rep(1, ncol(Y)))
#   #   RZY = Y - Y %*% Z %*% solve(t(Z) %*% Z) %*% t(Z)
#   #   W = svd(RZY[ctl, ])$v
#   #   W = W[, 1:k]
#   #   return(as.matrix(W, ncol = k))
#   # }
#   dset <- final
#   control <- matrix(NA, nrow = compounds, ncol = 6)
#   for (j in 1:compounds) {
#     control[j, 1] <- colnames(dset)[j]
#     control[j, 2] <- j
#     control[j, 3] <- mean(dset[, j], na.rm = TRUE)
#     control[j, 4] <- sd(dset[, j], na.rm = TRUE)
#     control[j, 5] <- sd(dset[, j], na.rm = TRUE)/mean(dset[, 
#                                                       j], na.rm = TRUE)
#     control[j, 6] <- sum(dset[, j] == 0, na.rm = TRUE)
#   }
#   control2 <- control[order(control[, 6], control[, 5]), ]
#   ctl <- as.numeric(control2[, 2][1:max(ncont, ncol(sva) + 
#                                         5)])
#   # if (length(controls) > 0) {
#   #   ctl <- controls
#   # }
#   # contout <- cbind(ctl, colnames(log_final)[ctl])
#   # ruv_raw <- as.data.frame(ruv(log_final, ncol(sva), ctl, pheno))
#   # for (i in 1:ncol(ruv_raw)) {
#   #   colnames(ruv_raw)[i] <- paste("f", i, sep = "")
#   # }
#   # genadj <- function(dset, factors) {
#   #   out <- sapply(1:compounds, function(j) lm(dset[, j] ~ 
#   #     as.matrix(factors, ncol = 1))$fitted.values)
#   #   colnames(out) <- colnames(dset)
#   #   rownames(out) <- rownames(dset)
#   #   return(out)
#   # }
#   # ruv_adj <- genadj(log_final, ruv_raw)
#   # sva_adj <- genadj(log_final, sva)
#   final_shift <- final
#   Y <- t(final_shift)
#   # Ya <- t(final_shift)[ctl, ]
#   # Ys <- t(final_shift)[-ctl, ]
#   j <- 1
#   isIS <- rep(FALSE, compounds)
#   ctlo <- ctl[order(ctl)]
#   for (i in 1:compounds) {
#     if (j <= 10) {
#       if (ctlo[j] == i) {
#         isIS[i] <- TRUE
#         j <- j + 1
#       }
#       else {
#         isIS[i] <- FALSE
#       }
#     }
#   }
#   # dset <- final_shift
#   # dset$mid <- rownames(dset)
#   # test <- merge(dset, clindat, by.x = "mid", by.y = link1)
#   # X <- subset(test, select = c(pheno))
#   # colnames(X)[1] <- "pheno"
#   G <- model.matrix(~-1 + rep(1, ncol(Y))) # Here, "G" comes from the Gold Stage phenotype.
#   # G <- model.matrix(~-1 + as.factor(apply(X, 1, paste0, collapse = "")))
#   normed.crmn <- normalize(Y, "crmn", factors = G, standards = isIS,
#                            ncomp = ncomp)
#   # lnormed.crmn <- convert(normed.crmn)
#   # normed.med <- normalize(Y, "median", factors = G, standards = isIS)
#   # lnormed.med <- convert(normed.med)
#   final_crmn <- as.data.frame(t(lnormed.crmn))
#   # final_med <- as.data.frame(t(lnormed.med))
#   # med_com <- combat(as.data.frame(t(normed.med)), pheno, batch)
#   list(data = final, data_combat = final_rc, quant = finaln, 
#        quant_combat = final_com)
# }
# 

# 
# # Additional imputation methods added by Harrison
# 
# 
# filterft_hpl <- function (metaf, filterpercent = 0.5) 
# {
#   print("Filtering")
#   count <- matrix(NA, nrow = ncol(metaf), ncol = 1)
#   rownames(count) <- colnames(metaf)
#   toss <- matrix(1, nrow = eval(ncol(metaf)), ncol = 1)
#   for (j in 1:eval(ncol(metaf))) {
#     k <- 0  
#     for (i in 1:dim(metaf)[1]) {
#       # sean jacobson added this so that missing data would be ignored
#       if(is.na(metaf[i, j])) next
#       if (metaf[i, j] == 0) {
#         k <- k + 1
#       }
#     }
#     count[j, 1] <- 1 - (k/nrow(metaf))
#     if (k >= eval(dim(metaf)[1] * (1 - filterpercent))) {
#       toss[j, 1] <- 0
#     }
#   }
#   colnames(toss) <- "toss"
#   tests <- cbind(toss, t(metaf))
#   metafin <- t(subset(tests, toss != 0)[, -1])
#   metaimp <- matrix(NA, nrow = nrow(metafin), ncol = ncol(metafin))
#   for (i in 1:nrow(metafin)) {
#     for (j in 1:ncol(metafin)) {
#       # sean jacobson added this so that missing data would be ignored
#       # if(is.na(metafin[i, j])) next
#       if (metafin[i, j] == 0) {
#         metaimp[i, j] <- NA
#       }
#       else {
#         metaimp[i, j] <- metafin[i, j]
#       }
#     }
#   }
#   present <- matrix(1, nrow = nrow(metaf), ncol = ncol(metaf))
#   colnames(present) <- colnames(metaf)
#   rownames(present) <- rownames(metaf)
#   for (j in 1:eval(ncol(metaf))) {
#     for (i in 1:dim(metaf)[1]) {
#       # sean jacobson added this for missing data
#       # if(is.na(metaf[i, j])){present[i,j] <- NA; next}
#       if (metaf[i, j] == 0) {
#         present[i, j] <- 0
#       }
#     }
#   }
#   list(metaf = metaf, metafin = metafin, metaimp = metaimp)
# }
# 
# impute.zero = function(metafin, metaimp) {
#   # Half minimum value imputation
#   print("Imputing zeros")
#   withzero <- metafin
#   colnames(withzero) <- colnames(metafin)
#   rownames(withzero) <- rownames(metafin)
#   for (i in 1:nrow(metafin)) {
#     for (j in 1:ncol(metafin)) {
#       # sean jacobson added this for missing data
#       # if(is.na(metafin[i, j])){withzero[i,j] <- NA; next}
#       if (metafin[i, j] == 0) {
#         withzero[i, j] <- 0.0001
#       }
#       else {
#         withzero[i, j] <- metafin[i, j]
#       }
#     }
#   }
#   as.data.frame(withzero)
# }
# 
# impute.min = function(metafin, metaimp) {
#   # Half minimum value imputation
#   print("Imputing half minimum")
#   minval <- metafin
#   colnames(minval) <- colnames(metafin)
#   rownames(minval) <- rownames(metafin)
#   for (i in 1:nrow(metafin)) {
#     for (j in 1:ncol(metafin)) {
#       # sean jacobson added this for missing data
#       # if(is.na(metafin[i, j])){minval[i,j] <- NA; next}
#       if (metafin[i, j] == 0) {
#         minval[i, j] <- summary(metaimp[, j])[1]/2
#       }
#       else {
#         minval[i, j] <- metafin[i, j]
#       }
#     }
#   }
#   as.data.frame(minval)
# }
# 
# impute.mean = function(metafin, metaimp) {
#   # Mean value imputation
#   print("Imputing mean")
#   meanval <- metafin
#   colnames(meanval) <- colnames(metafin)
#   rownames(meanval) <- rownames(metafin)
#   for (i in 1:nrow(metafin)) {
#     for (j in 1:ncol(metafin)) {
#       if (metafin[i, j] == 0) {
#         meanval[i, j] <- summary(metaimp[, j])[4]
#       }
#       else {
#         meanval[i, j] <- metafin[i, j]
#       }
#     }
#   }
#   as.data.frame(meanval)
# }
# 
# impute.median = function(metafin, metaimp) {
#   # Median value imputation
#   print("Imputing median")
#   medianval <- metafin
#   colnames(medianval) <- colnames(metafin)
#   rownames(medianval) <- rownames(metafin)
#   for (i in 1:nrow(metafin)) {
#     for (j in 1:ncol(metafin)) {
#       if (metafin[i, j] == 0) {
#         medianval[i, j] <- summary(metaimp[, j])[3]
#       }
#       else {
#         medianval[i, j] <- metafin[i, j]
#       }
#     }
#   }
#   as.data.frame(medianval)
# }
# 
# impute.bpca = function(metafin, metaimp, minval) {
#   # Bayesian PCA imputation
#   print("B PCA imputation")
#   metabpca <- pca_sj(metaimp, nPcs = 3, method = "bpca")
#   bpca <- completeObs(metabpca)
#   colnames(bpca) <- colnames(metafin)
#   rownames(bpca) <- rownames(metafin)
#   coldrop <- 1
#   for (i in 1:nrow(bpca)) {
#     for (j in 1:ncol(bpca)) {
#       if (bpca[i, j] < 0) {
#         bpca[i, j] <- minval[i, j]
#       }
#     }
#   }
#   as.data.frame(bpca)
# }
# 
# impute.knn = function(metafin, metaimp) {
#   #k-Nearest Neighbors imputation
#   print("KNN imputation")
#   knnimpute <- kNN(as.data.frame(metaimp), k=5)[1:dim(metafin)[1], 1:dim(metafin)[2]]
#   colnames(knnimpute) <- colnames(metafin)
#   rownames(knnimpute) <- rownames(metafin)
#   # coldrop <- 1
#   # for (i in 1:nrow(knnimpute)) {
#   #   for (j in 1:ncol(knnimpute)) {
#   #     if (knnimpute[i, j] < 0) {
#   #       knnimpute[i, j] <- minval[i, j]
#   #     }
#   #   }
#   # }
#   as.data.frame(knnimpute)
# }

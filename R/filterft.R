### Filterft to only include compounds that are found in specified percentage of subjects and perform imputation of missing data


filterft_sj <- function (metaf, filterpercent = 0.5) 
{
  count <- matrix(NA, nrow = ncol(metaf), ncol = 1)
  rownames(count) <- colnames(metaf)
  toss <- matrix(1, nrow = eval(ncol(metaf)), ncol = 1)
  for (j in 1:eval(ncol(metaf))) {
    k <- 0
    for (i in 1:dim(metaf)[1]) {
      # sean jacobson added this so that missing data would be ignored
      if(is.na(metaf[i, j])) next
      if (metaf[i, j] == 0) {
        k <- k + 1
      }
    }
    count[j, 1] <- 1 - (k/nrow(metaf))
    if (k >= eval(dim(metaf)[1] * (1 - filterpercent))) {
      toss[j, 1] <- 0
    }
  }
  colnames(toss) <- "toss"
  tests <- cbind(toss, t(metaf))
  metafin <- t(subset(tests, toss != 0)[, -1])
  metaimp <- matrix(NA, nrow = nrow(metafin), ncol = ncol(metafin))
  for (i in 1:nrow(metafin)) {
    for (j in 1:ncol(metafin)) {
      # sean jacobson added this so that missing data would be ignored
      # if(is.na(metafin[i, j])) next
      if (metafin[i, j] == 0) {
        metaimp[i, j] <- NA
      }
      else {
        metaimp[i, j] <- metafin[i, j]
      }
    }
  }
  present <- matrix(1, nrow = nrow(metaf), ncol = ncol(metaf))
  colnames(present) <- colnames(metaf)
  rownames(present) <- rownames(metaf)
  for (j in 1:eval(ncol(metaf))) {
    for (i in 1:dim(metaf)[1]) {
      # sean jacobson added this for missing data
      # if(is.na(metaf[i, j])){present[i,j] <- NA; next}
      if (metaf[i, j] == 0) {
        present[i, j] <- 0
      }
    }
  }
  minval <- metafin
  for (i in 1:nrow(metaf)) {
    for (j in 1:ncol(metafin)) {
      # sean jacobson added this for missing data
      # if(is.na(metafin[i, j])){minval[i,j] <- NA; next}
      if (metafin[i, j] == 0) {
        minval[i, j] <- summary(metaimp[, j])[1]/2
      }
      else {
        minval[i, j] <- metafin[i, j]
      }
    }
  }
  metabpca <- pca_sj(metaimp, nPcs = 3, method = "bpca")
  bpca <- completeObs(metabpca)
  colnames(bpca) <- colnames(metafin)
  rownames(bpca) <- rownames(metafin)
  coldrop <- 1
  for (i in 1:nrow(bpca)) {
    for (j in 1:ncol(bpca)) {
      if (bpca[i, j] < 0) {
        bpca[i, j] <- minval[i, j]
      }
    }
  }
  list(minval = as.data.frame(minval), withzero = as.data.frame(metafin), 
    bpca = as.data.frame(bpca), count = as.data.frame(count))
}



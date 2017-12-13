
#' Generates diagnostic plots for normalization routine
#' 
#' Generates a PCA plot and three boxplots for each normalization method.
#' 
#' @param testobj Object output from normdata() function
#' @param clindat Name of clinical dataset.
#' @param link1 Column in clinical dataset that links subject ID to row names
#' in quantification datasets.
#' @param batch Column name for batch identifier
#' @param pheno Column name for phenotype identifier
#' @param directory Directory to output the diagnostic plots to.
#' @param ylim1 Limits for median graphs.
#' @param ylim2 Limits for raw, quant, RUV, and SVA graphs.
#' @param ylim3 Limits for CRMN graph.
#' @return Placeholder
#' @examples
#'   load("test.Rdata")
#'   load("test3.Rdata")
#'   
#'   testobj   <- test3
#'   clindat   <- test$clinical
#'   link1     <- "SubjectID"
#'   pheno     <- "Spike"
#'   batch     <- "Operator"
#'   directory <- "C:/Users/Hazy/Dropbox/Metabolomics/Programs/MSProcess/"
#'   ylim2     <- c(10,28)
#'   ### For median
#'   ylim1     <- c(-15,15)
#'   ### for crmn
#'   ylim3     <- c(18,37)
#' 
#'   diagnosticgen(testobj, clindat, link1, batch, pheno, directory, ylim1,
#'                 ylim2, ylim3)
#'
#' @export 
diagnosticgen <- function (testobj,
                           clindat,
                           link1,
                           batch,
                           pheno,
                           directory,
                           ylim1,
                           ylim2,
                           ylim3)
{
  pca.point.pch.batch = function(model.id) {
    if (model.id == "null") {
      return("black")
    } else if (model.id == 1) {
      return("1")
    } else if (model.id == 2) {
      return("2")
    } else if (model.id == 3) {
      return("3")
    } else if (model.id == 4) {
      return("4")
    } else if (model.id == 5) {
      return("5")
    } else if (model.id == 6) {
      return("6")
    } else if (model.id == 7) {
      return("7")
    } else if (model.id == 8) {
      return("8")
    } else if (model.id == 9) {
      return("9")
    } else if (model.id == 10) {
      return("A")
    } else if (model.id == 11) {
      return("B")
    } else if (model.id == 13) {
      return("C")
    } else if (model.id == 14) {
      return("D")
    } else if (model.id == 15) {
      return("E")
    } else {
      return("99")
    }
  }

  pca.point.color.change = function(model.id) {
    if (is.na(model.id)) {
      return("black")
    } else if (model.id == 0) {
      return("green")
    } else if (model.id == 1) {
      return("blue")
    } else if (model.id == 2) {
      return("yellow")
    } else if (model.id == 3) {
      return("red")
    } else if (model.id == 4) {
      return("purple")
    } else {
      return("black")
    }
  }


  pcaplots <- function(dset, tlab, link1, batch, pheno) {
    dset$mid <- rownames(dset)
    test <- merge(dset, clindat, by.x = "mid", by.y = as.character(link1))
    pca1 <- prcomp(test[, 2:ncol(dset)])
    title = c(paste("PCA for ", tlab, sep = ""))
    title2 = c("Color Indicates Phenotype, Number Indicates Batch")
    xlab1 = paste("PCA 1, % Var = ", summary(pca1)$importance[2, 1], sep = "")
    ylab1 = paste("PCA 2, % Var = ", summary(pca1)$importance[2, 2], sep = "")
    temp2 <- subset(test, select = c(batch, pheno))
    colnames(temp2) <- c("batch", "pheno")
    plot(pca1$x[, 1], pca1$x[, 2],
         main = title,
         sub = title2,
         xlab = xlab1,
         ylab = ylab1,
         pch = sapply(temp2$batch, pca.point.pch.batch),
         col = sapply(temp2$pheno, pca.point.color.change))
    return()
  }
  stackdata <- function(dset, link1, pheno, batch) {
    temp <- matrix(NA, ncol = 5, nrow = ncol(dset) * nrow(dset))
    comps <- ncol(dset)
    subjects <- nrow(dset)
    dset$mid <- rownames(dset)
    dset <- merge(dset, clindat, by.x = "mid", by.y = link1)
    for (j in 0:eval(comps - 1)) {
      for (i in 1:subjects) {
        temp[(subjects * j + i), 1] <- as.character(dset[i, j + 2])
        temp[(subjects * j + i), 2] <- as.character(dset$mid[i])
        temp[(subjects * j + i), 3] <- as.character(subset(dset, select = pheno)[i, 1])
        temp[(subjects * j + i), 4] <- as.character(subset(dset, select = batch)[i, 1])
      }
    }
    return(temp)
  }

  data_r <- stackdata(testobj$log_data, link1, pheno, batch)
  data_n <- stackdata(testobj$log_quant, link1, pheno, batch)
  data_rc <- stackdata(as.data.frame(testobj$log_data_combat), link1, pheno, batch)
  data_nc <- stackdata(as.data.frame(testobj$log_quant_combat), link1, pheno, batch)
  data_med <- stackdata(testobj$med_adj, link1, pheno, batch)
  data_medc <- stackdata(as.data.frame(testobj$med_combat), link1, pheno, batch)
  data_ruv <- stackdata(as.data.frame(testobj$ruv_adj), link1, pheno, batch)
  data_sva <- stackdata(as.data.frame(testobj$sva_adj), link1, pheno, batch)
  data_crmn <- stackdata(testobj$crmn_adj, link1, pheno, batch)




  pdf(paste(directory, "Diagnostics", Sys.Date(), ".pdf", sep = ""))
  par(mfrow = c(2, 2))
  xlab1 <- link1
  xlab2 <- pheno
  xlab3 <- batch

  pcaplots(testobj$log_data, "Non-Normalized", link1, batch, pheno)
  temp2 <- data_r
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
          ylab = "Abundance", ylim = ylim2)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
          ylab = "Abundance", ylim = ylim2)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
          ylab = "Abundance", ylim = ylim2)

  pcaplots(as.data.frame(testobj$log_data_combat), "Raw with Combat", link1, batch, pheno)
  temp2 <- data_rc
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
          ylab = "Abundance", ylim = ylim2)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
          ylab = "Abundance", ylim = ylim2)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
          ylab = "Abundance", ylim = ylim2)

  pcaplots(testobj$log_quant, "Quantile Normalization", link1, batch, pheno)
  temp2 <- data_n
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
          ylab = "Abundance", ylim = ylim2)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
          ylab = "Abundance", ylim = ylim2)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
          ylab = "Abundance", ylim = ylim2)

  pcaplots(as.data.frame(testobj$log_quant_combat), "Quantile with Combat", 
           link1, batch, pheno)
  temp2 <- data_nc
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
          ylab = "Abundance", ylim = ylim2)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
          ylab = "Abundance", ylim = ylim2)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
          ylab = "Abundance", ylim = ylim2)
  pcaplots(testobj$med_adj, "Median", link1, batch, pheno)
  temp2 <- data_med
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
          ylab = "Abundance", ylim = ylim1)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
          ylab = "Abundance", ylim = ylim1)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
          ylab = "Abundance", ylim = ylim1)
  pcaplots(as.data.frame(testobj$med_combat), "Median with Combat", 
           link1, batch, pheno)
  temp2 <- data_medc
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
          ylab = "Abundance", ylim = ylim1)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
          ylab = "Abundance", ylim = ylim1)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
          ylab = "Abundance", ylim = ylim1)
  pcaplots(as.data.frame(testobj$sva_adj), "SVA Adjusted", 
           link1, batch, pheno)
  temp2 <- data_sva
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
          ylab = "Abundance", ylim = ylim2)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
          ylab = "Abundance", ylim = ylim2)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
          ylab = "Abundance", ylim = ylim2)
  pcaplots(as.data.frame(testobj$ruv_adj), "RUV Adjusted", 
           link1, batch, pheno)
  temp2 <- data_ruv
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
          ylab = "Abundance", ylim = ylim2)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
          ylab = "Abundance", ylim = ylim2)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
          ylab = "Abundance", ylim = ylim2)
  pcaplots(testobj$crmn_adj, "CRMN", link1, batch, pheno)
  temp2 <- data_crmn
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
          ylab = "Abundance", ylim = ylim3)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
          ylab = "Abundance", ylim = ylim3)
  boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
          ylab = "Abundance", ylim = ylim3)
  dev.off()
}


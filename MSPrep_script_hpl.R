setwd("C:/Users/Harrison/Google Drive/Kechris Lab/")
library(crmn)
library(affy)
library(preprocessCore)
library(sva)
library(psych)
library(Hmisc)
library(limma)
library(pcaMethods)
library(VIM)
source("Code/MSPrep_1.1_functions.R")

directory = "C:/Users/Harrison/Google Drive/Kechris Lab"
clinical.file = "./Data/Lipids_Clinical_29Nov2013.csv"
quantification.file = "/Data/Lipids_Raw_29Jan2014.csv"
link.file = "/Data/Lipids_SubjectLink_29Nov2013.csv"

cvmax = 0.5
missing = 1
linktxt = "LCMS_Run_Id"

test = readdata(directory, clinical.file, quantification.file, link.file, cvmax=cvmax, missing=missing, linktxt=linktxt)
rownames(test$sum_data1) = as.character(read.csv(clinical.file)$SubjectID)

# Read clinical data
clin.dat = read.table(paste(directory, clinical.file, sep=""), sep=",", header=TRUE)
clin.dat$Batch = as.factor(clin.dat$Batch)
link1 = "SubjectID"
pheno = c("Smoke","gender","BMI","race","Age_Enroll","ATS_PackYears","COPDExacTrt5",
          "pctEmph_Slicer","pctGasTrap_Slicer","Exacerbation_Frequency","Chronic_Bronchitis",
          "FEV1pp_utah","FEV1_FVC_utah","finalGold","Gold")
batch = "Batch"
ncont = 10
controls = c()
ncomp = 2

# Read annotation file
ref <- read.delim("Data/Lipids_Annotation_29Jan2014.REF")[,1:2]

ref$Compound <- as.character(ref$Compound)
for(anel in ref$DataLabel[duplicated(ref$DataLabel)]){
  ref$Compound[ref$DataLabel == anel] <-
    paste0(sapply(ref$Compound[ref$DataLabel == anel], as.character), collapse = ";")
}
ref <- unique(ref)
rownames(ref) <- ref$DataLabel

"%ni%" = Negate("%in%")

final_df <- test$sum_data1[, colnames(test$sum_data1) %ni% 
                       c("2378.7566_5.8669825", "2403.7534_5.869907", "1201.8774_5.869595",
                         "1200.8752_5.870417", "1200.8737_5.8697524", "1200.3741_5.8704977",
                         "1593.6711_5.8542724", "1198.372_5.8715267", "1200.3752_5.8703713",
                         "2401.75_5.8646865", "1199.3691_5.8723626", "1594.674_5.853273",
                         "1615.6342_5.869456", "2398.7422_5.874978")]


final_df <- rbind(colnames(final_df), ref[paste0("X", colnames(final_df)), "Compound"], final_df)
colnames(final_df) <- final_df[2, ]
rownames(final_df)[1:2] <- c("mz_rt", "Compound")
final_df = final_df[-1:-2, ]
class(final_df) = "numeric"

source("Code/MSPrep_1_1_functions_filterft_hpl.R")
filters = c(0.8, 0.5, 0.1)
for (filt in filters){
  imp = filterft_hpl(final_df, filt)
  metaf = imp$metaf
  metafin = imp$metafin
  metaimp = imp$metaimp
  
  numnas = c()
  for (j in 1:ncol(metafin)) {
    numnas = c(numnas, sum(is.na(metafin[, j]))/length(metaimp[, j]))
  }
  names(numnas) = colnames(metafin)
  write.csv(numnas, paste0("Data/Machine Learning Input/h131_met_num_nas_filt", filt, ".csv"))
  
  
  withzero = impute.zero(metafin, metaimp)
  minval = impute.min(metafin, metaimp)
  meanval = impute.mean(metafin, metaimp)
  medianval = impute.median(metafin, metaimp)
  bpca = impute.bpca(metafin, metaimp, minval)
  knn = impute.knn(metafin, metaimp)
  
  withzero = normdata(withzero, clin.dat, link1, pheno, batch, ncont=ncont, controls, ncomp, omit_combat = T)
  minval = normdata(minval, clin.dat, link1, pheno, batch, ncont=ncont, controls, ncomp, omit_combat = T)
  meanval = normdata(meanval, clin.dat, link1, pheno, batch, ncont=ncont, controls, ncomp, omit_combat = T)
  medianval = normdata(medianval, clin.dat, link1, pheno, batch, ncont=ncont, controls, ncomp, omit_combat = T)
  knn = normdata(knn, clin.dat, link1, pheno, batch, ncont=ncont, controls, ncomp, omit_combat = T)
  bpca = normdata(bpca, clin.dat, link1, pheno, batch, ncont=ncont, controls, ncomp, omit_combat = T)
  
  imp_norm_dfs = list(withzero=withzero, minval=minval, meanval=meanval, medianval=medianval, bpca=bpca, knn=knn)
  for (i in 1:length(imp_norm_dfs)){
    dfs = list(crmn=imp_norm_dfs[[i]]$crmn_adj, med=imp_norm_dfs[[i]]$med_adj, sva=imp_norm_dfs[[i]]$sva_adj, ruv=imp_norm_dfs[[i]]$ruv_adj)
    for (j in 1:length(dfs)){
      final_df = dfs[[j]]
        
      write.csv(final_df, paste0("Data/Machine Learning Input/h131_met_", names(imp_norm_dfs)[i], "_", names(dfs)[j], '_filt_', filt,  ".csv"), quote = T)
  
    }
  }
}  

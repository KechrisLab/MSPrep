\name{MSProcess-package}
\alias{MSProcess-package}
\alias{MSProcess}
\docType{package}
\title{
Package for summarizing, filtering, imputing, and normalizing metabolomics data.
}
\description{
Package performs summarization of replicates, filtering by frequency, three different options for missing data, and five options for normalizing dataset.  Also generates diagnostic plots of imputations and normalizations.
}
\details{
\tabular{ll}{
Package: \tab MSProcess\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-08-17\cr
License: \tab GPL-3 \cr
}
Package for processing of metabolomic quantification data.  Requires three input datasets.
The quantification dataset must have a row per compound, a column for retention time named rt, a column for mass named mz, and a column for each LC/MS run with column names corresponding to an identifier in the subject link dataset.  Each subject should have been run in triplicates. 
The subject link dataset must provide a subject link from each LC/MS run and the corresponding subject identifier from the clinical data file. 
The clinical dataset must contain the unique subject identifier and any clinically relevant outcomes, phenotypes, batch and/or covariates of interest.
These three datasets are utilized by the first function of the program readdata() which will summarize the quantification data as the mean or median of the observation depending on how many replicates the compound is found in and a specified cutoff level for the coefficient of variation between the replicates.
The filterft function takes the summarized data and filters to compounds that are found in a specified percentage of the subjects.  The filtered data can still contain missing observations.  Three datasets, one with the missing observations, one where missing are estimated by the BPCA algorithm, and one where missing are replaced with 1/2 the minimum observation for that compound, are generated from this function.
The graphimputations function will output a histogram of the densities of each of the three datasets generated in the filterft function.  Can restrict to a certain range to exclude extreme outliers or focus only a specific range of abundances.
The normdata function implements five normalization routines and one explicit batch effect adjustment algorithm. Median, Quantile, CRMN, RUV, and SVA are the five normalization alorithms.  Median, quantile, and non-normalized data all have the Combat batch removal algorithm performed.  CRMN and RUV require controls which can either be specified or data driven controls are selected.  The current implementation of the RUV algorithm selects the number of factors based on the number of significant factors found in the SVA algorithm.
The gendiagnostics function creates four diagnostic graphs for each normalization method. A PCA plot with number indicating batches and colors indicating primary categorical phenotype.  Boxplot of abundances by subject, batch, and phenotype are also included.


}
\author{
Grant Hughes

Maintainer: Grant Hughes <dydxhughes@gmail.com>
}
\references{
Bolstad, B.M.et al.(2003) A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. Bioinformatics, 19, 185-193
DeLivera, A.M.et al.(2012) Normalizing and Integrating Metabolomic Data. Anal. Chem, 84, 10768-10776.
Gagnon-Bartsh, J.A.et al.(2012) Using control genes to correct for unwanted variation in microarray data. Biostatistics, 13, 539-552.
Johnson, W.E.et al.(2007) Adjusting batch effects in microarray expression data using Empirical Bayes methods. Biostatistics, 8, 118-127.
Leek, J.T.et al.(2007) Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis. PLoS Genetics, 3(9), e161.
Oba, S.et al.(2003) A Bayesian missing value estimation for gene expression profile data. Bioinformatics, 19, 2088-2096
Redestig, H.et al.(2009) Compensation for Systematic Cross-Contribution Improves Normalization of Mass Spectrometry Based Metabolomics Data. Anal. Chem., 81, 7974-7980.
Stacklies, W.et al.(2007) pcaMethods—a bioconductor package providing PCA methods for incomplete data. Bioinformatics, 23, 1164-1167.
Wang, W.et al.(2003) Quantification of Proteins and Metabolites by Mass Spectrometry without Isotopic Labeling or Spiked Standards. Anal. Chem., 75, 4818-4826.
}
\keyword{ package }
\seealso{
}
\examples{
library(crmn)
library(preprocessCore)
library(sva)
library(psych)
library(Hmisc)
library(limma)
library(pcaMethods)
library(multcomp)


### Specify primary directory for input dataset
directory<-c("/home/grant/MSProcess/data/")

### Specify location of data files
clinicalfile<-c("Clinical.csv")
quantificationfile<-c("Quantification.csv")
linkfile<-c("SubjectLinks.csv")

### Set variables for program
cvmax<-0.5
missing<-1
linktxt<-"LCMS_Run_ID"

test<-readdata(directory, clinicalfile, quantificationfile, linkfile, cvmax=0.50, missing=1, linktxt)

test2<-filterft(test$sum_data1,0.80)

directory<-"/home/grant/"
minval<-test2$minval
withzero<-test2$withzero
bpca<-test2$bpca

graphimputations(directory, minval, withzero, bpca, meanval1=0, meanval2=200000, xmax1=400000, ymax1=800, xmax2=20, ymax2=600, xmax3=20, ymax3=175, nbreaks=200)

metafin<-test2$bpca
clindat<-test$clinical
link1<-"SubjectID"
pheno<-"Spike"
batch<-"Operator"
ncont<-10
controls<-c()
ncomp<-2

test3<-normdata(metafin, clindat, link1, pheno, batch, ncont=10, controls, ncomp)

testobj<-test3
clindat<-test$clinical
link1<-"SubjectID"
pheno<-"Spike"
batch<-"Operator"
directory<-"/home/grant/"
ylim2<-c(10,28)
### For median
ylim1<-c(-15,15)
### for crmn
ylim3<-c(18,37)

diagnosticgen(testobj, clindat, link1, batch, pheno, directory, ylim1, ylim2, ylim3)

}

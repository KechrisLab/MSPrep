
### Specify primary directory
directory<-c("C:/Users/Hazy/Dropbox/Metabolomics/Programs/MSProcess/data/")

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

directory<-"C:/Users/Hazy/Dropbox/Metabolomics/Programs/MSProcess/"
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
directory<-"C:/Users/Hazy/Dropbox/Metabolomics/Programs/MSProcess/"
ylim2<-c(10,28)
### For median
ylim1<-c(-15,15)
### for crmn
ylim3<-c(18,37)

diagnosticgen(testobj, clindat, link1, batch, pheno, directory, ylim1, ylim2, ylim3)

save(test, file=paste(directory,"test.Rdata",sep=""))
save(test2, file=paste(directory,"test2.Rdata",sep=""))
save(test3, file=paste(directory,"test3.Rdata",sep=""))

load("C:/Users/Hazy/Dropbox/Metabolomics/Programs/MSProcess/data/test.Rdata")
load("C:/Users/Hazy/Dropbox/Metabolomics/Programs/MSProcess/data/test2.Rdata")
load("C:/Users/Hazy/Dropbox/Metabolomics/Programs/MSProcess/data/test3.Rdata")
### Load required packages
### source("http://bioconductor.org/biocLite.R")
### biocLite("pcaMethods")

library(crmn)
library(preprocessCore)
library(sva)
library(psych)
library(Hmisc)
library(limma)
library(pcaMethods)
library(multcomp)





### Function readdata reads in all the datafiles and summarizes the replicates


readdata<-function(directory, clinicalfile, quantificationfile, linkfile, cvmax=0.50, missing=1){

### Read in data
quant<-read.table(paste(directory,quantificationfile,sep=""), header=TRUE, sep=",")
link<-read.table(paste(directory,linkfile,sep=""), header=TRUE, sep=",")
clinical<-read.table(paste(directory,clinicalfile,sep=""), header=TRUE, sep=",")

compounds<-nrow(quant)
subjects<-nrow(clinical)

metaf<-matrix(NA, ncol=compounds,nrow=subjects)
temp3<-matrix(NA, ncol=3, nrow=subjects*compounds)
temp4<-matrix(NA, ncol=2, nrow=subjects*compounds)
k<-1

rownames(quant)<-paste(quant$mz,"_",quant$rt,sep="")
metab<-subset(quant, select=-c(mz,rt))

for (i in 1:subjects){
	temp<-matrix(c(metab[,i*3-2],metab[,i*3-1],metab[,i*3]),ncol=3)
	temp2<-matrix(NA,ncol=4,nrow=compounds)
	for (j in 1:compounds){
			if (temp[j,1]==missing){temp[j,1]<-NA}
			if (temp[j,2]==missing){temp[j,2]<-NA}
			if (temp[j,3]==missing){temp[j,3]<-NA}
			temp2[j,1]<-mean(temp[j,],na.rm=TRUE)
			temp2[j,2]<-sd(temp[j,],na.rm=TRUE)
			temp2[j,3]<-temp2[j,2]/temp2[j,1]
			temp2[j,4]<-length(na.omit(temp[j,]))
			if (temp2[j,1]=="NaN"){temp2[j,1]<-0}
			### Compound must be present in at least two replicates
			### Compound must have CV less than CVMax
			if (temp2[j,4] < 2 || temp2[j,3] > cvmax){
				temp2[j,1]<-0
				if (temp2[j,4] == 3) {
					temp2[j,1]<-median(temp[j,],na.rm=TRUE)
					temp3[k,]<-temp[j,]
					temp4[k,1]<-rownames(metab)[j]
					temp4[k,2]<-colnames(metab)[i]
					k<-k+1
					}
				}
		}
		metaf[i,]<-t(temp2[,1])
		}
	med_comp<-temp3[1:eval((compounds*subjects)-summary(temp3[,1])[7]),]
	med_comp2<-temp4[1:eval((compounds*subjects)-summary(temp3[,1])[7]),]
	med_comp<-cbind(med_comp,med_comp2)

	colnames(metaf)<-rownames(metab)

	linktxt<-"LCMS_Run_ID"
	link2<-matrix(colnames(metab),ncol=1)
	colnames(link2)<-"LCMS_Run_ID"
	link3<-unique(merge(link2, link, by.x="LCMS_Run_ID", by.y=as.character(linktxt))[2])
	rownames(metaf)<-link3[,1]

	list(sum_data1=metaf, clinical=clinical, medians=med_comp)
}






### Filterft to only include compounds that are found in specified percentage of subjects and perform imputation of missing data

filterft<-function(metaf,filterpercent=0.50) {

	count<-matrix(NA,nrow=ncol(metaf),ncol=1)
	toss<-matrix(1,nrow=eval(ncol(metaf)),ncol=1)
	for (j in  1:eval(ncol(metaf))){
		k<-0
		for (i in 1:dim(metaf)[1]){
			if (metaf[i,j] == 0.0) {k<-k+1}
			}
		count[j,1]<-k
		if (k >= eval(dim(metaf)[1]*(1-filterpercent))) {toss[j,1]<-0}
		}
	colnames(toss)<-"toss"
	tests<-cbind(toss,t(metaf))		
	metafin<-t(subset(tests,toss != 0)[,-1])

	### Create new dataset with 0.00 replaced with NA.
	metaimp<-matrix(NA,nrow=nrow(metafin),ncol=ncol(metafin))
	for (i in 1:nrow(metafin)){
		for (j in 1:ncol(metafin)){
			if (metafin[i,j]==0.0){
				metaimp[i,j]<-NA
				}
			else {
				metaimp[i,j]<-metafin[i,j]
				}
			}	
		}

	### Create present/absent filter
	present<-matrix(1,nrow=nrow(metaf), ncol=ncol(metaf))
	colnames(present)<-colnames(metaf)
	rownames(present)<-rownames(metaf)
	for (j in  1:eval(ncol(metaf))){
		for (i in 1:dim(metaf)[1]){
			if (metaf[i,j] == 0.0) {present[i,j]<-0}
			}
		}

	### Create dataset with missing values replaced with 1/2 minimum value
	minval<-metafin
	for (i in 1:nrow(metaf)){
		for (j in 1:ncol(metafin)){
			if (metafin[i,j]==0.0){
				minval[i,j]<-summary(metaimp[,j])[1]/2
				}
			else {
				minval[i,j]<-metafin[i,j]
				}
			}	
		}


	metabpca <- pca(metaimp, nPcs=3, method="bpca")
	# metanlpca <- pca(metaimp, nPcs=3, method="nlpca")

	bpca<-completeObs(metabpca)
	colnames(bpca)<-colnames(metafin)
	rownames(bpca)<-rownames(metafin)
	### Replace compounds with negative abundance with minval
	coldrop<-1
	for (i in 1:nrow(bpca)){
		for (j in 1:ncol(bpca)) {
		if (bpca[i,j]<0) {
			bpca[i,j]<-minval[i,j]
		}}}



	list(minval=as.data.frame(minval),withzero=as.data.frame(metafin),bpca=as.data.frame(bpca),count=as.data.frame(count))
	}






### Create pdf of histogram of distribution of values

graphimputations<-function(directory, minval, withzero, bpca,fname, meanval1=0, meanval2=1000000, xmax1=500000, ymax1=25000, xmax2=20, ymax2=1000, xmax3=20, ymax3=1500, nbreaks=500) {

compounds<-ncol(bpca)
subjects<-nrow(bpca)

### Function for stacking data
stackdata<-function(dset){
	temp<-matrix(NA, ncol=5,nrow=ncol(dset)*nrow(dset))
	comps<-ncol(dset)
	subjects<-nrow(dset)
	for (j in 0:eval(comps-1)){
		for (i in 1:subjects) {
			temp[(subjects*j+i),1]<-as.character(dset[i,j+1])
		}
	}
	return(temp)
	}

keeper<-matrix(0,nrow=compounds,ncol=1)
log_zer<-matrix(NA,nrow=subjects, ncol=compounds)
log_min<-matrix(NA,nrow=subjects, ncol=compounds)
log_bpca<-matrix(NA,nrow=subjects, ncol=compounds)
for (i in 1:compounds) {
	if (mean(withzero[,i]) >= meanval1 && mean(withzero[,i]) <= meanval2) {  
		keeper[i,1]<-i
		}
	for (j in 1:subjects) {
		if (withzero[j,i]==0) {
			log_zer[j,i] <- log(1,2)
		} else {
			log_zer[j,i] <- log(withzero[j,i],2)
		}
		log_min[j,i] <- log(minval[j,i],2)
		log_bpca[j,i] <- log(bpca[j,i],2)
		}
	}

keep<-subset(keeper, keeper[,1]>0)

### Remove observations that were not within the range
minst<-stackdata(minval[,keep])
zerst<-stackdata(withzero[,keep])
bpcast<-stackdata(bpca[,keep])

### Stack remaining data
lminst<-stackdata(log_min[,keep])
lzerst<-stackdata(log_zer[,keep])
lbpcast<-stackdata(log_bpca[,keep])


pdf(paste(directory,"ImpMethods_CompsBetween_",meanval1,"-",as.character(meanval2),"_bk",nbreaks,"_",Sys.Date(),".pdf",sep=""))
par(mfrow = c(2,1))
ylim1<-c(0,ymax1)
xlim1<-c(0,xmax1)
ylim2<-c(0,ymax2)
xlim2<-c(0,xmax2)
ylim3<-c(0,ymax3)
xlim3<-c(0,xmax3)


x<-(as.numeric(as.character(zerst[,1])))
hist(x, breaks=nbreaks, col="red", xlab="With Zeros", main="Histogram of Abundances with Zeros", ylim=ylim1, xlim=xlim1) 

x<-(as.numeric(as.character(lzerst[,1])))
hist(x, breaks=nbreaks, col="red", xlab="With Zeros", main="Histogram of Log2 Abundances with Zeros", ylim=ylim2, xlim=xlim2) 

x<-(as.numeric(as.character(minst[,1])))
hist(x, breaks=nbreaks, col="red", xlab="With Min Val", main="Histogram of Abundances with 1/2 Minimum Value Replacement", ylim=ylim1, xlim=xlim1) 

x<-(as.numeric(as.character(lminst[,1])))
hist(x, breaks=nbreaks, col="red", xlab="With Min Val", main="Histogram of Log2 Abundances with 1/2 Minimum Value Replacement", ylim=ylim3, xlim=xlim2) 

x<-(as.numeric(as.character(bpcast[,1])))
hist(x, breaks=nbreaks, col="red", xlab="BPCA", main="Histogram of Abundances with BPCA Imputation", ylim=ylim1, xlim=xlim1) 

x<-(as.numeric(as.character(lbpcast[,1])))
hist(x, breaks=nbreaks, col="red", xlab="BPCA", main="Histogram of Log2 Abundances with BPCA Imputation", ylim=ylim3, xlim=xlim2) 

dev.off()

}












normdata<-function(metafin, clindat, link1, pheno, batch, ncont=10, controls=c(), ncomp=2) {

compounds<-ncol(metafin)
subjects<-nrow(metafin)

### Perform quantile normalization
metan<-t(normalize.quantiles(t(metafin)))
colnames(metan)<-colnames(metafin)
rownames(metan)<-rownames(metafin)
final<-metafin
finaln<-metan

convert<-function(dset){
	s_dset<-as.data.frame(matrix(NA, ncol=ncol(dset), nrow=nrow(dset)))
	for (j in 1:ncol(dset)) {
		s_dset[,j]<-log2(as.numeric(dset[,j]))
		}
	colnames(s_dset)<-colnames(dset)
	rownames(s_dset)<-rownames(dset)
	return(s_dset)
	}
log_final<-convert(final)
log_finaln<-convert(finaln)

combat<-function(dset, pheno, batch) {
	compounds<-ncol(dset)
	temp<-dset
	temp$mid<-rownames(temp)
	test<-merge(temp,clindat,by.x="mid",by.y=as.character(link1))
	### Create dataset 1
	d1<-t(dset)
	colnames(d1)<-dset[,1]

	### Create dataset 2
	d2<-subset(test,select=pheno)

	### Create dataset 3
	d3<-subset(test, select=batch)

	count<-matrix(0,nrow=compounds,ncol=1)
	d1mod<-matrix(0,nrow=compounds,ncol=subjects)
	for (j in 1:compounds) {
			d1mod[j,]<-log2(as.numeric(d1[j,]))
		}
	comadj<-t(ComBat(d1mod, d2, d3))
	colnames(comadj)<-colnames(dset)
	rownames(comadj)<-rownames(dset)

	return(comadj)
	}

final_rc<-combat(as.data.frame(final), pheno, batch)
final_com<-combat(as.data.frame(finaln), pheno, batch)


### Create surrogate variable estimates using SVA method
svaestimate<-function(dset, pheno, batch) {
	### Create dataset 1
	d1<-t(dset[,1:compounds])
	colnames(d1)<-rownames(dset)
	rownames(d1)<-colnames(dset)

	### Create dataset 2
	dset$mid<-rownames(dset)
	test<-merge(dset,clindat,by.x="mid",by.y=as.character(link1))
	d2<-subset(test, select=c(batch))

	### Create dataset 3
	d3<-subset(test, select=c(pheno))
	colnames(d3)[1]<-"pheno"
	d3m<-model.matrix(~as.factor(pheno),data=d3)

	count<-matrix(0,nrow=compounds,ncol=1)
	d1mod<-d1
	svaadj<-sva(d1mod,d3m,method="irw")
	return(as.matrix(svaadj$sv,ncol=svaadj$n.sv))
	}

sva<-svaestimate(log_final, pheno, batch)

model.func <- c()
for (i in 1:ncol(sva)) {
	colnames(sva)[i]<-paste("f",i,sep="")
	model.func[i]<-paste("test$f",i,sep="")
	}
sva_funct <- paste (model.func , collapse = " + ")

### Estimate surrogates using RUV method
ruv<-function(dset,k,ctl, pheno){
	d1<-t(dset[,1:compounds])
	count<-matrix(0,nrow=compounds,ncol=1)
	d1mod<-d1

	Y<-d1mod
	dset$mid<-rownames(dset)
	test<-merge(dset,clindat,by.x="mid",by.y=link1)
	X<-subset(test, select=pheno)
	Z=matrix(rep(1, ncol(Y)))

	# Project onto the orthogonal complement of Z
	RZY = Y - Y%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)

	# Perform SVD
	W = svd(RZY[ctl,])$v

	# Keep the first k factors
	W = W[,1:k]
	return(as.matrix(W,ncol=k))
}

### Control Genes
### Evaluation of potential control genes.  Top genes with the smallest CV and found in all subjects
dset<-final
control<-matrix(NA,nrow=compounds,ncol=6)
for (j in 1:compounds) {
	control[j, 1]<-colnames(dset)[j]
	control[j,2]<-j
	control[j,3]<-mean(dset[,j], na.rm=TRUE)
	control[j,4]<-sd(dset[,j], na.rm=TRUE)
	control[j,5]<-sd(dset[,j], na.rm=TRUE)/mean(dset[,j], na.rm=TRUE)
	control[j,6]<-sum(dset[,j]==0, na.rm=TRUE)
	}
control2<-control[order(control[,6],control[,5]), ]

### Select top ncont with complete data and minimum variation to use as controls
ctl<-as.numeric(control2[,2][1:ncont])
if (length(controls) > 0) {
	ctl<-controls
	}

contout<-cbind(ctl,colnames(log_final)[ctl])

ruv_raw<-ruv(log_final,ncol(sva), ctl, pheno)

model.func <- c()
for (i in 1:ncol(ruv_raw)) {
	colnames(ruv_raw)[i]<-paste("f",i,sep="")
	model.func[i]<-paste("as.numeric(as.character(test$f",i,"))",sep="")
	}
ruv_funct <- paste (model.func , collapse = " + ")

### Adjust RUV/SVA estimates
# adjust_factors<-function(dset, factors, funct) {

dset<-log_final
factors<-ruv_raw
funct<-ruv_funct

	out<-matrix(NA,ncol=ncol(dset),nrow=nrow(dset))
	rownames(out)<-rownames(dset)
	colnames(out)<-colnames(dset)
	dset$mid<-rownames(dset)
	test<-cbind(merge(dset,clindat,by.x="mid",by.y=link1),factors)
	for (k in 1:compounds) {
		funct2 <- c(paste("test[,k+1] ~ ",funct,sep=""))
		fit <- lm(funct2)
		out[,k]<-fit$fitted.values
		}

ruv_adj<-out

dset<-log_final
factors<-sva
funct<-sva_funct

	out<-matrix(NA,ncol=ncol(dset),nrow=nrow(dset))
	rownames(out)<-rownames(dset)
	colnames(out)<-colnames(dset)
	dset$mid<-rownames(dset)
	test<-cbind(merge(dset,clindat,by.x="mid",by.y=link1),factors)
	for (k in 1:compounds) {
		funct2 <- c(paste("test[,k+1] ~ ",funct,sep=""))
		fit <- lm(funct2)
		out[,k]<-fit$fitted.values
		}

sva_adj<-out

### CRMN, Median

final_shift<-final

Y<-t(final_shift)
Ya<-t(final_shift)[ctl,]
Ys<-t(final_shift)[-ctl,]


j<-1
isIS<-rep(FALSE, compounds)
ctlo<-ctl[order(ctl)]
for (i in 1:compounds) {
	if(j <= 10) {
	if (ctlo[j] == i) {
		isIS[i]<-TRUE
		j<-j+1
	} else {
		isIS[i]<-FALSE
		}
	}}

dset<-final_shift
dset$mid<-rownames(dset)
test<-merge(dset,clindat,by.x="mid",by.y=link1)
X<-subset(test, select=c(pheno))
colnames(X)[1]<-"pheno"
G <- model.matrix(~-1 + as.factor(X$pheno))

normed.crmn <- normalize(Y, "crmn", factors=G, standards = isIS, ncomp = ncomp)
lnormed.crmn<-convert(normed.crmn)
normed.med <- normalize(Y, "median", factors=G, standards = isIS)
lnormed.med<-convert(normed.med)

final_crmn<-as.data.frame(t(lnormed.crmn))
final_med<-as.data.frame(t(lnormed.med))

med_com<-combat(as.data.frame(t(normed.med)), pheno, batch)


list(log_data=log_final, log_data_combat=final_rc, log_quant=log_finaln, log_quant_combat=final_com, med_adj=final_med, med_combat=med_com,
 sva_factors=sva, sva_adj=sva_adj, ruv_factors=ruv_raw, ruv_adj=ruv_adj, crmn_adj=final_crmn, controls=contout)

}

















diagnosticgen<-function(testobj, link1, batch, pheno, directory, ylim1, ylim2


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


pcaplots<-function(dset,tlab, link1, batch, pheno) {
	dset$mid<-rownames(dset)
	test<-merge(dset,clindat,by.x="mid",by.y=as.character(link1))
	pca1<-prcomp(test[,2:ncol(dset)])
	title=c(paste("PCA for ",tlab,sep=""))
	title2=c("Color Indicates Phenotype, Number Indicates Batch")
	xlab1=paste("PCA 1, % Var = ",summary(pca1)$importance[2,1],sep="")
	ylab1=paste("PCA 2, % Var = ",summary(pca1)$importance[2,2],sep="")
	temp2<-subset(test, select=c(batch,pheno))
	colnames(temp2)<-c("batch","pheno")
	plot(pca1$x[,1], pca1$x[,2], 
		main=title, 
		sub=title2,
		xlab=xlab1,
		ylab=ylab1,
		pch=sapply(test$batch, pca.point.pch.batch),
		col=sapply(test$pheno, pca.point.color.change))
	return(pca1)
}

stackdata<-function(dset,link1,pheno,batch){
	temp<-matrix(NA, ncol=5,nrow=ncol(dset)*nrow(dset))
	comps<-ncol(dset)
	subjects<-nrow(dset)
	dset$mid<-rownames(dset)
	dset<-merge(dset,clindat,by.x="mid",by.y=link1)
	for (j in 0:eval(comps-1)){
		for (i in 1:subjects) {
			temp[(subjects*j+i),1]<-as.character(dset[i,j+2])
			temp[(subjects*j+i),2]<-as.character(dset$mid[i])
			temp[(subjects*j+i),3]<-as.character(subset(dset, select=pheno)[i,1])
			temp[(subjects*j+i),4]<-as.character(subset(dset, select=batch)[i,1])
		}
	}
	return(temp)
	}

data_r<-stackdata(testobj$log_data,link1,pheno,batch)
data_n<-stackdata(testobj$log_quant,link1,pheno,batch)
data_rc<-stackdata(as.data.frame(testobj$log_data_combat),link1,pheno,batch)
data_nc<-stackdata(as.data.frame(testobj$log_quant_combat),link1,pheno,batch)
data_med<-stackdata(testobj$med,link1,pheno,batch)
data_medc<-stackdata(as.data.frame(testobj$med_combat),link1,pheno,batch)
data_ruv<-stackdata(as.data.frame(testobj$ruv_adj),link1,pheno,batch)
data_sva<-stackdata(as.data.frame(testobj$sva_adj),link1,pheno,batch)
data_crmn<-stackdata(testobj$crmn_adj,link1,pheno,batch)



pdf(paste(directory,"Diagnostics",Sys.Date(),".pdf",sep=""))
par(mfrow = c(2,2))
xlab1<-link1
xlab2<-pheno
xlab3<-batch

pcaplots(testobj$log_data, "Non-Normalized", link1, batch, pheno)
temp2<-data_r
boxplot(as.numeric(temp2[,1])~temp2[,2],  main="", xlab=xlab1, ylab="Abundance",ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,3],  main="", xlab=xlab2, ylab="Abundance", ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,4],  main="", xlab=xlab3, ylab="Abundance", ylim=ylim2)

pcaplots(as.data.frame(testobj$log_data_combat), "Raw with Combat", link1, batch, pheno)
temp2<-data_rc
boxplot(as.numeric(temp2[,1])~temp2[,2],  main="", xlab=xlab1, ylab="Abundance",ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,3],  main="", xlab=xlab2, ylab="Abundance", ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,4],  main="", xlab=xlab3, ylab="Abundance", ylim=ylim2)

pcaplots(testobj$log_quant, "Quantile Normalization", link1, batch, pheno)
temp2<-data_n
boxplot(as.numeric(temp2[,1])~temp2[,2],  main="", xlab=xlab1, ylab="Abundance",ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,3],  main="", xlab=xlab2, ylab="Abundance", ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,4],  main="", xlab=xlab3, ylab="Abundance", ylim=ylim2)

pcaplots(as.data.frame(testobj$log_quant_combat), "Quantile with Combat", link1, batch, pheno)
temp2<-data_nc
boxplot(as.numeric(temp2[,1])~temp2[,2],  main="", xlab=xlab1, ylab="Abundance",ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,3],  main="", xlab=xlab2, ylab="Abundance", ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,4],  main="", xlab=xlab3, ylab="Abundance", ylim=ylim2)

pcaplots(testobj$med_adj, "Median", link1, batch, pheno)
temp2<-data_med
boxplot(as.numeric(temp2[,1])~temp2[,2],  main="", xlab=xlab1, ylab="Abundance",ylim=ylim1)
boxplot(as.numeric(temp2[,1])~temp2[,3],  main="", xlab=xlab2, ylab="Abundance", ylim=ylim1)
boxplot(as.numeric(temp2[,1])~temp2[,4],  main="", xlab=xlab3, ylab="Abundance", ylim=ylim1)

pcaplots(as.data.frame(testobj$med_combat), "Median with Combat", link1, batch, pheno)
temp2<-data_medc
boxplot(as.numeric(temp2[,1])~temp2[,2],  main="", xlab=xlab1, ylab="Abundance",ylim=ylim1)
boxplot(as.numeric(temp2[,1])~temp2[,3],  main="", xlab=xlab2, ylab="Abundance", ylim=ylim1)
boxplot(as.numeric(temp2[,1])~temp2[,4],  main="", xlab=xlab3, ylab="Abundance", ylim=ylim1)


pcaplots(as.data.frame(testobj$sva_adj), "SVA Adjusted", link1, batch, pheno)
temp2<-data_sva
boxplot(as.numeric(temp2[,1])~temp2[,2],  main="", xlab=xlab1, ylab="Abundance",ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,3],  main="", xlab=xlab2, ylab="Abundance", ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,4],  main="", xlab=xlab3, ylab="Abundance", ylim=ylim2)

pcaplots(as.data.frame(testobj$ruv_adj), "RUV Adjusted", link1, batch, pheno)
temp2<-data_ruv
boxplot(as.numeric(temp2[,1])~temp2[,2],  main="", xlab=xlab1, ylab="Abundance",ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,3],  main="", xlab=xlab2, ylab="Abundance", ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,4],  main="", xlab=xlab3, ylab="Abundance", ylim=ylim2)

pcaplots(testobj$crmn_adj, "CRMN", link1, batch, pheno)
temp2<-data_crmn
boxplot(as.numeric(temp2[,1])~temp2[,2],  main="", xlab=xlab1, ylab="Abundance",ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,3],  main="", xlab=xlab2, ylab="Abundance", ylim=ylim2)
boxplot(as.numeric(temp2[,1])~temp2[,4],  main="", xlab=xlab3, ylab="Abundance", ylim=ylim2)

dev.off()
}

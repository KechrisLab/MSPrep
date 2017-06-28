normdata <-
function(metafin, clindat, link1, pheno, batch, ncont=10, controls=c(), ncomp=2) {

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

for (i in 1:ncol(sva)) {
colnames(sva)[i]<-paste("f",i,sep="")
}

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

for (i in 1:ncol(ruv_raw)) {
colnames(ruv_raw)[i]<-paste("f",i,sep="")
}

### Adjust RUV/SVA estimates

genadj<-function(dset, factors) {
out <- sapply(1:compounds, function(j) lm(dset[,j] ~ as.matrix(factors,ncol=1))$fitted.values)
colnames(out)<-colnames(dset)
rownames(out)<-rownames(dset)
return(out)
}

ruv_adj<-genadj(log_final, ruv_raw)
sva_adj<-genadj(log_final, sva)

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

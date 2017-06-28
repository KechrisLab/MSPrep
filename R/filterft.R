filterft <- function(metaf, filterpercent = 0.50) {

    count<-matrix(NA,nrow=ncol(metaf),ncol=1)
    rownames(count)<-colnames(metaf)
    toss<-matrix(1,nrow=eval(ncol(metaf)),ncol=1)

    for (j in  1:eval(ncol(metaf))){
      k<-0
      for (i in 1:dim(metaf)[1]){
        if (metaf[i,j] == 0.0) {k<-k+1}
      }
      count[j,1]<- 1 - (k / nrow(metaf))
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
    colnames(bpca) <- colnames(metafin)
    rownames(bpca) <- rownames(metafin)
    ### Replace compounds with negative abundance with minval
    coldrop <- 1
    for (i in 1:nrow(bpca)) {
      for (j in 1:ncol(bpca)) {
        if (bpca[i, j] < 0) {
          bpca[i, j] <- minval[i, j]
    }}}



    list(minval   = as.data.frame(minval),
         withzero = as.data.frame(metafin),
         bpca     = as.data.frame(bpca),
         count    = as.data.frame(count))

  }

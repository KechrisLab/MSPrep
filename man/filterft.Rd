\name{filterft}
\alias{filterft}
\title{
filterft() function filters and imputes dataset
}
\description{
Filters compounds to those found in specified percentage of subjects and performs data imputation.
}
\usage{
filterft(metaf, filterpercent = 0.5)
}
\arguments{
  \item{metaf}{
Summarized dataset output as sum_data1 from readdata() function
}
  \item{filterpercent}{
Percent to filter the data
}
}
\details{
}
\value{
  \item{minval}{Filtered dataset with missing values replaced by 1/2 minimum observed value for that compound}
  \item{bpca}{Filtered dataset with missing values imputed by a Bayesian PCA from PCAMethods package}
  \item{withzero}{Filtered dataset with no imputation}
  \item{count}{List of all compounds and the percent present for each compound}
}
\references{
Oba, S.et al.(2003) A Bayesian missing value estimation for gene expression profile data. Bioinformatics, 19, 2088-2096
Stacklies, W.et al.(2007) pcaMethods—a bioconductor package providing PCA methods for incomplete data. Bioinformatics, 23, 1164-1167.
}
\author{
Grant Hughes
}
\note{

}


\seealso{
}
\examples{
### Load object generated from readdata() function
load("C:/Users/Hazy/Dropbox/Metabolomics/Programs/MSProcess/data/test.Rdata")

test2<-filterft(test$sum_data1,0.80)


## The function is currently defined as
function (metaf, filterpercent = 0.5) 
{
    count <- matrix(NA, nrow = ncol(metaf), ncol = 1)
    rownames(count)<-colnames(metaf)
    toss <- matrix(1, nrow = eval(ncol(metaf)), ncol = 1)
    for (j in 1:eval(ncol(metaf))) {
        k <- 0
        for (i in 1:dim(metaf)[1]) {
            if (metaf[i, j] == 0) {
                k <- k + 1
            }
        }
        count[j, 1] <- 1-(k / nrow(metaf))
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
            if (metaf[i, j] == 0) {
                present[i, j] <- 0
            }
        }
    }
    minval <- metafin
    for (i in 1:nrow(metaf)) {
        for (j in 1:ncol(metafin)) {
            if (metafin[i, j] == 0) {
                minval[i, j] <- summary(metaimp[, j])[1]/2
            }
            else {
                minval[i, j] <- metafin[i, j]
            }
        }
    }
    metabpca <- pca(metaimp, nPcs = 3, method = "bpca")
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

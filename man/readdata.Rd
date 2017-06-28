\name{readdata}
\alias{readdata}
\title{
ReadData() - Function for importing raw data and summarizing replicates
}
\description{
Function reads in raw data files and summarizes technical replicates as the mean of observations for compounds found in 2 or 3 replicates and with coefficient of variation below specied level, or median for those found in 3 replicates but excess CV.
}
\usage{
readdata(directory, clinicalfile, quantificationfile, linkfile, cvmax = 0.5, missing = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{directory}{
Directory where the required three data files are located.
}
  \item{clinicalfile}{
Name of the clinical data file.
}
  \item{quantificationfile}{
Name of the data file containing the quantification data.
}
  \item{linkfile}{
Name of Subject link file.
}
  \item{cvmax}{
Acceptable level of coefficient of variation between replicates.
}
  \item{missing}{
Value of missing data in the quantification data file.
}
  \item{linktxt}{
Column name for Run ID field in Subject Link dataset
}
}
\details{
}
\value{
  \item{sum_data1}{Matrix of summarized replicates, one obs per subject per compound}
  \item{clinical}{Clinical dataset}
  \item{medians}{List of compounds that had excess CV and utilized the median}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Grant Hughes
}
\note{
Must have unique compounds in raw data file.  Raw quantification data should have two columns titled mz and rt that are combined to make the column header. 
}

\seealso{
}
\examples{
### Specify primary directory
directory<-c("/MSProcess/data/")

### Specify location of data files
clinicalfile<-c("Clinical.csv")
quantificationfile<-c("Quantification.csv")
linkfile<-c("SubjectLinks.csv")

### Set variables for program
cvmax<-0.5
missing<-1
linktxt<-"LCMS_Run_ID"

test<-readdata(directory, clinicalfile, quantificationfile, linkfile, cvmax=0.50, missing=1, linktxt)
save(test, file=paste(directory,"test.Rdata",sep=""))

## The function is currently defined as
function (directory, clinicalfile, quantificationfile, linkfile, 
    cvmax = 0.5, missing = 1, linktxt) 
{
    quant <- read.table(paste(directory, quantificationfile, 
        sep = ""), header = TRUE, sep = ",")
    link <- read.table(paste(directory, linkfile, sep = ""), 
        header = TRUE, sep = ",")
    clinical <- read.table(paste(directory, clinicalfile, sep = ""), 
        header = TRUE, sep = ",")
    compounds <- nrow(quant)
    subjects <- nrow(clinical)
    metaf <- matrix(NA, ncol = compounds, nrow = subjects)
    temp3 <- matrix(NA, ncol = 3, nrow = subjects * compounds)
    temp4 <- matrix(NA, ncol = 2, nrow = subjects * compounds)
    k <- 1
    rownames(quant) <- paste(quant$mz, "_", quant$rt, sep = "")
    metab <- subset(quant, select = -c(mz, rt))
    for (i in 1:subjects) {
        temp <- matrix(c(metab[, i * 3 - 2], metab[, i * 3 - 
            1], metab[, i * 3]), ncol = 3)
        temp2 <- matrix(NA, ncol = 4, nrow = compounds)
        for (j in 1:compounds) {
            if (temp[j, 1] == missing) {
                temp[j, 1] <- NA
            }
            if (temp[j, 2] == missing) {
                temp[j, 2] <- NA
            }
            if (temp[j, 3] == missing) {
                temp[j, 3] <- NA
            }
            temp2[j, 1] <- mean(temp[j, ], na.rm = TRUE)
            temp2[j, 2] <- sd(temp[j, ], na.rm = TRUE)
            temp2[j, 3] <- temp2[j, 2]/temp2[j, 1]
            temp2[j, 4] <- length(na.omit(temp[j, ]))
            if (temp2[j, 1] == "NaN") {
                temp2[j, 1] <- 0
            }
            if (temp2[j, 4] < 2 || temp2[j, 3] > cvmax) {
                temp2[j, 1] <- 0
                if (temp2[j, 4] == 3) {
                  temp2[j, 1] <- median(temp[j, ], na.rm = TRUE)
                  temp3[k, ] <- temp[j, ]
                  temp4[k, 1] <- rownames(metab)[j]
                  temp4[k, 2] <- colnames(metab)[i]
                  k <- k + 1
                }
            }
        }
        metaf[i, ] <- t(temp2[, 1])
    }
    med_comp <- temp3[1:eval((compounds * subjects) - summary(temp3[, 
        1])[7]), ]
    med_comp2 <- temp4[1:eval((compounds * subjects) - summary(temp3[, 
        1])[7]), ]
    med_comp <- cbind(med_comp, med_comp2)
    colnames(metaf) <- rownames(metab)

    link2 <- matrix(colnames(metab), ncol = 1)
    colnames(link2) <- "LCMS_Run_ID"
    link3 <- unique(merge(link2, link, by.x = "LCMS_Run_ID", 
        by.y = as.character(linktxt))[2])
    rownames(metaf) <- link3[, 1]
    list(sum_data1 = metaf, clinical = clinical, medians = med_comp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

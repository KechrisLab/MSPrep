
#' Function for importing raw data and summarizing replicates
#'
#' Function reads in raw data files and summarizes technical replicates as the
#' mean of observations for compounds found in 2 or 3 replicates and with
#' coefficient of variation below specified level, or median for those found in
#' 3 replicates but excess CV.
#'
#' @param clinical_data Name of the clinical data file.
#' @param quantification_data Name of the data file containing the
#' quantification data.
#' @param link_data Name of Subject link file.
#' @param cvmax Acceptable level of coefficient of variation between replicates.
#' @param missing Value of missing data in the quantification data file.
#' @param linktxt Column name for Run ID field in Subject Link dataset
#' @return sum_data1 Matrix of summarized replicates, one obs per subject per
#' compound
#' @return clinical Clinical dataset
#' @return medians List of compounds that had excess CV and utilized the median
#' @note Must have unique compounds in raw data file.  Raw quantification data
#' should have two columns titled mz and rt that are combined to make the
#' column header. 
#' @examples
#'
#'   ### Read in data files
#'   clinical <- read.csv("./data-raw/Clinical.csv")
#'   quant    <- read.csv("./data-raw/Quantification.csv")
#'   link     <- read.csv("./data-raw/SubjectLinks.csv")
#'
#'
#'   ### Set variables for program
#'   cvmax   <- 0.5
#'   missing <- 1
#'   linktxt <- "LCMS_Run_ID"
#'
#'   test <- msprep(clinical_data       = clinical,
#'                  quantification_data = quant,
#'                  link_data           = link,
#'                  cvmax               = 0.50,
#'                  missing             = 1,
#'                  linktxt             = linktxt)
#'   #save(test, file = paste0("./data/test.Rdata"))
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @export

msprep <- function(directory,
                   clinicalfile,
                   quantificationfile,
                   linkfile,
                   cvmax   = 0.50,
                   missing = 1,
                   linktxt) {

  ### Read in data
  clinical <- clinical_data
  quant    <- quantification_data
  link     <- link_data

  compounds <- nrow(quant)
  subjects <- nrow(clinical)

  metaf <- matrix(NA, ncol = compounds,nrow = subjects)
  temp3 <- matrix(NA, ncol = 3, nrow = subjects * compounds)
  temp4 <- matrix(NA, ncol = 2, nrow = subjects * compounds)
  k <- 1

  rownames(quant) <- paste(quant$mz, quant$rt, sep = "_")
  metab <- subset(quant, select = -c(mz, rt))

  for (i in 1:subjects) {
    temp  <- matrix(c(metab[, (i * 3 - 2)], metab[, (i * 3 - 1)], metab[, (i * 3)]), ncol = 3)
    temp2 <- matrix(NA, ncol = 4, nrow = compounds)

    for (j in 1:compounds){
      if (temp[j, 1] == missing) temp[j, 1] <- NA 
      if (temp[j, 2] == missing) temp[j, 2] <- NA 
      if (temp[j, 3] == missing) temp[j, 3] <- NA 

      temp2[j, 1] <- mean(temp[j, ], na.rm = TRUE)
      temp2[j, 2] <- sd(temp[j, ], na.rm = TRUE)
      temp2[j, 3] <- temp2[j, 2]/temp2[j, 1]
      temp2[j, 4] <- length(na.omit(temp[j, ]))

      if (temp2[j, 1] == "NaN") temp2[j, 1] <- 0

      ### Compound must be present in at least two replicates
      ### Compound must have CV less than CVMax
      if (temp2[j, 4] < 2 || temp2[j, 3] > cvmax) {
        temp2[j, 1] <- 0
        if (temp2[j, 4]  ==  3) {
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

  med_comp  <- temp3[1:eval((compounds * subjects) - summary(temp3[, 1])[7]), ]
  med_comp2 <- temp4[1:eval((compounds * subjects) - summary(temp3[, 1])[7]), ]
  med_comp  <- cbind(med_comp,med_comp2)

  colnames(metaf) <- rownames(metab)


  link2 <- matrix(colnames(metab), ncol = 1)
  colnames(link2) <- "LCMS_Run_ID"
  link3 <- unique(merge(link2, link, by.x = "LCMS_Run_ID", 
                        by.y = as.character(linktxt))[2])
  rownames(metaf) <- link3[, 1]

  rtn <- structure(list(sum_data1 = metaf,
                        clinical  = clinical,
                        medians   = med_comp),
                   class = "msprep")
  return(rtn)

}


#' Function for importing raw data and summarizing replicates
#'
#' Function reads in raw data files and summarizes technical replicates as the
#' mean of observations for compounds found in 2 or 3 replicates and with
#' coefficient of variation below specied level, or median for those found in 3
#' replicates but excess CV.
#'
#' @param Directory where the required three data files are located.
#' @param clinicalfile Name of the clinical data file.
#' @param quantificationfile Name of the data file containing the
#' quantification data.
#' @param linkfile Name of Subject link file.
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
#'   ### Specify primary directory
#'   ##directory <- c("/MSProcess/data/")
#'
#'   ### Specify location of data files
#'   clinicalfile       <- c("./data-raw/Clinical.csv")
#'   quantificationfile <- c("./data-raw/Quantification.csv")
#'   linkfile <- c("./data-raw/SubjectLinks.csv")
#'
#'   ### Set variables for program
#'   cvmax   <- 0.5
#'   missing <- 1
#'   linktxt <- "LCMS_Run_ID"
#'
#'   test <- readdata_sj(clinicalfile, quantificationfile, linkfile,
#'                       cvmax = 0.50, missing = 1, linktxt)
#'   save(test, file = paste(directory, "test.Rdata", sep = ""))
#'
#' @export

readdata_sj <- function (clinicalfile,
                         quantificationfile,
                         linkfile,
                         cvmax = 0.5,
                         missing = 1,
                         linktxt,
                         sean_extra_test_col = F,
                         sean_extra_test_row = F,
                         sean_exclude = NULL,
                         name_list = NULL,
                         id_list = NULL,
                         id_replace_string = NULL)
{
  # sean_extra_test takes care of a situation in which there are duplicates in the "quant" file, but not replicates.
  # sean_exclude is the name of a metabolite to exclude.
  # obviously, sean_extra_test and sean_exclude functionality were added by Sean Jacobson
  # name_list is where we have a list of replicates - each element in the list is a vector of replicates. If it is NULL, the program 
  #   tries to figure out what the list is manually.
  # id_list is a list of ids that the program finds in the column labels and creates the name_list using ids (which are presumably found 
  #   in the column names)
  # id_replace_string is two elements, the string that should be replaced in the IDs, and what it should be replaced by.
  #   For example, in the Lipids, the ID is listed as 21732E and should be 21731E, so you would have
  #   id_replace_string = c("21732E", "21731E")
  #   it can also be a list of sets of 2 if more than one string needs to be fixed
  if(!is.null(name_list) & !is.null(id_list)) stop("only one of name_list and id_list can be specified")
  if(!is.null(id_replace_string) & typeof(id_replace_string) == "character" & length(id_replace_string) != 2) 
    stop("id_replace_string should have a length of 2, the string to be replaced, and the string with which to replace it")
  if(!is.null(id_replace_string) & typeof(id_replace_string) == "list" & any(sapply(id_replace_string, length) != 2))
    stop("elements in id_replace_string should have a length of 2, the strings to be replaced, and the strings with which to replace it")


  quant <- read.table(paste(quantificationfile, 
                            sep = ""), header = TRUE, sep = ",")
  if(!is.null(id_replace_string) & typeof(id_replace_string) == "character")
    names(quant) <- gsub(id_replace_string[1], id_replace_string[2], names(quant))
  if(!is.null(id_replace_string) & typeof(id_replace_string) == "list"){
    for(anel in id_replace_string)
      names(quant) <- gsub(anel[1], anel[2], names(quant))
  }
  if(length(sean_exclude) > 0){
    quant <- quant[,-which(names(quant) %in% sean_exclude)]
  }
  if(sean_extra_test_col){
    dupnames <- sapply(names(quant), function(x){unlist(strsplit(x, split = "_"))[1]})[
                                                                                       duplicated(sapply(names(quant), function(x){unlist(strsplit(x, split = "_"))[1]}))]
    dups <- which(sapply(names(quant), function(x){unlist(strsplit(x, split = "_"))[1]}) %in% dupnames)
    for(adupname in dupnames){
      aa <- quant[,grep(adupname, names(quant))]
      aa$comb <- apply(aa, 1, function(x){if(any(x == 1)) return(max(x)) else return(mean(x))})
      eval(parse(text = sprintf("quant$%s <- aa$comb", gsub("_B_Aq", "_AB_Aq", names(aa)[2]))))
    }
    quant <- quant[, -dups]
  }
  if(sean_extra_test_row){
    #rows
    mzdup_test <- function(){
      mzdup <- quant$mz[duplicated(quant$mz)]
      mzdupnum <- which(quant$mz %in% mzdup)
      quant_mzdup <- quant[mzdupnum, ]
      quant_mzdup <- quant_mzdup[order(quant_mzdup$mz), ]
    }
    mzrt <- apply(quant[,1:2], 1, function(x){paste0(as.character(x), collapse = "_")})
    duprows <- mzrt[duplicated(mzrt)]
    duprs <- which(mzrt %in% duprows)
    for(aduprow in duprows){
      aa <- quant[mzrt == aduprow, ]
      quant <- rbind(quant, apply(aa, 2, mean))
    }
    quant <- quant[-duprs, ]
  }
  # end of initial stuff added by Sean Jacobson
  link <- read.table(paste(linkfile, sep = ""), 
                     header = TRUE, sep = ",")
  if(!("COPDGeneID" %in% names(link))) stop("expected to have the string 'COPDGeneID' in the link file")

  if(!is.null(id_replace_string) & typeof(id_replace_string) == "character"){
    link$LCMS_Run_Id <- gsub(id_replace_string[1], id_replace_string[2], link$LCMS_Run_Id)
    link$COPDGeneID <- gsub(id_replace_string[1], id_replace_string[2], link$COPDGeneID)
  }
  if(!is.null(id_replace_string) & typeof(id_replace_string) == "list"){
    for(anel in id_replace_string){
      link$LCMS_Run_Id <- gsub(anel[1], anel[2], link$LCMS_Run_Id)
      link$COPDGeneID <- gsub(anel[1], anel[2], link$COPDGeneID)
    }
  }

  clinical <- read.table(paste(clinicalfile, sep = ""), 
                         header = TRUE, sep = ",")
  if(length(sean_exclude) > 0){
    excl_clin <- unique(sapply(sean_exclude, function(x){unlist(strsplit(unlist(strsplit(x, split = "X"))[2], split = "_"))[1]}))
    clinical <- clinical[-which(clinical$SubjectID %in% excl_clin), ]
  }
  compounds <- nrow(quant)
  subjects <- nrow(clinical)
  metaf <- matrix(NA, ncol = compounds, nrow = subjects)
  # replicates <- (ncol(quant) - 2)/subjects
  # if (replicates%%1 != 0) {
  #    stop("Check datafile for missing replicates")
  # }
  # temp3 <- matrix(NA, ncol = replicates, nrow = subjects * compounds)
  kk <- 1
  rownames(quant) <- paste(quant$mz, "_", quant$rt, sep = "")
  #metab <- subset(quant, select = -c(mz, rt))
  library(dplyr)
  metab <- dplyr::select(quant, -c(mz, rt))
  ## More stuff added: here's where we intelligently select replicates
  if(is.null(name_list)){
    if(is.null(id_list))
      unique_ids <- unique(sapply(names(metab), function(x) paste0(unlist(strsplit(unlist(strsplit(x, split = "_"))[1], split = ""))[-1], collapse = "")))
    if(!is.null(id_list))
      unique_ids <- unique(id_list)
    name_list <- lapply(unique_ids, function(x) names(metab)[grep(x, names(metab))])
  }

  temp3 <- matrix(NA, ncol = max(sapply(name_list, length)), nrow = subjects * compounds)
  temp4 <- matrix(NA, ncol = 2, nrow = subjects * compounds)

  for(ii in 1:length(name_list)){#1:subjects) {
    # bcall <- ""
    # for (m in 1:replicates) {
    #   bcall <- paste(bcall, "metab[,i*replicates-", eval(m - 
    #     1), "],", sep = "")
    # }
    # bcall2 <- paste("matrix(c(", substr(bcall, 1, nchar(bcall) - 
    #   1), "),ncol=", replicates, ")", sep = "")
    # temp <- eval(parse(text = bcall2))
    temp <- metab[,name_list[[ii]]]
    replicates <- length(name_list[[ii]])
    temp2 <- matrix(NA, ncol = 4, nrow = compounds)
    # thething <- numeric(compounds)
    for(jj in 1:compounds){
      if(replicates == 1){
        temp2[jj, 1] <- temp[jj]#, 1]
      }else{
        for(m in 1:replicates) {
          if(temp[jj, m] == missing) {
            temp[jj, m] <- NA
          }
        }
        temp2[jj, 1] <- mean(unlist(temp[jj, ]), na.rm = TRUE)
        temp2[jj, 2] <- sd(unlist(temp[jj, ]), na.rm = TRUE)
        temp2[jj, 3] <- temp2[jj, 2]/temp2[jj, 1]
        if (is.na(temp2[jj, 3])) {
          temp2[jj, 3] <- cvmax + 1
        }
        temp2[jj, 4] <- length(na.omit(unlist(temp[jj, ])))
        if (temp2[jj, 1] == "NaN") {
          temp2[jj, 1] <- 0
        }
        if(temp2[jj, 4] < ceiling(replicates/2) || temp2[jj, 3] > cvmax){
          # thething[jj] <- 1
          temp2[jj, 1] <- 0
          if(temp2[jj, 4] == replicates){
            temp2[jj, 1] <- median(unlist(temp[jj, ]), na.rm = TRUE)
            temp3[kk, ] <- c(unlist(temp[jj, ]), rep(NA, 3 - length(unlist(temp[jj, ]))))
            temp4[kk, 1] <- rownames(metab)[jj]
            temp4[kk, 2] <- name_list[[ii]][1]#colnames(metab)[i]
            kk <- kk + 1
          }
        }
      }
    }
    metaf[ii, ] <- t(temp2[, 1])
  }
  if(max(sapply(name_list, length)) != 1){
    med_comp <- temp3[1:eval((compounds * subjects) - sum(is.na(temp3[, 1]))), ]
    med_comp2 <- temp4[1:eval((compounds * subjects) - sum(is.na(temp3[, 1]))), ]
    med_comp <- cbind(med_comp, med_comp2)
  } else {
    med_comp <- 0
  }

  colnames(metaf) <- rownames(metab)

  link2 <- matrix(colnames(metab), ncol = 1)
  link2 <- matrix(sapply(link2, function(x){gsub("AB_A", "A_A", x)}))
  colnames(link2) <- "LCMS_Run_ID"

  if (all(sapply(link2[,1], function(x){unlist(strsplit(x, split = ""))[1] == "X"})))
    link2[,1] <- sapply(link2[,1], function(x){paste0(unlist(strsplit(x, split = ""))[-1], collapse = "")})
  # added by me:
  if(nrow(link) > 0 & ncol(link) > 0){
    link3 <- unique(merge(link2, link, by.x = "LCMS_Run_ID", 
                          by.y = as.character(linktxt))[2])
  } else link3 <- link2

  #link3 <- unique(merge(, link, by.x = "LCMS_Run_ID", 
  # by.y = as.character(linktxt))[2])
  rownames(metaf) <- link3[, 1]
  list(sum_data1 = metaf,
       clinical  = clinical,
       medians   = med_comp)
}



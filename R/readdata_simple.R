
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
#' @param rt Retention time.
#' @param mz Mass-to-charge ratio.
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
#'   # LCMS_Run_ID = operator/replicate (A-C), subject (01-03), concentration  (1x,2x,4x)
#'   # SubjectID   = subject (01-03), concentration  (1x,2x,4x)
#'
#'   ### Read in data files
#'   clinical <- read.csv("./data-raw/Clinical.csv")
#'   quant    <- read.csv("./data-raw/Quantification.csv")
#'   link     <- read.csv("./data-raw/SubjectLinks.csv")
#'
#'   tidyquant <- tidy_quant(quant, my_mz = "mz", my_rt = "rt")
#' 
#'   ### Set variables for program
#'   cvmax   <- 0.5
#'   missing <- 1
#'   linktxt <- "LCMS_Run_ID"
#'
#'   # Convert dataset to tidy-ish format
#'   #  need replicate, subject_id, mz, and rt as variables
#'
#'   test <- msprep(clinical_data       = clinical,
#'                  quantification_data = quant,
#'                  link_data           = link,
#'                  cvmax               = 0.50,
#'                  missing             = 1,
#'                  linktxt             = linktxt)
#'   #save(test, file = paste0("./data/test.Rdata"))
#'
#' str(quant)
#' str(test$sum_data)
#'
#'
#' save.image("readdata_workspace.RData")
#'
#'
#'
#' 
#' 
#' 
#' 
#' 
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @export


tidy_quant <- function(quantification_data, mz, rt) {

  quantification_data %>%
    tibble::as_data_frame(.) %>%
    tidyr::gather(key = LCMS_Run_ID, value = abundance, -mz, -rt) %>%
    dplyr::mutate(LCMS_Run_ID = stringr::str_replace_all(LCMS_Run_ID, "Neutral_Operator_Dif_Pos_", "")) %>%
    tidyr::separate(LCMS_Run_ID, sep = "_", into = c("spike", "subject_id", "replicate"))

}

replace_missing <- function(abundance, missing_val) {
  ifelse(abundance == missing, 0, abundance)
}

mean_or_median <- function(cv_abundance, cv_max, presence_count) {
  if_else(cv_abundance > cvmax & presence_count > 2, "median", "mean") %>%
    if_else(is.na(.) | is.nan(.), 0, .)
}

quant_tidy <- tidy_quant(quant, my_mz = "mz", my_rt = "rt")

#msprep <- function(clinical_data,
#                   quantification_data,
#                   link_data,
#                   subject_id,
#                   replicate,
#                   mz,
#                   rt,
#                   cvmax   = 0.50,
#                   missing = 1,
#                   linktxt) {

# Start temporary code #
clinical_data       <- clinical
quantification_data <- quant
link_data           <- link
cvmax               <- 0.50
missing             <- 1
linktxt             <- linktxt
# End temporary code #

  ### Read in data
  clinical <- clinical_data
  quant    <- quantification_data
  link     <- link_data

  # Replace and remove NAs
  quant_tidy <- 
    quant_tidy %>% 
      mutate(abundance = replace_missing(abundance, missing_val))

  replicate_count <- length(unique(quant_tidy$replicate))

  # Calculate initial summary measures
  quant_summary <-
    quant_tidy %>%
      group_by(subject_id, spike, mz, rt) %>%
      arrange(mz, rt, subject_id, spike, replicate) %>%
      summarise(presence_count   = n(),
                mean_abundance   = mean(abundance),
                sd_abundance     = sd(abundance),
                median_abundance = median(abundance)) 

  # Identify and select summary measure
  quant_summary <- 
    quant_summary %>%
    mutate(cv_abundance      = sd_abundance / mean_abundance) %>%
    mutate(summary_measure   =
             mean_or_median(cv_abundance, cv_max, presence_count)) %>%
    mutate(abundance_summary =
             if_else(summary_measure == "median", median_abundance,
                     mean_abundance)) %>%
    ungroup

  compounds <- quant_summary %>% select(mz, rt) %>% distinct %>% nrow
  # actually count of subject,spike pairs
  subjects  <- quant_summary %>% select(subject_id, spike) %>% distinct %>% nrow

  test$sum_data1[, 1:10]
  sum_data <-
    quant_summary %>% select(subject_id, spike, mz, rt, abundance_summary) %>%
    unite(id, spike, subject_id, sep = "_") %>% 
    unite(metabolite, mz, rt, sep = "_")

  new_sum_data <- sum_data %>% spread(key = id, value = abundance_summary)

  old_sum_data <- 
    test$sum_data1 %>% as.data.frame %>% 
      rownames_to_column(var = "id") %>%
      gather(key = metabolite, value = abundance_summary, -id) %>%
      as_data_frame  %>%
      select(id, metabolite, abundance_summary) %>%
      separate(metabolite, into = c("mz", "rt"), sep = "_") %>%
      mutate(mz = as.numeric(mz)) %>%
      arrange(mz) %>%
      unite(metabolite, mz, rt, sep = "_")

  # NOTE: difference is missings were set to 0 at some point?
  anti_join(old_sum_data, sum_data)

    





  ########################################
  # Old version
  ########################################
  metaf <- matrix(NA, ncol = compounds, nrow = subjects)
  temp3 <- matrix(NA, ncol = 3, nrow = subjects * compounds)
  temp4 <- matrix(NA, ncol = 2, nrow = subjects * compounds)
  k <- 1

  rownames(quant) <- paste(quant$mz, quant$rt, sep = "_")
  metab <- subset(quant, select = -c(mz, rt))

  for (i in 1:subjects) {
    # selects subject and all operators
    temp  <- matrix(c(metab[, (i * 3 - 2)], metab[, (i * 3 - 1)], metab[, (i * 3)]), ncol = 3)
    temp2 <- matrix(NA, ncol = 4, nrow = compounds)

    # Replace missing values
    temp <- ifelse(temp == missing, NA, temp) 
    # Calculate mean, sd, cv across replicates, and count of reps where metabolite found
    temp2[, 1] <- rowMeans(temp, na.rm = TRUE)
    temp2[, 2] <- matrixStats::rowSds(temp, na.rm = TRUE)
    temp2[, 3] <- temp2[, 2] / temp2[, 1]
    temp2[, 4] <- apply(temp, 1, function(x) length(na.omit(x)))
    # if cv > cvmax and all replicates, use median instead of mean
    temp2[, 1] <- sapply(1:nrow(temp2), 
                         function(x) {
                           ifelse((temp2[x, 3] > cvmax) & (temp2[x, 4] == 3),
                                  median(temp[x, ]),
                                  temp2[x, 1])
                         })
    # Replace NaNs w/ NA
    temp2 <- ifelse(is.nan(temp2), NA, temp2)

    # Not currently re-implemented, but provide data frame of metabolites
    # summarized by median -- including all repl values, mz, rt, id
    for (j in 1:compounds) {
      ### Compound must be present in at least two replicates
      ### Compound must have CV less than CVMax
      if (temp2[j, 4] < 2 || temp2[j, 3] > cvmax) {
        temp2[j, 1] <- 0
        if (temp2[j, 4]  ==  3) {

          # if cv too big and all replicates
          temp3[k, ]  <- temp[j, ]
          temp4[k, 1] <- rownames(metab)[j]
          temp4[k, 2] <- colnames(metab)[i]
          k           <- k + 1
        }
      }
    }
    metaf[i, ] <- t(temp2[, 1])
  }

  # result of this for loop is:
  #   metaf - rows filled in by subject using transposed column 1 (means) of temp2
  #   temp2 - columns are mean (1), sd (2), coefficient of variation CV (3),
  #           count of replications where metabolite is present (4)
  #   temp  - 

  # Metabolites summarized by medians (cv too high and have 3 replicates)
  med_comp  <- temp3[1:eval((compounds * subjects) - summary(temp3[, 1])[7]), ]
  med_comp2 <- temp4[1:eval((compounds * subjects) - summary(temp3[, 1])[7]), ]
  med_comp  <- cbind(med_comp,med_comp2)

  colnames(metaf) <- rownames(metab)


  link2 <- matrix(colnames(metab), ncol = 1)
  colnames(link2) <- "LCMS_Run_ID"
  link3 <- unique(merge(link2, link, by.x = "LCMS_Run_ID", 
                        by.y = as.character(linktxt))[2])
  rownames(metaf) <- link3[, 1]

  rtn <- structure(list(sum_data  = metaf,
                        clinical  = clinical,
                        medians   = med_comp),
                   class = "msprep")
#  return(rtn)
#
#}

#' 
#' Prepare a mass spec quantification data frame for filtering, imputation,
#' and normalization. Also provides summaries of data structure (replicates,
#' subjects, spike, etc.) 
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
#'   tidyquant  <- tidy_quant(quant, mz = "mz", rt = "rt")
#'   msprep_obj <- msprep(quant, 
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
#' ########################
#' # Start temporary code #
#'
#' # source 3 fns at bottom of script,
#' clinical      <- read.csv("./data-raw/Clinical.csv")
#' quant         <- read.csv("./data-raw/Quantification.csv")
#' link          <- read.csv("./data-raw/SubjectLinks.csv")
#' .data         <- tidy_quant(quant, mz = "mz", rt = "rt")
#' clinical_data <- clinical
#' link_data     <- link
#' cvmax         <- 0.50
#' missing       <- 1
#' linktxt       <- "LCMS_Run_ID"
#' # End temporary code #
#' ########################
#'
#'
#'
#'
#'
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @export
msprep <- function(.data,
                   clinical_data,
                   link_data,
                   subject_id,
                   replicate,
                   mz,
                   rt,
                   cvmax   = 0.50,
                   missing = 1,
                   linktxt,
                   
                   min_proportion_present = 1/3
                   ) {

  # Check args
  stopifnot(is.data.frame(.data))

  # Replace miss val with NAs 
  .data <- .data %>% mutate(abundance = replace_missing(abundance, missing_val))

  # Get replicate count for each mz/rt/spike/subject combo
  # NOTE: verify that replicate is a scalar for a given dataset
  replicate_count <- length(unique(.data[[replicate]]))
  #replicate_count <-
  #  .data %>%
  #  select(mz, rt, spike, subject_id, replicate) %>%
  #  group_by(mz, rt, spike, subject_id) %>%
  #  distinct %>%
  #  summarise(replicate_count = n())

  # Calculate initial summary measures
  quant_summary <-
    .data %>%
      group_by(subject_id, spike, mz, rt) %>%
      arrange(mz, rt, subject_id, spike, replicate) %>%
      summarise(n_present        = sum(!is.na(abundance)),
                prop_present     = n_present / replicate_count,
                mean_abundance   = mean(abundance, na.rm = T),
                sd_abundance     = sd(abundance, na.rm = T),
                median_abundance = median(abundance, na.rm = T)) #%>%
      #full_join(replicate_count)

  # Identify and select summary measure
  quant_summary <-
    quant_summary %>%
    mutate(cv_abundance      = sd_abundance / mean_abundance,
           summary_measure   = select_summary_measure(min_proportion_present,
                                                      cv_abundance, cvmax,
                                                      replicate_count, n_present), 
           abundance_summary = ifelse(summary_measure == "median", median_abundance, mean_abundance)) %>% # else if NA, NA
    ungroup %>%
    mutate(abundance_summary = abundance_summary %>% ifelse(is.na(.) | is.nan(.), 0, .))

  compounds <- quant_summary %>% select(mz, rt) %>% distinct %>% nrow
  # actually count of subject,spike pairs (not subjects only)
  subjects  <- quant_summary %>% select(subject_id, spike) %>% distinct %>% nrow

############################################################
# Rules from manuscript
############################################################

  cv = sd/mean
  prop_present <- n_present / n_replicates

  #  - **Only abundances that are found in at least two of three replicates are kept. **
  # TODO: check with Katerina on when to set this to 0 -- should it be only if 1
  # rep is present or some proportion of total replicates
  if (prop_present <= min_proportion_present) summary_measure <- NA # or <- 0 # default min_proportion_present = 1/3 

  #  - If CV is below the user-specified level, the average of the replicates is used. 
  # base summary stat, if not needed
  # if (cv <= cv_max) mean()

  #  - If the CV is above the specified level and found in exactly two of three
  #    replicates, the summarization is not used and the observation is left
  #    blank. 
  # TODO: check with Katerina -- should this be 2/3rds for any replicate number
  # or just when 3 replicates and 2 only present (and cv > cvmax obvi)?
  if (cv > cv_max & (n_replicates == 3 & n_present == 2)) summary_measure <- NA # or <- 0

  #  - If the compound was found in all three replicates but with unacceptable CV,
  #    the median is used as the summarization measure. 
  # TODO: check with Katerina -- should this be median if at least 3 present of 3+ replicates and cv >
  # cvmax, or strictly all present?
  if (cv > cv_max & n_present == n_replicates) summary_measure <- median()


# NOTE: text from manuscript
#  The first processing step is summarization of technical replicates, three
#  replicates required per subject/sample. MSPrep provides options to remove
#  erroneous data and to reduce the effect of extreme observations. 
#
#  - The user specifies a cutoff for the coefficient of variation (CV),
#    calculated by dividing the standard deviation of the replicates by the
#    average, yielding a measure for magnitude of the variation between
#    replicates.
#
#  The summarization routine summarizes each compound by subject (or sample) and
#  returns a single observation per compound per subject. 
#
#  - **Only abundances that are found in at least two of three replicates are kept. **
#  - If CV is below the user-specified level, the average of the replicates is used. 
#  - If the CV is above the specified level and found in exactly two of three
#    replicates, the summarization is not used and the observation is left
#    blank. 
#  - If the compound was found in all three replicates but with unacceptable CV,
#    the median is used as the summarization measure. 
# 
#  This approach removes potential erroneous data. We have found that
#  most compounds with high CV have two consistent and one extreme observation.
#  Using the median reduces the effect of the extreme observation.

############################################################

# COMPARE to old version
  test$sum_data1[, 1:10]

  sum_data <-
    quant_summary %>% select(subject_id, spike, mz, rt, abundance_summary, summary_measure) %>%
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
      #mutate(rt = as.numeric(rt)) %>%
      arrange(mz) %>%
      unite(metabolite, mz, rt, sep = "_")

  # 
  new <- sum_data %>% separate(metabolite, into = c("mz", "rt"), sep = "_") %>% arrange(id, mz, rt)
  old <- old_sum_data %>% separate(metabolite, into = c("mz", "rt"), sep = "_") %>% arrange(id, mz, rt)
  diffrows <- !(select(new, -summary_measure) == old)[, 4]
  newdiffs <- new[diffrows, ] %>% rename(new_summary = abundance_summary)
  olddiffs <- old[diffrows, ] %>% rename(old_summary = abundance_summary)

  anti_join(old_sum_data, sum_data)
  comparing <- 
    .data %>% mutate(mz = as.character(mz), rt = as.character(rt)) %>% unite(id, spike, subject_id) %>% 
    right_join(., newdiffs) %>%
    left_join(., olddiffs)








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



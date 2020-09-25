# ms_summarize <- function(data,
#                          abundance     = "abundance",
#                          met_id        = NULL,
#                          mz            = NULL,
#                          rt            = NULL,
#                          subject_id    = "subject_id",
#                          replicate     = NULL,
#                          batch         = NULL,
#                          groupingvars  = NULL,
#                          cvmax         = 0.50,
#                          missing_val   = 1,
#                          min_proportion_present = 1/3) {
# 
#     ## Check args
#     stopifnot(is.data.frame(data))
#     ## my_args  <- mget(names(formals()), sys.frame(sys.nframe()))
#     stopifnot(is.null(batch) | length(batch) == 1)
#     stopifnot(is.null(replicate) | length(replicate) == 1)
#   
#     if (is.null(met_id) & (is.null(mz) | is.null(rt))) {
#         stop("Must include 'met_id' or both 'mz' and 'rt'")
#     }
#   
#     ## Create vector of metabolite id column names
#     met_vars <- c()
#     if (!is.null(mz) & !is.null(rt)) {
#         met_vars <- c(met_vars, "mz", "rt")
#     }
#     if (!is.null(met_id)) {
#         met_vars <- c(met_vars, "met_id")
#     }
#   
#     ## rlang magic 
#     grouping_quo = syms(groupingvars)
#     met_syms <- syms(met_vars)
# 
#     ## Convert to tibble
#     data <- as_tibble(data)
# 
#     ## Replace provided variable names with standardized ones
#     data <- standardize_dataset(data, subject_id, replicate, abundance,
#                                 grouping_quo, batch, met_id, mz, rt)
#   
#     ## Create ordered vector of column names to later store as object attribute
#     col_order <- colnames(data)[!(colnames(data) %in% c("met_id", "mz", "rt",
#                                                         "abundance",
#                                                         "replicate"))]
# 
#     ## Replace miss val with NAs 
#     data <- mutate_at(data, vars("abundance"), replace_missing, missing_val)
# 
#     ## Get replicate count for each mz/rt/grouping_quo/subject combo
#     replicate_count <- length(unique(data[["replicate"]]))
# 
#     ## Roughly check if all compounds are present in each replicate
#     stopifnot(replicate_count == 0 | nrow(data) %% replicate_count == 0)
#   
#     ## Skips all summarization for data sets with no replicates
#     if (is.null(replicate)) {
#         
#         if (!is.null(groupingvars)) {
#             summary_data <- ms_arrange(data, batch, groupingvars)
#     }
#         else {
#             summary_data <- ms_arrange(data)
#         }
#     
#         ## Rename "abundance" to "abundance_summary" to comply with remainder of
#         ## package
#         summary_data <- summary_data %>% 
#             rename("abundance_summary" = "abundance")
#         
#         ## Replace NAs w/ 0
#         summary_data$abundance_summary <- 
#           replace_na(summary_data$abundance_summary, 0)
#         
#         ## Order columns
#         summary_data  <- select_at(summary_data,
#                                    vars(subject_id, batch, `!!!`(grouping_quo),
#                                         `!!!`(met_syms), "abundance_summary"))
#     
#         ## Create return object & return
#         return(structure(list("data" = summary_data,
#                               replicate_info  = NULL,
#                               medians         = NULL),
#                          replicate_count        = NULL,
#                          cvmax                  = cvmax,
#                          min_proportion_present = min_proportion_present,
#                          met_vars               = met_vars,
#                          groupingvars           = groupingvars,
#                          batch_var              = batch,
#                          replicate_var          = replicate,
#                          col_order              = col_order,
#                          stage = "summarized",
#                          class = "msprep"))
#     }
#   
#     ## Arrange and group data according to present function args
#     if (is.null(batch) & is.null(groupingvars)) {
#         quant_summary <- ms_arrange(data, replicate)
#         quant_summary <- group_by(quant_summary, subject_id, `!!!`(met_syms))
#     } else if(is.null(batch)) {
#         quant_summary <- ms_arrange(data, groupingvars, replicate)
#         quant_summary <- group_by(quant_summary, subject_id, `!!!`(met_syms),
#                                   `!!!`(grouping_quo))
#     } else if (is.null(groupingvars)) {
#         quant_summary <- ms_arrange(data, batch, replicate)
#         quant_summary <- group_by(quant_summary, subject_id, batch,
#                                   `!!!`(met_syms))
#     } else {
#         quant_summary <- ms_arrange(data, batch, groupingvars, replicate)
#         quant_summary <- group_by(quant_summary, subject_id, batch, 
#                                   `!!!`(met_syms),
#                                   `!!!`(grouping_quo))
#     }
# 
#     # Calculate remaining summary measures
#     quant_summary <- summarise(quant_summary,
#                                n_present = sum(!is.na(.data$abundance)),
#                                prop_present = UQ(sym("n_present")) / 
#                                    replicate_count,
#                                mean_abundance = mean(.data$abundance,
#                                                      na.rm = TRUE),
#                                sd_abundance = sd(.data$abundance, na.rm = TRUE),
#                                median_abundance = median(.data$abundance, 
#                                                          na.rm = TRUE))
#     quant_summary <- mutate(quant_summary, 
#                             cv_abundance = .data$sd_abundance / 
#                                 .data$mean_abundance)
#     quant_summary <- ungroup(quant_summary)
# 
#   ## Identify and select summary measure
#   quant_summary <-
#     mutate(quant_summary,
#            summary_measure = 
#              select_summary_measure(.data$n_present, .data$cv_abundance, 
#                                     replicate_count,
#                                     min_proportion_present, cvmax))
#   quant_summary <-
#     mutate(quant_summary,
#            abundance_summary = 
#              case_when(.data$summary_measure == "median" ~ 
#                          .data$median_abundance,
#                        .data$summary_measure == "mean"   ~ 
#                          .data$mean_abundance,
#                        TRUE                              ~ 0))
#   quant_summary <- mutate_at(quant_summary,
#                              vars(subject_id, "summary_measure", batch, 
#                                   `!!!`(grouping_quo)),
#                              factor)
# 
#   # Extract summarized dataset
#   summary_data  <- select_at(quant_summary, 
#                              vars(subject_id, batch, `!!!`(grouping_quo), 
#                                   `!!!`(met_syms), "abundance_summary"))
#   if(!is.null(batch) & !is.null(groupingvars)){
#     summary_data <- ms_arrange(summary_data, batch, groupingvars)
#   }
#   else {
#     summary_data <- ms_arrange(summary_data)
#   }
#   
# 
#   # Additional info extracted in summarizing replicates
#   replicate_info <- select_at(quant_summary, 
#                               vars(subject_id, batch, `!!!`(grouping_quo), 
#                                    `!!!`(met_syms), "n_present",
#                                    "cv_abundance", "summary_measure"))
# 
#   # Summaries that used medians
#   medians        <- filter(quant_summary, .data$summary_measure == "median")
#   medians        <- select_at(medians, vars(subject_id, batch, 
#                                             `!!!`(grouping_quo), 
#                                             `!!!`(met_syms),
#                                             "abundance_summary"))
#   
#   # Replace medians with NULL if no medians used
#   if(nrow(medians) == 0){
#     medians = NULL
#   }
# 
#   # Total number of compounds identified
#   n_compounds    <- nrow(distinct(select(quant_summary, `!!!`(met_syms))))
# 
#   # Create return object & return
#   rtn <- structure(list("data" = summary_data,
#                            replicate_info  = replicate_info,
#                            medians         = medians),
#                       replicate_count        = replicate_count,
#                       cvmax                  = cvmax,
#                       min_proportion_present = min_proportion_present,
#                       met_vars               = met_vars,
#                       groupingvars           = groupingvars,
#                       batch_var              = batch,
#                       replicate_var          = replicate,
#                       col_order              = col_order,
#                       stage = "summarized",
#                       class = "msprep")
#   
#   attr(rtn, "patient_count") <- length(unique(rtn$data$subject_id))
#   
#   return(rtn)
# 
# }

# #' @export
# print.msprep <- function(x, ...) {
# 
#   stage <- stage(x)
# 
#   # Print object info
#   cat("msprep dataset\n")
#   cat(paste("    Stage:", stage, "\n"))
#   cat("    Patient count: ", attr(x, "patient_count"), "\n")
#   
#   # Print grouping vars (may be NULL)
#   if (!is.null(grouping_vars(x))) {
#     cat("    Grouping vars:", paste(grouping_vars(x), collapse = ", "), "\n")
#   }
#   else {
#     cat("    Grouping vars: none present\n")
#   }
#   
#   # Print batch var (may be NULL)
#   if (!is.null(batch_var(x))) {
#     cat("    Batch var: ", batch_var(x), "\n")
#   }
#   else {
#     cat("    Batch var: none present\n")
#   }
#   
#   # Print met_vars (never NULL)
#   cat("    Met vars: ", paste(attr(x, "met_vars"), collapse = ", "), "\n")
#   
#   # Print summary info
#   cat("    Summarization:\n")
#   cat("        Replicate count: ", attr(x, "replicate_count"), "\n")
#   cat("        User-defined parameters: \n")
#   cat("          cvmax = ", attr(x, "cvmax"), "\n")
#   cat("          min_proportion_present = ", 
#       round(attr(x, "min_proportion_present"), digits = 3), "\n")
#   cat("        Resulting stats:\n")
#   # Print median summarized count (may be NULL)
#   if (!is.null(x$medians)) {
#     cat("          Count of patient-compounds summarized by median: ", 
#         nrow(x$medians), "\n")
#   }
#   else {
#     cat("          Count of patient-compounds summarized by median: 0\n")
#   }
#   
#   # Print filter info
#   if (stage %in% msprep_stages()[2:4]) {
#     cat("    Filtering:\n")
#     cat("        User-defined parameters: \n")
#     cat("          Filter percent = ", attr(x, "filter_percent"), "\n")
#     cat("        Resulting stats: \n")
#     cat("          Count of compounds present in >= ", 
#         round(attr(x, "filter_percent")*100, digits = 3), "% of samples = ", 
#         sum(attr(x, "filter_status")$keep), "\n", sep = "")
#   }
#   
#   # Print imputation info
#   if (stage %in% msprep_stages()[3:4]) {
#     cat("    Imputation:\n")
#     cat("        Imputation Method: ", attr(x, "impute_method"), "\n")
#     cat("        Number of imputed values: ", attr(x, "missing_count"), "\n")
#   }
#   
#   # Print normalization info
#   if (stage %in% msprep_stages()[4]) {
#     cat("    Normalization:\n")
#     cat("        Normalization Method: ", attr(x, "normalize_method"), "\n")
#     cat("        Transformation: ", attr(x, "transformation"), "\n")
#   }
# 
#   # Print dataset (head)
#   cat("    Dataset:\n")
#   print(x$data, n = 6)
# 
# }



# #
# # Internal ms_summarize functions
# #
# 
# # #' @importFrom dplyr case_when
# select_summary_measure <- function(n_present,
#                                    cv_abundance,
#                                    n_replicates,
#                                    min_proportion_present,
#                                    cv_max) {
# 
#   case_when((n_present / n_replicates) <= min_proportion_present
#               ~ "none: proportion present <= min_proportion_present",
#             cv_abundance > cv_max & (n_replicates == 3 & n_present == 2) 
#               ~ "none: cv > cvmax & 2 present",
#             cv_abundance > cv_max & (n_present == n_replicates) 
#               ~ "median",
#             TRUE ~ "mean")
# 
# }

# #' @importFrom dplyr rename
# #' @importFrom rlang sym
# #' @importFrom rlang UQ
# standardize_dataset <- function(data, subject_id, replicate, abundance, 
#                                 groupingvars, batch, met_id, mz, rt) {
# 
#   # Rename required variables
#   subject_id = sym(subject_id)
#   abundance  = sym(abundance)
#   
#   data <- data %>% 
#     rename("subject_id" = UQ(subject_id),
#            "abundance"  = UQ(abundance))
#   
#   # Rename mz/rt if used
#   if (!is.null(mz) & !is.null(rt)) {
#     mz         = sym(mz)
#     rt         = sym(rt)
#   
#     data <- data %>% 
#       rename("mz"         = UQ(mz),
#              "rt"         = UQ(rt))
#   }
#   
#   # Rename met_id if used
#   if(!is.null(met_id)){
#     met_id         = sym(met_id)
#     
#     data <- data %>% 
#       rename("met_id" = UQ(met_id))
#   }
# 
#   data <- standardize_datatypes(data, groupingvars, batch)
# 
#   # Rename optional variables if present
#   if (!is.null(replicate)) {
#     replicate <- sym(replicate)
#     data      <- data %>% rename("replicate" = UQ(replicate))
#   }
# 
#   if (!is.null(batch)) {
#     batch <- sym(batch)
#     data      <- data %>% rename("batch" = UQ(batch))
#   }
# 
#   return(data)

# }

# #' @importFrom dplyr mutate_at
# standardize_datatypes <- function(data, groupingvars, batch) {
# 
#   count_var   <- c("abundance_summary", "abundance")
#   count_var   <- count_var[count_var %in% colnames(data)]
#   factor_vars <- c("subject_id", as.character(groupingvars), batch)
#   data <- mutate_at(data, factor_vars, as.factor)
#   
#   if ("rt" %in% colnames(data) & "mz" %in% colnames(data)) {
#     data <- mutate_at(data, c("mz", "rt", count_var), as.numeric) 
#   }
#   if ("met_id" %in% colnames(data)) {
#     data <- mutate_at(data, "met_id", as.factor)
#   }
#   
#   return(data)
# }

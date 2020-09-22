# context("check that new equals old")
# 
# # Function for comparing new and old data
# compareOldNew <- function(old, new) {
#     
#     if(any(is.null(old$data[["1x_O1_01"]]), is.null(old$data[["2x_O1_01"]]),
#            is.null(old$data[["4x_O1_01"]]), is.null(old$data[["1x_O2_01"]]),
#            is.null(old$data[["2x_O2_01"]]), is.null(old$data[["4x_O2_01"]]),
#            is.null(old$data[["1x_O3_01"]]), is.null(old$data[["2x_O3_01"]]),
#            is.null(old$data[["4x_O3_01"]]), is.null(old$data[["1x_O1_02"]]),
#            is.null(old$data[["2x_O1_02"]]), is.null(old$data[["4x_O1_02"]]),
#            is.null(old$data[["1x_O2_02"]]), is.null(old$data[["2x_O2_02"]]),
#            is.null(old$data[["4x_O2_02"]]), is.null(old$data[["1x_O3_02"]]),
#            is.null(old$data[["2x_O3_02"]]), is.null(old$data[["4x_O3_02"]]),
#            is.null(new[["1x_O1_01"]]), is.null(new[["2x_O1_01"]]),
#            is.null(new[["4x_O1_01"]]), is.null(new[["1x_O2_01"]]),
#            is.null(new[["2x_O2_01"]]), is.null(new[["4x_O2_01"]]),
#            is.null(new[["1x_O3_01"]]), is.null(new[["2x_O3_01"]]),
#            is.null(new[["4x_O3_01"]]), is.null(new[["1x_O1_02"]]),
#            is.null(new[["2x_O1_02"]]), is.null(new[["4x_O1_02"]]),
#            is.null(new[["1x_O2_02"]]), is.null(new[["2x_O2_02"]]),
#            is.null(new[["4x_O2_02"]]), is.null(new[["1x_O3_02"]]),
#            is.null(new[["2x_O3_02"]]), is.null(new[["4x_O3_02"]]))) {
#         return(FALSE)
#     }
#         
#     all(old$data[["1x_O1_01"]] == new[["1x_O1_01"]],
#         old$data[["2x_O1_01"]] == new[["2x_O1_01"]],
#         old$data[["4x_O1_01"]] == new[["4x_O1_01"]],
#         old$data[["1x_O2_01"]] == new[["1x_O2_01"]],
#         old$data[["2x_O2_01"]] == new[["2x_O2_01"]],
#         old$data[["4x_O2_01"]] == new[["4x_O2_01"]],
#         old$data[["1x_O3_01"]] == new[["1x_O3_01"]],
#         old$data[["2x_O3_01"]] == new[["2x_O3_01"]],
#         old$data[["4x_O3_01"]] == new[["4x_O3_01"]],
#         old$data[["1x_O1_02"]] == new[["1x_O1_02"]],
#         old$data[["2x_O1_02"]] == new[["2x_O1_02"]],
#         old$data[["4x_O1_02"]] == new[["4x_O1_02"]],
#         old$data[["1x_O2_02"]] == new[["1x_O2_02"]],
#         old$data[["2x_O2_02"]] == new[["2x_O2_02"]],
#         old$data[["4x_O2_02"]] == new[["4x_O2_02"]],
#         old$data[["1x_O3_02"]] == new[["1x_O3_02"]],
#         old$data[["2x_O3_02"]] == new[["2x_O3_02"]],
#         old$data[["4x_O3_02"]] == new[["4x_O3_02"]])
# }
# 
# set.seed(9997)
# 
# data("msquant")
# 
# tidied_data <- ms_tidy(msquant, mz = "mz", rt = "rt",
#                        col_extra_txt = "Neutral_Operator_Dif_Pos_",
#                        separator = "_", 
#                        col_names = c("spike", "batch", "replicate", 
#                                      "subject_id"))
# 
# summarized_data <- ms_summarize(tidied_data, mz = "mz", rt = "rt", 
#                                 replicate = "replicate", 
#                                 batch = "batch", groupingvars = "spike", 
#                                 subject_id = "subject_id", 
#                                 cvmax = 0.50, min_proportion_present = 1/3, 
#                                 missing_val = 1)
# 
# filtered_data <- ms_filter(summarized_data, filter_percent = .8)
# 
# hm_imputed_data <- ms_impute(filtered_data, imputeMethod = "halfmin")
# bpca_imputed_data <- ms_impute(filtered_data, imputeMethod = "bpca")
# knn_imputed_data <- ms_impute(filtered_data, imputeMethod = "knn")
# 
# med_normalized_data <- ms_normalize(hm_imputed_data, 
#                                     normalizeMethod ="median", 
#                                     transform = "log10")
# quant_normalized_data <- ms_normalize(hm_imputed_data, 
#                                       normalizeMethod ="quantile",
#                                       transform = "log10")
# suppressWarnings(combat_normalized_data <- ms_normalize(hm_imputed_data, 
#                                        normalizeMethod ="ComBat",
#                                        transform = "log10"))
# suppressWarnings(sva_normalized_data <- ms_normalize(hm_imputed_data, 
#                                                      normalizeMethod ="SVA",
#                                     transform = "log10"))
# suppressWarnings(quant_comb_normalized_data <- 
#                      ms_normalize(hm_imputed_data, 
#                                   normalizeMethod = "quantile + ComBat",
#                                   transform = "log10"))
# suppressWarnings(med_comb_normalized_data <- ms_normalize(hm_imputed_data,
#                                          normalizeMethod = 
#                                              "median + ComBat",
#                                          transform = "log10"))
# suppressWarnings(crmn_normalized_data <- ms_normalize(hm_imputed_data,
#                                      normalizeMethod = "CRMN",
#                                      transform = "log10"))
# suppressWarnings(ruv_normalized_data <- ms_normalize(hm_imputed_data, 
#                                                      normalizeMethod ="RUV",
#                                     controls = NULL,  n_control = 5, k_ruv = 5, 
#                                     transform = "log10"))
# 
# summarized_data <- summarized_data %>% ms_return()
# 
# filtered_data <- filtered_data %>% ms_return()
# 
# hm_imputed_data <- hm_imputed_data %>% ms_return()
# bpca_imputed_data <- bpca_imputed_data %>% ms_return()
# knn_imputed_data <- knn_imputed_data %>% ms_return()
# 
# med_normalized_data <- med_normalized_data %>% ms_return()
# quant_normalized_data <- quant_normalized_data %>% ms_return()
# combat_normalized_data <- combat_normalized_data %>% ms_return()
# sva_normalized_data <- sva_normalized_data %>% ms_return()
# quant_comb_normalized_data <- quant_comb_normalized_data %>% ms_return()
# med_comb_normalized_data <- med_comb_normalized_data %>% ms_return()
# crmn_normalized_data <- crmn_normalized_data %>% ms_return()
# ruv_normalized_data <- ruv_normalized_data %>% ms_return()
# 
# summarizedDF <- msSummarize(msquant,
#                             compVars = c("mz", "rt"),
#                             sampleVars = c("spike", "batch", "replicate", 
#                                          "subject_id"),
#                             cvMax = 0.50,
#                             minPropPresent = 1/3,
#                             returnSummaryDetails = FALSE,
#                             colExtraText = "Neutral_Operator_Dif_Pos_",
#                             separator = "_",
#                             missingValue = 1)
# 
# filteredDF <- msFilter(summarizedDF,
#                        filterPercent = 0.8,
#                        compVars = c("mz", "rt"),
#                        sampleVars = c("spike", "batch", "subject_id"),
#                        separator = "_",
#                        returnToSE = FALSE)
# 
# hmImputedDF <- msImpute(filteredDF, imputeMethod = "halfmin",
#                         compVars = c("mz", "rt"),
#                         sampleVars = c("spike", "batch", "subject_id"),
#                         separator = "_",
#                         returnToSE = FALSE,
#                         missingValue = 0)
# bpcaImputedDF <- msImpute(filteredDF, imputeMethod = "bpca",
#                         compVars = c("mz", "rt"),
#                         sampleVars = c("spike", "batch", "subject_id"),
#                         separator = "_",
#                         returnToSE = FALSE,
#                         missingValue = 0)
# knnImputedDF <- msImpute(filteredDF, imputeMethod = "knn",
#                          compVars = c("mz", "rt"),
#                          sampleVars = c("spike", "batch", "subject_id"),
#                          separator = "_",
#                          returnToSE = FALSE,
#                          missingValue = 0)
# 
# medianNormalized <- msNormalize(hmImputedDF, normalizeMethod = "median",
#                                 transform = "log10",
#                                 compVars = c("mz", "rt"),
#                                 sampleVars = c("spike", "batch", "subject_id"),
#                                 separator = "_",
#                                 returnToSE = FALSE)
# quantNormalized <- msNormalize(hmImputedDF, normalizeMethod = "quantile",
#                                transform = "log10",
#                                compVars = c("mz", "rt"),
#                                sampleVars = c("spike", "batch", "subject_id"),
#                                separator = "_",
#                                returnToSE = FALSE)
# combatNormalized <- msNormalize(hmImputedDF, normalizeMethod = "ComBat", 
#                                 compVars = c("mz", "rt"),
#                                 sampleVars = c("spike", "batch", "subject_id"),
#                                 covariatesOfInterest = c("spike"),
#                                 separator = "_",
#                                 returnToSE = FALSE)
# svaNormalized <- msNormalize(hmImputedDF, normalizeMethod = "SVA", 
#                              compVars = c("mz", "rt"),
#                              sampleVars = c("spike", "batch", "subject_id"),
#                              covariatesOfInterest = c("spike"),
#                              separator = "_",
#                              returnToSE = FALSE)
# quantCombatNormalized <- msNormalize(hmImputedDF, 
#                                      normalizeMethod = "quantile + ComBat",
#                                      compVars = c("mz", "rt"),
#                                      sampleVars = c("spike", "batch", 
#                                                   "subject_id"),
#                                      covariatesOfInterest = c("spike"),
#                                      separator = "_",
#                                      returnToSE = FALSE)
# medCombatNormalized <- msNormalize(hmImputedDF,
#                                    normalizeMethod = "median + ComBat",
#                                    compVars = c("mz", "rt"),
#                                    sampleVars = c("spike", "batch", 
#                                                   "subject_id"),
#                                    covariatesOfInterest = c("spike"),
#                                    separator = "_",
#                                    returnToSE = FALSE)
# crmnNormalized <- msNormalize(hmImputedDF, 
#                               normalizeMethod = "CRMN",
#                               compVars = c("mz", "rt"),
#                               sampleVars = c("spike", "batch", "subject_id"),
#                               covariatesOfInterest = c("spike"),
#                               separator = "_",
#                               returnToSE = FALSE)
# ruvNormalized <- msNormalize(hmImputedDF, 
#                              normalizeMethod = "RUV",
#                              compVars = c("mz", "rt"),
#                              sampleVars = c("spike", "batch", "subject_id"),
#                              covariatesOfInterest = c("spike"),
#                              controls = NULL,  nControl = 5, kRUV = 5,
#                              separator = "_",
#                              returnToSE = FALSE)
# 
# 
# test_that("msSummarize() == ms_summarize()", {
#     expect_true(all(summarized_data$data[["1x_O1_01"]] == 
#                         summarizedDF[["1x_O1_01"]],
#                     summarized_data$data[["2x_O1_01"]] == 
#                         summarizedDF[["2x_O1_01"]],
#                     summarized_data$data[["4x_O1_01"]] == 
#                         summarizedDF[["4x_O1_01"]],
#                     summarized_data$data[["1x_O2_01"]] == 
#                         summarizedDF[["1x_O2_01"]],
#                     summarized_data$data[["2x_O2_01"]] == 
#                         summarizedDF[["2x_O2_01"]],
#                     summarized_data$data[["4x_O2_01"]] == 
#                         summarizedDF[["4x_O2_01"]],
#                     summarized_data$data[["1x_O3_01"]] == 
#                         summarizedDF[["1x_O3_01"]],
#                     summarized_data$data[["2x_O3_01"]] == 
#                         summarizedDF[["2x_O3_01"]],
#                     summarized_data$data[["4x_O3_01"]] == 
#                         summarizedDF[["4x_O3_01"]],
#                     summarized_data$data[["1x_O1_02"]] == 
#                         summarizedDF[["1x_O1_02"]],
#                     summarized_data$data[["2x_O1_02"]] == 
#                         summarizedDF[["2x_O1_02"]],
#                     summarized_data$data[["4x_O1_02"]] == 
#                         summarizedDF[["4x_O1_02"]],
#                     summarized_data$data[["1x_O2_02"]] == 
#                         summarizedDF[["1x_O2_02"]],
#                     summarized_data$data[["2x_O2_02"]] == 
#                         summarizedDF[["2x_O2_02"]],
#                     summarized_data$data[["4x_O2_02"]] == 
#                         summarizedDF[["4x_O2_02"]],
#                     summarized_data$data[["1x_O3_02"]] == 
#                         summarizedDF[["1x_O3_02"]],
#                     summarized_data$data[["2x_O3_02"]] == 
#                         summarizedDF[["2x_O3_02"]],
#                     summarized_data$data[["4x_O3_02"]] == 
#                         summarizedDF[["4x_O3_02"]]))
# })
# 
# test_that("2msSummarize() == ms_summarize()", {
#     expect_true(compareOldNew(summarized_data, summarizedDF))
# })
# 
# test_that("msFilter() == ms_filter()", {
#     expect_true(all(filtered_data$data[["1x_O1_01"]] == 
#                         filteredDF[["1x_O1_01"]],
#                     filtered_data$data[["2x_O1_01"]] == 
#                         filteredDF[["2x_O1_01"]],
#                     filtered_data$data[["4x_O1_01"]] == 
#                         filteredDF[["4x_O1_01"]],
#                     filtered_data$data[["1x_O2_01"]] == 
#                         filteredDF[["1x_O2_01"]],
#                     filtered_data$data[["2x_O2_01"]] == 
#                         filteredDF[["2x_O2_01"]],
#                     filtered_data$data[["4x_O2_01"]] == 
#                         filteredDF[["4x_O2_01"]],
#                     filtered_data$data[["1x_O3_01"]] == 
#                         filteredDF[["1x_O3_01"]],
#                     filtered_data$data[["2x_O3_01"]] == 
#                         filteredDF[["2x_O3_01"]],
#                     filtered_data$data[["4x_O3_01"]] == 
#                         filteredDF[["4x_O3_01"]],
#                     filtered_data$data[["1x_O1_02"]] == 
#                         filteredDF[["1x_O1_02"]],
#                     filtered_data$data[["2x_O1_02"]] == 
#                         filteredDF[["2x_O1_02"]],
#                     filtered_data$data[["4x_O1_02"]] == 
#                         filteredDF[["4x_O1_02"]],
#                     filtered_data$data[["1x_O2_02"]] == 
#                         filteredDF[["1x_O2_02"]],
#                     filtered_data$data[["2x_O2_02"]] == 
#                         filteredDF[["2x_O2_02"]],
#                     filtered_data$data[["4x_O2_02"]] == 
#                         filteredDF[["4x_O2_02"]],
#                     filtered_data$data[["1x_O3_02"]] == 
#                         filteredDF[["1x_O3_02"]],
#                     filtered_data$data[["2x_O3_02"]] == 
#                         filteredDF[["2x_O3_02"]],
#                     filtered_data$data[["4x_O3_02"]] == 
#                         filteredDF[["4x_O3_02"]]))
# })
# 
# test_that("2msFilter() == ms_filter()", {
#     expect_true(compareOldNew(filtered_data, filteredDF))
# })
# 
# test_that("msFilter(halfmin) == ms_filter(halfmin)", {
#     expect_true(all(hm_imputed_data$data[["1x_O1_01"]] == 
#                         hmImputedDF[["1x_O1_01"]],
#                     hm_imputed_data$data[["2x_O1_01"]] == 
#                         hmImputedDF[["2x_O1_01"]],
#                     hm_imputed_data$data[["4x_O1_01"]] == 
#                         hmImputedDF[["4x_O1_01"]],
#                     hm_imputed_data$data[["1x_O2_01"]] == 
#                         hmImputedDF[["1x_O2_01"]],
#                     hm_imputed_data$data[["2x_O2_01"]] == 
#                         hmImputedDF[["2x_O2_01"]],
#                     hm_imputed_data$data[["4x_O2_01"]] == 
#                         hmImputedDF[["4x_O2_01"]],
#                     hm_imputed_data$data[["1x_O3_01"]] == 
#                         hmImputedDF[["1x_O3_01"]],
#                     hm_imputed_data$data[["2x_O3_01"]] == 
#                         hmImputedDF[["2x_O3_01"]],
#                     hm_imputed_data$data[["4x_O3_01"]] == 
#                         hmImputedDF[["4x_O3_01"]],
#                     hm_imputed_data$data[["1x_O1_02"]] == 
#                         hmImputedDF[["1x_O1_02"]],
#                     hm_imputed_data$data[["2x_O1_02"]] == 
#                         hmImputedDF[["2x_O1_02"]],
#                     hm_imputed_data$data[["4x_O1_02"]] == 
#                         hmImputedDF[["4x_O1_02"]],
#                     hm_imputed_data$data[["1x_O2_02"]] == 
#                         hmImputedDF[["1x_O2_02"]],
#                     hm_imputed_data$data[["2x_O2_02"]] == 
#                         hmImputedDF[["2x_O2_02"]],
#                     hm_imputed_data$data[["4x_O2_02"]] == 
#                         hmImputedDF[["4x_O2_02"]],
#                     hm_imputed_data$data[["1x_O3_02"]] == 
#                         hmImputedDF[["1x_O3_02"]],
#                     hm_imputed_data$data[["2x_O3_02"]] == 
#                         hmImputedDF[["2x_O3_02"]],
#                     hm_imputed_data$data[["4x_O3_02"]] == 
#                         hmImputedDF[["4x_O3_02"]]))
# })
# 
# test_that("2msImpute(hm) == ms_impute(hm)", {
#     expect_true(compareOldNew(hm_imputed_data, hmImputedDF))
# })
# 
# ## Note: bpca needs to be rounded
# bpca_imputed_data$data <- bpca_imputed_data$data[, 3:20] %>%
#     mutate_all(round, digits = 5)
# bpcaImputedDF <- bpcaImputedDF[, 3:20] %>%
#     mutate_all(round, digits = 5)
# 
# test_that("msImpute(bpca) == ms_impute(bpca)", {
#     expect_true(all(bpca_imputed_data$data[["1x_O1_01"]] == 
#                         bpcaImputedDF[["1x_O1_01"]],
#                     bpca_imputed_data$data[["2x_O1_01"]] == 
#                         bpcaImputedDF[["2x_O1_01"]],
#                     bpca_imputed_data$data[["4x_O1_01"]] == 
#                         bpcaImputedDF[["4x_O1_01"]],
#                     bpca_imputed_data$data[["1x_O2_01"]] == 
#                         bpcaImputedDF[["1x_O2_01"]],
#                     bpca_imputed_data$data[["2x_O2_01"]] == 
#                         bpcaImputedDF[["2x_O2_01"]],
#                     bpca_imputed_data$data[["4x_O2_01"]] == 
#                         bpcaImputedDF[["4x_O2_01"]],
#                     bpca_imputed_data$data[["1x_O3_01"]] == 
#                         bpcaImputedDF[["1x_O3_01"]],
#                     bpca_imputed_data$data[["2x_O3_01"]] == 
#                         bpcaImputedDF[["2x_O3_01"]],
#                     bpca_imputed_data$data[["4x_O3_01"]] == 
#                         bpcaImputedDF[["4x_O3_01"]],
#                     bpca_imputed_data$data[["1x_O1_02"]] == 
#                         bpcaImputedDF[["1x_O1_02"]],
#                     bpca_imputed_data$data[["2x_O1_02"]] == 
#                         bpcaImputedDF[["2x_O1_02"]],
#                     bpca_imputed_data$data[["4x_O1_02"]] == 
#                         bpcaImputedDF[["4x_O1_02"]],
#                     bpca_imputed_data$data[["1x_O2_02"]] == 
#                         bpcaImputedDF[["1x_O2_02"]],
#                     bpca_imputed_data$data[["2x_O2_02"]] == 
#                         bpcaImputedDF[["2x_O2_02"]],
#                     bpca_imputed_data$data[["4x_O2_02"]] == 
#                         bpcaImputedDF[["4x_O2_02"]],
#                     bpca_imputed_data$data[["1x_O3_02"]] == 
#                         bpcaImputedDF[["1x_O3_02"]],
#                     bpca_imputed_data$data[["2x_O3_02"]] == 
#                         bpcaImputedDF[["2x_O3_02"]],
#                     bpca_imputed_data$data[["4x_O3_02"]] == 
#                         bpcaImputedDF[["4x_O3_02"]]))
# })
# 
# test_that("2msImpute(quant) == ms_impute(quant)", {
#     expect_true(compareOldNew(bpca_imputed_data, bpcaImputedDF))
# })
# 
# test_that("msImpute(knn) == ms_impute(knn)", {
#     expect_true(all(knn_imputed_data$data[["1x_O1_01"]] == 
#                         knnImputedDF[["1x_O1_01"]],
#                     knn_imputed_data$data[["2x_O1_01"]] == 
#                         knnImputedDF[["2x_O1_01"]],
#                     knn_imputed_data$data[["4x_O1_01"]] == 
#                         knnImputedDF[["4x_O1_01"]],
#                     knn_imputed_data$data[["1x_O2_01"]] == 
#                         knnImputedDF[["1x_O2_01"]],
#                     knn_imputed_data$data[["2x_O2_01"]] == 
#                         knnImputedDF[["2x_O2_01"]],
#                     knn_imputed_data$data[["4x_O2_01"]] == 
#                         knnImputedDF[["4x_O2_01"]],
#                     knn_imputed_data$data[["1x_O3_01"]] == 
#                         knnImputedDF[["1x_O3_01"]],
#                     knn_imputed_data$data[["2x_O3_01"]] == 
#                         knnImputedDF[["2x_O3_01"]],
#                     knn_imputed_data$data[["4x_O3_01"]] == 
#                         knnImputedDF[["4x_O3_01"]],
#                     knn_imputed_data$data[["1x_O1_02"]] == 
#                         knnImputedDF[["1x_O1_02"]],
#                     knn_imputed_data$data[["2x_O1_02"]] == 
#                         knnImputedDF[["2x_O1_02"]],
#                     knn_imputed_data$data[["4x_O1_02"]] == 
#                         knnImputedDF[["4x_O1_02"]],
#                     knn_imputed_data$data[["1x_O2_02"]] == 
#                         knnImputedDF[["1x_O2_02"]],
#                     knn_imputed_data$data[["2x_O2_02"]] == 
#                         knnImputedDF[["2x_O2_02"]],
#                     knn_imputed_data$data[["4x_O2_02"]] == 
#                         knnImputedDF[["4x_O2_02"]],
#                     knn_imputed_data$data[["1x_O3_02"]] == 
#                         knnImputedDF[["1x_O3_02"]],
#                     knn_imputed_data$data[["2x_O3_02"]] == 
#                         knnImputedDF[["2x_O3_02"]],
#                     knn_imputed_data$data[["4x_O3_02"]] == 
#                         knnImputedDF[["4x_O3_02"]]))
# })
# 
# test_that("2msImpute(knn) == ms_impute(knn)", {
#     expect_true(compareOldNew(knn_imputed_data, knnImputedDF))
# })
# 
# test_that("1msNorm(median) == ms_norm(median)", {
#     expect_true(all(med_normalized_data$data[["1x_O1_01"]] == 
#                         medianNormalized[["1x_O1_01"]],
#                     med_normalized_data$data[["2x_O1_01"]] == 
#                         medianNormalized[["2x_O1_01"]],
#                     med_normalized_data$data[["4x_O1_01"]] == 
#                         medianNormalized[["4x_O1_01"]],
#                     med_normalized_data$data[["1x_O2_01"]] == 
#                         medianNormalized[["1x_O2_01"]],
#                     med_normalized_data$data[["2x_O2_01"]] == 
#                         medianNormalized[["2x_O2_01"]],
#                     med_normalized_data$data[["4x_O2_01"]] == 
#                         medianNormalized[["4x_O2_01"]],
#                     med_normalized_data$data[["1x_O3_01"]] == 
#                         medianNormalized[["1x_O3_01"]],
#                     med_normalized_data$data[["2x_O3_01"]] == 
#                         medianNormalized[["2x_O3_01"]],
#                     med_normalized_data$data[["4x_O3_01"]] == 
#                         medianNormalized[["4x_O3_01"]],
#                     med_normalized_data$data[["1x_O1_02"]] == 
#                         medianNormalized[["1x_O1_02"]],
#                     med_normalized_data$data[["2x_O1_02"]] == 
#                         medianNormalized[["2x_O1_02"]],
#                     med_normalized_data$data[["4x_O1_02"]] == 
#                         medianNormalized[["4x_O1_02"]],
#                     med_normalized_data$data[["1x_O2_02"]] == 
#                         medianNormalized[["1x_O2_02"]],
#                     med_normalized_data$data[["2x_O2_02"]] == 
#                         medianNormalized[["2x_O2_02"]],
#                     med_normalized_data$data[["4x_O2_02"]] == 
#                         medianNormalized[["4x_O2_02"]],
#                     med_normalized_data$data[["1x_O3_02"]] == 
#                         medianNormalized[["1x_O3_02"]],
#                     med_normalized_data$data[["2x_O3_02"]] == 
#                         medianNormalized[["2x_O3_02"]],
#                     med_normalized_data$data[["4x_O3_02"]] == 
#                         medianNormalized[["4x_O3_02"]]))
# })
# 
# test_that("2msNorm(median) == ms_norm(median)", {
#           expect_true(compareOldNew(med_normalized_data, medianNormalized))
# })
# 
# test_that("1msNorm(median) == ms_norm(median)", {
#     expect_true(all(quant_normalized_data$data[["1x_O1_01"]] == 
#                         quantNormalized[["1x_O1_01"]],
#                     quant_normalized_data$data[["2x_O1_01"]] == 
#                         quantNormalized[["2x_O1_01"]],
#                     quant_normalized_data$data[["4x_O1_01"]] == 
#                         quantNormalized[["4x_O1_01"]],
#                     quant_normalized_data$data[["1x_O2_01"]] == 
#                         quantNormalized[["1x_O2_01"]],
#                     quant_normalized_data$data[["2x_O2_01"]] == 
#                         quantNormalized[["2x_O2_01"]],
#                     quant_normalized_data$data[["4x_O2_01"]] == 
#                         quantNormalized[["4x_O2_01"]],
#                     quant_normalized_data$data[["1x_O3_01"]] == 
#                         quantNormalized[["1x_O3_01"]],
#                     quant_normalized_data$data[["2x_O3_01"]] == 
#                         quantNormalized[["2x_O3_01"]],
#                     quant_normalized_data$data[["4x_O3_01"]] == 
#                         quantNormalized[["4x_O3_01"]],
#                     quant_normalized_data$data[["1x_O1_02"]] == 
#                         quantNormalized[["1x_O1_02"]],
#                     quant_normalized_data$data[["2x_O1_02"]] == 
#                         quantNormalized[["2x_O1_02"]],
#                     quant_normalized_data$data[["4x_O1_02"]] == 
#                         quantNormalized[["4x_O1_02"]],
#                     quant_normalized_data$data[["1x_O2_02"]] == 
#                         quantNormalized[["1x_O2_02"]],
#                     quant_normalized_data$data[["2x_O2_02"]] == 
#                         quantNormalized[["2x_O2_02"]],
#                     quant_normalized_data$data[["4x_O2_02"]] == 
#                         quantNormalized[["4x_O2_02"]],
#                     quant_normalized_data$data[["1x_O3_02"]] == 
#                         quantNormalized[["1x_O3_02"]],
#                     quant_normalized_data$data[["2x_O3_02"]] == 
#                         quantNormalized[["2x_O3_02"]],
#                     quant_normalized_data$data[["4x_O3_02"]] == 
#                         quantNormalized[["4x_O3_02"]]))
# })
# 
# test_that("2msNorm(quant) == ms_norm(quant)", {
#     expect_true(compareOldNew(quant_normalized_data, quantNormalized))
# })
# 
# ## Note: ComBat needs to be rounded
# combat_normalized_data$data <- combat_normalized_data$data[, 3:20] %>%
#     mutate_all(round, digits = 5)
# combatNormalized <- combatNormalized[, 3:20] %>%
#     mutate_all(round, digits = 5)
# 
# test_that("1msNorm(combat) == ms_norm(combat)", {
#     expect_true(all(combat_normalized_data$data[["1x_O1_01"]] == 
#                         combatNormalized[["1x_O1_01"]],
#                     combat_normalized_data$data[["2x_O1_01"]] == 
#                         combatNormalized[["2x_O1_01"]],
#                     combat_normalized_data$data[["4x_O1_01"]] == 
#                         combatNormalized[["4x_O1_01"]],
#                     combat_normalized_data$data[["1x_O2_01"]] == 
#                         combatNormalized[["1x_O2_01"]],
#                     combat_normalized_data$data[["2x_O2_01"]] == 
#                         combatNormalized[["2x_O2_01"]],
#                     combat_normalized_data$data[["4x_O2_01"]] == 
#                         combatNormalized[["4x_O2_01"]],
#                     combat_normalized_data$data[["1x_O3_01"]] == 
#                         combatNormalized[["1x_O3_01"]],
#                     combat_normalized_data$data[["2x_O3_01"]] == 
#                         combatNormalized[["2x_O3_01"]],
#                     combat_normalized_data$data[["4x_O3_01"]] == 
#                         combatNormalized[["4x_O3_01"]],
#                     combat_normalized_data$data[["1x_O1_02"]] == 
#                         combatNormalized[["1x_O1_02"]],
#                     combat_normalized_data$data[["2x_O1_02"]] == 
#                         combatNormalized[["2x_O1_02"]],
#                     combat_normalized_data$data[["4x_O1_02"]] == 
#                         combatNormalized[["4x_O1_02"]],
#                     combat_normalized_data$data[["1x_O2_02"]] == 
#                         combatNormalized[["1x_O2_02"]],
#                     combat_normalized_data$data[["2x_O2_02"]] == 
#                         combatNormalized[["2x_O2_02"]],
#                     combat_normalized_data$data[["4x_O2_02"]] == 
#                         combatNormalized[["4x_O2_02"]],
#                     combat_normalized_data$data[["1x_O3_02"]] == 
#                         combatNormalized[["1x_O3_02"]],
#                     combat_normalized_data$data[["2x_O3_02"]] == 
#                         combatNormalized[["2x_O3_02"]],
#                     combat_normalized_data$data[["4x_O3_02"]] == 
#                         combatNormalized[["4x_O3_02"]]))
# })
# 
# test_that("2msNorm(combat) == ms_norm(combat)", {
#     expect_true(compareOldNew(combat_normalized_data, combatNormalized))
# })
# 
# ## Note: SVA needs to be rounded
# sva_normalized_data$data <- sva_normalized_data$data[, 3:20] %>%
#     mutate_all(round, digits = 5)
# svaNormalized <- svaNormalized[, 3:20] %>%
#     mutate_all(round, digits = 5)
# 
# test_that("1msNorm(sva) == ms_norm(sva)", {
#     expect_true(all(sva_normalized_data$data[["1x_O1_01"]] == 
#                         svaNormalized[["1x_O1_01"]],
#                     sva_normalized_data$data[["2x_O1_01"]] == 
#                         svaNormalized[["2x_O1_01"]],
#                     sva_normalized_data$data[["4x_O1_01"]] == 
#                         svaNormalized[["4x_O1_01"]],
#                     sva_normalized_data$data[["1x_O2_01"]] == 
#                         svaNormalized[["1x_O2_01"]],
#                     sva_normalized_data$data[["2x_O2_01"]] == 
#                         svaNormalized[["2x_O2_01"]],
#                     sva_normalized_data$data[["4x_O2_01"]] == 
#                         svaNormalized[["4x_O2_01"]],
#                     sva_normalized_data$data[["1x_O3_01"]] == 
#                         svaNormalized[["1x_O3_01"]],
#                     sva_normalized_data$data[["2x_O3_01"]] == 
#                         svaNormalized[["2x_O3_01"]],
#                     sva_normalized_data$data[["4x_O3_01"]] == 
#                         svaNormalized[["4x_O3_01"]],
#                     sva_normalized_data$data[["1x_O1_02"]] == 
#                         svaNormalized[["1x_O1_02"]],
#                     sva_normalized_data$data[["2x_O1_02"]] == 
#                         svaNormalized[["2x_O1_02"]],
#                     sva_normalized_data$data[["4x_O1_02"]] == 
#                         svaNormalized[["4x_O1_02"]],
#                     sva_normalized_data$data[["1x_O2_02"]] == 
#                         svaNormalized[["1x_O2_02"]],
#                     sva_normalized_data$data[["2x_O2_02"]] == 
#                         svaNormalized[["2x_O2_02"]],
#                     sva_normalized_data$data[["4x_O2_02"]] == 
#                         svaNormalized[["4x_O2_02"]],
#                     sva_normalized_data$data[["1x_O3_02"]] == 
#                         svaNormalized[["1x_O3_02"]],
#                     sva_normalized_data$data[["2x_O3_02"]] == 
#                         svaNormalized[["2x_O3_02"]],
#                     sva_normalized_data$data[["4x_O3_02"]] == 
#                         svaNormalized[["4x_O3_02"]]))
# })
# 
# test_that("2msNorm(sva) == ms_norm(sva)", {
#     expect_true(compareOldNew(sva_normalized_data, svaNormalized))
# })
# 
# ## Note: ComBat needs to be rounded
# quant_comb_normalized_data$data <- quant_comb_normalized_data$data[, 3:20] %>%
#     mutate_all(round, digits = 5)
# quantCombatNormalized <- quantCombatNormalized[, 3:20] %>%
#     mutate_all(round, digits = 5)
# 
# test_that("1msNorm(quantComb) == ms_norm(quantComb)", {
#     expect_true(all(quant_comb_normalized_data$data[["1x_O1_01"]] == 
#                         quantCombatNormalized[["1x_O1_01"]],
#                     quant_comb_normalized_data$data[["2x_O1_01"]] == 
#                         quantCombatNormalized[["2x_O1_01"]],
#                     quant_comb_normalized_data$data[["4x_O1_01"]] == 
#                         quantCombatNormalized[["4x_O1_01"]],
#                     quant_comb_normalized_data$data[["1x_O2_01"]] == 
#                         quantCombatNormalized[["1x_O2_01"]],
#                     quant_comb_normalized_data$data[["2x_O2_01"]] == 
#                         quantCombatNormalized[["2x_O2_01"]],
#                     quant_comb_normalized_data$data[["4x_O2_01"]] == 
#                         quantCombatNormalized[["4x_O2_01"]],
#                     quant_comb_normalized_data$data[["1x_O3_01"]] == 
#                         quantCombatNormalized[["1x_O3_01"]],
#                     quant_comb_normalized_data$data[["2x_O3_01"]] == 
#                         quantCombatNormalized[["2x_O3_01"]],
#                     quant_comb_normalized_data$data[["4x_O3_01"]] == 
#                         quantCombatNormalized[["4x_O3_01"]],
#                     quant_comb_normalized_data$data[["1x_O1_02"]] == 
#                         quantCombatNormalized[["1x_O1_02"]],
#                     quant_comb_normalized_data$data[["2x_O1_02"]] == 
#                         quantCombatNormalized[["2x_O1_02"]],
#                     quant_comb_normalized_data$data[["4x_O1_02"]] == 
#                         quantCombatNormalized[["4x_O1_02"]],
#                     quant_comb_normalized_data$data[["1x_O2_02"]] == 
#                         quantCombatNormalized[["1x_O2_02"]],
#                     quant_comb_normalized_data$data[["2x_O2_02"]] == 
#                         quantCombatNormalized[["2x_O2_02"]],
#                     quant_comb_normalized_data$data[["4x_O2_02"]] == 
#                         quantCombatNormalized[["4x_O2_02"]],
#                     quant_comb_normalized_data$data[["1x_O3_02"]] == 
#                         quantCombatNormalized[["1x_O3_02"]],
#                     quant_comb_normalized_data$data[["2x_O3_02"]] == 
#                         quantCombatNormalized[["2x_O3_02"]],
#                     quant_comb_normalized_data$data[["4x_O3_02"]] == 
#                         quantCombatNormalized[["4x_O3_02"]]))
# })
# 
# test_that("2msNorm(quantComb) == ms_norm(quantComb)", {
#     expect_true(compareOldNew(quant_comb_normalized_data, 
#                               quantCombatNormalized))
# })
# 
# ## Note: ComBat needs to be rounded
# med_comb_normalized_data$data <- med_comb_normalized_data$data[, 3:20] %>%
#     mutate_all(round, digits = 5)
# medCombatNormalized <- medCombatNormalized[, 3:20] %>%
#     mutate_all(round, digits = 5)
# 
# test_that("1msNorm(medComb) == ms_norm(medComb)", {
#     expect_true(all(med_comb_normalized_data$data[["1x_O1_01"]] == 
#                         medCombatNormalized[["1x_O1_01"]],
#                     med_comb_normalized_data$data[["2x_O1_01"]] == 
#                         medCombatNormalized[["2x_O1_01"]],
#                     med_comb_normalized_data$data[["4x_O1_01"]] == 
#                         medCombatNormalized[["4x_O1_01"]],
#                     med_comb_normalized_data$data[["1x_O2_01"]] == 
#                         medCombatNormalized[["1x_O2_01"]],
#                     med_comb_normalized_data$data[["2x_O2_01"]] == 
#                         medCombatNormalized[["2x_O2_01"]],
#                     med_comb_normalized_data$data[["4x_O2_01"]] == 
#                         medCombatNormalized[["4x_O2_01"]],
#                     med_comb_normalized_data$data[["1x_O3_01"]] == 
#                         medCombatNormalized[["1x_O3_01"]],
#                     med_comb_normalized_data$data[["2x_O3_01"]] == 
#                         medCombatNormalized[["2x_O3_01"]],
#                     med_comb_normalized_data$data[["4x_O3_01"]] == 
#                         medCombatNormalized[["4x_O3_01"]],
#                     med_comb_normalized_data$data[["1x_O1_02"]] == 
#                         medCombatNormalized[["1x_O1_02"]],
#                     med_comb_normalized_data$data[["2x_O1_02"]] == 
#                         medCombatNormalized[["2x_O1_02"]],
#                     med_comb_normalized_data$data[["4x_O1_02"]] == 
#                         medCombatNormalized[["4x_O1_02"]],
#                     med_comb_normalized_data$data[["1x_O2_02"]] == 
#                         medCombatNormalized[["1x_O2_02"]],
#                     med_comb_normalized_data$data[["2x_O2_02"]] == 
#                         medCombatNormalized[["2x_O2_02"]],
#                     med_comb_normalized_data$data[["4x_O2_02"]] == 
#                         medCombatNormalized[["4x_O2_02"]],
#                     med_comb_normalized_data$data[["1x_O3_02"]] == 
#                         medCombatNormalized[["1x_O3_02"]],
#                     med_comb_normalized_data$data[["2x_O3_02"]] == 
#                         medCombatNormalized[["2x_O3_02"]],
#                     med_comb_normalized_data$data[["4x_O3_02"]] == 
#                         medCombatNormalized[["4x_O3_02"]]))
# })
# 
# test_that("2msNorm(medComb) == ms_norm(medComb)", {
#     expect_true(compareOldNew(med_comb_normalized_data, 
#                               medCombatNormalized))
# })
# 
# test_that("1msNorm(RUV) == ms_norm(RUV)", {
#     expect_true(all(ruv_normalized_data$data[["1x_O1_01"]] == 
#                         ruvNormalized[["1x_O1_01"]],
#                     ruv_normalized_data$data[["2x_O1_01"]] == 
#                         ruvNormalized[["2x_O1_01"]],
#                     ruv_normalized_data$data[["4x_O1_01"]] == 
#                         ruvNormalized[["4x_O1_01"]],
#                     ruv_normalized_data$data[["1x_O2_01"]] == 
#                         ruvNormalized[["1x_O2_01"]],
#                     ruv_normalized_data$data[["2x_O2_01"]] == 
#                         ruvNormalized[["2x_O2_01"]],
#                     ruv_normalized_data$data[["4x_O2_01"]] == 
#                         ruvNormalized[["4x_O2_01"]],
#                     ruv_normalized_data$data[["1x_O3_01"]] == 
#                         ruvNormalized[["1x_O3_01"]],
#                     ruv_normalized_data$data[["2x_O3_01"]] == 
#                         ruvNormalized[["2x_O3_01"]],
#                     ruv_normalized_data$data[["4x_O3_01"]] == 
#                         ruvNormalized[["4x_O3_01"]],
#                     ruv_normalized_data$data[["1x_O1_02"]] == 
#                         ruvNormalized[["1x_O1_02"]],
#                     ruv_normalized_data$data[["2x_O1_02"]] == 
#                         ruvNormalized[["2x_O1_02"]],
#                     ruv_normalized_data$data[["4x_O1_02"]] == 
#                         ruvNormalized[["4x_O1_02"]],
#                     ruv_normalized_data$data[["1x_O2_02"]] == 
#                         ruvNormalized[["1x_O2_02"]],
#                     ruv_normalized_data$data[["2x_O2_02"]] == 
#                         ruvNormalized[["2x_O2_02"]],
#                     ruv_normalized_data$data[["4x_O2_02"]] == 
#                         ruvNormalized[["4x_O2_02"]],
#                     ruv_normalized_data$data[["1x_O3_02"]] == 
#                         ruvNormalized[["1x_O3_02"]],
#                     ruv_normalized_data$data[["2x_O3_02"]] == 
#                         ruvNormalized[["2x_O3_02"]],
#                     ruv_normalized_data$data[["4x_O3_02"]] == 
#                         ruvNormalized[["4x_O3_02"]]))
# })
# 
# test_that("2msNorm(RUV) == ms_norm(RUV)", {
#     expect_true(compareOldNew(ruv_normalized_data, 
#                               ruvNormalized))
# })
# 
# test_that("1msNorm(CRMN) == ms_norm(CRMN)", {
#     expect_true(all(crmn_normalized_data$data[["1x_O1_01"]] == 
#                         crmnNormalized[["1x_O1_01"]],
#                     crmn_normalized_data$data[["2x_O1_01"]] == 
#                         crmnNormalized[["2x_O1_01"]],
#                     crmn_normalized_data$data[["4x_O1_01"]] == 
#                         crmnNormalized[["4x_O1_01"]],
#                     crmn_normalized_data$data[["1x_O2_01"]] == 
#                         crmnNormalized[["1x_O2_01"]],
#                     crmn_normalized_data$data[["2x_O2_01"]] == 
#                         crmnNormalized[["2x_O2_01"]],
#                     crmn_normalized_data$data[["4x_O2_01"]] == 
#                         crmnNormalized[["4x_O2_01"]],
#                     crmn_normalized_data$data[["1x_O3_01"]] == 
#                         crmnNormalized[["1x_O3_01"]],
#                     crmn_normalized_data$data[["2x_O3_01"]] == 
#                         crmnNormalized[["2x_O3_01"]],
#                     crmn_normalized_data$data[["4x_O3_01"]] == 
#                         crmnNormalized[["4x_O3_01"]],
#                     crmn_normalized_data$data[["1x_O1_02"]] == 
#                         crmnNormalized[["1x_O1_02"]],
#                     crmn_normalized_data$data[["2x_O1_02"]] == 
#                         crmnNormalized[["2x_O1_02"]],
#                     crmn_normalized_data$data[["4x_O1_02"]] == 
#                         crmnNormalized[["4x_O1_02"]],
#                     crmn_normalized_data$data[["1x_O2_02"]] == 
#                         crmnNormalized[["1x_O2_02"]],
#                     crmn_normalized_data$data[["2x_O2_02"]] == 
#                         crmnNormalized[["2x_O2_02"]],
#                     crmn_normalized_data$data[["4x_O2_02"]] == 
#                         crmnNormalized[["4x_O2_02"]],
#                     crmn_normalized_data$data[["1x_O3_02"]] == 
#                         crmnNormalized[["1x_O3_02"]],
#                     crmn_normalized_data$data[["2x_O3_02"]] == 
#                         crmnNormalized[["2x_O3_02"]],
#                     crmn_normalized_data$data[["4x_O3_02"]] == 
#                         crmnNormalized[["4x_O3_02"]]))
# })
# 
# test_that("2msNorm(CRMN) == ms_norm(CRMN)", {
#     expect_true(compareOldNew(crmn_normalized_data, 
#                               crmnNormalized))
# })
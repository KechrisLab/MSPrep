context("msImpute()")

set.seed(9999)
data(msquant)

preppedData <- msSummarize(msquant,
                           compVars = c("mz", "rt"),
                           sampleVars = c("spike", "batch", "replicate",
                                          "subject_id"),
                           colExtraText = "Neutral_Operator_Dif_Pos_",
                           separator = "_",
                           missingValue = 1) %>%
    msFilter(compVars = c("mz", "rt"),
             sampleVars = c("spike", "batch", "subject_id"),
             separator = "_")

hmImputedDF <- msImpute(preppedData, imputeMethod = "halfmin",
                        compVars = c("mz", "rt"),
                        sampleVars = c("spike", "batch", "subject_id"),
                        separator = "_",
                        returnToSE = FALSE,
                        missingValue = 0)

bpcaImputedDF <- msImpute(preppedData, imputeMethod = "bpca",
                          compVars = c("mz", "rt"),
                          sampleVars = c("spike", "batch", "subject_id"),
                          separator = "_",
                          returnToSE = FALSE,
                          missingValue = 0)

knnImputedDF <- msImpute(preppedData, imputeMethod = "knn",
                         compVars = c("mz", "rt"),
                         sampleVars = c("spike", "batch", "subject_id"),
                         separator = "_",
                         returnToSE = FALSE,
                         missingValue = 0)

# test_that("Check dataset columns", {
# 
#     halfmin_colnames <- sort(colnames(halfmin_imputed$data))
#     knn_colnames     <- sort(colnames(knn_imputed$data))
#     bpca_colnames    <- sort(colnames(bpca_imputed$data))
#     halfmin_extravars <- c(grouping_vars(halfmin_imputed),
#                            batch_var(halfmin_imputed))
#     knn_extravars     <- c(grouping_vars(knn_imputed), 
#                            batch_var(knn_imputed))
#     bpca_extravars    <- c(grouping_vars(bpca_imputed), 
#                            batch_var(bpca_imputed))
# 
#     expect_equal(halfmin_colnames, valid_cols(halfmin_extravars))
#     expect_equal(knn_colnames, valid_cols(knn_extravars))
#     expect_equal(bpca_colnames, valid_cols(bpca_extravars))
# 
# })

compVars <- c("mz", "rt")

# Convert abundance data to vector of values
abundanceVec <- select(preppedData, -all_of(compVars)) %>%
    as.matrix() %>%
    c()
hmImputedVec <- select(hmImputedDF, -all_of(compVars)) %>%
    as.matrix() %>%
    c()
bpcaImputedVec <- select(bpcaImputedDF, -all_of(compVars)) %>%
    as.matrix() %>%
    c()
knnImputedVec <- select(knnImputedDF, -all_of(compVars)) %>%
    as.matrix() %>%
    c()

# Get values which were originally missing
missingVals <- which(abundanceVec == 0)

test_that("Check imputations",{
    expect_true(all(hmImputedVec != 0))
    expect_true(all(bpcaImputedVec != 0))
    expect_true(all(knnImputedVec != 0))
    expect_true(sum(hmImputedVec[missingVals] == 
                        bpcaImputedVec[missingVals]) == 8)
})

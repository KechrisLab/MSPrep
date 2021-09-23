context("equality between DF and SE functions")

set.seed(9997)

data(msquant)

se <- .dfToSE(msquant,
              compVars = c("mz", "rt"),
              sampleVars = c("spike", "batch", "replicate", "subject_id"),
              colExtraText = "Neutral_Operator_Dif_Pos_",
              separator = "_")

metadata(se) <- list(test1 = "test",
                     test2 = matrix(c(1, 2, 3, 4, 5, 6), nrow = 3),
                     test3 = data.frame(test1 = c(1, 2, 3),
                                        test2 = c(4, 5, 6)),
                     test4 = TRUE)

## Run msTidy with SE and data frame
tidySE <- .msTidy(se, missingValue = 1)

tidyDF <- .msTidy(msquant,
                 compVars = c("mz", "rt"),
                 sampleVars = c("spike", "batch", "replicate", "subject_id"),
                 colExtraText = "Neutral_Operator_Dif_Pos_",
                 separator = "_",
                 missingValue = 1)

## Run msSummarize with SE and data frame
summarizedDF <- msSummarize(msquant,
                            compVars = c("mz", "rt"),
                            sampleVars = c("spike", "batch", "replicate", 
                                           "subject_id"),
                            cvMax = 0.50,
                            minPropPresent = 1/3,
                            returnSummaryDetails = TRUE,
                            colExtraText = "Neutral_Operator_Dif_Pos_",
                            separator = "_",
                            missingValue = 1)

summarizedSE <- msSummarize(se,
                            cvMax = 0.50,
                            minPropPresent = 1/3,
                            missingValue = 1,
                            returnSummaryDetails = TRUE)

## Run msFilter with SE and data frame
filteredDF <- msFilter(summarizedDF$data,
                       filterPercent = 0.8,
                       compVars = c("mz", "rt"),
                       sampleVars = c("spike", "batch", "subject_id"),
                       separator = "_",
                       returnToSE = FALSE)

filteredSE <- msFilter(summarizedSE,
                       filterPercent = .8)

## Run msImpute with DF and SE
hmImputedDF <- msImpute(filteredDF, imputeMethod = "halfmin",
                        compVars = c("mz", "rt"),
                        sampleVars = c("spike", "batch", "subject_id"),
                        separator = "_",
                        returnToSE = FALSE,
                        missingValue = 0)
bpcaImputedDF <- msImpute(filteredDF, imputeMethod = "bpca",
                          compVars = c("mz", "rt"),
                          sampleVars = c("spike", "batch", "subject_id"),
                          separator = "_",
                          returnToSE = FALSE,
                          missingValue = 0)
knnImputedDF <- msImpute(filteredDF, imputeMethod = "knn",
                         compVars = c("mz", "rt"),
                         sampleVars = c("spike", "batch", "subject_id"),
                         separator = "_",
                         returnToSE = FALSE,
                         missingValue = 0)
set.seed(123)
rfImputedDF <- msImpute(filteredDF, imputeMethod = "rf",
                        compVars = c("mz", "rt"),
                        sampleVars = c("spike", "batch", "subject_id"),
                        separator = "_",
                        returnToSE = FALSE,
                        missingValue = 0)

hmImputedSE <- msImpute(filteredSE, imputeMethod = "halfmin", 
                        returnToSE = TRUE,
                        missingValue = 0)
bpcaImputedSE <- msImpute(filteredSE, imputeMethod = "bpca",
                          returnToSE = TRUE,
                          missingValue = 0)
knnImputedSE <- msImpute(filteredSE, imputeMethod = "knn",
                         returnToSE = TRUE,
                         missingValue = 0)
set.seed(123)
rfImputedSE <- msImpute(filteredSE, imputeMethod = "rf",
                        returnToSE = TRUE,
                        missingValue = 0)

## Run msNormalize with DF and SE
svaNormalizedDF <- msNormalize(hmImputedDF,
                                 normalizeMethod = "SVA",
                                 compVars = c("mz", "rt"),
                                 sampleVars = c("spike", "batch", "subject_id"),
                                 covariatesOfInterest = c("spike"),
                                 separator = "_",
                                 returnToSE = FALSE)
combatNormalizedDF <- msNormalize(hmImputedDF, normalizeMethod = "ComBat", 
                                compVars = c("mz", "rt"),
                                sampleVars = c("spike", "batch", "subject_id"),
                                covariatesOfInterest = c("spike"),
                                separator = "_",
                                returnToSE = FALSE)
quantNormalizedDF <- msNormalize(hmImputedDF, normalizeMethod = "quantile",
                                 compVars = c("mz", "rt"),
                                 sampleVars = c("spike", "batch", "subject_id"),
                                 separator = "_",
                                 returnToSE = FALSE)
medianNormalizedDF <- msNormalize(hmImputedDF, normalizeMethod = "median",
                                compVars = c("mz", "rt"),
                                sampleVars = c("spike", "batch", "subject_id"),
                                separator = "_",
                                returnToSE = FALSE)
quantCombatNormalizedDF <- msNormalize(hmImputedDF,
                                       normalizeMethod = "quantile + ComBat",
                                       compVars = c("mz", "rt"),
                                       sampleVars = c("spike", "batch", 
                                                      "subject_id"),
                                       covariatesOfInterest = c("spike"),
                                       separator = "_",
                                       returnToSE = FALSE)
medCombatNormalizedDF <- msNormalize(hmImputedDF,
                                   normalizeMethod = "median + ComBat",
                                   compVars = c("mz", "rt"),
                                   sampleVars = c("spike", "batch", 
                                                  "subject_id"),
                                   covariatesOfInterest = c("spike"),
                                   separator = "_",
                                   returnToSE = FALSE)
crmnNormalizedDF <- msNormalize(hmImputedDF, 
                                normalizeMethod = "CRMN",
                                compVars = c("mz", "rt"),
                                sampleVars = c("spike", "batch", "subject_id"),
                                covariatesOfInterest = c("spike"),
                                separator = "_",
                                returnToSE = FALSE)
ruvNormalizedDF <- msNormalize(hmImputedDF, 
                               normalizeMethod = "RUV",
                               compVars = c("mz", "rt"),
                               sampleVars = c("spike", "batch", "subject_id"),
                               covariatesOfInterest = c("spike"),
                               controls = NULL,  nControl = 5, kRUV = 5,
                               separator = "_",
                               returnToSE = FALSE)

svaNormalizedSE <- msNormalize(hmImputedSE,
                               normalizeMethod = "SVA",
                               covariatesOfInterest = c("spike"))
combatNormalizedSE <- msNormalize(hmImputedSE, normalizeMethod = "ComBat",
                                  covariatesOfInterest = c("spike"),
                                  returnToSE = TRUE)
quantileNormalizedSE <- msNormalize(hmImputedSE, normalizeMethod = "quantile",
                                    returnToSE = TRUE)
medNormalizedSE <- msNormalize(hmImputedSE, normalizeMethod = "median",
                               returnToSE = TRUE)
quantCombatNormalizedSE <- msNormalize(hmImputedSE,
                                       normalizeMethod = "quantile + ComBat",
                                       covariatesOfInterest = c("spike"),
                                       returnToSE = TRUE)
medCombatNormalizedSE <- msNormalize(hmImputedSE,
                                     normalizeMethod = "median + ComBat",
                                     covariatesOfInterest = c("spike"),
                                     returnToSE = TRUE)
crmnNormalizedSE <- msNormalize(hmImputedSE, 
                                normalizeMethod = "CRMN",
                                covariatesOfInterest = c("spike"),
                                returnToSE = TRUE)
ruvNormalizedSE <- msNormalize(hmImputedSE, 
                               normalizeMethod = "RUV",
                               controls = NULL,  nControl = 5, kRUV = 5,
                               separator = "_",
                               returnToSE = TRUE)


## Check that results are equal
## msTidy()
test_that("Check .msTidy()", {
    expect_true(all(tidySE == tidyDF | is.na(tidySE == tidyDF)))
})

## msSummarize()
test_that("Check msSummarize()", {
    expect_true(identical(assay(summarizedSE), 
                          as.matrix(summarizedDF$data[, 3:20])))
    expect_true(all(all(rowData(summarizedSE) == 
                            as.data.frame(summarizedDF$data[, 1:2]))))
    expect_true(identical(S4Vectors::metadata(summarizedSE)$summaryDetails,
                          summarizedDF$summaryDetails))
})

## msFilter()
test_that("Check msFilter()", {
    expect_true(identical(assay(filteredSE), 
                          as.matrix(filteredDF[, 3:20])))
    expect_true(all(all(rowData(filteredSE) == 
                            as.data.frame(filteredDF[, 1:2]))))
    expect_true(all(rownames(colData(filteredSE)) == 
                        colnames(filteredDF[3:20])))
})

## msImpute()
test_that("Check msImpute(hm)", {
    expect_true(identical(assay(hmImputedSE), 
                          as.matrix(hmImputedDF[, 3:20])))
    expect_true(all(all(rowData(hmImputedSE) == 
                            as.data.frame(hmImputedDF[, 1:2]))))
    expect_true(all(rownames(colData(hmImputedSE)) == 
                        colnames(hmImputedDF[3:20])))
})
test_that("Check msImpute(bpca)", {
    expect_true(identical(assay(bpcaImputedSE), 
                          as.matrix(bpcaImputedDF[, 3:20])))
    expect_true(all(all(rowData(bpcaImputedSE) == 
                            as.data.frame(bpcaImputedDF[, 1:2]))))
    expect_true(all(rownames(colData(bpcaImputedSE)) == 
                        colnames(bpcaImputedDF[3:20])))
})
test_that("Check msImpute(knn)", {
    expect_true(identical(assay(knnImputedSE), 
                          as.matrix(knnImputedDF[, 3:20])))
    expect_true(all(all(rowData(knnImputedSE) == 
                            as.data.frame(knnImputedDF[, 1:2]))))
    expect_true(all(rownames(colData(knnImputedSE)) == 
                        colnames(knnImputedDF[3:20])))
})
test_that("Check msImpute(rf)", {
    expect_true(identical(assay(rfImputedSE), 
                          as.matrix(rfImputedDF[, 3:20])))
    expect_true(all(all(rowData(rfImputedSE) == 
                            as.data.frame(rfImputedDF[, 1:2]))))
    expect_true(all(rownames(colData(rfImputedSE)) == 
                        colnames(rfImputedDF[3:20])))
})

## msNormalize()
test_that("Check msNormalize(sva)", {
    expect_true(identical(assay(svaNormalizedSE), 
                          as.matrix(svaNormalizedDF[, 3:20])))
    expect_true(all(all(rowData(svaNormalizedSE) == 
                            as.data.frame(svaNormalizedDF[, 1:2]))))
    expect_true(all(rownames(colData(svaNormalizedSE)) == 
                        colnames(svaNormalizedDF[3:20])))
})
test_that("Check msNormalize(combat)", {
    expect_true(identical(assay(combatNormalizedSE), 
                          as.matrix(combatNormalizedDF[, 3:20])))
    expect_true(all(all(rowData(combatNormalizedSE) == 
                            as.data.frame(combatNormalizedDF[, 1:2]))))
    expect_true(all(rownames(colData(combatNormalizedSE)) == 
                        colnames(combatNormalizedDF[3:20])))
})
test_that("Check msNormalize(quantile)", {
    expect_true(identical(assay(quantileNormalizedSE), 
                          as.matrix(quantNormalizedDF[, 3:20])))
    expect_true(all(all(rowData(quantileNormalizedSE) == 
                            as.data.frame(quantNormalizedDF[, 1:2]))))
    expect_true(all(rownames(colData(quantileNormalizedSE)) == 
                        colnames(quantNormalizedDF[3:20])))
})
test_that("Check msNormalize(median)", {
    expect_true(identical(assay(medNormalizedSE), 
                          as.matrix(medianNormalizedDF[, 3:20])))
    expect_true(all(all(rowData(medNormalizedSE) == 
                            as.data.frame(medianNormalizedDF[, 1:2]))))
    expect_true(all(rownames(colData(medNormalizedSE)) == 
                        colnames(medianNormalizedDF[3:20])))
})
test_that("Check msNormalize(quantile + ComBat)", {
    expect_true(identical(assay(quantCombatNormalizedSE), 
                          as.matrix(quantCombatNormalizedDF[, 3:20])))
    expect_true(all(all(rowData(quantCombatNormalizedSE) == 
                            as.data.frame(quantCombatNormalizedDF[, 1:2]))))
    expect_true(all(rownames(colData(quantCombatNormalizedSE)) == 
                        colnames(quantCombatNormalizedDF[3:20])))
})
test_that("Check msNormalize(median + ComBat)", {
    expect_true(identical(assay(medCombatNormalizedSE), 
                          as.matrix(medCombatNormalizedDF[, 3:20])))
    expect_true(all(all(rowData(medCombatNormalizedSE) == 
                            as.data.frame(medCombatNormalizedDF[, 1:2]))))
    expect_true(all(rownames(colData(medCombatNormalizedSE)) == 
                        colnames(medCombatNormalizedDF[3:20])))
})
test_that("Check msNormalize(crmn)", {
    expect_true(all(assay(crmnNormalizedSE) ==
                          as.matrix(crmnNormalizedDF[, 3:20])))
    expect_true(all(all(rowData(crmnNormalizedSE) == 
                            as.data.frame(crmnNormalizedDF[, 1:2]))))
    expect_true(all(rownames(colData(crmnNormalizedSE)) == 
                        colnames(crmnNormalizedDF[3:20])))
})
test_that("Check msNormalize(ruv)", {
    expect_true(identical(assay(ruvNormalizedSE), 
                          as.matrix(ruvNormalizedDF[, 3:20])))
    expect_true(all(all(rowData(ruvNormalizedSE) == 
                            as.data.frame(ruvNormalizedDF[, 1:2]))))
    expect_true(all(rownames(colData(ruvNormalizedSE)) == 
                        colnames(ruvNormalizedDF[3:20])))
})


context("msNormalize()")

summarizedDF <- msSummarize(msquant,
                            compVars = c("mz", "rt"),
                            sampleVars = c("spike", "batch", "replicate", 
                                           "subject_id"),
                            cvMax = 0.50,
                            minPropPresent = 1/3,
                            returnSummaryDetails = FALSE,
                            colExtraText = "Neutral_Operator_Dif_Pos_",
                            separator = "_",
                            missingValue = 1,
                            returnToSE = FALSE)

filteredDF <- msFilter(summarizedDF,
                       filterPercent = 0.8,
                       compVars = c("mz", "rt"),
                       sampleVars = c("spike", "batch", "subject_id"),
                       separator = "_",
                       returnToSE = FALSE)

hmImputedDF <- msImpute(filteredDF, imputeMethod = "halfmin",
                        compVars = c("mz", "rt"),
                        sampleVars = c("spike", "batch", "subject_id"),
                        separator = "_",
                        returnToSE = FALSE,
                        missingValue = 0)

medianNormalizedDF <- msNormalize(hmImputedDF, normalizeMethod = "median",
                                  compVars = c("mz", "rt"),
                                  sampleVars = c("spike", "batch", 
                                                 "subject_id"),
                                  separator = "_",
                                  returnToSE = FALSE)
quantNormalizedDF <- msNormalize(hmImputedDF, normalizeMethod = "quantile",
                                 compVars = c("mz", "rt"),
                                 sampleVars = c("spike", "batch", 
                                                "subject_id"),
                                 separator = "_",
                                 returnToSE = FALSE)

quantiles <- c(0, .25, .5, .75, 1)

quantilesDF <- summarise_all(quantNormalizedDF[3:20],
                             quantile, probs = quantiles)

test_that("msNormalize(quantile)", {
    expect_true(length(unique(c(quantilesDF[1, ]))) == 1)
    expect_true(length(unique(c(quantilesDF[2, ]))) == 1)
    expect_true(length(unique(c(quantilesDF[3, ]))) == 1)
    expect_true(length(unique(c(quantilesDF[4, ]))) == 1)
    expect_true(length(unique(c(quantilesDF[5, ]))) == 1)
})

mediansDF <- summarise_all(medianNormalizedDF[3:20], median)

test_that("msNormalize(median)", {
    expect_true(all(as.numeric(mediansDF[1, ]) == 0))
})


context("Parameter Validation")

data <- msquant

SE <- .dfToSE(msquant,
              compVars = c("mz", "rt"),
              sampleVars = c("spike", "batch", "replicate", "subject_id"),
              colExtraText = "Neutral_Operator_Dif_Pos_",
              separator = "_")

test_that(".dfParamValidation", {
    
    colnames(data)[3] <- "Neutral_Operator_Dif_Pos_1x_O1_A_01_extra"
    
    expect_error(.dfParamValidation(data, compVars = c("mz", "rt"),
                                    sampleVars = c("spike", "batch", 
                                                   "replicate", "subject_id"),
                                    colExtraText = "Neutral_Operator_Dif_Pos_",
                                    separator = "_"))
    
    colnames(data)[3] <- "Neutral_Operator_Dif_Pos_1x_O1_A_01"
    
    expect_error(.dfParamValidation(data, compVars = c("mz", "rt"),
                                    sampleVars = c("spike", "batch", 
                                                   "replicate"),
                                    colExtraText = "Neutral_Operator_Dif_Pos_",
                                    separator = "_"))
    
    expect_error(.dfParamValidation(data, compVars = c("mz", "rt", "test"),
                                    sampleVars = c("spike", "batch", 
                                                   "replicate", "subject_id"),
                                    colExtraText = "Neutral_Operator_Dif_Pos_",
                                    separator = "_"))
})

test_that(".summarizeParamValidation", {
    
    colnames(colData(SE))[3] <- "test"
    
    expect_error(msSummarize(SE, cvMax = 0.50, minPropPresent = 1/3, 
                             missingValue = 1))
    
    colnames(colData(SE))[3] <- "replicate"
    
    expect_error(msSummarize(SE, cvMax = 0.50, minPropPresent = 5, 
                             missingValue = 1))
    
    expect_error(msSummarize(SE, cvMax = 0.50, minPropPresent = 5, 
                             missingValue = 1, returnToDF = TRUE,
                             returnToSE = TRUE))
})

test_that(".filterParamValidation", {
    
    expect_error(msFilter(summarizedDF,
                          filterPercent = 2,
                          compVars = c("mz", "rt"),
                          sampleVars = c("spike", "batch", "subject_id"),
                          separator = "_",
                          returnToSE = FALSE))
    
    expect_error(msFilter(summarizedSE,
                          filterPercent = -2))
    
})


context("msFilter()")

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

filterPercent1 <- 0.8
filterPercent2 <- 0.5
filterPercent3 <- 0.3

filteredDF1 <- msFilter(summarizedDF,
                       filterPercent = filterPercent1,
                       compVars = c("mz", "rt"),
                       sampleVars = c("spike", "batch", "subject_id"),
                       separator = "_",
                       returnToSE = FALSE)
filteredDF2 <- msFilter(summarizedDF,
                        filterPercent = filterPercent2,
                        compVars = c("mz", "rt"),
                        sampleVars = c("spike", "batch", "subject_id"),
                        separator = "_",
                        returnToSE = FALSE)
filteredDF3 <- msFilter(summarizedDF,
                        filterPercent = filterPercent3,
                        compVars = c("mz", "rt"),
                        sampleVars = c("spike", "batch", "subject_id"),
                        separator = "_",
                        returnToSE = FALSE)

tidyData1 <- msTidy(filteredDF1,
                 compVars = c("mz", "rt"),
                 sampleVars = c("spike", "batch", "subject_id"),
                 separator = "_",
                 missingValue = 0,
                 setMissing = 0)
tidyData2 <- msTidy(filteredDF2,
                    compVars = c("mz", "rt"),
                    sampleVars = c("spike", "batch", "subject_id"),
                    separator = "_",
                    missingValue = 0,
                    setMissing = 0)
tidyData3 <- msTidy(filteredDF3,
                    compVars = c("mz", "rt"),
                    sampleVars = c("spike", "batch", "subject_id"),
                    separator = "_",
                    missingValue = 0,
                    setMissing = 0)

compVarsSyms <- syms(c("mz", "rt"))

filterStatus1 <- group_by(tidyData1, `!!!`(compVarsSyms)) %>%
    summarise(percPresent = sum(.data$abundance != 0) / n())
filterStatus2 <- group_by(tidyData2, `!!!`(compVarsSyms)) %>%
    summarise(percPresent = sum(.data$abundance != 0) / n())
filterStatus3 <- group_by(tidyData3, `!!!`(compVarsSyms)) %>%
    summarise(percPresent = sum(.data$abundance != 0) / n())

test_that("msFilter()", {
    expect_true(all(filterStatus1$percPresent > filterPercent1))
    expect_true(all(filterStatus2$percPresent > filterPercent2))
    expect_true(all(filterStatus3$percPresent > filterPercent3))
})


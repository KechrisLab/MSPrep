context(".msTidy() and return functions")

data(msquant)

SE <- .dfToSE(msquant,
              compVars = c("mz", "rt"),
              sampleVars = c("spike", "batch", "replicate", "subject_id"),
              colExtraText = "Neutral_Operator_Dif_Pos_",
              separator = "_")

metadata(SE) <- list(test1 = "test",
                     test2 = matrix(c(1, 2, 3, 4, 5, 6), nrow = 3),
                     test3 = data.frame(test1 = c(1, 2, 3),
                                        test2 = c(4, 5, 6)),
                     test4 = TRUE)

tidySE <- .msTidy(SE, missingValue = 1, setMissing = 1)

returnedSE <- .tidyReturn(tidySE,
                        compVars = c("mz", "rt"),
                        sampleVars = c("spike", "batch", "replicate", 
                                     "subject_id"),
                        metaData = metadata(SE),
                        toSE = TRUE)

test_that("check equality before and after tidy/return", {
    expect_true(all(assay(SE) == assay(returnedSE)))
    expect_true(identical(rowData(SE), rowData(returnedSE)))
    expect_true(identical(colData(SE)@listData, colData(returnedSE)@listData))  
    expect_true(identical(metadata(SE), 
                          metadata(returnedSE)))
    })


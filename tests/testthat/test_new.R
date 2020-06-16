context("equality between DF and SE functions")

data(msquant)

## Create SE from MSQuant
colnames <- colnames(msquant)[3:56] %>% 
    str_replace_all("Neutral_Operator_Dif_Pos_", "")

spike <- vector("character", 54)
replicate <- vector("character", 54)
batch <- vector("character", 54)
subject_id <- vector("character", 54)

for(i in 1:54) {
    spike[i] <- stringr::str_split(colnames, "_")[[i]][1]
    batch[i] <- stringr::str_split(colnames, "_")[[i]][2]
    replicate[i] <- stringr::str_split(colnames, "_")[[i]][3]
    subject_id[i] <- stringr::str_split(colnames, "_")[[i]][4]
}

colnamesDF <- data.frame(spike, batch, replicate, subject_id)
rownamesDF <- msquant[1:2]
abundanceMatrix <- as.matrix(msquant[3:56])
rownames(abundanceMatrix) <- c()

se <- SummarizedExperiment(abundanceMatrix, colData = colnamesDF, 
                           rowData = rownamesDF, 
                           metadata = "Test Metadata")

## Run msTidy with SE and data frame
tidySE <- msTidy(se, missingValue = 1)

tidyDF <- msTidy(msquant,
                 compVars = c("mz", "rt"),
                 techVars = c("spike", "batch", "replicate", "subject_id"),
                 colExtraText = "Neutral_Operator_Dif_Pos_",
                 separator = "_",
                 missingValue = 1)

## Run msSummarize with SE and data frame
summarizedDF <- msSummarize(msquant,
                            compVars = c("mz", "rt"),
                            techVars = c("spike", "batch", "replicate", "subject_id"),
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

# tidyMatrix <- msTidy(as.matrix(msquant),
#                      compVars = c("mz", "rt"),
#                      techVars = c("spike", "batch", "replicate", 
#                                   "subject_id"),
#                      colExtraText = "Neutral_Operator_Dif_Pos_",
#                      separator = "_",
#                      missingValue = 1)

## Check that results are equal
test_that("Check msTidy()", {
    expect_true(all(tidySE == tidyDF | is.na(tidySE == tidyDF)))
})

# test_that("Check msTidy", {expect_true(all(tidySE == tidyMatrix | 
#                                                is.na(tidySE == tidyMatrix)))})

test_that("Check msSummarize()", {
    expect_true(identical(assay(summarizedSE), 
                          as.matrix(summarizedDF$data[, 3:20])))
    expect_true(all(all(rowData(summarizedSE) == 
                            as.data.frame(summarizedDF$data[, 1:2]))))
    expect_true(identical(S4Vectors::metadata(summarizedSE)$summaryDetails,
                          summarizedDF$summaryDetails))
})


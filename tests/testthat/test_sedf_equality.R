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

## Run msImpute with SE and data frame
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

hmImputedSE <- msImpute(filteredSE, imputeMethod = "halfmin", 
                        returnToSE = TRUE,
                        missingValue = 0)

bpcaImputedSE <- msImpute(filteredSE, imputeMethod = "bpca",
                          returnToSE = TRUE,
                          missingValue = 0)

knnImputedSE <- msImpute(filteredSE, imputeMethod = "knn",
                         returnToSE = TRUE,
                         missingValue = 0)

## Check that results are equal
test_that("Check msTidy()", {
    expect_true(all(tidySE == tidyDF | is.na(tidySE == tidyDF)))
})

test_that("Check msSummarize()", {
    expect_true(identical(assay(summarizedSE), 
                          as.matrix(summarizedDF$data[, 3:20])))
    expect_true(all(all(rowData(summarizedSE) == 
                            as.data.frame(summarizedDF$data[, 1:2]))))
    expect_true(identical(S4Vectors::metadata(summarizedSE)$summaryDetails,
                          summarizedDF$summaryDetails))
})

test_that("Check msFilter()", {
    expect_true(identical(assay(filteredSE), 
                          as.matrix(filteredDF[, 3:20])))
    expect_true(all(all(rowData(filteredSE) == 
                            as.data.frame(filteredDF[, 1:2]))))
    expect_true(all(rownames(colData(filteredSE)) == 
                        colnames(filteredDF[3:20])))
})

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

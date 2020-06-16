context("msTidy() and msReturn()")

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
metaData <- data.frame(test = c(1, 2, 3), test2 = c(2, 3, 4))

se <- SummarizedExperiment(abundanceMatrix, colData = colnamesDF, 
                           rowData = rownamesDF, 
                           metadata = metaData)

tidySE <- msTidy(se, missingValue = 1)
returnedSE <- .msReturn(tidySE,
                        compVars = c("mz", "rt"),
                        techVars = c("spike", "batch", "replicate", 
                                     "subject_id"),
                        metaData = S4Vectors::metadata(se),
                        toSE = TRUE)

test_that("check equality before and after tidy/return", {
    expect_true(all(assay(se) == assay(returnedSE) | is.na(assay(se) == 
                                                        assay(returnedSE))))
    
    expect_true(identical(rowData(se), rowData(returnedSE)))
    expect_true(identical(colData(se)@listData, colData(returnedSE)@listData))
    expect_true(identical(S4Vectors::metadata(se), 
                          S4Vectors::metadata(returnedSE)))
    })
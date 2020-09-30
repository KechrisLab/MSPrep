context(".msTidy() and return functions")

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

SE <- SummarizedExperiment(abundanceMatrix, colData = colnamesDF, 
                           rowData = rownamesDF, 
                           metadata = metaData)

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


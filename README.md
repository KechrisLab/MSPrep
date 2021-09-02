MSPrep
======

### Introduction

`MSPrep` provides a convenient set of functionalities used in the pre-analytic
processing pipeline for mass spectrometry based metabolomics data. Functions are
included for the following processes commonly performed prior to analysis of
such data:

1. Summarization of technical replicates (if available)
2. Filtering of metabolites
3. Imputation of missing values
4. Transformation, normalization, and batch correction

Original manuscript published in
[Bioinformatics](https://academic.oup.com/bioinformatics/article/30/1/133/236721),
and package is hosted by [Bioconductor](https://bioconductor.org/packages/release/bioc/html/MSPrep.html).

Additional helpful links:
1. [Vignette providing detailed instructions with examples](https://bioconductor.org/packages/release/bioc/vignettes/MSPrep/inst/doc/using_MSPrep.html)
2. [Reference Manual describing function usage](https://bioconductor.org/packages/release/bioc/manuals/MSPrep/man/MSPrep.pdf)


### Installation

Install via Bioconductor:

    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")

    BiocManager::install("MSPrep")

Install via Github:

    if (!require("devtools")) install.packages("devtools")
    devtools::install_github("KechrisLab/MSPrep")

### Examples

Two examples are provided below. For more detailed information see the
package Vignette which can be accessed [via Bioconductor](https://bioconductor.org/packages/release/bioc/vignettes/MSPrep/inst/doc/using_MSPrep.html)
or by using the following R command following package installation:

```s
vignette("using_MSPrep", package = "MSPrep")
```

The following code loads the example data set, `MSQuant`, summarizes its
technical replicates, filters metabolites by only keeping those which are
present in 80% of samples, imputes missing values using k-nearest neighbors,
applies a log base ten transformation, and finally normalizes and batch corrects
the data set using quantile normalization and ComBat batch correction. Data is
then returned as a `data.frame`.

```s
library(MSPrep)
data(msquant)

preparedDF <- msPrepare(msquant,
                        minPropPresent = 1/3,
                        missingValue = 1,
                        filterPercent = 0.8,
                        imputeMethod = "knn",
                        transform = "log10",
                        normalizeMethod = "quantile + ComBat",
                        covariatesOfInterest = c("spike"),
                        compVars = c("mz", "rt"),
                        sampleVars = c("spike", "batch", "replicate",
                                       "subject_id"),
                        colExtraText = "Neutral_Operator_Dif_Pos_",
                        separator = "_")
```

The second example uses the data set `COPD_131`. The raw data set can be found [here, at Metabolomics Workbench.](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR000438). The code loads the data set,
summarizes its
technical replicates, filters metabolites by only keeping those which are
present in 80% of samples, imputes missing values using BPCA imputation,
and finally normalizes the data set using median normalization. Data is then
returned as a `SummarizedExperiment` by setting the argument
`returnToSE = TRUE`.

```s
library(MSPrep)
data(COPD_131)

preparedSE <- msPrepare(COPD_131,
                        minPropPresent = 1/3,
                        filterPercent = 0.8,
                        missingValue = 0,
                        imputeMethod = "bpca",
                        nPcs = 3,
                        normalizeMethod = "median",
                        transform = "none",
                        compVars = c("Mass", "Retention.Time",
                                     "Compound.Name"),
                        sampleVars = c("subject_id", "replicate"),
                        colExtraText = "X",
                        separator = "_",
                        returnToSE = TRUE)
```
### Bug Reports

Report bugs as issues on the [GitHub repository new
issue](https://github.com/KechrisLab/MSPrep/issues/new)

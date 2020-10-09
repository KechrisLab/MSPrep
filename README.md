MSPrep
======

A processing pipeline for summarizing, filtering, imputing, and normalizing 
LC/MS metabolomics data.

Original manuscript published in
[Bioinformatics](https://academic.oup.com/bioinformatics/article/30/1/133/236721)

### Installation

Install via Bioconductor:

    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    
    BiocManager::install("MSPrep")

Install via Github:

    if (!require("devtools")) install.packages("devtools")
    devtools::install_github("KechrisLab/MSPrep")

### Bug Reports

Report bugs as issues on the [GitHub repository new
issue](https://github.com/KechrisLab/MSPrep/issues/new)
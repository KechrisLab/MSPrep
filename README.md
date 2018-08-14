
[![Build Status](https://travis-ci.org/KechrisLab/MSPrep.svg?branch=master)](https://travis-ci.org/KechrisLab/MSPrep)
[![codecov](https://codecov.io/gh/KechrisLab/MSPrep/branch/master/graph/badge.svg)](https://codecov.io/gh/KechrisLab/MSPrep)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)<Paste>
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/MSPrep)](https://cran.r-project.org/package=MSPrep)

# MSPrep 

A processing pipeline for the summarization, normalization and diagnostics of
mass spectrometry–based metabolomics data.

Original manuscript published in
[Biometrics](https://academic.oup.com/bioinformatics/article/30/1/133/236721)

This refactoring of the original package seeks to extend the package to any number
of replicates and allow for more normalization and imputation methods.

To install and use the package follow the examples below.  Note that you will
may need to tidy your own dataset, replacing the ms_tidy function in this
pipeline.  See `?ms_tidy` or the dataset resulting from `ms_tidy(msquant)` for
more detail.

```{r, install-and-example}

# Install package
if (!require("devtools")) install.packages("devtools")
devtools::install_github("KechrisLab/MSPrep")

# Load package
library(MSPrep)
library(tidyverse)

# Load example quantification dataset and view 
data(msquant)
as_data_frame(msquant)
names(msquant)

# Use pipe to tidy and create summarized dataset
dat <- 
  msquant %>% ms_tidy %>% ms_prepare %>% ms_filter(0.8)

# Use MSPrep functions one at a time
dat <- ms_tidy(msquant)
dat <- ms_prepare(dat)
dat <- ms_filter(dat)
dat <- ms_impute(method = "knn") 
# dat <- ms_normalize(method = "knn") # not yet implemented

```

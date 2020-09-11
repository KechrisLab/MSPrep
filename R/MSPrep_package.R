#' Package for summarizing, filtering, imputing, and normalizing metabolomics 
#' data.
#' 
#' This package performs summarization of replicates, filtering by frequency,
#' three different options for handling/imputing missing data, and five options 
#' for normalizing data. 
#'
#' @author Max McGrath
#' @author Matt Mulvahill
#' @author Grant Hughes
#' @author Sean Jacobson
#' @author Harrison Pielke-Lombardo
#' @author Katerina Kechris
#' @docType package
#' @name MSPrep 
#' @details
#' Package for pre-analytic processing of mass spectrometry quantification data.
#' Six functions are provided and are intended to be used in sequence (as a 
#' pipeline) to produce cleaned and normalized data. These are ms_tidy(), 
#' ms_summarize(), ms_filter(), ms_impute(), ms_normalize(), and ms_return(). 
#' The function ms_prepare() is also provided as a wrapper function combining 
#' the six previously mentioned functions.
#' 
#' @references
#' Bolstad, B.M.et al.(2003) A comparison of normalization methods for high
#' density oligonucleotide array data based on variance and bias.
#' Bioinformatics, 19, 185-193
#' 
#' DeLivera, A.M.et al.(2012) Normalizing and Integrating Metabolomic Data.
#' Anal. Chem, 84, 10768-10776.
#' 
#' Gagnon-Bartsh, J.A.et al.(2012) Using control genes to correct for unwanted
#' variation in microarray data. Biostatistics, 13, 539-552.
#' 
#' Johnson, W.E.et al.(2007) Adjusting batch effects in microarray expression
#' data using Empirical Bayes methods. Biostatistics, 8, 118-127.
#' 
#' Leek, J.T.et al.(2007) Capturing Heterogeneity in Gene Expression Studies by
#' Surrogate Variable Analysis. PLoS Genetics, 3(9), e161.
#' 
#' Oba, S.et al.(2003) A Bayesian missing value estimation for gene expression
#' profile data. Bioinformatics, 19, 2088-2096
#' 
#' Redestig, H.et al.(2009) Compensation for Systematic Cross-Contribution
#' Improves Normalization of Mass Spectrometry Based Metabolomics Data. Anal.
#' Chem., 81, 7974-7980.
#' 
#' Stacklies, W.et al.(2007) pcaMethods: A bioconductor package providing PCA
#' methods for incomplete data. Bioinformatics, 23, 1164-1167.
#' 
#' Wang, W.et al.(2003) Quantification of Proteins and Metabolites by Mass
#' Spectrometry without Isotopic Labeling or Spiked Standards. Anal. Chem., 75,
#' 4818-4826.
#' 
#' @examples
#' # Load example data
#' data(msquant)
#' 
#' # Call function to tidy, summarize, filter, impute, and normalize data
#' preparedDF <- msPrepare(msquant,
#'                         minPropPresent = 1/3,
#'                         missingValue = 1,
#'                         filterPercent = 0.8,
#'                         imputeMethod = "knn",
#'                         normalizeMethod = "quantile + ComBat",
#'                         transform = "log10",
#'                         covariatesOfInterest = c("spike"),
#'                         compVars = c("mz", "rt"),
#'                         sampleVars = c("spike", "batch", "replicate", 
#'                                        "subject_id"),
#'                         colExtraText = "Neutral_Operator_Dif_Pos_",
#'                         separator = "_")
#' 
NULL

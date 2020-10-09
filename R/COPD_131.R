#' Example mass spectrometry dataset
#' 
#' Data contains LC-MS metabolite analysis for samples from 131 subjects 
#' with 3 technical replicates per subject. The first three columns indicate 
#' "Mass" (mass-to-charge ratio), 
#' "Retention.Time", and "Compound.Name" for the 662 unique metabolites observed
#' in the samples. The remaining columns indicate abundance for each of the 662
#' mass/retention-time combination for each subject/replicate combination.
#' 
#' @docType data
#' @format Data frame containing 662 observations of 396 samples
#' \describe{
#'   \item{Mass}{Mass-to-charge ratio}
#'   \item{Retention.Time}{Retention-time}
#'   \item{Compound.Name}{Compound name for each mass/retention time 
#'   combination}
#'   \item{X10062C_1}{The remaining columns indicate metabolite abundances 
#'   found in each Subject/Replicate combination. Each column begins with an 
#'   'X', followed by the subject ID, and then the replicate (1, 2, or 3), 
#'   each separated by '_'.}
#' }
#' @source \url{https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR000438}
#' @references 
#' Nichole Reisdorph (NaN). Untargeted LC-MS metabolomics analysis of human 
#' COPD plasma, HILIC & C18, metabolomics_workbench, V1.
#' 
#' This data is available at the NIH Common Fund's National Metabolomics Data 
#' Repository (NMDR) website, the Metabolomics Workbench, 
#' \url{https://www.metabolomicsworkbench.org}, where it has been assigned 
#' Project ID PR000438. The data can be accessed directly via its Project DOI: 
#' 10.21228/M8FC7C. This work is supported by NIH grant, U2C-DK119886.
#' @keywords datasets
#' @usage data(COPD_131)
#' @examples
#' data(COPD_131)
"COPD_131"
#' Example mass spectrometry dataset. 
#'
#' Data contains LC-MS samples for 2 subjects, each with three features: 
#' spike-in (1x, 2x, 4x), batch (01, 02, 03), and technical replicate (A, B, C).
#' The first two columns indicate mass-to-charge ratio and retention-time for
#' the 2644 unique metabolites observed in the samples. The remaining 54 
#' columns indicate metabolite abundance for each subject/spike-in/
#' batch/replicate combination.
#'
#' @docType data
#' @format Data frame containing 2644 observations of 56 variables
#' \describe{
#'   \item{mz}{Mass-to-charge ratio}
#'   \item{rt}{Retention-time}
#'   \item{Neutral_Operator_Dif_Pos_1x_O1_A_01}{The remaining columns specify
#'   metabolite abundances found in each subject/spike-in/batch/replicate 
#'   combination. Each columns begins with
#'   'Neutral_Operator_Dif_Pos' followed by the spike-in (1x, 2x, or 4x), then
#'   the batch (01, 02, or 03), the replicate (A, B, or C), and finally the
#'   subject ID (01 or 02), each seperated by '_'.}
#' }
#' @keywords datasets
#' @references
#'     Hughes, G., Cruickshank-Quinn, C., Reisdorph, R., Lutz, S., Petrache, I.,
#'     Reisdorph, N., … Kechris, K. (2014). MSPrep--summarization, normalization
#'     and diagnostics for processing of mass spectrometry-based metabolomic 
#'     data. 
#'     Bioinformatics (Oxford, England), 30(1), 133–134. 
#'     doi:10.1093/bioinformatics/btt589
#' @examples
#'   data(msquant)
"msquant"




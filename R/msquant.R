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
#' @format A data frame with 2644 observations of 29 variables
#' \describe{
#'   \item{mz}{Mass-to-charge ratio}
#'   \item{rt}{Retention-time}
#'   \item{Neutral_Operator_Dif_Pos_1x_O1_A_01}{The remaining columns specify
#'   metabolite abundances found in each subject/spike-in/batch/replicate 
#'   combinations. Each columns begins with
#'   'Neutral_Operator_Dif_Pos' followed by the spike-in (1x, 2x, or 4x), then
#'   the batch (01, 02, or 03), the replicate (A, B, or C), and finally the
#'   subject ID (01 or 02), each seperated by '_'.}
#' }
#' @keywords datasets
#' @examples
#'   data(msquant)
#'   str(msquant)
"msquant"




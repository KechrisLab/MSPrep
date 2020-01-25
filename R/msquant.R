#' Example mass spectrometry dataset. 
#' 
#' Size: 2644 obs. of 56 variables.
#'
#' Data contains samples for 2 subjects. First two variables are:
#' mz = Mass/Charge ratio
#' rt = Retention time
#' 
#' The other 54 variables are all combinations of spike-in spike-in (1x, 2x, 4x), 
#' batch (O1, O2, O3), technical replicate (A, B, C), and subject ID (01, 02). These
#' are indicated in the column names, separated by '_'.
#'
#' @docType data
#' @format 
#' \code{
#' 'data.frame':   2654 obs. of  29 variables:
#'  $ mz                              : num  577 539 724 525 483 ...
#'  $ rt                              : Factor w/ 2510 levels "0.481","0.58933336",..: 3 2 6 15 10 8 22 18 17 13 ...
#'  $ Neutral_Operator_Dif_Pos_1x_O1_A_01: num  1 1 36270 105305 1 ...
#'  $ Neutral_Operator_Dif_Pos_1x_O1_B_01 num  1 1 35271 92489 1 ...
#'  $ Neutral_Operator_Dif_Pos_1x_O1_C_01 num  1 1 39355 105189 1 ...
#'  $ Neutral_Operator_Dif_Pos_1x_O2_A_01 num  24182 20686 57021 141999 67760 ...
#'  $ Neutral_Operator_Dif_Pos_1x_O2_B_01 num  25857 24500 60641 151657 73151 ...
#'  $ Neutral_Operator_Dif_Pos_1x_O2_C_01 num  24616 22820 61023 144739 71789 ...
#' }
#' @keywords datasets
#' @examples
#'   data(msquant)
#'   str(msquant)
"msquant"




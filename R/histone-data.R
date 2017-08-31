#' Data on histone modifications.
#'
#' Data for the prediction of steady-state transcript abundance in developing secondary xylem from histone modification ChIP-seq data in E. grandis.
#'
#' @docType data
#'
#' @usage data(histone)
#'
#' @format A \code{data.frame} with five variables:
#' \describe{
#' \item{Gene_ID}{Factor containing 35711 unique gene identifiers.}
#' \item{Raw_FPKM}{A numeric value containing the raw FPKM signal.}
#' \item{H3K4me3_signal_bin25}{A numeric value containing the H3K4me3 signal that has maximum positive correlation with gene expression level, at bin 25.}
#' \item{H3K27me3_signal_bin21}{A numeric value containing the H3K27me3 signal that has maximum negative correlation with gene expression level, at bin 21.}
#' \item{Train}{A logical indicator giving a 60\% split for validation purposes.}
#' }
#'
#' @keywords datasets
#' @source Hussey, S.G., Loots, M.T., van der Merwe, K., Mizrachi, E. and Myburg, A.A., 2017. Integrated analysis and transcript abundance modelling of H3K4me3 and H3K27me3 in developing secondary xylem. Scientific Reports, 7.
#'
#' @examples
#' library(alR)
#' str(histone)
#' summary(histone$Raw_FPKM)
#' plot(Raw_FPKM~H3K4me3_signal_bin25, data=histone)
#' plot(Raw_FPKM~H3K27me3_signal_bin21, data=histone)
'histone'
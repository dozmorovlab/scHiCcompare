#' scHi-C table from microglia (MG) cell type - chromosome 22 at 1 MB resolution
#'
#' A table that gathers 10 single-cell Hi-C data modified into a sparse upper
#'  triangular matrix. It contains each cell's interacting region coordinates
#'  (cell, chr, start1, start2) and the corresponding Interaction Frequencies
#'  (IF) for each single cell
#'  (IF_1, IF_2, ..., IF_10).
#'
#' @format A data frame with 335 rows and 5 columns:
#' \describe{
#'   \item{cell}{Name of the target cell type}
#'   \item{chr}{The chromosome of the interacting region's starting coordinate}
#'   \item{region1}{The starting coordinate of the first interacting region}
#'   \item{region2}{The starting coordinate of the second interacting region}
#'   \item{IF_1}{Interaction frequency values for single cell 1}
#'   \item{IF_2}{Interaction frequency values for single cell 2}
#'   \item{IF_3}{Interaction frequency values for single cell 3}
#'   \item{IF_4}{Interaction frequency values for single cell 4}
#'   \item{IF_5}{Interaction frequency values for single cell 5}
#'   \item{IF_6}{Interaction frequency values for single cell 6}
#'   \item{IF_7}{Interaction frequency values for single cell 7}
#'   \item{IF_8}{Interaction frequency values for single cell 8}
#'   \item{IF_9}{Interaction frequency values for single cell 9}
#'   \item{IF_10}{Interaction frequency values for single cell 10}
#' }
#'
#' @source Single-cell data of MG cell type downloaded using `download_schic()`
#'  from `Bandnorm`. See their website at
#'  \url{https://sshen82.github.io/BandNorm/articles/BandNorm-tutorial.
#'  html#download-existing-single-cell-hi-c-data}.
#'
#' @return A data frame
#' @usage data("scHiC.table_MG_chr22")
#' @examples
#' data("scHiC.table_MG_chr22")
#' head(scHiC.table_MG_chr22)
"scHiC.table_MG_chr22"

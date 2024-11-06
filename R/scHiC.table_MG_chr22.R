#' scHi-C table from microglia (MG) cell type - chromosome 22 at 1 MB resolution
#'
#' A table that gathers 10 single-cell Hi-C data modified into a sparse upper triangular matrix. 
#' It contains each cell's interacting region coordinates (cell, chr, start1, start2) 
#' and the corresponding Interaction Frequencies (IF) for each single cell (IF_1, IF_2, ..., IF_10).
#'
#' @format A data frame with 335 rows and 5 columns:
#' \describe{
#'   \item{cell}{Name of the target cell type}
#'   \item{chr}{The chromosome of the interacting region's starting coordinate}
#'   \item{start1}{The starting coordinate of the first interacting region}
#'   \item{start2}{The starting coordinate of the second interacting region}
#'   \item{IF_1, IF_2, ..., IF_10}{Interaction frequency values for each single cell}
#' }
#'
#' @source Single-cell data of MG cell type downloaded using `download_schic()` from `Bandnorm`. 
#' See their website at \url{https://sshen82.github.io/BandNorm/articles/BandNorm-tutorial.html#download-existing-single-cell-hi-c-data}.
#'
#' @return A data frame
#' @usage data("scHiC.table_MG_chr22")
#' @examples
#' data("scHiC.table_MG_chr22")
#' head(scHiC.table_MG_chr22)
"scHiC.table_MG_chr22"

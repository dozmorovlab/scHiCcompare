#' scHi-C table from oligodendrocyte (ODC) cell type - chromosome 22 at 1 MB resolution
#'
#' A table gather 10 single cells Hi-C data modified in sparse upper triangular matrix containing 
#' the cell's interacting regions cordination (cell, chr, start1, start2) 
#' and the corresponding Interaction Frequencies for each single cell (IF_1, IF2,... IF10).
#'
#'
#'
#' @format A dataframe with 335 rows and 5 columns:
#'     \describe{
#'     \item{cell}{Name of the target cell type }
#'     \item{chr}{The interacting region starting cordination chromosome}
#'     \item{start1}{The first interacting region starting cordination}
#'     \item{start2}{The second interacting region starting cordination}
#'     \item{IF_1, IF_2, ..., IF_10}{Interaction frequency values for each single cell}
#' }
#'
#' @source Single-cells data of ODC cell type download by download_schic() by `Bandnorm`. See their website at
#'     \url{https://sshen82.github.io/BandNorm/articles/BandNorm-tutorial.html#download-existing-single-cell-hi-c-data}
#'
#' @return A data frame

#'
"scHiC.table_ODC_chr22"

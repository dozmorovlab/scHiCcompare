#' scHi-C table from oligodendrocyte (ODC) cell type - chromosome 22 at 1 MB resolution
#'
#' A table gather 10 single cells Hi-C data modified in sparse upper triangular matrix containing
#' the cell's interacting regions cordination (cell, chr, start1, start2)
#' and the corresponding Interaction Frequencies for each single cell (IF_1, IF2,... IF10).
#'
#'
#'
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
#'
#' @source Single-cells data of ODC cell type download by download_schic() by `Bandnorm`. See their website at
#'     \url{https://sshen82.github.io/BandNorm/articles/BandNorm-tutorial.html#download-existing-single-cell-hi-c-data}
#'
#' @return A data frame

#'
"scHiC.table_ODC_chr22"

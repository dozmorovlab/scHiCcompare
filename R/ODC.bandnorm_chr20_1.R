#' scHi-C data from oligodendrocyte (ODC) cell type - chromosome 20 at 1 MB resolution
#'
#' A Modified sparse upper triangular matrix containing the interacting regions cordination (chr1, start1, chr2, start2)
#'  and the corresponding Interaction Frequency (IF).
#'
#'
#'
#' @format A dataframe with 335 rows and 5 columns:
#'     \describe{
#'     \item{chr}{The first interacting region chromosome}
#'     \item{start1}{The first interacting region starting cordination}
#'     \item{chr}{The second interacting region chromosome}
#'     \item{start2}{The second interacting region starting cordination}
#'     \item{IF}{The interaction frequency value between two regions}
#'     }
#'
#' @source Data downloaded by download_schic() by `Bandnorm`. See their website at
#'     \url{https://sshen82.github.io/BandNorm/articles/BandNorm-tutorial.html#download-existing-single-cell-hi-c-data}
#'     The data is a single-cell Hi-C of ODC cell type of Lee et al. public dataset. The GEO link to download the data
#'     \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130711}
#'
#' @return A data frame

#'
"ODC.bandnorm_chr20_1"

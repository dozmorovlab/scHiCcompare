#' Generate Pseudo-bulk Hi-C Data
#'
#' This function generates pseudo-bulk Hi-C data from a single-cell Hi-C
#'  interaction frequency table. It returns the data in either sparse or full
#'  matrix format.
#'
#' @param scHiC.table A data frame containing interaction frequencies across
#'  cells that are created by the \code{scHiC_table()} function.
#'  The first four columns should represent 'cell', 'chr', 'region1', and
#'  'region2',followed by columns representing interaction frequencies ('IF')
#'  for individual cells.
#' @param out A character string specifying the output format. It must be either
#'  'sparse' or 'full'. Default is 'sparse'.
#'
#' @return A data frame representing the pseudo-bulk Hi-C data.
#'  If `out` is 'sparse'; returns a data frame containing a pseudo-bulk matrix
#'  in the sparse upper triangular Hi-C matrix format with five columns:
#'   chromosome 1 (chr1), start position 1 (start1), chromosome 2 (chr2),
#'   start position 2 (start2), and interaction frequency (IF).
#'   If `out` is 'full';
#'         returns a full matrix representation of the pseudo-bulk data.
#'
#' @references
#'
#' Stansfield JC, Cresswell KG, Vladimirov VI  et al (2018).
#' Hiccompare: an R-package for joint normalization and comparison of HI-C
#' datasets. BMC Bioinformatics  2018;19:279.
#'
#' @examples
#' data("scHiC.table_MG_chr22")
#' data("scHiC.table_ODC_chr22")
#' pseudo_bulk_MG <- pseudo_bulkHic(scHiC.table_MG_chr22)
#' pseudo_bulk_ODC <- pseudo_bulkHic(scHiC.table_ODC_chr22)
#' head(pseudo_bulk_MG)
#' head(pseudo_bulk_ODC)
#'
#' @import HiCcompare
#'
#' @export

## Create PseudoBulk Sparse ----
pseudo_bulkHic <- function(scHiC.table, out = "sparse") {
  # Function to generate pseudo-bulk Hi-C data
  # Input: scHiC.table (table of interaction frequencies across single cells)
  #        out ('sparse' or 'full' output format)
  # Output: Pseudo-bulk sparse data or full matrix

  # Extract IF columns and calculate pseudo-bulk IF as row sums
  if_scs <- scHiC.table[, -c(1, 2, 3, 4)]
  IF <- rowSums(if_scs)
  bulk.table <- cbind(scHiC.table[, c(1, 2, 3, 4)], IF)

  # Remove regions where pseudo-bulk IF is 0
  bulk.table <- bulk.table[bulk.table$IF != 0, ]

  # Remove 'cell' and 'chr' columns
  bulk.table <- bulk.table[, !names(bulk.table) %in% c("cell", "chr")]

  # Return either sparse or full matrix
  if (out == "sparse") {
    message("\nTransfering into pseudo-bulk sparse matrix.")
    output.table <- bulk.table
  } else {
    message("\nTransfering pseudo-bulk full matrix.")
    output.table <- HiCcompare::sparse2full(sparse.mat = bulk.table)
  }

  return(output.table)
}

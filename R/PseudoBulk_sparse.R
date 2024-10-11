#' Generate Pseudo-bulk Hi-C Data
#'
#' This function generates pseudo-bulk Hi-C data from a single-cell Hi-C interaction frequency table. 
#' It returns the data in either sparse or full matrix format.
#'
#' @param scHiC.table A data frame containing interaction frequencies across single cells that is created by the \code{scHiC_table()} function. 
#'                     The first four columns should represent 'cell', 'chr', 'region1', and 'region2', 
#'                     followed by columns representing interaction frequencies ('IF') for individual cells.
#' @param out A character string specifying the output format. It must be either 'sparse' or 'full'.
#'            Default is 'sparse'.
#'
#' @return A data frame representing the pseudo-bulk Hi-C data. If `out` is 'sparse', 
#'         it returns a data frame with pseudo-bulk matrix in the format of a sparse upper triangular Hi-C matrix with five columns:
#' chromosome 1 (chr1), start position 1 (start1), chromosome 2 (chr2), start position 2 (start2), and interaction frequency (IF).  If `out` is 'full', 
#'         it returns a full matrix representation of the pseudo-bulk data.
#'
#' @references
#' 
#' Stansfield JC, Cresswell KG, Vladimirov VI  et al (2018). Hiccompare: an R-package for joint normalization and comparison of HI-C datasets. BMC Bioinformatics  2018;19:279.
#' 
#' @examples
#' \dontrun{
#' # Load MG datasets folder example
#' load_example_MGfolder()
#' # Generate scHiC table for single cell Hi-C data of MG cell type
#' scHiC.table <- scHiC_table(file.path = "MGs_example",  cell.type = "MG", position.dataset = 1:50, type= "txt", select.chromosome = 22)
#' 
#' # Transform into the pseudo bulk HiC
#' pseudo_bulk_data <- pseudo_bulkHic(scHiC.table, out = 'sparse')
#' full_matrix_data <- pseudo_bulkHic(scHiC.table, out = 'full')
#' }
#'
#' @export

## Create PseudoBulk Sparse ----
pseudo_bulkHic <- function(scHiC.table, out = 'sparse'){ 
  # Function to generate pseudo-bulk Hi-C data
  # Input: scHiC.table (table of interaction frequencies across single cells)
  #        out ('sparse' or 'full' output format)
  # Output: Pseudo-bulk sparse data or full matrix
  
  # Check for valid 'out' input
  if(!out %in% c('sparse', 'full')){
    stop("Error: 'out' must be either 'sparse' or 'full'.")
  }
  
  # Extract IF columns and calculate pseudo-bulk IF as row sums
  if_scs <- scHiC.table[,-c(1,2,3,4)] # Remove 'cell', 'chr', 'region1', 'region2' columns
  IF <- rowSums(if_scs) # Sum interaction frequencies across cells
  bulk.table <- cbind(scHiC.table[,c(1,2,3,4)], IF) # Add pseudo-bulk IF to table
  
  # Remove regions where pseudo-bulk IF is 0
  bulk.table <- bulk.table[bulk.table$IF != 0, ]
  
  # Remove 'cell' and 'chr' columns
  bulk.table <- bulk.table[, !names(bulk.table) %in% c('cell', 'chr')]
  
  # Return either sparse or full matrix
  if(out == 'sparse'){
    message("Returning pseudo-bulk sparse matrix.")
    output.table <- bulk.table
  } else {
    message("Returning pseudo-bulk full matrix.")
    output.table <- HiCcompare::sparse2full(sparse.mat = bulk.table)
  }
  
  return(output.table)
}

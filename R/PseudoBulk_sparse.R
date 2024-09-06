## Create PseudoBulk Sparse ----
pseudo_bulkHic <- function(scHiC_table, out = 'sparse'){ 
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

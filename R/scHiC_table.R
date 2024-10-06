
#' Create scHiC Interaction Frequency Table
#'
#' This function generates a single-cell Hi-C interaction frequency (IF) table for all single cells at a selected chromosome. The resulting table is usable for `Pooling_RF_impute` and `pseudo_bulkHic` functions. It reads the input files, extracts the relevant data, and outputs a table of interaction frequencies between genomic regions for each single cell.
#'
#' @param file.list The list object where elements are processed scHi-C data file.
#' @param cell.type The cell type used in the analysis (e.g., 'NSN', 'SN').
#' @param position.dataset A vector of indices specifying the file positions to read from the directory.
#' @param select.chromosome The chromosome name to be studied (e.g., 'chr1' or 'chrX').
#' @return A data frame containing the interaction frequency table with genomic loci (cell, chromosome, start1, end1) and interaction frequencies (IF) of each single cell. This table can be used with `Pooling_RF_impute` and `pseudo_bulkHic` functions.
#' @details
#' This function processes single-cell Hi-C data in a list object. Then the function transforms them into a single 'scHiC table' data frame. Each element in the list should be in the form of a sparse upper triangular Hi-C matrix with foyr tab-separated columns (chr, start1,  start2, IF) and no row or column names or quotes around character strings. 
#' @examples
#' \dontrun{
#' # Load MG data folder example
#' load_example_MGfolder()
#' MGs_example = read_files(file.path = "MGs_example", cell = "MG", position.dataset = c(1, 2, 3), type = "txt", txt.sparse.heads.position = c(1,2,4,5))
#' # Create scHiC table to be used in scHiCcompare
#' IF_table <- scHiC_table(file.list = MGs_example, cell.type = 'MG', position.dataset = 1:3, select.chromosome = 'chr22')
#' }
#' @import dplyr
#' @export

## Create scHiC.table -----
scHiC_table <- function(file.list = NULL, cell.type, position.dataset, select.chromosome){
  # Input validation
  if(missing(cell.type) || missing(position.dataset) || missing(select.chromosome)){
    stop("Error: Missing one or more required arguments.")
  }
  
  # Read data 
  datasets = file.list[position.dataset]
  n_sc <- length(datasets)
  
  if (n_sc == 0) {
    stop("Error: No datasets available for processing.")
  }
  
  
  
  regions <- NULL
  options(scipen = 999)
  
  for(i in 1:n_sc){
    data <- datasets[[i]]  # Get the dataset from the list

    # Select chromosome
    data <- data[data[,1] == select.chromosome,]
    
    # Transform dataset into sparse
    data <- data[, -c(1)]
    names(data) <- c('region1', 'region2', 'IF')
    
    # Remove rows with 0 values in any column
    data <- data[data[,1] != 0 & data[,2] != 0,]
    
    # Extract single cell regions
    region1sc <- unique(c(data$region1, data$region2))
    regions <- unique(c(regions, region1sc))
  }
  
  if (is.null(regions) || length(regions) == 0) {
    stop("Error: No valid regions found in the datasets.")
  }
  
  start.region <- min(regions, na.rm = TRUE)
  end.region <- max(regions, na.rm = TRUE)
  bin <- min(abs(diff(as.numeric(regions))), na.rm = TRUE)
  
  if (bin == 0 || is.na(bin)) {
    stop("Error: Unable to calculate bin size for regions.")
  }
  
  regions.seq <- seq(start.region, end.region, by = bin)
  
  # Coordination of pair of bins
  grid1 <- data.frame(X1 = 1:length(regions.seq), X2 = 1:length(regions.seq))  # diagonal line
  grid2 <- data.frame(t(combn(1:length(regions.seq), 2)))  # off diagonal line value
  grid <- rbind(grid1, grid2)
  region1 <- regions.seq[grid[,1]]
  region2 <- regions.seq[grid[,2]]
  cordination <- cbind(region1, region2)
  
  # Preallocate memory for the table
  table <- data.frame(
    region1 = cordination[,'region1'],
    region2 = cordination[,'region2']
  )
  
  for (i in 1:n_sc) {
    data <- datasets[[i]]  # Get the dataset from the list
    
    # Filter rows based on chromosome
    data <- data[data[,1] == select.chromosome & data[,3], ]
    if (nrow(data) == 0) {
      warning(paste("Warning: No data found for chromosome", select.chromosome, "in dataset", i))
      next
    }
    
    data <- data[, -c(1)]
    names(data)[c(1, 2)] <- c('region1', 'region2')
    names(data)[3] <- paste0('IF_', i)
    data <- data[rowSums(data[,c(1,2)] == 0) == 0, ]
    
    table <- tryCatch(
      dplyr::full_join(table, data, by = c('region1', 'region2')),
      error = function(e) stop("Error in joining tables: ", e$message)
    )
  }
  
  table[is.na(table)] <- 0
  
  up.table <- data.frame(
    cell = rep(cell.type, nrow(table)),
    chr = rep(select.chromosome, nrow(table)),
    table
  )
  
  return(up.table)
}

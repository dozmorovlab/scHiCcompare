
#' Create scHiC Interaction Frequency Table
#'
#' This function generates a single-cell Hi-C interaction frequency (IF) table for all single cells for a selected chromosome. The resulting table is usable for the \code{Pooling_RF_impute()} and \code{pseudo_bulkHic()} functions. It reads the input files, extracts the relevant data, and outputs a table of interaction frequencies between genomic regions for each single cell dataset.
#'
#' @param file.path Character string specifying the directory containing scHi-C data for condition (a cell-type group). The folder should contain '.txt' scHi-C files in modified sparse upper triangular format (chr1, start1, chr2, start2, IF)
#' @param cell.type The cell type name used in the analysis (e.g., 'NSN', 'SN').
#' @param select.chromosome The chromosome name to be studied (e.g., 'chr1' or 'chrX').
#' @return A data frame containing the interaction frequency table with genomic loci (cell, chromosome, start1, end1) and interaction frequencies (IF) of each single cell. This table can be used with the \code{Pooling_RF_impute()} and \code{pseudo_bulkHic()} functions.
#' @details
#' This function processes single-cell Hi-C data in a folder directory, then transforms them into a single 'scHiC table' data frame. Each element in the list should be in the form of a sparse upper triangular Hi-C matrix with five tab-separated columns (chr1, start1, chr2, start2, IF) with no row or column names and no quotes around character strings.
#' @examples
#' \dontrun{
#' # Load MG data folder example
#' load_example_MGfolder()
#'
#' # Create scHiC table to be used in scHiCcompare
#' IF_table <- scHiC_table(file.path = "MGs_example", cell.type = 'MG', position.dataset = 1:3, select.chromosome = 'chr22')
#' }
#' @import dplyr
#' @export

## Create scHiC.table -----
scHiC_table <- function(file.path, cell.type, select.chromosome) {
  # Input validation
  if (missing(cell.type) || missing(select.chromosome)) {
    stop("Error: Missing one or more required arguments.")
  }
  
  library(scHiCcompare)
  
  # Read data
  datasets <- read_files(file.path = file.path, type='txt',
                         txt.sparse.heads.position = c(1,2,3,4,5), out = 'sparse')
  n_sc <- length(datasets)
  
  if (n_sc == 0) {
    stop("Error: No datasets available for processing.")
  }
  
  regions <- NULL
  options(scipen = 999)
  
  # Process datasets and extract unique regions
  regions <- unique(unlist(lapply(datasets, function(data) {
    data <- data[data$chr1 == select.chromosome, ]  # Select chromosome
    region1sc <- unique(c(data$start1, data$start2))  # Extract regions
    region1sc <- region1sc[region1sc != 0]  # Remove zero regions
    return(region1sc)
  })))
  regions <- as.numeric(regions)
  
  if (is.null(regions) || length(regions) == 0) {
    stop("Error: No valid regions found in the datasets.")
  }
  
  # Define region sequences
  start.region <- min(regions, na.rm = TRUE)
  end.region <- max(regions, na.rm = TRUE)
  bin <- min(diff(sort(unique(regions))), na.rm = TRUE)
  
  if (bin <= 0 || is.na(bin)) {
    stop("Error: Unable to calculate bin size for regions.")
  }
  
  regions.seq <- seq(start.region, end.region, by = bin)
  
  # Generate coordinate pairs (combinations of regions)
  grid <- expand.grid(region1 = regions.seq, region2 = regions.seq)
  grid <- grid[grid$region1 <= grid$region2, ]  # Only unique pairs and diagonal
  
  # Preallocate memory for the table
  table <- grid
  
  # Join datasets
  for (i in 1:n_sc) {
    data <- datasets[[i]]  # Get the dataset
    data <- data[data$chr1 == select.chromosome & data$IF != 0, -c(1,3)]  # Filter rows
    
    if (nrow(data) == 0) {
      warning(paste("Warning: No data found for chromosome", select.chromosome, "in dataset", i))
      next
    }
    
    names(data) <- c('region1', 'region2', paste0('IF_', i))
    
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


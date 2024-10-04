#' Read Files for scHi-C Dataset
#'
#' This function reads single-cell Hi-C data from a specified file path.
#' It supports two file formats: 'txt' and 'cool'. For 'txt', it reads tab-delimited files and assumes the format contains five columns: 
#' chromosome 1 (chr1), start position 1 (start1), chromosome 2 (chr2), start position 2 (start2), and interaction frequency (IF).
#' For 'cool', it uses the `HiCcompare` package to transform cooler files to BEDPE format.
#'
#' @param file.path The directory path where the data files are stored. 
#' @param cell A character string specifying a cell type (eg., 'MG', etc).
#' @param position.dataset A vector of indices specifying the file positions to read from the directory.
#' @param type The file type, either 'txt' or 'cool'. Default is 'txt'. Each 'txt; file should be in the format of a sparse upper triangular Hi-C matrix with five columns: 
#' chromosome 1 (chr1), start position 1 (start1), chromosome 2 (chr2), start position 2 (start2), and interaction frequency (IF). 
#' If the 'cool' type is selected, the .cool files are HDF5 containers that store Hi-C data, which will be read using `cooler2bedpe()` from the `HiCcompare` package.
#' @return A list of datasets, where each element corresponds to a dataset from the selected files. If 'txt' type is selected, each element contains a data frame of a sparse upper triangular Hi-C matrix. 
#' If 'cool' is selected, each element contains a list with two items: 
#' \itemize{
#'   \item "cis": Contains the intra-chromosomal contact matrices, one per chromosome.
#'   \item "trans": Contains the inter-chromosomal contact matrix. 
#' }
#' @details
#' This function input single cell Hi-C data in '.txt' or 'cool' format. Each 'txt; file should be in the form of a sparse upper triangular Hi-C matrix with five columns: 
#' chromosome 1 (chr1), start genomic position 1 (start1), chromosome 2 (chr2), start genomic position 2 (start2), and interaction frequency (IF), with format of tab-separated columns without row names, column names, or quotes around character strings.
#' If the 'cool' type is selected, the .cool files are HDF5 containers that store Hi-C data, which will be read using `cooler2bedpe()` from the `HiCcompare` package.
#'  
#' 
#' @examples
#' \dontrun{
#' datasets <- read_files("scHiCcompare/Example/Data/MG", "NSN", c(1, 2, 3), type = "txt")
#' }
#' @import HiCcompare
#' @export



read_files <- function(file.path, cell, position.dataset, type='txt'){
  # Check if the file path exists
  if(!dir.exists(file.path)){
    stop("Error: The specified file path does not exist.")
  }
  
  data_names <- list.files(path = file.path, full.names = TRUE, recursive = TRUE)[position.dataset] 
  
  # Check if the selected files are valid
  if (length(data_names) == 0) {
    stop("Error: No files found for the specified position.dataset.")
  }
  
  datasets <- list()  # Initialize a list to store the datasets
  
  if (type == 'txt') {
    for(i in 1:length(position.dataset)){
      dataset <- tryCatch(
        read.delim(data_names[i]),
        error = function(e) stop(paste("Error reading file:", data_names[i]))
      )
      datasets[[i]] <- dataset  # Store each dataset in the list
    }
  } else if (type == 'cool') {
    if (!requireNamespace("HiCcompare", quietly = TRUE)) {
      stop("Error: HiCcompare package not installed. Please install the HiCcompare package.")
    }
    
    library(HiCcompare)
    for(i in 1:length(position.dataset)){
      dataset <- tryCatch(
        cooler2bedpe(data_names[i]),
        error = function(e) stop(paste("Error reading .cool file:", data_names[i]))
      )
      datasets[[i]] <- dataset  # Store each dataset in the list
    }
  } else {
    stop("Error: Unsupported file type. Use 'txt' or 'cool'.")
  }
  
  return(datasets)
}





#' Create scHiC Interaction Frequency Table
#'
#' This function generates a single-cell Hi-C interaction frequency (IF) table for all single cells at a selected chromosome. The resulting table is usable for `Pooling_RF_impute` and `pseudo_bulkHic` functions. It reads the input files, extracts the relevant data, and outputs a table of interaction frequencies between genomic regions for each single cell.
#'
#' @param file.path The directory path where the data files are stored.
#' @param cell.type The cell type used in the analysis (e.g., 'NSN', 'SN').
#' @param position.dataset A vector of indices specifying the file positions to read from the directory.
#' @param type The file type, either 'txt' or 'cool'. Default is 'txt'. If 'txt', the files should be in the format of a sparse upper triangular Hi-C matrix with five columns: chromosome 1 (chr1), start position 1 (start1), chromosome 2 (chr2), start position 2 (start2), and interaction frequency (IF). If 'cool', the .cool files will be read using `cooler2bedpe()` from the `HiCcompare` package.
#' @param select.chromosome The chromosome name to be studied (e.g., 'chr1' or 'chrX').
#' @return A data frame containing the interaction frequency table with genomic loci (cell, chromosome, start1, end1) and interaction frequencies (IF). This table can be used with `Pooling_RF_impute` and `pseudo_bulkHic` functions.
#' @details
#' This function processes single-cell Hi-C data in '.txt' or 'cool' format, transforming them into a single 'scHiC table' data frame that is . Each input '.txt' file in the folder should be in the form of a sparse upper triangular Hi-C matrix with five tab-separated columns (chr1, start1, chr2, start2, IF) and no row or column names or quotes around character strings. If the 'cool' type is selected, the function uses `cooler2bedpe()` from the `HiCcompare` package to read the data.
#' @examples
#' \dontrun{
#' # Load MG data folder example
#' load_example_MGfolder()
#' # Create scHiC table to be used in scHiCcompare
#' IF_table <- scHiC_table(file.path = "MGs_example", cell.type = 'MG', position.dataset = 1:50, type = 'txt', select.chromosome = 'chr22')
#' }
#' @import dplyr
#' @export

## Create scHiC.table -----
scHiC_table <- function(file.path, cell.type, position.dataset, type='txt', select.chromosome){
  # Input validation
  if(missing(file.path) || missing(cell.type) || missing(position.dataset) || missing(select.chromosome)){
    stop("Error: Missing one or more required arguments.")
  }
  
  # Read data with error handling
  datasets <- tryCatch(
    read_files(file.path = file.path, cell = cell.type, position.dataset = position.dataset, type = type),
    error = function(e) stop("Error in reading files: ", e$message)
  )
  
  n_sc <- length(datasets)
  
  if (n_sc == 0) {
    stop("Error: No datasets available for processing.")
  }
  
  regions <- NULL
  options(scipen = 999)
  
  for(i in 1:n_sc){
    data <- datasets[[i]]  # Get the dataset from the list

    # Select chromosome
    data <- data[data[,1] == select.chromosome & data[,3] == select.chromosome,]
    
    # Transform dataset into sparse
    data <- data[, c(2, 4, 5)]
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
    data <- data[data[,1] == select.chromosome & data[,3] == select.chromosome, ]
    if (nrow(data) == 0) {
      warning(paste("Warning: No data found for chromosome", select.chromosome, "in dataset", i))
      next
    }
    
    data <- data[, c(2, 4, 5)]
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

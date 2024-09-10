.read_files <- function(file_path, cell, position_dataset, type='txt'){
  # Check if the file path exists
  if(!dir.exists(file_path)){
    stop("Error: The specified file path does not exist.")
  }
  
  data_names <- list.files(path = file_path, full.names = TRUE, recursive = TRUE)[position_dataset] 
  
  # Check if the selected files are valid
  if (length(data_names) == 0) {
    stop("Error: No files found for the specified position_dataset.")
  }
  
  datasets <- list()  # Initialize a list to store the datasets
  
  if (type == 'txt') {
    for(i in 1:length(position_dataset)){
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
    for(i in 1:length(position_dataset)){
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



## Create scHiC.table -----
scHiC.table <- function(file_path, cell_type, position_dataset, type='txt', select_chromosome){
  # Input validation
  if(missing(file_path) || missing(cell_type) || missing(position_dataset) || missing(select_chromosome)){
    stop("Error: Missing one or more required arguments.")
  }
  
  # Read data with error handling
  datasets <- tryCatch(
    .read_files(file_path = file_path, cell = cell_type, position_dataset = position_dataset, type = type),
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
    data <- data[data[,1] == select_chromosome & data[,3] == select_chromosome,]
    
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
    data <- data[data[,1] == select_chromosome & data[,3] == select_chromosome, ]
    if (nrow(data) == 0) {
      warning(paste("Warning: No data found for chromosome", select_chromosome, "in dataset", i))
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
    cell = rep(cell_type, nrow(table)),
    chr = rep(select_chromosome, nrow(table)),
    table
  )
  
  return(up.table)
}

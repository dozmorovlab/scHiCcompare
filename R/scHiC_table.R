scHiC_table <- function(cell_type, n_sc, select_chromosome) {
  # Input: single-cell Hi-C data with columns: chr1, start1, chr2, start2, IF
  # Arguments:
  # - cell_type: Prefix for single-cell data (assumed to be named as "cell_type_i" for i=1 to n_sc)
  # - n_sc: Number of single cells
  # - select_chromosome: Chromosome of interest
  # Output: A table with columns for interaction frequencies (IF) across single cells
  
  # Check for valid inputs
  if (!is.character(cell_type) || length(cell_type) != 1) {
    stop("Error: 'cell_type' must be a single string.")
  }
  
  if (!is.numeric(n_sc) || n_sc <= 0) {
    stop("Error: 'n_sc' must be a positive number.")
  }
  
  if (!is.character(select_chromosome) || length(select_chromosome) != 1) {
    stop("Error: 'select_chromosome' must be a single string.")
  }
  
  # Initialize regions and prevent scientific notation in output
  regions <- NULL
  options(scipen = 999)
  
  # Loop through each single-cell dataset to extract regions
  for (i in 1:n_sc) {
    dataset_name <- paste0(cell_type, "_", i) # Construct dataset name
    
    # Check if dataset exists
    if (!exists(dataset_name)) {
      warning(paste("Warning: Dataset", dataset_name, "not found. Skipping this dataset."))
      next
    }
    
    data <- get(dataset_name) # Retrieve dataset using name
    
    # Check if dataset has the correct columns
    if (ncol(data) < 5 || !all(c('chr1', 'start1', 'chr2', 'start2', 'IF') %in% colnames(data))) {
      warning(paste("Warning: Dataset", dataset_name, "does not have the expected columns. Skipping this dataset."))
      next
    }
    
    # Filter by selected chromosome
    data <- data[data[,1] == select_chromosome, c(2, 4, 5)]
    
    # Check if there's data for the selected chromosome
    if (nrow(data) == 0) {
      warning(paste("Warning: No data for chromosome", select_chromosome, "in dataset", dataset_name, ". Skipping this dataset."))
      next
    }
    
    names(data) <- c('region1', 'region2', 'IF')
    
    # Remove rows with 0 values in region columns
    data <- data[data$region1 != 0 & data$region2 != 0,]
    
    # If there's no valid data after filtering, show a warning
    if (nrow(data) == 0) {
      warning(paste("Warning: No valid data (non-zero regions) in dataset", dataset_name, ". Skipping this dataset."))
      next
    }
    
    # Collect unique regions from both region1 and region2
    region1sc <- unique(c(data$region1, data$region2))
    regions <- unique(c(regions, region1sc))
  }
  
  # Check if there are any regions to process
  if (is.null(regions) || length(regions) == 0) {
    stop("Error: No valid regions found across all single cells. Check your input data.")
  }
  
  # Determine start, end regions, and bin size
  start.region <- min(regions, na.rm = TRUE)
  end.region <- max(regions, na.rm = TRUE)
  bin <- min(abs(diff(sort(regions))), na.rm = TRUE)
  
  # Create a sequence of regions with the determined bin size
  regions.seq <- seq(start.region, end.region, by = bin)
  
  # Create all possible region pairs (combinations of regions)
  grid1 <- data.frame(X1 = 1:length(regions.seq), X2 = 1:length(regions.seq)) # diagonal pairs
  grid2 <- data.frame(t(combn(1:length(regions.seq), 2))) # off-diagonal pairs
  grid <- rbind(grid1, grid2)
  region1 <- regions.seq[grid[, 1]]
  region2 <- regions.seq[grid[, 2]]
  coordinates <- data.frame(region1, region2)
  
  # Initialize the table with region pairs
  table <- data.frame(
    region1 = coordinates$region1,
    region2 = coordinates$region2
  )
  
  # Loop through each single-cell dataset to add IF columns
  for (i in 1:n_sc) {
    dataset_name <- paste0(cell_type, "_", i)
    
    # Check if dataset exists before proceeding
    if (!exists(dataset_name)) {
      next
    }
    
    data <- get(dataset_name)
    
    # Filter by selected chromosome and rename columns
    data <- data[data[,1] == select_chromosome, c(2, 4, 5)]
    
    # Check for data existence after filtering
    if (nrow(data) == 0) {
      next
    }
    
    names(data) <- c('region1', 'region2', paste0('IF_', i))
    
    # Remove rows with 0 values in region columns
    data <- data[data$region1 != 0 & data$region2 != 0, ]
    
    # Merge the current dataset's IF with the main table
    table <- dplyr::full_join(table, data, by = c('region1', 'region2'))
  }
  
  # If no valid data is merged, stop with an error
  if (ncol(table) == 2) {
    stop("Error: No valid interaction frequency data found for any single cells.")
  }
  
  # Replace NA values with 0 in the final table
  table[is.na(table)] <- 0
  
  # Add cell and chromosome columns
  up.table <- data.frame(
    cell = rep(cell_type, nrow(table)),
    chr = rep(select_chromosome, nrow(table)),
    table
  )
  
  message("scHiC.table completed successfully.")
  return(up.table)
}

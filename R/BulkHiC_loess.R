
BulkHiC_loess <- function(bulkHiC_sparse_1, bulkHiC_sparse_2,
                            selected_chromosome, scale = FALSE, include.zeros = FALSE, subset.dist = NA, subset.index = NA,
                            exclude.regions = NA, exclude.overlap = 0.2,
                            smooth.degree = 1, smooth.span = NA, loess.criterion = "gcv",
                            Plot = FALSE, Plot.smooth = TRUE, parallel = FALSE,
                            BP_param = bpparam()) {
  # Function to apply loess normalization to pseudo-bulk Hi-C data
  # Input: bulkHiC_sparse_1, bulkHiC_sparse_2 (sparse matrices for two conditions)
  #        selected_chromosome (chromosome for analysis)
  #        scale, include.zeros, subset.dist, etc. are additional parameters for the normalization
  # Output: Loess-normalized Hi-C interaction table
  
  # Check for valid input types
  if (!is.data.frame(bulkHiC_sparse_1) || !is.data.frame(bulkHiC_sparse_2)) {
    stop("Error: bulkHiC_sparse_1 and bulkHiC_sparse_2 must be data frames.")
  }
  
  if (!is.character(selected_chromosome) || length(selected_chromosome) != 1) {
    stop("Error: 'selected_chromosome' must be a single string representing the chromosome.")
  }
  
  # Check that smoothing degree is valid
  if (!is.numeric(smooth.degree) || smooth.degree < 0) {
    stop("Error: 'smooth.degree' must be a non-negative numeric value.")
  }
  
  # Check that the exclude.overlap is between 0 and 1
  if (!is.numeric(exclude.overlap) || exclude.overlap < 0 || exclude.overlap > 1) {
    stop("Error: 'exclude.overlap' must be a numeric value between 0 and 1.")
  }
  
  # Create the Hi-C table
  hic.table <- create.hic.table(bulkHiC_sparse_1, bulkHiC_sparse_2, chr = selected_chromosome, scale = scale,
                                include.zeros = include.zeros, subset.dist = subset.dist, subset.index = subset.index,
                                exclude.regions = exclude.regions, exclude.overlap = exclude.overlap)
  
  # Perform LOESS normalization
  hic.table_cc <- hic_loess(hic.table, degree = smooth.degree, span = smooth.span, loess.criterion = loess.criterion,
                            Plot = Plot, Plot.smooth = Plot.smooth, parallel = parallel,
                            BP_param = BP_param)
  
  # Return the loess-normalized Hi-C interaction table
  message("LOESS normalization complete.")
  return(hic.table_cc)
}


  

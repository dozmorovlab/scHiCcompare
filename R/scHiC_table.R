#' Create scHiC Interaction Frequency Table
#'
#' This function generates a single-cell Hi-C interaction frequency (IF) table for all single cells for a selected chromosome. The resulting table is usable for the \code{Pooling_RF_impute()} and \code{pseudo_bulkHic()} functions. It reads the input files, extracts the relevant data, and outputs a table of interaction frequencies between genomic regions for each single cell dataset.
#'
#' @param file.path Character string specifying the directory containing scHi-C data for condition (a cell-type group). The folder should contain '.txt' scHi-C files in a modified sparse upper triangular format (chr1, start1, chr2, start2, IF)
#' @param cell.type The cell type name used in the analysis (e.g., 'NSN', 'SN').
#' @param select.chromosome The chromosome name to be studied (e.g., 'chr1' or 'chrX').
#' @return A data frame containing the interaction frequency table with genomic loci (cell, chromosome, start1, end1) and interaction frequencies (IF) of each single cell. This table can be used with the \code{Pooling_RF_impute()} and \code{pseudo_bulkHic()} functions.
#' @details
#' This function processes single-cell Hi-C data in a folder directory, then transforms them into a single 'scHiC table' data frame. Each element in the list should be in the form of a sparse upper triangular Hi-C matrix with five tab-separated columns (chr1, start1, chr2, start2, IF) with no row or column names and no quotes around character strings.
#' @examples
#' # Load MG data folder example
#' MGs_example <- system.file("MGs_example", package = "scHiCcompare")
#'
#' # Create scHiC table to be used in scHiCcompare
#' IF_table <- scHiC_table(file.path = MGs_example, cell.type = "MG", select.chromosome = "chr20")
#'
#' @import dplyr
#' @export

## Create scHiC.table -----
scHiC_table <- function(file.path, cell.type, select.chromosome) {
  # Input validation
  if (missing(cell.type) || missing(select.chromosome)) {
    stop("Missing one or more required arguments: 'cell.type' and 'select.chromosome' are needed.")
  }

  # Read data
  datasets <- read_files(
    file.path = file.path, type = "txt",
    txt.sparse.heads.position = c(1, 2, 3, 4, 5), out = "sparse"
  )
  n_sc <- length(datasets)

  if (n_sc == 0) {
    stop("No datasets available for processing at the specified 'file.path'.")
  }

  # Extract unique regions from specified chromosome
  regions <- unique(unlist(lapply(datasets, function(data) {
    data <- data[data$chr1 == select.chromosome, ]
    region1sc <- unique(c(data$start1, data$start2))
    region1sc <- region1sc[region1sc != 0]
    return(region1sc)
  })))

  if (length(regions) == 0) {
    stop("No valid regions found for the specified chromosome in the datasets.")
  }

  # Calculate region sequences
  regions <- as.numeric(regions)
  start.region <- min(regions, na.rm = TRUE)
  end.region <- max(regions, na.rm = TRUE)
  bin <- min(diff(sort(unique(regions))), na.rm = TRUE)

  if (bin <= 0 || is.na(bin)) {
    stop("Unable to calculate a valid bin size for the regions.")
  }

  regions.seq <- seq(start.region, end.region, by = bin)

  # Generate coordinate pairs (combinations of regions)
  grid <- expand.grid(region1 = regions.seq, region2 = regions.seq)
  grid <- grid[grid$region1 <= grid$region2, ]

  # Preallocate memory for the table
  table <- grid

  # Join datasets by merging interaction frequencies
  for (i in seq_len(n_sc)) {
    data <- datasets[[i]]
    data <- data[data$chr1 == select.chromosome & data$IF != 0, -c(1, 3)]

    if (nrow(data) == 0) {
      message(sprintf("No data found for chromosome %s in dataset %d.", select.chromosome, i))
      next
    }

    names(data) <- c("region1", "region2", paste0("IF_", i))

    table <- tryCatch(
      dplyr::full_join(table, data, by = c("region1", "region2")),
      error = function(e) stop("Error in joining tables: ", e$message)
    )
  }

  table[is.na(table)] <- 0

  # Final output table with metadata
  up.table <- data.frame(
    cell = rep(cell.type, nrow(table)),
    chr = rep(select.chromosome, nrow(table)),
    table
  )

  return(up.table)
}

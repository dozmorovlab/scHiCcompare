#' Read Files for scHi-C Dataset
#'
#' This function reads single-cell Hi-C data from a specified file path.
#' It supports two file formats: 'txt' and 'cool'. For 'txt', it reads
#'  tab-delimited files and assumes the format contains five columns:
#'  chromosome 1 (chr1), start position 1 (start1), chromosome 2 (chr2),
#'  start position 2 (start2), and interaction frequency (IF).
#'   For 'cool', it uses the `HiCcompare` package to transform cooler files
#'    to BEDPE format.
#'
#' @param file.path The directory path where the data files are stored.
#' @param position.dataset A vector of indices specifying the file positions to
#'  read from the directory. These indices help select specific files,
#'  determining which files will be included and the sequence in which they are
#'  processed. If all single cell data should be included, set this to NULL.
#' @param type The file type, either 'txt'. The default is 'txt'. Each 'txt'
#' file should be in the format of a sparse upper triangular Hi-C matrix, where
#' each row contains the interaction frequency value (IFs) of two interacting
#' regions.
#' @param txt.sparse.heads.position A vector of four integers specifying the
#' column positions of chromosome, start1, start2, and IF in the 'txt' file.
#' @param out Character string specifying the output format. Options are
#' "sparse" (for a sparse matrix format) or "original" (to retain the original
#' structure of each single-cell Hi-C dataset). Default is "sparse".
#' @return A list of datasets, where each element corresponds to a dataset from
#' the selected files.
#' If `out` is "sparse", each dataset element is transformed into a sparse
#'  matrix format (chr, start1, start2, IF). If `out` is "original", the
#'  original structure of each single-cell Hi-C dataset is preserved.
#'
#' @details
#' This function reads single-cell Hi-C data in 'txt', with output options of
#' 'sparse' and 'original'.
#' Each input 'txt' file should be in the form of a sparse upper triangular
#'  Hi-C matrix, storing pair-wise interaction frequencies of loci pairs.
#'  The 'txt' dataset should have one column indicating the interaction
#'  frequency (IF) of each pair of interacting regions, with tab-separated
#'  columns and no row names, column names, or quotes around character strings.
#'
#' @examples
#' # Load MG data folder example
#' MGs_example <- system.file("extdata/MGs_example", package = "scHiCcompare")
#' datasets <- read_files(
#'   file.path = MGs_example, position.dataset = c(1, 2, 3, 4, 5),
#'   txt.sparse.heads.position = c(1, 2, 3, 4, 5)
#' )
#' @import HiCcompare
#' @importFrom utils read.delim

#' @export

read_files <- function(file.path, position.dataset = NULL, type = "txt",
                       txt.sparse.heads.position = NULL, out = "sparse") {
  # List all files in the directory with full names
  all_files <- list.files(path = file.path, full.names = TRUE)

  # Filter for specified file type (e.g., .txt)
  txt_files <- all_files[grepl(sprintf("\\.%s$", type), all_files)]

  # If position.dataset is NULL, use all files
  if (is.null(position.dataset)) {
    position.dataset <- seq_along(txt_files)
  }

  # Subset txt_files based on position.dataset
  txt_files <- txt_files[position.dataset]

  # Initialize a list to store the datasets
  datasets <- list()

  # Read the files if type is 'txt'
  for (i in seq_along(txt_files)) {
    dataset <- read.table(txt_files[i], header = TRUE)

    # Process the data for sparse matrix format if 'out' is 'sparse'
    if (out == "sparse") {
      # Identify head positions for sparse matrix format
      chr1.position <- txt.sparse.heads.position[1]
      start1.position <- txt.sparse.heads.position[2]
      chr2.position <- txt.sparse.heads.position[3]
      start2.position <- txt.sparse.heads.position[4]
      if.position <- txt.sparse.heads.position[5]

      dataset <- dataset[, c(
        chr1.position, start1.position,
        chr2.position, start2.position, if.position
      )]
      names(dataset) <- c("chr1", "start1", "chr2", "start2", "IF")
    }

    datasets[[i]] <- dataset
  }

  return(datasets)
}

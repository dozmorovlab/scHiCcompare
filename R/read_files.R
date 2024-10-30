#' Read Files for scHi-C Dataset
#'
#' This function reads single-cell Hi-C data from a specified file path.
#' It supports two file formats: 'txt' and 'cool'. For 'txt', it reads tab-delimited files and assumes the format contains five columns:
#' chromosome 1 (chr1), start position 1 (start1), chromosome 2 (chr2), start position 2 (start2), and interaction frequency (IF).
#' For 'cool', it uses the `HiCcompare` package to transform cooler files to BEDPE format.
#'
#' @param file.path The directory path where the data files are stored.
#' @param position.dataset A vector of indices specifying the file positions to read from the directory. These indices help select specific files, determining which files will be included and the sequence in which they are processed. If all single cell data should be include, choose NULL.
#' @param type The file type, either 'txt'. Default is 'txt'. Each 'txt' file should be in the format of a sparse upper triangular Hi-C matrix, where each row contains the interaction frequency value (IFs) of two interacting regions.
#' @param txt.sparse.heads.position A vector of four integers specifying the column positions of chromosome, start1, start2, and IF in the 'txt' file.
#' @return A list of datasets, where each element corresponds to a dataset from the selected files. If `out` is 'sparse', each dataset element is transformed into sparse matrix format (chr, start1, start2, IF).
#' If `out` is 'original', the original structure of each single-cell Hi-C dataset is preserved.
#'
#' @details
#' This function reads single-cell Hi-C data in 'txt', with output options of 'sparse' and 'original'. Each input 'txt' file should be in the form of a sparse upper triangular Hi-C matrix,
#' storing pair-wise interaction frequencies of loci pairs. The 'txt' dataset should have one column indicating the interaction frequency (IF) of each pair of interacting regions, with tab-separated columns and no row names, column names, or quotes around character strings.
#'
#' @examples
#' \dontrun{
#' # Load MG data folder example
#' load_example_MGfolder()
#' datasets <- read_files(file.path = "MGs_example", cell = "MG", position.dataset = c(1, 2, 3), type = "txt", txt.sparse.heads.position = c(1, 2, 4, 5))
#' }
#' @import HiCcompare
#' @export


read_files <- function(file.path, position.dataset = NULL, type = "txt",
                       txt.sparse.heads.position = NULL, out = "sparse") {
  # Check if the file path exists
  if (!dir.exists(file.path)) {
    stop("Error: The specified file path does not exist.")
  }

  # List all files in the directory with full names
  all_files <- list.files(path = file.path, full.names = TRUE)

  # Filter for .txt files
  txt_files <- all_files[grepl(paste0("\\.", type, "$"), all_files)]

  # Check if there are any .txt files in the directory
  if (length(txt_files) == 0) {
    stop("Error: No .txt files found in the specified directory.")
  }

  # If position.dataset is NULL, use all files
  if (is.null(position.dataset)) {
    position.dataset <- 1:length(txt_files)
  }

  # Ensure the selected indices are valid
  if (any(position.dataset > length(txt_files))) {
    stop("Error: Invalid positions specified in position.dataset.")
  }

  # Subset txt_files based on position.dataset
  txt_files <- txt_files[position.dataset]

  # Initialize a list to store the datasets
  datasets <- list()

  # Read the files if type is 'txt'
  for (i in seq_along(txt_files)) {
    dataset <- tryCatch(
      {
        read.delim(txt_files[i])
      },
      error = function(e) stop(paste("Error reading file:", txt_files[i]))
    )

    # Process the data for sparse matrix format if 'out' is 'sparse'
    if (out == "sparse") {
      ## Identify head positions for sparse matrix format
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

    datasets[[i]] <- dataset # Store each dataset in the list
  }

  return(datasets)
}

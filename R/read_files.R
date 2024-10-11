#' Read Files for scHi-C Dataset
#'
#' This function reads single-cell Hi-C data from a specified file path.
#' It supports two file formats: 'txt' and 'cool'. For 'txt', it reads tab-delimited files and assumes the format contains five columns:
#' chromosome 1 (chr1), start position 1 (start1), chromosome 2 (chr2), start position 2 (start2), and interaction frequency (IF).
#' For 'cool', it uses the `HiCcompare` package to transform cooler files to BEDPE format.
#'
#' @param file.path The directory path where the data files are stored.
#' @param cell A character string specifying a cell type (e.g., 'MG', etc). [??? Where to find it?]
#' @param position.dataset A vector of indices specifying the file positions to read from the directory. [??? What is a file index here? What does this do?]
#' @param type The file type, either 'txt' or 'cool'. Default is 'txt'. [??? Is position.dataset needed for txt?] Each 'txt' file should be in the format of a sparse upper triangular Hi-C matrix, where each row contains the interaction frequency value (IFs) of two interacting regions.
#' If the 'cool' type is selected, the 'cool' files are HDF5 containers that store Hi-C data, which will be read using `cooler2bedpe()` [??? Why not just use this function? writing a wrapper is reinventing the weel] from the `HiCcompare` package. [??? This parameter can be automated by reading the file extension from the 'file.path' argument]
#' @param txt.sparse.heads.position A vector of four integers specifying the column positions of chromosome, start1, start2, and IF in the 'txt' file. [??? This could have a default order; 'out' is sparse, so (chr, start1, start2, IF)?] [??? why it is needed at all? very confusing]
#' @param out Output format with options 'sparse' and 'original'. If 'sparse', the sparse upper triangular matrix format is returned. If 'original', the dataset will retain its original structure. Default is 'sparse'. [??? Why even care if the package uses sparse? HiCcompare has funcions to confert to full]
#'
#' @return A list of datasets, where each element corresponds to a dataset from the selected files. If `out` is 'sparse', each dataset element is transformed into sparse matrix format (chr, start1, start2, IF).
#' If `out` is 'original', the original structure of each single-cell Hi-C dataset is preserved. If the 'cool' type is selected and `out` is 'original', each dataset will retain the structure from `cooler2bedpe()`, containing a list with two items:
#' \itemize{
#'   \item "cis": Contains the intra-chromosomal contact matrices, one per chromosome.
#'   \item "trans": Contains the inter-chromosomal contact matrix.
#' }
#'
#' @details
#' This function reads single-cell Hi-C data in 'txt' or 'cool' format, with output options of 'sparse' and 'original'. Each input 'txt' file should be in the form of a sparse upper triangular Hi-C matrix, 
#' storing pair-wise interaction frequencies of loci pairs. The 'txt' dataset should have one column indicating the interaction frequency (IF) of each pair of interacting regions, with tab-separated columns and no row names, column names, or quotes around character strings.
#' If the 'cool' type is selected, the input 'cool' files are HDF5 containers that store Hi-C data, which will be read using `cooler2bedpe()` from the `HiCcompare` package. The function can provide output in the dataset's 'original' structure or in 'sparse' matrix format.
#'
#' @examples
#' \dontrun{
#' # Load MG data folder example
#' load_example_MGfolder()
#' datasets <- read_files(file.path = "MGs_example", cell = "MG", position.dataset = c(1, 2, 3), type = "txt", txt.sparse.heads.position = c(1,2,4,5))
#' }
#' @import HiCcompare
#' @export


read_files <- function(file.path, cell, position.dataset, type='txt', txt.sparse.heads.position = NULL, out = 'sparse'){
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
  
  if (type == 'txt') { #[??? may want to capture anything not 'cool' or '.cool' and treat it like txt]
    for(i in 1:length(position.dataset)){ #[??? perhaps just loop through seq_along(data_names)?]
      dataset <- tryCatch(
        
        read.delim(data_names[i]),
        error = function(e) stop(paste("Error reading file:", data_names[i]))
      )
      if(out == 'sparse'){
        ## Identify head position in sparse matrix format
        chr.position = txt.sparse.heads.position[1]
        start1.position = txt.sparse.heads.position[2]
        start2.position = txt.sparse.heads.position[3]
        if.position = txt.sparse.heads.position[4]
        dataset <- dataset[ ,c(chr.position, start1.position, start2.position, if.position)]
        names(dataset) = c('chr', 'start1','start2','IF')
      }
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
        error = function(e) stop(paste("Error reading cool file:", data_names[i]))
      )
      if(out == 'sparse'){
        ## Identify head position in sparse matrix format
        data_all <- do.call(rbind,dataset$cis)
        dataset <- data_all[,c('chr1', 'start1', 'start2','IF')]
        names(dataset) <- c('chr', 'region1', 'region2','IF')
        row.names(dataset) = NULL
      }
      datasets[[i]] <- dataset  # Store each dataset in the list
    }
  } else {
    stop("Error: Unsupported file type. Use 'txt' or 'cool'.")
  }
  
  return(datasets)
}




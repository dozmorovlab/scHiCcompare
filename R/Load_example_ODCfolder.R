#' Load ODC Datasets Folder Example 
#'
#' This function loads example of ODC datasets and saves them in a newly created folder named 'ODCs_example' 
#' within the specified directory. This is the example data folder to be used as example for `scHiC_table` function. If no directory is specified, it defaults to 
#' the current working directory.
#'
#' @param save_dir A character string specifying the directory where the 'ODCs_example' 
#' folder will be created. If NULL (default), the function uses the current working directory.
#'
#' @return Saves the datasets as '.txt' files in the 'ODCs_example' folder.
#'
#' @examples
#' load_example_ODCfolder() # Saves in 'ODCs_example' in the current working directory
#' load_example_ODCfolder("/path/to/save/directory") # Saves in specified directory
#'
#' @export
#' 
load_example_ODCfolder <- function(save_dir = NULL) {
  # Get the list of all .rda files in the package's data directory that contain 'ODC' in their name
  data_dir <- system.file("data", package = "scHiCcompare")
  dataset_files <- list.files(data_dir, pattern = "ODC.bandnorm.*\\.rda$", full.names = TRUE)  # Pattern includes 'ODC'
  
  # Set save_dir to the current working directory if NULL
  if (is.null(save_dir)) {
    save_dir <- getwd()  # Get the current working directory
  }
  
  # Create the 'ODCs_example' folder within save_dir
  save_dir <- file.path(save_dir, "ODCs_example")
  
  # Ensure the save directory exists
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Load each dataset and save it to the specified directory
  for (dataset_file in dataset_files) {
    dataset_name <- tools::file_path_sans_ext(basename(dataset_file))  # Get dataset name without extension
    load(dataset_file)  # Load into the current environment
    
    # Save the dataset into the 'ODCs_example' folder as an RData file
    write.table( get(dataset_name), file = file.path(save_dir, paste0(dataset_name, ".txt")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

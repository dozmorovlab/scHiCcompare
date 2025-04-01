#' Print summary for ScHiCcompare object
#'
#' This function prints a summary of the internal state of a ScHiCcompare 
#' object, includingthe number of cells in each condition, the selected 
#' chromosome, and the processes involved.
#'
#' @param obj An object of class 'ScHiCcompare'.
#' @export
print.scHiCcompare <- function(obj) {
  message("-----------------------------------------------------------------")
  message("   ScHiCcompare - Differential Analysis for single cell Hi-C")
  message("-----------------------------------------------------------------")
  
  cond1_cells <- if (!is.null(obj$Intermediate$Imputation$condition1)) {
    ncol(obj$Intermediate$Imputation$condition1) - 4 
  } else {
    0
  }
  
  cond2_cells <- if (!is.null(obj$Intermediate$Imputation$condition2)) {
    ncol(obj$Intermediate$Imputation$condition2) - 4 
  } else {
    0
  }
  
  selected_chromosome <- if (!is.null(obj$Intermediate$Bulk.Normalization$chr1)) {
    paste(unique(obj$Intermediate$Bulk.Normalization$chr1), collapse = ", ")
  } else {
    "Unknown"
  }
  
  # total_differences <- if (!is.null(obj$Differential_Analysis$Difference.cluster)) {
  #   sum(obj$Differential_Analysis$Difference.cluster, na.rm = TRUE)
  # } else {
  #   NA
  # }
  
  impute <- if (!is.null(obj$Intermediate$Imputation$condition1)) "Imputation" else NULL
  norm <- if (!is.null(obj$Intermediate$Bulk.Normalization)) "Bulk Normalization" else NULL
  
  processes <- c(impute, norm)
  processes <- processes[!is.null(processes)]
  
  message(sprintf(
    "ScHiCcompare analyzes %d cells of condition 1 group and %d cells of condition 2 group at chromosome %s.",
    cond1_cells, cond2_cells, selected_chromosome
  ))
  
  if (length(processes) > 0) {
    message("\nThe process includes:")
    message(sprintf("%s, Differential Analysis", paste(processes, collapse = ", ")))
  } else {
    message("\nNo processes available for analysis.")
  }
  
  message("\nNote: See full differential result in $Differential_Analysis. Intermediate results can be accessed with $Intermediate.\n")
}
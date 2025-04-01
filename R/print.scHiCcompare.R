#' Print method for `scHiCcompare` objects
#'
#' This method is used to print a summary of the `scHiCcompare` object.
#' It provides an overview of the number of cells in each condition, the
#' chromosome analyzed, and the processes involved in the analysis.
#'
#' @param x An object of class `scHiCcompare`. This object should contain
#'   intermediate results from the differential analysis of single-cell
#'   Hi-C data.
#' @param ... Further arguments to be passed to or from other methods.
#'
#' @return This function does not return any value. It is called for its
#'  side effect, which is printing a summary of the `scHiCcompare` object.
#'
#' @examples
#' ## Load example data for ODC and MG file paths
#' ODCs_example <- system.file("extdata/ODCs_example",
#'   package = "scHiCcompare"
#' )
#' MGs_example <- system.file("extdata/MGs_example",
#'   package = "scHiCcompare"
#' )
#'
#' ## Run scHiCcompare on example data
#' result <- scHiCcompare(
#'   file.path.1 = MGs_example,
#'   file.path.2 = ODCs_example,
#'   select.chromosome = "chr20",
#'   Plot = FALSE
#' )
#' print(result)
#'
#' @export print.scHiCcompare
#' @export
#'
print.scHiCcompare <- function(x, ...) {
  message("-----------------------------------------------------------------")
  message("   ScHiCcompare - Differential Analysis for single cell Hi-C")
  message("-----------------------------------------------------------------")

  cond1_cells <- if (!is.null(x$Intermediate$Imputation$condition1)) {
    ncol(x$Intermediate$Imputation$condition1) - 4
  } else {
    0
  }

  cond2_cells <- if (!is.null(x$Intermediate$Imputation$condition2)) {
    ncol(x$Intermediate$Imputation$condition2) - 4
  } else {
    0
  }

  selected_chromosome <- if (!is.null(x$Intermediate$Bulk.Normalization$chr1)) {
    paste(unique(x$Intermediate$Bulk.Normalization$chr1), collapse = ", ")
  } else {
    "Unknown"
  }

  impute <- if (!is.null(x$Intermediate$Imputation$condition1)) "Imputation" else NULL
  norm <- if (!is.null(x$Intermediate$Bulk.Normalization)) "Bulk Normalization" else NULL

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

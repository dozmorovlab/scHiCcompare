#' Plot Hi-C Interaction Matrix Heatmap
#'
#' This function generates a heatmap of Hi-C interaction frequencies from a given sparse matrix,
#' allowing users to visualize either original or imputed Hi-C data.
#'
#' @param scHiC.sparse A modified sparse matrix of Hi-C interaction frequencies
#' in the format (chr1, start1, chr2, start2, IF).
#' @param zlim A numeric vector of length 2 specifying the limits of the color
#' scale. If `zlim` is not specified, it will include the minimum and maximum
#' values of the matrix. For example, `zlim = c(0, 10)` limits the color scale to values between 0 and 10.
#' @param color_low A character string specifying the color for the lowest values. Default
#' is "white". Other options include colors such as "lightblue", "yellow", or colors defined using hex codes (e.g., "#FFFFFF" for white).
#'  Users can refer to the R color names documentation by running `colors()` in R.
#' @param color_high A character string specifying the color for the highest values.
#' Default is "red". Other options include colors such as "lightblue", "yellow", or colors defined using hex codes (e.g., "#FFFFFF" for white).
#'  Users can refer to the R color names documentation by running `colors()` in R.
#' @param main A character string for the main title of the plot. The default is NULL.
#' @param figure_name A character string for additional figure labeling. The default is NULL.
#'
#' @return A heatmap plot visualizing the Hi-C interaction matrix.
#'
#' @examples
#' data("ODC.bandnorm_chr20_1")
#'
#' # Call the function with this sparse matrix to generate the heatmap
#' plot_HiCmatrix_heatmap(
#'   scHiC.sparse = ODC.bandnorm_chr20_1,
#'   zlim = c(0, 7), # Log scale color limits
#'   color_low = "white", # Color for low values
#'   color_high = "red", # Color for high values
#'   main = "Single-Cell Hi-C Heatmap", # Title of the plot
#'   figure_name = "Example Heatmap" # Subtitle for the plot
#' )
#' @importFrom HiCcompare sparse2full
#' @importFrom lattice levelplot
#'
#' @export

plot_HiCmatrix_heatmap <- function(scHiC.sparse, zlim = NULL, color_low = "white",
                                   color_high = "red", main = NULL, figure_name = NULL) {
  
  # Transform sparse matrix to a full matrix with only required columns
  scHiC.sparse <- scHiC.sparse[, c(2, 4, 5)]
  org_sc_full <- HiCcompare::sparse2full(scHiC.sparse)
  
  # Validate if sparse2full returned a matrix
  if (!is.matrix(org_sc_full)) {
    stop("Conversion to full matrix failed. Please check 'sparse2full' function.")
  }
  
  # Define zlim if not provided
  if (is.null(zlim)) {
    zlim <- range(org_sc_full, na.rm = TRUE)
  }
  
  # Define the color palette
  color_scale <- colorRampPalette(c(color_low, color_high))
  colors <- color_scale(16)
  
  # Rotate matrix for correct plotting orientation
  rotated_org_sc_full <- t(apply(org_sc_full, 2, rev))
  
  # Plot with log transformation, adding a small constant to avoid log(0)
  lattice::levelplot(
    log(rotated_org_sc_full + 1e-6), # Adding small constant to avoid log(0)
    pretty = TRUE,
    xlab = "",
    ylab = "",
    scales = list(x = list(at = NULL), y = list(at = NULL)),
    col.regions = colors,
    at = seq(log(zlim[1] + 1e-6), log(zlim[2] + 1e-6), length.out = length(colors) + 1),
    main = main,
    sub = figure_name,
    aspect = 1 # Controls the plot's height/width ratio
  )
}

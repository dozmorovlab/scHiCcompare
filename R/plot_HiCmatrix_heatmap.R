#' Plot Hi-C Interaction Matrix Heatmap
#'
#' This function generates a heatmap of Hi-C interaction frequencies from a given sparse matrix, 
#' allowing users to visualize either original or imputed Hi-C data.
#'
#' @param scHiC.sparse A modified sparse matrix of Hi-C interaction frequencies in format (chr1, start1, chr2, start2, IF).
#' @param zlim A numeric vector of length 2 specifying the limits of the color scale. If the zlim is not specified, the zlim vector include minimum and maximum values of the matrix.
#' @param color_low A character specifying the color for the lowest values. Default is "white".
#' @param color_high A character specifying the color for the highest values. Default is "red".
#' @param main A character string for the main title of the plot. Default is NULL.
#' @param figure_name A character string for additional figure labeling. Default is NULL.
#' 
#' @return A heatmap plot visualizing the Hi-C interaction matrix.
#' 
#' @examples
#' # Create a 36x36 matrix with random integer values between 0 and 9
#' matrix_36x36 <- matrix(sample(0:9, 36*36, replace = TRUE), nrow = 36, ncol = 36)
#' diag(matrix_36x36) <- diag(matrix_36x36) * 6
#' 
#' # Transform full matrix into sparse
#' library(HiCcompare)
#' sparse = full2sparse(matrix_36x36)
#' 
#' # Call the function with this sparse matrix to generate the heatmap
#'plot_imputation_heatmap(
#'  scHiC.sparse = sparse, 
#'  zlim = c(0, 7),              # Log scale color limits
#'  color_low = "white",          # Color for low values
#'  color_high = "red",           # Color for high values
#'  main = "Single-Cell Hi-C Heatmap",  # Title of the plot
#'  figure_name = "Example Heatmap"     # Subtitle for the plot
#')
#' 
#' @export

plot_HiCmatrix_heatmap <- function(scHiC.sparse, zlim = NULL, color_low = "white", color_high = "red", main= NULL, figure_name = NULL) {
  library(lattice)
  
  # Transform sparse matrices to full matrices
  scHiC.sparse <- scHiC.sparse[ ,c(2,4,5)]
  org_sc_full <- sparse2full(scHiC.sparse)
  
  if(is.null(zlim)){
    z.max = max(org_sc_full); z.min = min(org_sc_full)
    zlim = c(z.min, z.max)
  }

  # Define the color palette from low to high
  color_scale <- colorRampPalette(c(color_low, color_high))
  # Generate 16 colors for the scale (corresponding to values 0 to 12)
  colors <- color_scale(16)
  
  # Rotate the original matrix for correct plotting orientation
  rotated_org_sc_full <- t(apply(org_sc_full, 2, rev))
  
  # Plot single-cell matrix heatmap with specified aspect ratio
  
  levelplot(
    log(rotated_org_sc_full), 
    pretty = TRUE, 
    xlab = "", 
    ylab = "", 
    scales = list(x = list(at = NULL), y = list(at = NULL)), 
    col.regions = colors, 
    at = seq(zlim[1], zlim[2], length.out = length(colors) + 1), 
    main = main,
    sub = figure_name,
    aspect = 1  # This controls the plot's height/width ratio
  )
}

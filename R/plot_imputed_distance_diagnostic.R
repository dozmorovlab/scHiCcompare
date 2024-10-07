#' Plot Imputed Distance Diagnostic
#'
#' This function plots the density curves of interaction frequencies 
#' from original and imputed single-cell Hi-C data for a specified distance.
#'
#' @param org_sc_data A data frame containing the original single-cell Hi-C data in scHiC Table object.
#'   The data frame should include columns for `region1`,`region2`,`Cell`, `Chr`, and `IF_i` for i single cells.
#' @param imp_sc_data A data frame containing the imputed single-cell Hi-C datain scHiC Table object.
#'   The data frame should include columns for `region1`,`region2`,`Cell`, `Chr`, and `IF_i` for i single cells.
#' @param D An integer specifying the genomic distance for which to plot the density 
#'   curves of interaction frequencies. 
#'
#' @return A ggplot2 object representing the density plot comparing original 
#'   and imputed interaction frequencies.
#' 
#' @examples
#' # Assuming scHiC.table_ODC and imp.cc.table.1 are your data frames:
#' plot_imputed_distance_diagnostic(scHiC.table_ODC, imp.cc.table.1, D = 3)
#' 
#' @import ggplot2
#' 
#' @export
plot_imputed_distance_diagnostic <- function(org_sc_data, imp_sc_data, D) {
  # Calculate the resolution
  res <- min(abs(diff(unique(org_sc_data$region1))))
  
  # Process original single-cell data
  org_sc_data$D <- (org_sc_data$region2 - org_sc_data$region1) / res
  org_D_data <- org_sc_data[org_sc_data$D == D, -c(1:4)]
  org_D_data <- unlist(org_D_data)
  org_D_data[org_D_data == 0] <- NA  # Replace 0 with NA
  
  # Process imputed single-cell data
  imp_sc_data$D <- (imp_sc_data$region2 - imp_sc_data$region1) / res
  imp_D_data <- imp_sc_data[imp_sc_data$D == D, -c(1:4)]
  imp_D_data <- unlist(imp_D_data)
  
  # Create a data frame for plotting
  df <- data.frame(
    IF = c(org_D_data, imp_D_data), # Interaction Frequencies (IF)
    Group = factor(rep(c("Original", "Imputed"), c(length(org_D_data), length(imp_D_data))))
  )
  
  # Plot density curves
  ggplot(df, aes(x = IF, fill = Group)) +
    geom_density(alpha = 0.5) +  # Use alpha for transparency to overlap the curves
    labs(title = paste("Density plot for D =", D), x = "Interaction Frequency", y = "Density") +
    scale_fill_manual(values = c("Original" = "blue", "Imputed" = "red")) +  # Customize colors
    theme_classic()  # Apply a clean theme
}

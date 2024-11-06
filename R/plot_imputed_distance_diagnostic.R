#' Plot Imputed Distance Diagnostic
#'
#' This function plots the distribution density curves of interaction frequencies
#' from original and imputed single-cell Hi-C data at a given genomic distance.
#'
#' @param raw_sc_data A data frame containing the original single-cell Hi-C data in scHiC Table format.
#'   The data frame should include columns for `region1`, `region2`, `Cell`, `Chr`, and `IF_i` for each of the "i" cells.
#' @param imp_sc_data A data frame containing the imputed single-cell Hi-C data in scHiC Table format.
#'   The data frame should include columns for `region1`, `region2`, `Cell`, `Chr`, and `IF_i` for each of the "i" cells.
#' @param D An integer specifying the genomic distance for which to plot the density curves of interaction frequencies.
#'   Genomic distance refers to the distance between two regions in the genome, expressed in units of resolution bins
#'   (e.g., D = (start2 - start1)/resolution). For example, if the genomic regions are at positions 16,000,000 and 17,000,000
#'   with a resolution of 1,000,000, then D = (17,000,000 - 16,000,000)/1,000,000 = 1.
#' @details
#' The distance D represents how far apart two genomic loci are. Observations indicate that chromatin interaction frequency (IF) in Hi-C data tends to decrease as the genomic distance between two loci increases. It is assumed that scHi-C interaction frequencies at the same distance share similar statistical properties. To diagnose whether the imputed IF values retain these statistical properties at each distance, the
#' `plot_imputed_distance_diagnostic()` function plots the density curves of interaction frequencies from original and imputed single-cell Hi-C data for a specified genomic distance.
#'
#' @return A ggplot2 object representing the density plot comparing the original
#'   and imputed interaction frequencies.
#'
#' @examples
#' data("scHiC.table_MG_chr22")
#' ## Impute data above
#' scHiC.table_MG_imp <- scHiCcompare_impute(scHiC.table = scHiC.table_MG_chr22)
#'
#' # Call the function with this sparse matrix to generate the plot
#' plot_imputed_distance_diagnostic(
#'   raw_sc_data = scHiC.table_MG_chr22,
#'   imp_sc_data = scHiC.table_MG_imp, D = 1
#' )
#'
#' @import ggplot2
#'
#' @export
plot_imputed_distance_diagnostic <- function(raw_sc_data, imp_sc_data, D) {
  # Calculate the resolution
  res <- min(abs(diff(unique(raw_sc_data$region1))))

  # Process original single-cell data
  raw_sc_data$D <- (raw_sc_data$region2 - raw_sc_data$region1) / res
  org_D_data <- raw_sc_data[raw_sc_data$D == D, -c(1:4)]
  org_D_data <- unlist(org_D_data)
  org_D_data[org_D_data == 0] <- NA # Replace 0 with NA

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
  plot <- ggplot2::ggplot(df, aes(x = IF, fill = Group)) +
    geom_density(alpha = 0.5) + # Use alpha for transparency to overlap the curves
    labs(title = paste0("Diagnostic Density for imputed D = ", D), x = " Interaction Frequency", y = " Density") +
    scale_fill_manual(values = c("Original" = "blue", "Imputed" = "red")) + # Customize colors
    theme_classic() # Apply a clean theme

  return(plot)
}

##### Identify best A threshold in HiCcompare
best_A <- function (hic.table, SD = 2, numChanges = 35, FC = 3, alpha = 0.05) 
{
  # Check if hic.table is not a list
  if (is(hic.table, "list")) {
    stop("Error: Enter a single hic.table object, not a list of hic.tables.")
  }
  
  # Check for valid numChanges and FC
  if (!is.numeric(numChanges) || numChanges <= 0) {
    stop("Error: 'numChanges' must be a positive numeric value.")
  }
  
  if (!is.numeric(FC) || FC <= 0) {
    stop("Error: 'FC' (fold change) must be a positive numeric value.")
  }
  
  if (!is.numeric(alpha) || alpha <= 0 || alpha > 1) {
    stop("Error: 'alpha' must be a numeric value between 0 and 1.")
  }
  
  
  # Remove rows where abs(M) >= SD
  new.table <- new.table[abs(M) < SD, ]
  
  # Define the sample space and exclude low A values
  sample_space <- 1:nrow(new.table)
  tmp_A <- (new.table$IF1 + new.table$IF2) / 2
  low_A <- which(tmp_A < quantile(tmp_A, 0.1))
  changes <- sample(sample_space[-low_A], numChanges)
  
  # Apply mean interaction frequencies for selected changes
  meanIF <- ((new.table[changes, ]$IF1 + new.table[changes, ]$IF2) / 2) %>% round() %>% as.integer()
  suppressWarnings(new.table[changes, `:=`(IF1, meanIF)])
  suppressWarnings(new.table[changes, `:=`(IF2, meanIF)])
  
  # Split the changes for applying the fold change (FC)
  midpoint <- floor(numChanges / 2)
  newIF1 <- new.table[changes[1:midpoint], ]$IF1 * FC %>% as.integer()
  newIF2 <- new.table[changes[(midpoint + 1):numChanges], ]$IF2 * FC %>% as.integer()
  
  # Update the table with new interaction frequencies
  new.table[changes[1:midpoint], `:=`(IF1, newIF1)]
  new.table[changes[(midpoint + 1):numChanges], `:=`(IF2, newIF2)]
  new.table = new.table[, `:=`(M, log2(IF2 / IF1))]
  
  # Add truth column to track changes
  truth <- rep(0, nrow(new.table))
  truth[changes] <- 1
  new.table[, `:=`(truth, truth)]
  
  # Normalize the table with hic_loess
  new.table <- hic_loess(new.table, Plot = Plot)
  new.table <- suppressMessages(hic_compare(new.table, Plot = Plot))
  
  # Initialize vectors for performance metrics
  TP <- vector(length = 50)
  FP <- vector(length = 50)
  FN <- vector(length = 50)
  TN <- vector(length = 50)
  A_seq <- seq(1, 50, by = 1)
  
  # Loop through A values and compute metrics
  for (i in seq_along(A_seq)) {
    tmp.table <- suppressMessages(hic_compare(new.table, A.min = A_seq[i], adjust.dist = TRUE, p.method = "fdr", Plot = FALSE))
    
    TP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 1)
    FP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 0)
    FN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 1)
    TN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 0)
  }
  
  # Calculate performance metrics
  MCC <- ((TP * TN) - (FP * FN)) / (sqrt((TP + FP)) * sqrt((TP + FN)) * sqrt((TN + FP)) * sqrt((TN + FN)))
  FPR <- FP / (FP + TP)
  FNR <- FN / (FN + TN)
  TPR <- TP / (TP + FP)
  
  # Identify the best A values
  MCC.A <- which(MCC == max(MCC))
  TPR.A <- which(TPR == max(TPR))
  FPR.A <- which(FPR == min(FPR))
  intersect.MCC_TPR.A <- intersect(MCC.A, TPR.A)
  best_A <- intersect(intersect.MCC_TPR.A, FPR.A)[1]
  
  message("Best A threshold filtered: ", best_A)
  
  return(best_A)
}



#### GMM cluster layer
GMM_layer <- function(hic_table, D.interval = 1:10, threshold = 0.8) {
  hic_result <- NULL
  if(D.interval == 'full'){
    D.interval = unique(hic_table$D)
    D.interval = D.interval[!D.interval == 0]
  }
  library(mclust)
  for (d in D.interval) {
    # Subset the table for the current distance 'd'
    hic_d <- hic_table[hic_table$D == d, ]
    
    # Initialize the 'p.value_final' column
    hic_d$p.value_final <- 1
    
    # Extract the 'adj.M' values
    x <- hic_d$adj.M
    
    # Pre-existing significant values (p < 0.05)
    hiccompare_sig <- hic_d$adj.M[hic_d$p.value < 0.05]
    
    # Perform Shapiro-Wilk normality test
    if(length(x) >2){
      norm_test <- shapiro.test(x)
      p_norm.test <- norm_test$p.value
    } else {p_norm.test <- 1}
    
    # If data is not normally distributed, fit GMM
    if (p_norm.test < 0.05) {
      # Fit a Gaussian Mixture Model with three components
      gmm_model <- Mclust(x, G = 3)
      
      # Extract the means of the components and rank them
      means <- gmm_model$parameters$mean
      mean_ranks <- rank(means)
      
      # Classify data points based on the GMM classification
      classification <- gmm_model$classification
      labels <- ifelse(classification == which(mean_ranks == 3), 1,  # Highest mean
                       ifelse(classification == which(mean_ranks == 2), 2, 0))  # Middle and lowest means
      
      # Create a dataframe with 'x' values and assigned labels
      data_with_labels <- data.frame(x = x, label = as.factor(labels))
      
      # Identify significant values (based on thresholds)
      sig <- data_with_labels$x[data_with_labels$label %in% c(0, 1)]
      sig <- sig[sig > threshold | sig < -threshold]
      
      # Combine with pre-existing significant values
      final_sig <- unique(c(sig, hiccompare_sig))
      
    } else {
      # If data is normally distributed, use pre-existing significant values
      final_sig <- hiccompare_sig
    }
    
    # Find positions of significant values in 'x'
    final_sig_pos <- which(x %in% final_sig)
    
    # Update 'p.value_final' for the significant positions
    hic_d$p.value_final[final_sig_pos] <- 0
    
    # Append the results to 'hic_result'
    hic_result <- rbind(hic_result, hic_d)
  }
  
  return(hic_result)
}



#### Differential plot 
require(ggplot2)
differential_result_plot <- function(hic.table.result){
  ggplot(hic.table.result, aes(y = adj.M, x = D, fill = factor(Difference.cluster))) +
    geom_point(shape = 21, size = 1.5) +  # Shape 21 allows fill color
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at 0
    scale_fill_manual(values = c("red", "black")) +  # Set custom fill colors for clusters
    theme_classic() +
    theme(legend.position = "none") + 
    ggtitle('MD Plot of Differential Result')  # Use ggtitle() for plot title
}









#' Compare Bulk Hi-C Data
#'
#' This function compares single-cell Hi-C data between two groups using the `scHicCompare` differential analysis workflow. It detects chromatin interaction differences between the single-cell Hi-C data of two cell types or conditions.
#'
#' @param norm.hic.table A data frame representing a jointly normalized pseudo-bulk Hi-C table output from two conditions, generated by the `hic_loess` function. 
#'                  This should be a pre-processed table that has been jointly normalized before differential analysis.
#' @param D.interval A numeric vector defining the distance intervals to consider in the analysis, or a character string 'full' indicating the inclusion of all genomic distances in the analysis. The distance should be scaled by dividing by the data resolution, D = (start2 - start1)/resolution (e.g., D = (16,000,000 - 17,000,000)/1,000,000 -> D = 1).
#' @param fprControl.logfc A numeric value controlling the false positive rate by setting the threshold for the log fold change. Detected differences should have a value larger than this threshold. 
#'                       Default is 0.8.
#' @param alpha A numeric value for the significance level of outlier detection. Default is 0.05.
#' @param A.min A numeric value or NA, specifying the A-value quantile cutoff to filter lower average expression in the `hic_compare` function from HiCcompare. 
#'               If not provided, an optimized minimum A threshold that maximizes MCC and TPR while minimizing FPR is used.
#' @param SD A numeric value specifying the standard deviation threshold for fuzzing, used to produce a Hi-C matrix from data with few true differences. This value is used to identify the optimal 'A.min' quantile cutoff in the `filter_params` and `hic_compare` functions for detecting significant outliers. Default is 2.
#' @param numChanges An integer indicating the number of changes to add to the Hi-C matrix. This value is used to identify the optimal 'A.min' quantile cutoff in the `filter_params` and `hic_compare` functions for detecting significant outliers. The number of changes should be proportional to the resolution of the data. High-resolution data should use more changes (e.g., 1MB resolution - 300 changes, 100KB resolution - 1000 changes, etc.). Default is 300.
#' @param FC A numeric value representing the fold change threshold added to the Hi-C matrix. This value is used to identify the optimal 'A.min' quantile cutoff in the `filter_params` and `hic_compare` functions for detecting significant outliers. Default is 3.
#' @param Plot A logical value indicating whether to plot the differential results in an MD plot. Default is TRUE.
#' @param parallel A logical value indicating whether to use parallel processing for computations. 
#'                 Default is FALSE. This option only works on Unix-based operating systems and is useful when analyzing a list of Hi-C tables.
#' @param BP_param Parameters for `BiocParallel`, to be passed to the `bpparam()` function.
#'
#' @return A data frame containing the results of the differential Hi-C analysis, including a normalized report and a 
#'         'Difference.cluster' column indicating the clusters of differences identified in the analysis.
#'
#' @details 
#' The `BulkHiC_compare` function performs differential chromatin interaction comparisons between single-cell Hi-C data from two groups. The workflow includes clustering normalized log fold changes between interaction frequencies into "difference" and "non-difference" groups. The non-difference group is assumed to follow a normal distribution centered around 0 and is clustered by a Gaussian Mixture Model. The difference cluster consists of points belonging to other distributions. In cases where the size of the differences is not large enough to form distinct distributions, these differences are likely outliers of the normal distribution, which are identified by the `hic_compare` function.
#' 
#' @references
#' 
#' Stansfield JC, Cresswell KG, Vladimirov VI  et al (2018). Hiccompare: an R-package for joint normalization and comparison of HI-C datasets. BMC Bioinformatics  2018;19:279.
#' 
#' Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models, The R Journal, 8/1, pp. 289-317.
#' 
#' C. Fraley and A. E. Raftery (2007) Bayesian regularization for normal mixture estimation and model-based clustering. Journal of Classification, 24, 155-181.
#' 
#' Fraley C. and Raftery A. E. (2002) Model-based clustering, discriminant analysis and density estimation, Journal of the American Statistical Association, 97/458, pp. 611-631.
#' 
#' Patrick Royston (1982). Algorithm AS 181: The W test for Normality. Applied Statistics, 31, 176–180. doi:10.2307/2347986.
#' 
#' Patrick Royston (1982). An extension of Shapiro and Wilk's W test for normality to large samples. Applied Statistics, 31, 115–124. doi:10.2307/2347973.
#' 
#' 
#' @examples
#' \dontrun{
#' # Load data folder example to current working directory
#' load_example_MGfolder()
#' load_example_ODCfolder()
#' # Input single-cell Hi-C in sparse format (.txt) from a path
#' scHiC.table_ODC <- scHiC_table(file.path = "ODCs_example", 
#'                                cell.type = 'ODC', position.dataset =  1:50, type = 'txt', 
#'                                select.chromosome = 'chr22')
#' scHiC.table_MG <- scHiC_table(file.path = "MGs_example", 
#'                               cell.type = 'MG', position.dataset =  1:50, type = 'txt', 
#'                               select.chromosome = 'chr22')
#' # Bulk matrix in sparse format
#' bulk.sparse.1 <- na.omit(pseudo_bulkHic(scHiC.table = scHiC.table_ODC, out = 'sparse'))
#' bulk.sparse.2 <- na.omit(pseudo_bulkHic(scHiC.table = scHiC.table_MG, out = 'sparse'))
#' # Create the `hic.table` object
#' bulk.hic.table <- create.hic.table(bulk.sparse.1, bulk.sparse.2, chr = 'chr22', scale = FALSE)
#' # Jointly normalize data for a single chromosome
#' hic.table_normalize <- hic_loess(bulk.hic.table, Plot = TRUE, Plot.smooth = FALSE)
#' # Example usage of the BulkHiC_compare function
#' result <- BulkHiC_compare(hic.table_normalize, D.interval = c(1, 100), fprControl.logfc = 0.8)
#' }
#'
#' @export


scHiC_bulk_compare <- function(norm.hic.table, D.interval, fprControl.logfc = 0.8,  alpha = 0.05,
                               SD = 2, numChanges = 300, FC = 3, A.min = NA,
                               Plot = T,  parallel = FALSE, BP_param = bpparam()){
  
  if(!is.na(A.min)){
    A.min <- best_A(norm.hic.table, SD = SD, numChanges = numChanges, FC = FC, alpha = alpha)
  }
  
  hic.table_result <- hic_compare(norm.hic.table, A.min = A.min,Plot = F, Plot.smooth = F,
                                  parallel = parallel, BP_param =BP_param)
  hic.table.GMM_result <- GMM_layer(hic_table = hic.table_result, D.interval = D.interval, threshold = fprControl.logfc)
  hic.table.GMM_result <- hic.table.GMM_result[,-c(17,18)]
  names(hic.table.GMM_result)[ncol(hic.table.GMM_result)] <- 'Difference.cluster'
  
  if(Plot == T){
    plot = differential_result_plot(hic.table.GMM_result)
    print(plot)
  }
  
  return(hic.table.GMM_result)
}





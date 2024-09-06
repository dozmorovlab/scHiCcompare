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



BulkHiC_compare <- function(hic.table, SD = 2, numChanges = 35, FC = 3, alpha = 0.05,
                            A.min = NA, adjust.dist = TRUE,
                            p.method = "fdr", Plot = FALSE, Plot.smooth = TRUE,
                            parallel = FALSE, BP_param = bpparam()){
  if(!is.na(A.min)){
    A.min <- best_A(hic.table, SD = SD, numChanges = numChanges, FC = FC, alpha = alpha)
  }
  
  hic.table_result <- hic_compare(hic.table_cc, A.min = A.min, adjust.dist = adjust.dist,
                            p.method = p.method, Plot = Plot, Plot.smooth = Plot.smooth,
                            parallel = parallel, BP_param =BP_param)
  return(hic.table_result)
}





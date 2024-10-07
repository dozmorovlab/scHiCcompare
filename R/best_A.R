##### Identify best A threshold in HiCcompare
# Randomize_IFs - function to add noise to IFs of one matrix to create a similar 2 matrix
randomize_IFs <- function(hic.table, SD) {
  # copy first IF vector
  newIF2 <- hic.table$IF1
  # add constant offset
  newIF2 <- newIF2 + 5
  # add random noise
  newIF2 <- newIF2 + rnorm(length(newIF2), 0, SD)
  # check for 0's and negatives
  newIF2[newIF2 <= 0] <- 1
  # create new hic.table with new IF vectors
  sparse1 <- cbind(hic.table$start1, hic.table$start2, hic.table$IF1)
  sparse2 <- cbind(hic.table$start1, hic.table$start2, newIF2)
  temp.table <- create.hic.table(sparse1, sparse2, chr = hic.table$chr1[1])
  return(temp.table)
}


#' Find the Best Quantile for A quantile cutoff in `filter_params` {HiCcompare}
#'
#' This function identifies the best quantile for the parameter A quantile cutoff in `filter_params` {HiCcompare} by evaluating 
#' different quantile levels based on performance metrics such as 
#' Matthews Correlation Coefficient (MCC), True Positive Rate (TPR), and 
#' False Positive Rate (FPR).
#'
#' @param hic.table A hic.table object
#' @param SD Numeric. The standard deviation of the fuzzing used to produce a Hi-C matrix from your data with few true differences.
#' @param numChanges Integer. The number of changes to add into the Hi-C matrix created. This should be proportional to the resolution of the data. High resolution data should use more changes i.e. 1MB resolution - 300 changes, 100KB resolution - 1000 changes, etc.
#' @param FC Numeric. The fold change of the changes added to the Hi-C matrix.
#' @param alpha Numeric. Significance level for adjusting p-values (default is 0.05).
#' @param Plot Logical. If TRUE, plots will be generated during the processing.
#'
#' @return A data frame containing the best quantile for A, the corresponding best A value,
#'         and performance metrics (MCC, TPR, FPR).
#'
#' @details The function randomizes interaction frequencies, calculates mean interaction 
#' frequencies, and performs a series of evaluations to determine the best quantile.
#' The results include the quantile value, best A value, and performance metrics for 
#' each evaluated quantile.
#'
#' @examples
#' # Assuming 'hic.table' is a valid Hi-C interaction frequency table
#' best_result <- best_quantile_A(hic.table, SD = 2, numChanges = 35, FC = 3, alpha = 0.05)
#' print(best_result)
#'
#' @export
best_A <- function(hic.table, SD = 2, numChanges = 35, FC = 3, alpha = 0.05,
                            Plot = FALSE) {
  
  if (is(hic.table, "list")) {
    stop("Enter a single hic.table object, not a list of hic.tables.")
  }
  
  new.table <- randomize_IFs(hic.table, SD)
  new.table <- new.table[abs(new.table$M) < 2, ]
  
  # Initialize vectors for performance metrics
  results <- data.frame(quantile = numeric(), best_A = integer(), MCC = numeric(), TPR = numeric(), FPR = numeric())
  
  # Loop over each quantile
  for (q in quantile_seq) {
    
    tmp_A <- (new.table$IF1 + new.table$IF2) / 2
    low_A <- which(tmp_A < quantile(tmp_A, q))
    sample_space <- 1:nrow(new.table)
    changes <- sample(sample_space[-low_A], numChanges)
    
    meanIF <- ((new.table[changes, ]$IF1 + new.table[changes, ]$IF2) / 2) %>%
      round() %>%
      as.integer()
    
    suppressWarnings(new.table[changes, `:=`(IF1, meanIF)])
    suppressWarnings(new.table[changes, `:=`(IF2, meanIF)])
    
    midpoint <- floor(numChanges / 2)
    newIF1 <- new.table[changes[1:midpoint], ]$IF1 * FC %>% as.integer()
    newIF2 <- new.table[changes[(midpoint + 1):numChanges], ]$IF2 * FC %>% as.integer()
    
    new.table[changes[1:midpoint], `:=`(IF1, newIF1)]
    new.table[changes[(midpoint + 1):numChanges], `:=`(IF2, newIF2)]
    new.table[, `:=`(M, log2(IF2 / IF1))]
    
    truth <- rep(0, nrow(new.table))
    truth[changes] <- 1
    new.table[, `:=`(truth, truth)]
    
    new.table <- hic_loess(new.table, Plot = Plot)
    new.table <- suppressMessages(hic_compare(new.table, Plot = Plot))
    
    TP <- vector(length = 50)
    FP <- vector(length = 50)
    FN <- vector(length = 50)
    TN <- vector(length = 50)
    A_seq <- seq(1, 50, by = 1)
    
    for (i in seq_along(A_seq)) {
      tmp.table <- suppressMessages(hic_compare(new.table, 
                                                A.min = A_seq[i], adjust.dist = TRUE, p.method = "fdr", 
                                                Plot = FALSE))
      TP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 1)
      FP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 0)
      FN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 1)
      TN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 0)
    }
    
    # Calculate performance metrics
    MCC <- ((TP * TN) - (FP * FN)) / (sqrt((TP + FP)) * sqrt((TP + FN)) * sqrt((TN + FP)) * sqrt((TN + FN)))
    FPR <- FP / (FP + TP)
    TPR <- TP / (TP + FP)
    
    # Best A in terms of MCC, TPR, FPR
    MCC.A = which(MCC == max(MCC, na.rm = TRUE))
    TPR.A = which(TPR == max(TPR, na.rm = TRUE))
    FPR.A = which(FPR == min(FPR, na.rm = TRUE))
    
    intersect.MCC_TPR.A <- intersect(MCC.A, TPR.A)
    best_A <- intersect(intersect.MCC_TPR.A, FPR.A)[1]
    
    # Store results for current quantile
    results <- rbind(results, data.frame(quantile = q, best_A = best_A, MCC = max(MCC, na.rm = TRUE),
                                         TPR = max(TPR, na.rm = TRUE), FPR = min(FPR, na.rm = TRUE)))
  }
  
  # Find the overall best quantile based on MCC or other metrics
  best_result <- results[which.max(results$MCC), ]
  message("Best quantile for A: ", best_result$quantile, " with best A: ", best_result$best_A)
  
  return(best_result)
}

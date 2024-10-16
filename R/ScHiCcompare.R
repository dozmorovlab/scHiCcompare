## Not-Normalize schic - hic table ----
withoutNorm_hicTable <- function(hic.table){
  adj.IF1 <- hic.table$IF1
  adj.IF2 <- hic.table$IF2
  adj.M <- hic.table$M
  mc <- rep(0, nrow(hic.table))
  table <- cbind(hic.table, adj.IF1, adj.IF2, adj.M, mc)
  table[, `:=`(A, (adj.IF1 + adj.IF2)/2)]
  return(table)
}




#' ScHiCcompare: Differential Analysis of Single-Cell Hi-C Data
#'
#' This function performs a differential analysis between two single-cell Hi-C data groups. It includes the 
#' steps of imputation, normalization, and detection of differential chromatin interactions (DCIs).
#'
#' @param file.path.1 Character string specifying the directory containing scHi-C data for the first condition (first cell-type group). The folder should contain '.txt' scHi-C files in sparse upper triangular format (chr, start1, end1, IF)
#' @param file.path.2 Character string specifying the directory containing Hi-C data for the second condition (second cell-type group). The folder should contain '.txt' scHi-C files in sparse upper triangular format (chr, start1, end1, IF)
#' @param imputation Character string or NULL of indicating the imputation method. Default is 'RF' for Random Forest imputation.
#' @param normalization Character string or NULL indicating the normalization method. Default is 'Loess'.
#' @param differential.detect Character string indicating the differential detection method. Default is 'MD.cluster'.
#' @param select.chromosome Integer or character indicating the chromosome to be analyzed (e.g., 'chr1' or 'chrX'.)
#' @param main.Distances Numeric vector or character 'full' indicating the interacting range (or full range) of genomic distances (in bp) for the method to focus.
#'  Genomic distances (in bp) is the number of base pairs between two regions in the genome (e.g., loci or bins). Default is 1 to 10,000,000.
#' @param save.output.path Character string specifying the directory to save outputs, including the imputed cells in the form of a sparse upper triangular format,
#'  normalization result table, and differential analysis result table. If NULL, no files are saved.
#' @param Plot Logical. If TRUE, a plot of the differential results will be generated. Default is TRUE.
#' @param pool.style Character string specifying the pooling style for `imputation`. Options are 'progressive' or 'Fibonacci'. Default is 'progressive'. If the 'imputation
#'  is skipped as NULL, the `pool.style` also should be NULL
#' @param n.imputation Integer specifying the number of imputations for the `imputation` step. Default is 5.
#' @param maxit Integer specifying the maximum number of iterations for internal refinement process within a single `imputation` cycle. Default is 1.
#' @param outlier.rm Logical. If TRUE, outliers are removed during `imputation`. Default is TRUE.
#' @param missPerc.threshold Numeric value specifying the maximum allowable percentage of missing data in each pool band.
#'  Only pool bands within the `main.Distances` range and having a missing percentage below this threshold will be imputed
#'  by the `imputation` method. Default is 95%.
#' @param fprControl.logfc A numeric value controlling the false positive rate of `differential.detect` step by setting the threshold for the log fold change in the 'difference' cluster. 
#'  Detected differences identified by Gaussian Mixed Model (GMM) clusters only include values with log fold change that are larger than this threshold. Default is 0.8.
#' @param alpha A numeric value for the significance level of outlier detection of `differential.detect` step by the analysis of `hic_compare()` function from HiCcompare. Default is 0.05. 
#' @param A.min A numeric value or NULL, specifying the A-value quantile cutoff to filter lower average expression in `differential.detect` step of `hic_compare()` function (from HiCcompare). 
#'  `hic_compare()` is used to detect outliers, which is assumed to be 'differences' bins in case of its number is too small (or none) to be cluster by GMM method.
#'  If not provided, an optimized minimum A threshold that maximizes MCC and TPR while minimizing FPR in the simulated Hi-C matrix.
#' @param SD A numeric value specifying the standard deviation threshold for fuzzing, used to produce a simulated Hi-C matrix. This value is used to modify the process finding optimal
#'  'A.min' quantile cutoff for detecting significant outliers in `differential.detect` step. Users can select the value based on their assumption of the scHi-C data. Default is 2.
#' @param numChanges An integer or NULL, indicating the number of changes to add to the simulated Hi-C matrix. This value is used to modify the process finding the optimal 'A.min' quantile cutoff
#'  for detecting significant outliers. Based on the users assumption about possible number of difference, they can set the number of changes that should be proportional
#'  to the resolution of the data. High-resolution data should be assumed more changes. If `numChanges` = NULL, the function is set number of changes (or simulated difference) is scaled by a factor of 30
#'  (e.g., 1MB resolution - 30 changes, 500KB resolution - 60 changes, etc.) Default is NULL. 
#' @param FC A numeric value representing the fold change threshold added to the simulated Hi-C matrix. This value is used to identify the optimal 'A.min' quantile cutoff for detecting significant outliers
#'  in `differential.detect` step. Users can select the FC value based on their assumption of difference fold change in their data. Default is 3. 
#' @param Plot A logical value indicating whether to plot the `differential.dect` results in an MD plot. Default is TRUE.
#' @param Plot.normalize A logical value indicating whether to plot the `normalization` results in an MD plot. Default is FALSE.
#' @param BP_param Parameters for `BiocParallel`, to be passed to the `bpparam()` function. See `?bpparam()` for more info.
#' 
#' @details
#' 
#' This function implements the ScHiCcompare workflow. It first reads sparse Hi-C data from two conditions 
#' and imputes missing interaction frequencies (if specified) using a random forest model (RF) along with a chosen 
#' pooling method (either progressive or Fibonacci). In progressive pooling, distances are combined consecutively 
#' to form larger sets, while Fibonacci pooling uses a Fibonacci sequence for combination.
#' 
#' Next, pseudo-bulk Hi-C matrices are generated, followed by joint normalization using Loess regression (from HiCcompare) 
#' before detecting differential chromatin interactions via a Gaussian Mixture Model (GMM) clustering approach. 
#' 
#' The differential analysis clusters normalized log fold changes in interaction frequencies between the two cell types 
#' at each genomic distance into "difference" and "non-difference" groups. The non-difference group is assumed to 
#' follow a normal distribution centered around 0 and is clustered using a Gaussian Mixture Model. 
#' The difference cluster comprises points that belong to other distributions. If the size of the differences is 
#' insufficient to form distinct distributions, these differences are considered outliers of the normal distribution, 
#' identified by the `hic_compare()` function.

#' 
#'
#' @return A list containing the differential analysis results and intermediate results (imputation, pseudo-bulk, normalization).
#' If `save.output.path` is provided , the imputed results for both conditions are saved in sparse format in the given firectory. Normalization and differential analysis results 
#' are also saved if `save.output.path` is provided. 
#'
#' @examples
#' \dontrun{
#' Load_example_MGFolder()
#' Load_example_ODCFolder()
#' result <- ScHiCcompare(
#'   file.path.1 = "MGs_example", 
#'   file.path.2 = "ODCs_example", 
#'   select.chromosome = "chr22", 
#' )
#' print(result)
#' }
#'
#' @export


ScHiCcompare <- function(file.path.1, file.path.2, imputation = 'RF', normalization = 'Loess', differential.detect = 'MD.cluster',
                         select.chromosome, main.Distances = 1:10000000, save.output.path =  NULL, Plot = T, Plot.normalize = F,
                         pool.style = 'progressive' ,n.imputation = 5,  maxit = 1, outlier.rm = TRUE, missPerc.threshold = 95,
                         A.min = NULL, fprControl.logfc = 0.8, alpha = 0.05, SD = 2, numChanges = 30, FC = 3, 
                         BP_param = bpparam()){
  
  # Read file 'txt' from 2 folder path
  # cond1_list <- read_files(file.path = file.path.1, type='txt',
  #                        txt.sparse.heads.position = c(1,2,3,4), out = 'sparse')
  # cond2_list <- read_files(file.path = file.path.2, cell = 'condition2', type='txt',
  #                          txt.sparse.heads.position  = c(1,2,3,4), out = 'sparse')
    
  
  # Step 0 : Transfer into scHiC table oject
  scHiC.table_cond1 <- scHiC_table(file.path = file.path.1, cell.type = 'condition1',
                                   select.chromosome = select.chromosome)
  scHiC.table_cond2 <- scHiC_table(file.path = file.path.2, cell.type = 'condition2',
                                   select.chromosome = select.chromosome)
  ## main distance of scHiC table - transfer to D
  res = min( abs( diff(unique(scHiC.table_cond1$region1)) ) )
  D.max = abs(max(main.Distances))/res 
  D.min = max(1, abs(min(main.Distances))/res)
  # Check if D.min and D.max are whole numbers
  if (D.max %% 1 != 0) {
    message("Warning: The end of main.Distances is not proportional to the resolution. Adjust the input.")
  }
  
  if (D.min %% 1 != 0) {
    message("Warning: The start of main.Distances is not proportional to the resolution. Adjust the input.")
  }
  D.interval = c(D.min : D.max)

  
  
  # Step 1: Imputation 
  if (!is.null(imputation)){
    cat('Imputing Condition 1 group cells in: ')
    scHiC.table_cond1 <- Pooling_RF_impute(scHiC.table = scHiC.table_cond1, n.imputation = n.imputation,  maxit = maxit, outlier.rm = outlier.rm, 
                             main.Distances = main.Distances, pool.style = pool.style, missPerc.threshold = missPerc.threshold)
    cat('\nImputing Condition 2 group cells in: ')
    scHiC.table_cond2 <- Pooling_RF_impute(scHiC.table = scHiC.table_cond2, n.imputation = n.imputation,  maxit = maxit, outlier.rm = outlier.rm, 
                                   main.Distances = main.Distances, pool.style =pool.style, missPerc.threshold = missPerc.threshold)
    impute1_result = scHiC.table_cond1; impute2_result = scHiC.table_cond2
  } else {
    impute1_result = NULL; impute2_result = NULL
  }
  
  
  
  # Step 2: Pseudobulk
  bulk_sparse_cond1 = pseudo_bulkHic(scHiC.table = scHiC.table_cond1, out = 'sparse')
  bulk_sparse_cond2 =  pseudo_bulkHic(scHiC.table = scHiC.table_cond2, out = 'sparse')
 
  
  # Step 2: Normalization
  if (is.null(normalization)){
    library(HiCcompare)
    bulk.hic.table <- create.hic.table(bulk_sparse_cond1, bulk_sparse_cond2, chr = select.chromosome, scale = F)
    #jointly normalize data for a single chromosome
    norm.hic.table <- withoutNorm_hicTable(hic.table = bulk.hic.table)
    norm.result <- NULL
  } else {
    bulk.hic.table <- create.hic.table(bulk_sparse_cond1, bulk_sparse_cond2, chr = select.chromosome, scale = F)
    #jointly normalize data for a single chromosome
    cat('\nJointly normalizing pseudo bulk matrices ')
    norm.hic.table <- suppressWarnings(hic_loess(bulk.hic.table, Plot = Plot.normalize, Plot.smooth = F))
    norm.result <- norm.hic.table
    names(norm.result)[c(7,8,11,12)] <- c("bulk.IF1", "bulk.IF2" ,"adj.bulk.IF1" ,"bulk.adj.IF2")
  }
  
  

  library(HiCcompare)
  cat('\nProcessing detect differential chromotin interaction ')
  if(is.null(A.min)){
    A.min <-  suppressWarnings(best_A(hic.table = norm.hic.table, SD = SD, numChanges = numChanges, FC = FC, alpha = alpha))
  }
  # HiCcompare
  hic.table_result <- suppressWarnings(hic_compare(norm.hic.table, A.min = A.min,Plot = F, Plot.smooth = F,
                                  BP_param =BP_param))

  ## add GMM layer
  hic.table.GMM_result <- suppressWarnings(GMM_layer(hic_table = hic.table_result, D.interval = D.interval, threshold = fprControl.logfc))
  hic.table.GMM_result <- hic.table.GMM_result[,-c(17,18)]
  names(hic.table.GMM_result)[ncol(hic.table.GMM_result)] <- 'Difference.cluster'
  names(hic.table.GMM_result)[c(7,8,11,12)] <- c("bulk.IF1", "bulk.IF2" ,"adj.bulk.IF1" ,"bulk.adj.IF2")
  
  
  
  ## Save output option
  if(!is.null(save.output.path)){
    ##### Imputed cell #####
    library(tidyverse)
    library(data.table)
    
    
    if(!is.null(impute1_result)){
     
      ##### Group 1 ##### 
      df = impute1_result
      folder_name <- "Imputed_Condition1_cells"
      full_output_path <- file.path(save.output.path, folder_name)
      dir.create(full_output_path, recursive = TRUE)
      cat("\nImputed cells in condition1 saved", "to:", full_output_path, "\n")
      
      # Transform to long format, rename cells, and select necessary columns in one step
      sparse_df <- df %>%
        pivot_longer(cols = starts_with("IF_"),
                     names_to = "cell_id",
                     values_to = "IF") %>%
        mutate(cell_id = sub("IF_", "cell", cell_id)) %>%  # Change IF_1, IF_2 to cell1, cell2, etc.
        select(cell_id,chr, region1, region2, IF)  # Changed 'IF' to 'IF_value'
      
      
      # Split the data into a list of data frames by cell
      cell_data_list <- split(sparse_df, sparse_df$cell_id)
      
      # Loop through each cell data frame and save to individual .txt files using fwrite
      lapply(names(cell_data_list), function(cell) {
        # Define the output file path for the current cell
        output_file_path <- file.path(full_output_path, paste0('grp1.imp_',cell, ".txt"))
        
        # Save the data frame to a .txt file
        #fwrite(cell_data_list[[cell]], file = output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)
        write.table(cell_data_list[[cell]][,-1], output_file_path,  row.names = FALSE, quote = FALSE)

      })
    }
    
    ##### Group 2 ##### 
    df = impute2_result
    folder_name <- "Imputed_Condition2_cells"
    full_output_path <- file.path(save.output.path, folder_name)
    dir.create(full_output_path, recursive = TRUE)
    cat("\nImputed cells in condition2 saved", "to:", full_output_path, "\n")
    
    # Transform to long format, rename cells, and select necessary columns in one step
    sparse_df <- df %>%
      pivot_longer(cols = starts_with("IF_"),
                   names_to = "cell_id",
                   values_to = "IF") %>%
      mutate(cell_id = sub("IF_", "cell", cell_id)) %>%  # Change IF_1, IF_2 to cell1, cell2, etc.
      select(cell_id, chr, region1, region2, IF)  # Changed 'IF' to 'IF_value'
    
    # Split the data into a list of data frames by cell
    cell_data_list <- split(sparse_df, sparse_df$cell_id)
    
    # Loop through each cell data frame and save to individual .txt files using fwrite
    lapply(names(cell_data_list), function(cell) {
      # Define the output file path for the current cell
      output_file_path <- file.path(full_output_path, paste0('grp2.imp_',cell, ".txt"))
      
      # Save the data frame to a .txt file
      write.table(cell_data_list[[cell]][,-1], output_file_path,  row.names = FALSE, quote = FALSE)
    
    })
    
    ##### Normalized result #####
    df = norm.result
    cat("\nNormalization result table saved", "to:", save.output.path, "\n")
    # Define the output file path for the current cell
    output_file_path <- file.path(save.output.path, paste0("Bulk_normalization_table.txt"))
    # Save the data frame to a .txt file
    write.table(df, output_file_path, quote = FALSE)
    
    ##### Differential result #####
    df = hic.table.GMM_result
    cat("\nDifferential analysis result table saved", "to:", save.output.path, "\n")
    # Define the output file path for the current cell
    output_file_path <- file.path(save.output.path, paste0("Differential_analysis_table.txt"))
    # Save the data frame to a .txt file
    write.table(df, output_file_path, quote = FALSE)
  }
  
 
  
  
  
  # Print result
  if(Plot == T){
    plot = differential_result_plot(hic.table.GMM_result)
    print(plot)
  }
  
  result <- list(
    Differential_Analysis = hic.table.GMM_result,
    Intermediate = list(
      Imputation = list(condition1 = impute1_result, condition2 = impute2_result),
      PseudoBulk = list(condition1 = bulk_sparse_cond1, condition2 = bulk_sparse_cond2),
      Bulk.Normalization = norm.result
    )
  )
  
  # Assign a custom class to the result
  class(result) <- "checkNumbers"
  
  return(result)
}



print.checkNumbers <- function(obj) {
  # Print a header for the output
  cat("---------------------------------------------------------------------------------\n")
  cat("           ScHiCcompare - Differential Analysis for single cell Hi-C\n")
  cat("---------------------------------------------------------------------------------\n")
  
  # Access the necessary information from the object
  cond1_cells <- ncol(obj$Intermediate$Imputation$condition1) - 4  # Adjust if necessary
  cond2_cells <- ncol(obj$Intermediate$Imputation$condition2) - 4  # Adjust if necessary
  selected_chromosome <- unique(obj$Intermediate$Bulk.Normalization$chr1)  # Assuming chromosome info is stored here
  total_differences <- sum(obj$Differential_Analysis$Difference.cluster, na.rm = TRUE)  # Handle NA values
  
  # Determine which processes are included, avoiding NULL values
  impute <- if (!is.null(obj$Intermediate$Imputation$condition1)) {
    'Imputation'
  } else {
    NULL
  }
  
  norm <- if (!is.null(obj$Intermediate$BulkNormalization)) {
    'sc Bulk Normalization'
  } else {
    NULL
  }
  
  # Combine the processes and ensure only non-null values are printed
  processes <- na.omit(c(impute, norm))
  
  # Print the summary message
  cat(paste("ScHiCcompare analyzes", cond1_cells, "cells of condition 1 group and",  
            cond2_cells, "cells of condition 2 group at chromosome", selected_chromosome, "\n"))
  
  # Check if there are any processes to print
  if (length(processes) > 0) {
    cat("\nThe process includes:\n")
    cat(paste(processes, collapse = ", "), "Differential Analysis\n")
  } else {
    cat("\nNo processes available for analysis.\n")
  }
  
  
  # Note about accessing intermediate results
  cat("\nNote: See full differential result in $Differential_Analysis. Intermediate results can be accessed with $Intermediate\n")
}


 # result <- ScHiCcompare(file.path.1 = "MGs_example", file.path.2 = "ODCs_example",
 #                      select.chromosome = "chr22", save.output.path = 'Result')



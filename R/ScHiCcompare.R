utils::globalVariables(c("cell_id", "region1", "region2"))


## Not-Normalize schic - hic table ----
withoutNorm_hicTable <- function(hic.table) {
  adj.IF1 <- hic.table$IF1
  adj.IF2 <- hic.table$IF2
  adj.M <- hic.table$M
  mc <- rep(0, nrow(hic.table))
  table <- cbind(hic.table, adj.IF1, adj.IF2, adj.M, mc)
  A <- NULL
  table[, `:=`(A, (adj.IF1 + adj.IF2) / 2)]
  return(table)
}


#' ScHiCcompare: Differential Analysis of Single-Cell Hi-C Data
#'
#' This function performs a differential analysis between two single-cell Hi-C
#'  data groups. It includes the steps of imputation, normalization, and
#'  detection of differential chromatin interactions (DCIs).
#'
#' @param file.path.1 Required character string specifying the directory
#'  containing scHi-C data for the first condition (first cell-type group).
#'  The folder should contain '.txt' scHi-C files in modified sparse upper
#'  triangular format with 5 columns (chr1, start1, chr2, start2, IF).
#' @param file.path.2 Required character string specifying the directory
#'  containing Hi-C data for the second condition (second cell-type group).
#'  The folder should contain '.txt' scHi-C files in modified sparse upper
#'  triangular format with 5 columns (chr1, start1, chr2, start2, IF).
#' @param select.chromosome Required integer or character indicating the
#'  chromosome to be analyzed (e.g., 'chr1' or 'chr10').
#' @param imputation Character string 'RF' or NULL indicating the
#'  imputation method. Default is 'RF' for Random Forest imputation.
#' @param normalization Character string 'LOESS' or NULL indicating the
#'  normalization method. Default is 'LOESS'.
#' @param differential.detect Character string 'MD.cluster' indicating the
#'  differential detection method. The default is 'MD.cluster'.
#' @param main.Distances Numeric vector indicating the range of interacting
#'  genomic distances (in base pairs) between two regions (e.g., loci or bins)
#'  to focus on (e.g., 1:100000, Inf, etc). The `main.Distance` vector needs to
#'  be proportional to the data's resolution (e.g., for 10kb - 1:10000,
#'  1:50000, 1:100000, Inf, etc). Selecting a large distance range at higher
#'  resolution (e.g., below 200kb) can make the function take longer to run due
#'  to extreme sparsity. The default is 1:10000000.
#' @param pool.style Character string specifying the pooling style for
#'  `imputation`. Options are 'none', 'progressive', or 'Fibonacci'.
#'  The default is 'progressive'. If `imputation` is NULL, then `pool.style`
#'  should also be NULL.
#' @param n.imputation Integer specifying the number of multiple imputations
#'  for the imputation step, with final imputed values computed as the average
#'  of these multiple imputation values. Increasing the number of imputations
#'  enhances the accuracy of imputed values, though it may increase the
#'  imputation runtime. The default is 5.
#' @param maxit Integer specifying the maximum number of iterations for the
#'  internal refinement process within a single `imputation` cycle. Increasing
#'  `maxit` can help stabilize imputed values, though it may increase the
#'  imputation runtime. Default is 1.
#' @param outlier.rm Logical. If TRUE, outliers are removed during `imputation`.
#'  The default is TRUE.
#' @param missPerc.threshold Numeric value specifying the maximum allowable
#'  percentage of missing data in pool bands outside the `main.Distances` to be
#'  imputed by the `imputation` method. A higher threshold includes more sparse
#'  distances for imputation (e.g., above 95 percent), increasing memory and
#'  runtime, while a lower threshold (e.g., below 50 percent) might reduce the
#'  number of distances imputed. The default is 95.
#' @param A.min Numeric value or NULL that sets the A-value quantile cutoff
#'  (e.g., 7, 10, etc) for filtering low average interaction frequencies in
#'  outlier detection during the differential step of `hic_compare()`
#'  from `HiCcompare`. If not provided (NULL), A is auto-detected.
#' @param fprControl.logfc Numeric value controlling the false positive rate for
#'  GMM difference clusters (`differential.detect`) (e.g., 0.5, 0.8, 1, 1.5,
#'  etc). Increasing `fprControl.logfc` may reduce the false positive rate but
#'  can also reduce the number of chromatin interaction differences detected.
#'  Default is 0.8, equivalent to a 2-fold change.
#' @param alpha Numeric value for the significance level of outlier detection
#'  during the `differential.detect` step by `hic_compare()` from HiCcompare.
#'  The default is 0.05.
#' @param Plot Logical value indicates whether to plot the
#' `differential.detect` results in an MD plot. The default is TRUE.
#' @param Plot.normalize Logical value indicates whether to plot the
#'  `normalization` results in an MD plot. The default is FALSE.
#' @param save.output.path Character string specifying the directory to save
#'  outputs, including the imputed cells in a modified sparse upper triangular
#'  format, a normalization result table, and a differential analysis result
#'  table. If NULL, no files are saved. The default is NULL.
#' @param BPPARAM Parameters for `BiocParallel`, to be passed to the `bpparam()`
#'  function. See `?bpparam()` for more info.
#'
#' @details
#'
#' This function implements the ScHiCcompare workflow. It first reads sparse
#'  Hi-C data from two conditions and, by default, imputes missing interaction
#'  frequencies using a random forest model (RF) with the choice of
#'  `pool.style` (either progressive or Fibonacci). With the progressive
#'  pooling of interaction frequencies, genomic distance ranges increase
#'  linearly to form subsequent pooled bands, while Fibonacci pooling uses the
#'  Fibonacci sequence to increase the size of genomic distance ranges.
#'  Then, the random forest method is applied to individual genomic distance
#'  ('none' pooling) or pooled bands (when `pool.style` is selected).
#'
#' Next, pseudo-bulk Hi-C matrices are generated, followed by joint
#'  normalization using Loess regression (from `HiCcompare`) before detecting
#'  differential chromatin interactions via a Gaussian Mixture Model (GMM)
#'  clustering approach. GMM clusters normalized log fold changes in
#'  interaction frequencies between the two cell types at each genomic distance
#'  into "difference" and "non-difference" groups. The non-difference group is
#'  assumed to follow a normal distribution centered around 0.
#'  The difference cluster comprises points that belong to other distributions.
#'  If the size of the differences is insufficient to form distinct
#'  distributions, these differences are identified by the
#'  `HiCcompare::hic_compare()` function.
#'
#'
#' @return A custom class object ("checkNumbers") summarize information of
#'  workflow's steps. At the same time, the object result also contain the
#'  differential analysis results and intermediate results (imputation,
#'  pseudo-bulk, normalization). If `save.output.path` is provided, the imputed
#'  results for both conditions are saved in a sparse format in the given
#'  directory. Normalization and differential analysis results are also saved
#'  if `save.output.path` is provided. See the vignette for more details.
#'
#' @examples
#' ## Load example data for ODC and MG file paths
#' ODCs_example <- system.file("extdata/ODCs_example", package = "scHiCcompare")
#' MGs_example <- system.file("extdata/MGs_example", package = "scHiCcompare")
#'
#' ## Run scHiCcompare on example data
#' result <- scHiCcompare(
#'   file.path.1 = MGs_example,
#'   file.path.2 = ODCs_example,
#'   select.chromosome = "chr20"
#' )
#' print(result)
#'
#' @import HiCcompare
#' @import gtools
#' @import mclust
#' @importFrom rlang :=
#' @importFrom stats na.omit
#' @importFrom BiocParallel bpparam
#'
#' @export print.checkNumbers



scHiCcompare <- function(file.path.1, file.path.2, select.chromosome,
                         imputation = "RF", normalization = "LOESS",
                         differential.detect = "MD.cluster",
                         main.Distances = seq(10000000),
                         pool.style = "progressive", n.imputation = 5,
                         maxit = 1, outlier.rm = TRUE, missPerc.threshold = 95,
                         A.min = NULL, fprControl.logfc = 0.8, alpha = 0.05,
                         Plot = TRUE, Plot.normalize = FALSE,
                         save.output.path = NULL, 
                         BPPARAM = BiocParallel::bpparam()) {
  # Read file 'txt' from 2 folder path
  file_name1 <- list.files(file.path.1)
  file_name2 <- list.files(file.path.2)
  # Step 0 : Transfer into scHiC table oject
  scHiC.table_cond1 <- scHiC_table(
    file.path = file.path.1, cell.type = "condition1",
    select.chromosome = select.chromosome
  )
  scHiC.table_cond2 <- scHiC_table(
    file.path = file.path.2, cell.type = "condition2",
    select.chromosome = select.chromosome
  )
  ## main distance of scHiC table - transfer to D
  if (is.infinite(max(main.Distances))) {
    D.interval <- Inf
  } else {
    res <- min(abs(diff(unique(scHiC.table_cond1$region1))))
    D.max <- abs(max(main.Distances)) / res
    D.min <- max(1, abs(min(main.Distances)) / res)
    # Check if D.min and D.max are whole numbers
    D.interval <- c(D.min:D.max)
  }




  # Step 1: Imputation
  if (!is.null(imputation)) {
    message("Imputing Condition 1 group cells in: ")
    scHiC.table_cond1 <- scHiCcompare_impute(
      scHiC.table = scHiC.table_cond1, n.imputation = n.imputation,
      maxit = maxit, outlier.rm = outlier.rm,
      main.Distances = main.Distances, pool.style = pool.style,
      missPerc.threshold = missPerc.threshold
    )
    message("\nImputing Condition 2 group cells in: ")
    scHiC.table_cond2 <- scHiCcompare_impute(
      scHiC.table = scHiC.table_cond2, n.imputation = n.imputation,
      maxit = maxit, outlier.rm = outlier.rm,
      main.Distances = main.Distances, pool.style = pool.style,
      missPerc.threshold = missPerc.threshold
    )
    impute1_result <- scHiC.table_cond1
    impute2_result <- scHiC.table_cond2

    ### Change cell name
    imp_cell_name1 <- sub("\\.txt$", "", basename(file_name1))
    sorted_file_names1 <- gtools::mixedsort(imp_cell_name1)
    names(impute1_result)[5:ncol(impute1_result)] <- sprintf(
      "imp.IF_%s",
      sorted_file_names1
    )


    imp_cell_name2 <- sub("\\.txt$", "", basename(file_name2))
    sorted_file_names2 <- gtools::mixedsort(imp_cell_name2)
    names(impute2_result)[5:ncol(impute2_result)] <- sprintf(
      "imp.IF_%s",
      sorted_file_names2
    )
  } else {
    impute1_result <- NULL
    impute2_result <- NULL
  }



  # Step 2: Pseudobulk
  bulk_sparse_cond1 <- pseudo_bulkHic(
    scHiC.table = scHiC.table_cond1,
    out = "sparse"
  )
  bulk_sparse_cond2 <- pseudo_bulkHic(
    scHiC.table = scHiC.table_cond2,
    out = "sparse"
  )


  # Step 2: Normalization
  if (is.null(normalization)) {
    bulk.hic.table <- HiCcompare::create.hic.table(bulk_sparse_cond1,
      bulk_sparse_cond2,
      chr = select.chromosome, scale = FALSE
    )
    # jointly normalize data for a single chromosome
    norm.hic.table <- withoutNorm_hicTable(hic.table = bulk.hic.table)
    norm.result <- NULL
  } else {
    bulk.hic.table <- HiCcompare::create.hic.table(bulk_sparse_cond1,
      bulk_sparse_cond2,
      chr = select.chromosome, scale = FALSE
    )
    # jointly normalize data for a single chromosome
    message("\nJointly normalizing pseudo bulk matrices ")
    norm.hic.table <- HiCcompare::hic_loess(bulk.hic.table,
      Plot = Plot.normalize,
      Plot.smooth = FALSE, BP_param = BPPARAM
    )
    norm.result <- norm.hic.table
    names(norm.result)[c(7, 8, 11, 12)] <- c(
      "bulk.IF1", "bulk.IF2", "adj.bulk.IF1",
      "bulk.adj.IF2"
    )
  }


  message("\nProcessing detect differential chromotin interaction ")
  if (is.null(A.min)) {
    SD <- 2
    FC <- 3
    # numChanges proprtion with 30
    numChanges <- (1000000 / res) * 30
    A.min <- suppressMessages(suppressWarnings(best_A(
      hic.table = norm.hic.table, SD = SD, numChanges = numChanges,
      FC = FC, alpha = alpha
      )
    ))
  }
  # HiCcompare
  hic.table_result <- HiCcompare::hic_compare(norm.hic.table,
    A.min = A.min, Plot = FALSE, Plot.smooth = FALSE,
    BP_param = BPPARAM
  )

  ## add GMM layer
  hic.table.GMM_result <- GMM_layer(
    hic_table = hic.table_result,
    D.interval = D.interval, threshold = fprControl.logfc
  )
  hic.table.GMM_result <- hic.table.GMM_result[, -c(17, 18)]
  names(hic.table.GMM_result)[ncol(hic.table.GMM_result)] <-
    "Difference.cluster"
  names(hic.table.GMM_result)[c(7, 8, 11, 12)] <- c(
    "bulk.IF1", "bulk.IF2",
    "adj.bulk.IF1", "bulk.adj.IF2"
  )



  ## Save output option
  if (!is.null(save.output.path)) {
    ##### Imputed cell #####

    if (!is.null(impute1_result)) {
      ##### Group 1 #####
      df <- scHiC.table_cond1
      # folder_name <- paste0("imp_", basename(file.path.1))
      folder_name <- sprintf("imp_%s", basename(file.path.1))
      full_output_path <- file.path(save.output.path, folder_name)
      dir.create(full_output_path, recursive = TRUE)
      message(
        "\nImputed cells in condition1 saved to:",
        full_output_path, "\n"
      )

      # Transform to long format and rename cells
      sparse_df <- tidyr::pivot_longer(
        df,
        cols = starts_with("IF_"),
        names_to = "cell_id",
        values_to = "IF"
      ) %>%
        dplyr::mutate(cell_id = sub("IF_", "cell", cell_id)) %>%
        dplyr::select(cell_id, chr, region1, region2, IF)

      # Split the data into a list of data frames by cell
      cell_data_list <- split(sparse_df, sparse_df$cell_id)

      # Loop through each cell data frame and save to individual .txt files
      lapply(names(cell_data_list), function(cell) {
        cell_index <- as.numeric(gsub("[^0-9]", "", cell))
        org_name <- file_name1[cell_index]
        # Define the output file path for the current cell
        output_file_path <- file.path(
          full_output_path,
          sprintf("imp_%s", org_name)
        )
        save.data.format <- data.frame(
          chr1 = cell_data_list[[cell]]$chr,
          start1 = cell_data_list[[cell]]$region1,
          chr2 = cell_data_list[[cell]]$chr,
          start2 = cell_data_list[[cell]]$region2,
          IF = cell_data_list[[cell]]$IF
        )
        # Save the data frame to a .txt file
        write.table(save.data.format, output_file_path,
          row.names = FALSE, quote = FALSE
        )
      })
    }

    ##### Group 2 #####
    df <- scHiC.table_cond2
    # folder_name <- paste0("imp_", basename(file.path.2))
    folder_name <- sprintf("imp_%s", basename(file.path.2))
    full_output_path <- file.path(save.output.path, folder_name)
    dir.create(full_output_path, recursive = TRUE)
    message(
      "\nImputed cells in condition2 saved to:",
      full_output_path, "\n"
    )

    # Transform to long format and rename cells
    sparse_df <- tidyr::pivot_longer(
      df,
      cols = starts_with("IF_"),
      names_to = "cell_id",
      values_to = "IF"
    ) %>%
      dplyr::mutate(cell_id = sub("IF_", "cell", cell_id)) %>%
      dplyr::select(cell_id, chr, region1, region2, IF)

    # Split the data into a list of data frames by cell
    cell_data_list <- split(sparse_df, sparse_df$cell_id)

    # Loop through each cell data frame and save to individual .txt files
    lapply(names(cell_data_list), function(cell) {
      cell_index <- as.numeric(gsub("[^0-9]", "", cell))
      org_name <- file_name2[cell_index]
      # Define the output file path for the current cell
      output_file_path <- file.path(
        full_output_path,
        sprintf("imp_%s", org_name)
      )
      save.data.format <- data.frame(
        chr1 = cell_data_list[[cell]]$chr,
        start1 = cell_data_list[[cell]]$region1,
        chr2 = cell_data_list[[cell]]$chr,
        start2 = cell_data_list[[cell]]$region2,
        IF = cell_data_list[[cell]]$IF
      )
      # Save the data frame to a .txt file
      write.table(save.data.format, output_file_path,
        row.names = FALSE, quote = FALSE
      )
    })

    ##### Normalized result #####
    df <- norm.result
    message("\nNormalization result table saved to:", save.output.path, "\n")
    # Define the output file path for the normalization table
    output_file_path <- file.path(
      save.output.path,
      "Bulk_normalization_table.txt"
    )
    write.table(df, output_file_path, quote = FALSE)

    ##### Differential result #####
    df <- hic.table.GMM_result
    message(
      "\nDifferential analysis result table saved to:",
      save.output.path, "\n"
    )
    # Define the output file path for the differential analysis table
    output_file_path <- file.path(
      save.output.path,
      "Differential_analysis_table.txt"
    )
    write.table(df, output_file_path, quote = FALSE)
  }

  # Print result
  if (Plot == TRUE) {
    plot <- differential_result_plot(hic.table.GMM_result)
    print(plot)
  }

  result <- list(
    Differential_Analysis = hic.table.GMM_result,
    Intermediate = list(
      Imputation = list(condition1 = impute1_result, 
                        condition2 = impute2_result),
      PseudoBulk = list(condition1 = bulk_sparse_cond1, 
                        condition2 = bulk_sparse_cond2),
      Bulk.Normalization = norm.result
    )
  )
  
  class(result) <- "checkNumbers"
  return(result)
}


print.checkNumbers <- function(obj) {
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
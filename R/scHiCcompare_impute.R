# Declare global variables to avoid R CMD check NOTE
utils::globalVariables(c("Single_cell", "IF"))


################## Checking collinear ##################
find.collinear <- function(x, threshold = 0.999, ...) {
  nvar <- ncol(x)
  x <- data.matrix(x)
  r <- !is.na(x)
  nr <- apply(r, 2, sum, na.rm = TRUE)
  ord <- order(nr, decreasing = TRUE)
  xo <- x[, ord, drop = FALSE] ## SvB 24mar2011
  varnames <- dimnames(xo)[[2]]
  z <- stats::cor(xo, use = "pairwise.complete.obs")
  hit <- outer(seq_len(nvar), seq_len(nvar), "<") &
    (abs(z) >= threshold)
  out <- apply(hit, 2, any, na.rm = TRUE)
  return(varnames[out])
}





##################### Pooling + MICE RF imputation  #####################

################# Progressive pooling ########################
.all_progressive_pooling <- function(vector_distance) {
  # Exclude distance 0
  vector_distance <- vector_distance[-1]

  # Calculate the number of elements required for each pool
  pool_sizes <- seq_len(ceiling(
    (sqrt(8 * length(vector_distance) + 1) - 1) / 2
  ))

  # Initialize list to store each pool
  poolings_list <- list()

  # Track the current index in vector_distance
  current_index <- 1

  for (size in pool_sizes) {
    # Determine the end index for the current pool
    end_index <- min(current_index + size - 1, length(vector_distance))
    current_pool <- vector_distance[current_index:end_index]

    # Store the current pool
    poolings_list[[length(poolings_list) + 1]] <- current_pool

    # Update current_index to move to the next position in vector_distance
    current_index <- end_index + 1
  }

  # Check if the last pool needs additional elements to meet its target size
  last_pool <- utils::tail(poolings_list, 1)[[1]]
  last_pool_size <- length(last_pool)
  required_size <- utils::tail(pool_sizes, 1)

  # If the last pool is too small, add backward elements to complete it
  if (last_pool_size < required_size) {
    last_element <- min(utils::tail(poolings_list, 1)[[1]])
    additional_elements <- seq(
      last_element - 1,
      by = -1,
      length.out = required_size - last_pool_size
    )
    poolings_list[[length(poolings_list)]] <- c(last_pool, additional_elements)
  }

  # Prepend distance 0 as the first element
  poolings_list <- c(list(0), poolings_list)

  return(poolings_list)
}


################# Fibonacci pooling ########################
.fibonacciR <- function(x) {
  if (x < 2) {
    out <- x
  } else {
    out <- .fibonacciR(x - 1) + .fibonacciR(x - 2)
  }
  return(out)
}

.all_fib_pooling <- function(vector_distance) {
  ## input:
  ##    + vector_distance = vector of distance D (can include D=0)
  ## output: all poolings list

  # Length of distance vector
  l <- length(vector_distance)

  # Initialize variables
  allD_pooled <- NULL
  poolings_list <- list() # Pre-allocate the list to improve performance

  for (i in seq_len(l)) { # Safe iteration with seq_len
    # Extract number of elements to pool for each iteration (Fibonacci number)
    nElement_apool <- .fibonacciR(i)

    # Identify which distances will be pooled in iteration i
    left.D.i <- vector_distance[!(vector_distance %in% allD_pooled)]
    whichD_apool <- left.D.i[seq_len(nElement_apool)]

    # Update all pooled distances
    allD_pooled <- c(allD_pooled, whichD_apool)

    # Store pooled elements in list
    poolings_list[[i]] <- whichD_apool

    # Break if the last distance is pooled
    if (max(vector_distance) %in% whichD_apool) {
      break
    }
  }

  # Check and replace any NA values in the pools
  checkNA.pool <- vapply(
    poolings_list, function(x) any(is.na(x)),
    FUN.VALUE = logical(1)
  )

  if (any(checkNA.pool)) {
    for (i in which(checkNA.pool)) {
      pool.contain.NA <- poolings_list[[i]]
      n_NA <- sum(is.na(pool.contain.NA))

      # Replace NA values with the backward distances
      backforward_Ds <- seq(
        min(pool.contain.NA, na.rm = TRUE) - n_NA,
        min(pool.contain.NA, na.rm = TRUE) - 1
      )

      pool.contain.NA[is.na(pool.contain.NA)] <- backforward_Ds
      poolings_list[[i]] <- pool.contain.NA
    }
  }

  return(poolings_list)
}



################# pooling RF imputation process #################

############ Function to design data matrix
predictorMatrixNP_sc_D <- function(scHiC.table, distance) {
  ## input: schic.table and assigned distance
  ## output: table of single cell name, region1, region2, bin, and if

  ## Add D - distance column into scHiC table
  res <- min(abs(diff(unique(scHiC.table$region1))))
  D <- abs(scHiC.table$region1 - scHiC.table$region2) / res
  scHiC.table.up <- cbind(D, scHiC.table)

  ## Extract IF vector at corresponding distance
  entire.if <- scHiC.table[, -c(1, 2, 3, 4)]
  row.index.D <- which(D == distance)
  IF <- as.vector(unlist(c(entire.if[row.index.D, ])))
  IF[IF == 0] <- NA

  ## Extract name of level 1: single cells for each bin
  matrix_position <- ncol(scHiC.table.up) - 5
  single_cell <- rep(seq_len(matrix_position), each = length(row.index.D))

  ## Extract name of level 2: bins
  region1 <- region2 <- bin <- NULL
  for (i in seq_along(row.index.D)) {
    rg1 <- scHiC.table$region1[row.index.D[i]]
    rg2 <- scHiC.table$region2[row.index.D[i]]
    b <- sprintf("(%s,%s)", rg1, rg2)
    # get vector name of bin
    region1 <- c(region1, rg1)
    region2 <- c(region2, rg2)
    bin <- c(bin, b)
  }

  ## Combining to data.frame
  report <- data.frame(
    Single_cell = single_cell, region1, region2,
    bin, IF = IF
  )

  return(report)
}



######## Function to set up RF imputation
mice.rf_impute <- function(data_input, n.imputation = 5, maxit = 1,
                           outlier.rm = TRUE, seed = 123) {
  # Temporarily remove outlier
  if (outlier.rm == TRUE) {
    outlier.threshold <- stats::quantile(data_input$IF, 0.85,
      na.rm = TRUE
    ) +
      3 * stats::IQR(data_input$IF, na.rm = TRUE)
    # Ensure a minimum threshold value
    outlier.threshold <- max(2, outlier.threshold)
    outlier <- unique(stats::na.omit(data_input$IF)[
      stats::na.omit(data_input$IF) > outlier.threshold
    ])
    outlier.pos <- which(data_input$IF %in% stats::na.omit(outlier))
    outlier.value <- data_input$IF[outlier.pos]
    # Replace outlier with NA for now
    data_input$IF[outlier.pos] <- NA
  }

  ## Classify variable class
  data_input$Single_cell <- as.numeric(data_input$Single_cell)

  ## Classify method and predictor matrix
  # Get initial default imputation setting
  ini <- mice::mice(data_input, maxit = 0)
  # Set up method
  meth <- ini$meth
  meth["IF"] <- "rf"
  # Set up predictor matrix
  pred <- ini$pred

  ## Imputation
  # Perform multiple imputations with specified iterations
  imp <- mice::mice(data_input,
    method = meth, predictorMatrix = pred,
    print = FALSE, maxit = maxit, seed = seed,
    m = n.imputation
  )
  # Retrieve completed data in long format
  imp_data <- mice::complete(imp, action = "long", include = FALSE)
  # Aggregate mean of all imputed complete data
  agg_new_if2 <- round(stats::aggregate(imp_data$IF,
    by = list(imp_data$.id),
    FUN = mean
  )$x)

  # Add back outliers if they were removed
  if (outlier.rm == TRUE) {
    agg_new_if2[outlier.pos] <- outlier.value
  }

  return(agg_new_if2)
}





#### Impute process for each pool band
pools_impute <- function(scHiC.table, n.imputation = 5, outlier.rm = TRUE,
                         pool.style = "progressive", maxit = 1) {
  # Input: scHiC.table, n.imputation, option for outlier removal,
  #        option for impute at main closer distance
  # Output: schic table (all single cells in columns) with imputed IF

  ## Distance of scHiC table
  res <- min(abs(diff(unique(scHiC.table$region1))))
  D <- unique(abs(scHiC.table$region1 - scHiC.table$region2) / res)

  ## List of all Distance Pools
  if (pool.style == "progressive") {
    Dpool.list.full <- .all_progressive_pooling(
      vector_distance = c(min(D):max(D))
    )
  } else if (pool.style == "Fibonancci") {
    Dpool.list.full <- .all_fib_pooling(
      vector_distance = c(min(D):max(D))
    )
  }

  ## Only pool list from the input in correct original order
  Dpool.list.pos <- which(vapply(Dpool.list.full, function(x) {
    all(x %in% D)
  }, logical(1)))
  Dpool.list <- Dpool.list.full[Dpool.list.pos]
  length.Dpool.list <- length(Dpool.list)

  ## Impute all pools by RF
  process_pool <- function(i) {
    message(sprintf("pooled band %s, ", Dpool.list.pos[i]))

    ## Pooling data
    data.pool <- do.call(rbind, lapply(Dpool.list[[i]], function(x) {
      predictorMatrixNP_sc_D(scHiC.table, distance = x)
    }))

    ## Impute by pools
    data <- unique(data.pool)
    data <- data[order(data$region1), ]
    data_input <- data[, -c(2, 3, 4)]
    rm.na.data.input <- stats::na.omit(data_input)

    if (all(rm.na.data.input$IF == rm.na.data.input$IF[1]) ||
      anyDuplicated(rm.na.data.input$Single_cell) == 0 ||
      length(find.collinear(data_input)) > 0) {
      agg_new_if2 <- data_input$IF
      agg_new_if2[is.na(agg_new_if2)] <- 0
      agg_new_if2[agg_new_if2 == 0] <-
        round(mean(agg_new_if2, na.rm = TRUE))
    } else if (nrow(rm.na.data.input) == 0) {
      agg_new_if2 <- data_input$IF
      agg_new_if2[is.na(agg_new_if2)] <- rep(0, nrow(data_input))
    } else {
      agg_new_if2 <- mice.rf_impute(
        data_input = data_input, n.imputation = n.imputation,
        outlier.rm = outlier.rm, maxit = maxit
      )
    }

    ## Create new scHiC table
    new_data <- data.frame(
      Single_cell = data$Single_cell,
      region1 = data$region1,
      region2 = data$region2,
      IF = agg_new_if2
    )
    new_data <- unique(new_data)

    ## Transform new_data into wide format
    update_new_table <- tidyr::pivot_wider(
      new_data,
      names_from = c(Single_cell), values_from = IF
    )
    update_new_table <- data.frame(
      region1 = update_new_table$region1,
      region2 = update_new_table$region2,
      cell = rep(unique(scHiC.table$cell), nrow(update_new_table)),
      chr = rep(unique(scHiC.table$chr), nrow(update_new_table)),
      update_new_table[, -c(1, 2)]
    )
    colnames(update_new_table)[5:ncol(update_new_table)] <-
      sprintf("IF_%d", seq_len(ncol(update_new_table) - 4))

    return(update_new_table)
  }

  ## Run lapply for faster processing
  new_table_list <- lapply(seq_len(length.Dpool.list), process_pool)

  ## Handle backward D (NA replacement)
  lastD.poolPosition <- which(Dpool.list[[length.Dpool.list]] == max(D))
  lastPool <- Dpool.list[[length.Dpool.list]]

  if (lastD.poolPosition < length(lastPool)) {
    backward_D <- lastPool[(lastD.poolPosition + 1):length(lastPool)]
    new_table_lastPool <- new_table_list[[length.Dpool.list]]
    backward_D_position <- which(
      abs((new_table_lastPool$region1 - new_table_lastPool$region2) / res) %in%
        backward_D
    )
    new_table_list[[length.Dpool.list]] <-
      new_table_lastPool[-backward_D_position, ]
  }

  ## Merge all elements of new_table_list into one data frame
  new_table <- do.call(rbind, new_table_list)

  return(new_table)
}







#####################  RF without Pooling imputation  #####################
## Function to design data matrix
predictorMatrix_sc_D <- function(scHiC.table, distance) {
  ## Add D - distance column into scHiC table
  res <- min(abs(diff(unique(scHiC.table$region1))))
  D <- abs(scHiC.table$region1 - scHiC.table$region2) / res
  scHiC.table.up <- cbind(D, scHiC.table)

  ## Extract IF vector at corresponding distance
  entire.if <- scHiC.table[, -c(1, 2, 3, 4)]
  row.index.D <- which(D == distance)
  IF <- as.vector(unlist(c(entire.if[row.index.D, ])))
  IF[IF == 0] <- NA

  ## Extract name of level 1: single cells for each bin
  matrix_position <- ncol(scHiC.table.up) - 5
  single_cell <- rep(seq_len(matrix_position), each = length(row.index.D))

  ## Extract name of level 2: bins
  region1 <- region2 <- NULL
  for (i in seq_along(row.index.D)) {
    rg1 <- scHiC.table$region1[row.index.D[i]]
    rg2 <- scHiC.table$region2[row.index.D[i]]
    # get vector name of bin
    region1 <- c(region1, rg1)
    region2 <- c(region2, rg2)
  }

  ## Combining to data.frame
  report <- data.frame(Single_cell = single_cell, region1, region2, IF = IF)

  return(report)
}




#### Implement full range RF
RF_impute.outrm.schic <- function(scHiC.table, n_imputation = 5,
                                  outlier.rm = TRUE, maxit = 1) {
  res <- min(abs(diff(unique(scHiC.table$region1))))
  D <- unique(abs(scHiC.table$region1 - scHiC.table$region2) / res)
  n_distance <- length(D)
  new_table_list <- vector("list", n_distance)

  for (i in seq_len(n_distance)) {
    D_pos <- seq_len(n_distance - 1) - 1
    message(sprintf("band %s, ", D_pos[i]))
    data <- predictorMatrix_sc_D(scHiC.table, distance = D_pos[i])
    data <- data[order(data$region1), ]
    data_input <- data[, -c(2, 3)]

    if (outlier.rm) {
      outlier <- unique(grDevices::boxplot.stats(data_input$IF)$out)
      outlier.pos <- which(data_input$IF %in% outlier)
      outlier.value <- data_input$IF[outlier.pos]
      data_input$IF[outlier.pos] <- NA
    }

    data_input$Single_cell <- as.character(data_input$Single_cell)
    rm.na.data.input <- stats::na.omit(data_input)

    if (all(rm.na.data.input$IF == rm.na.data.input$IF[1]) ||
      anyDuplicated(rm.na.data.input$Single_cell) == 0 ||
      all(table(stats::na.omit(data)$bin) == 1) ||
      length(find.collinear(data_input)) > 0) {
      agg_new_if2 <- data_input$IF
      agg_new_if2[is.na(agg_new_if2)] <- round(
        mean(rm.na.data.input$IF, na.rm = TRUE)
      )
    } else if (nrow(rm.na.data.input) == 0) {
      agg_new_if2 <- rep(1, nrow(data_input))
    } else {
      ini <- mice::mice(data_input, maxit = 0)
      ini$meth["IF"] <- "rf"
      imp <- mice::mice(data_input,
        method = ini$meth,
        predictorMatrix = ini$pred,
        print = FALSE, maxit = maxit,
        m = n_imputation
      )
      imp_data <- mice::complete(imp,
        action = "long",
        include = FALSE
      )
      agg_new_if2 <- round(stats::aggregate(
        imp_data[, 4],
        by = list(imp_data$.id), FUN = mean
      )$x)
    }

    if (outlier.rm) agg_new_if2[outlier.pos] <- outlier.value

    new_data <- data.frame(
      Single_cell = data$Single_cell,
      region1 = data$region1,
      region2 = data$region2, IF = agg_new_if2
    )
    update_new_table <- tidyr::pivot_wider(
      new_data,
      names_from = Single_cell, values_from = IF
    )
    update_new_table <- data.frame(
      region1 = update_new_table$region1,
      region2 = update_new_table$region2,
      cell = rep(
        unique(scHiC.table$cell),
        nrow(update_new_table)
      ),
      chr = rep(
        unique(scHiC.table$chr),
        nrow(update_new_table)
      ),
      update_new_table[, -c(1, 2)]
    )
    colnames(update_new_table)[5:ncol(update_new_table)] <- sprintf(
      "IF_%d", seq_len(ncol(update_new_table) - 4)
    )

    new_table_list[[i]] <- update_new_table
  }

  return(dplyr::bind_rows(new_table_list))
}






###### Entire of Random Forest Process
RF_process <- function(scHiC.table, n_imputation = 5, outlier.rm = TRUE,
                       maxit = 1, main_Distances,
                       missPerc.threshold = 95) {
  ############################# Identify important Pools and Parameters

  ## Calculate resolution and unique distances
  res <- min(abs(diff(unique(scHiC.table$region1))))
  D <- unique(abs(scHiC.table$region1 - scHiC.table$region2) / res)
  n_distance <- length(D)

  ## Define main distance range and NA in each distance
  maxD_res <- max(main_Distances) / res
  main_D_range <- seq_len(maxD_res) - 1

  na_perc_all <- vapply(D, function(d) {
    data_input <- predictorMatrixNP_sc_D(scHiC.table, distance = d)$IF
    mean(is.na(data_input)) * 100
  }, numeric(1)) # Specify the expected output type

  ### Impute if all distances have < 95% NA
  if (all(na_perc_all <= missPerc.threshold)) {
    new_table <- RF_impute.outrm.schic(
      scHiC.table = scHiC.table, n_imputation = n_imputation,
      outlier.rm = outlier.rm, maxit = maxit
    )
  } else {
    ## Identify distances with NA > 95%
    D_outside.mainD <- D[!(D %in% main_D_range)]
    D_above_threshold <- D[na_perc_all > missPerc.threshold]
    D_outside_main_above_threshold <-
      D_outside.mainD[D_outside.mainD %in% D_above_threshold]

    ## Process each distance for imputation
    new_table_list <- vector("list", n_distance)
    for (i in seq_len(n_distance)) {
      D_pos <- seq_len(n_distance - 1) - 1
      message(sprintf("band %s, ", D_pos[i]))
      data <- predictorMatrix_sc_D(scHiC.table, distance = D_pos[i])
      data <- data[order(data$region1), ]
      data_input <- data[, -c(2, 3)]

      # Temporarily remove outliers if needed
      if (outlier.rm) {
        outlier <- unique(grDevices::boxplot.stats(data_input$IF)$out)
        outlier.pos <- which(data_input$IF %in% outlier)
        outlier.value <- data_input$IF[outlier.pos]
        data_input$IF[outlier.pos] <- NA
      }

      ## Imputation conditions
      rm.na.data.input <- stats::na.omit(data_input)
      if (D_pos[i] %in% D_outside_main_above_threshold ||
        all(rm.na.data.input$IF == rm.na.data.input$IF[1]) ||
        anyDuplicated(rm.na.data.input$Single_cell) == 0 ||
        all(table(stats::na.omit(data)$bin) == 1) ||
        length(find.collinear(data_input)) > 0) {
        # Use mean or default imputation for high NA distances
        agg_new_if2 <- if (nrow(rm.na.data.input) == 0) {
          rep(1, nrow(data_input))
        } else {
          data_input$IF
          data_input$IF[is.na(data_input$IF)] <-
            round(mean(rm.na.data.input$IF, na.rm = TRUE))
        }
      } else {
        # Multiple imputation using `mice`
        ini <- mice::mice(data_input, maxit = 0)
        ini$meth["IF"] <- "rf"
        imp <- mice::mice(
          data_input,
          method = ini$meth,
          predictorMatrix = ini$pred, print = FALSE,
          maxit = maxit, m = n_imputation
        )
        imp_data <- mice::complete(imp, action = "long", include = FALSE)
        agg_new_if2 <- round(
          stats::aggregate(imp_data[, 4],
            by = list(imp_data$.id), FUN = mean
          )$x
        )
      }

      # Restore outliers if removed
      if (outlier.rm) {
        agg_new_if2[outlier.pos] <- outlier.value
      }

      # Construct new data frame for this distance
      new_data <- data.frame(
        Single_cell = data$Single_cell, region1 = data$region1,
        region2 = data$region2, IF = agg_new_if2
      )

      # Transform to wide format
      update_new_table <- tidyr::pivot_wider(
        new_data,
        names_from = Single_cell, values_from = IF
      )
      update_new_table <- data.frame(
        region1 = update_new_table$region1,
        region2 = update_new_table$region2,
        cell = rep(unique(scHiC.table$cell), nrow(update_new_table)),
        chr = rep(unique(scHiC.table$chr), nrow(update_new_table)),
        update_new_table[, -c(1, 2)]
      )
      colnames(update_new_table)[5:ncol(update_new_table)] <-
        sprintf("IF_%d", seq_len(ncol(update_new_table) - 4))

      new_table_list[[i]] <- update_new_table
    }

    # Combine all distance tables
    new_table <- dplyr::bind_rows(new_table_list)
  }

  return(new_table)
}










#' Random Forest Imputation with Pooling options for scHi-C Data
#'
#' This function performs imputation of single-cell Hi-C (scHi-C) interaction
#'  frequencies (IF) using Random Forest imputation methods with different
#'  options for distance-based pooling strategies.
#'
#' @param scHiC.table A data frame containing interaction frequencies across
#'  single cells, created by the `scHiC_table` function. The first four columns
#'  should represent 'cell', 'chr', 'region1', and 'region2', followed by
#'  columns representing interaction frequencies ('IF') for individual cells.
#' @param n.imputation An integer specifying the number of imputations to be
#'  performed. The default is 5.
#' @param maxit An integer specifying the number of iterations for the internal
#'  refinement process within a single imputation cycle. The default is 1.
#' @param outlier.rm A logical value indicate whether to remove outliers during
#'  the imputation process. The default is TRUE.
#' @param main.Distances A vector of integers or 'full' representing the scHiC
#'  data in the main distance range to focus the imputation on, in bp units
#'  (e.g., 1:1,000,000). Genomic distances (in bp) are the number of base pairs
#'  between two regions in the genome (e.g., loci or bins).
#'  The default is from 1 to 10,000,000.
#' @param pool.style A string specifying the pooling technique to use. Options
#'  are 'none', 'progressive', or 'Fibonacci'. Default is 'progressive'.
#' @param missPerc.threshold An integer specifying the missing value percentage
#'  threshold in each pool band.
#'
#' @return A table in the format of a scHiC table (same structure as the output
#'  of the `scHiC_table` function) with imputed interaction frequencies (IF)
#'  across all single cells. The output table is formatted with regions and
#'  single cells in wide format, with one column per single cell containing
#'  imputed IF values.
#'
#' @details
#' The function first identifies important pools based on the given scHi-C
#'  genomic distance effect by pooling distance data according to the chosen
#'  method. For progressive pooling, pools of distances are consecutively
#'  combined to form larger sets, while Fibonacci pooling uses a Fibonacci
#'  sequence to combine distances. If the pooling style `none` is selected, the
#'  band contains individual 1 genomic distance. During the imputation process,
#'  the function imputes all missing values (NAs) within each pool within the
#'  main distance range. For distances outside this main focus range, if any
#'  pool contains more than `missPerc.threshold` missing values, it triggers an
#'  alternative imputation method, filling in missing values based on the mean
#'  for distances.
#'
#' @references
#'
#' Doove, L.L., van Buuren, S., Dusseldorp, E. (2014), Recursive partitioning
#' for missing data imputation in the presence of interaction Effects.
#' Computational Statistics & Data Analysis, 72, 92-104.
#'
#' Shah, A.D., Bartlett, J.W., Carpenter, J., Nicholas, O., Hemingway, H.
#' (2014), Comparison of random forest and parametric imputation models for
#' imputing missing data using MICE: A CALIBER study. American Journal of
#' Epidemiology, doi:10.1093/aje/kwt312.
#'
#' Van Buuren, S. (2018). Flexible Imputation of Missing Data. Second Edition.
#' Chapman & Hall/CRC. Boca Raton, FL.
#'
#' @examples
#' # Load MG data folder example
#' MGs_example <- system.file("extdata/MGs_example", package = "scHiCcompare")
#' # Create scHicCompare table to be used in scHicCompare
#' IF_table <- scHiC_table(
#'   file.path = MGs_example, cell.type = "MG",
#'   select.chromosome = "chr20"
#' )
#' # Example usage of Pooling_RF_impute
#' library(tidyr)
#' imputed_table <- scHiCcompare_impute(IF_table,
#'   n.imputation = 5, outlier.rm = TRUE,
#'   main.Distances = 1:10000000, pool.style = "progressive"
#' )
#'
#' @import tidyr
#' @import mice
#' @import lattice
#' @import miceadds
#' @import rstatix
#' @import ranger
#' @import utils
#' @importFrom stats na.omit aggregate quantile IQR cor
#' @importFrom graphics boxplot
#' @importFrom grDevices boxplot.stats
#'
#' @export



scHiCcompare_impute <- function(scHiC.table, n.imputation = 5, maxit = 1,
                                outlier.rm = TRUE,
                                main.Distances = seq(10000000),
                                pool.style = "progressive",
                                missPerc.threshold = 95) {
  #############################################################################
  ############## Identify important Pools and Parameter #######################

  ## distance of scHiC table
  res <- min(abs(diff(unique(scHiC.table$region1))))
  D <- unique(abs(scHiC.table$region2 - scHiC.table$region1) / res)
  n_all.distance <- length(D)

  ## Identify main Distances
  if (is.infinite(max(main.Distances))) {
    main_D_range <- seq_len(max(D))
  } else {
    maxD <- max(main.Distances)
    maxD_res <- maxD / res
    main_D_range <- seq_len(maxD_res)
  }

  #############################################################################
  ############################# Imputation process ############################

  ## If Pooling styles is applied
  if (pool.style == "none") { # If Single pooling is selected
    new_table <- RF_process(
      scHiC.table = scHiC.table, n_imputation = n.imputation,
      outlier.rm = outlier.rm, main_Distances = main.Distances,
      missPerc.threshold = missPerc.threshold
    )
  } else { # If Pooling is selected

    ## List of all Distance Pool
    if (pool.style == "progressive") {
      Dpool.list <- .all_progressive_pooling(vector_distance = D)
    } else if (pool.style == "Fibonancci") {
      Dpool.list <- .all_fib_pooling(vector_distance = D)
    }
    length.Dpool.list <- length(Dpool.list)

    new_table <- NULL
    new_table_list <- rep(list(NULL), (length.Dpool.list))
    na_perc_all <- NULL
    for (i in seq_len(length.Dpool.list)) {
      ## Test na percent
      data.pool <- do.call(rbind, lapply(Dpool.list[[i]], function(x) {
        data <- predictorMatrixNP_sc_D(scHiC.table, distance = x)
        return(data)
      }))
      data_input <- data.pool$IF
      na_perc <- (sum(is.na(data_input)) / length(data_input)) * 100
      na_perc_all <- c(na_perc_all, na_perc)
      ## Test if any element in pools D within main Distances range
      which.in.mainD <- Dpool.list[[i]] %in% main_D_range
      mainD.range <- any(which.in.mainD)
    }

    ### If all Distance have NA < `missPerc.threshold`
    if (length(which(na_perc_all > missPerc.threshold)) == 0) {
      new_table <- pools_impute(
        scHiC.table = scHiC.table, n.imputation = n.imputation,
        outlier.rm = outlier.rm, pool.style = pool.style
      )
    } else { ### if any pool have NA > `missPerc.threshold`
      ## Check which pool have above `missPerc.threshold` -> impute mean only
      which_pool_aboveNA <- which(na_perc_all > missPerc.threshold)
      pool_aboveNA_Dlist <- lapply(which_pool_aboveNA, function(x) {
        return(Dpool.list[[x]])
      }) # list of detail D in which_pool_aboveNA
      pool_aboveNA_D <- unlist(pool_aboveNA_Dlist)

      ## check which pool above NA level and in main distance
      pool_aboveNA_mainD <- which_pool_aboveNA[
        unlist(lapply(which_pool_aboveNA, function(x) {
          any(main_D_range %in% Dpool.list[[x]])
        }))
      ]
      pool_aboveNA_NOTmainD <- which_pool_aboveNA[!which_pool_aboveNA %in%
        pool_aboveNA_mainD]
      pool_mainD_imp <- c(seq_len(length.Dpool.list))[
        !c(seq_len(length.Dpool.list)) %in% pool_aboveNA_NOTmainD
      ]
      D_aboveNA_NOTmainD <- unlist(lapply(pool_aboveNA_NOTmainD, function(x) {
        return(Dpool.list[[x]])
      })) # which D above NA level and NOT mainD
      D_aboveNA_mainD <- unlist(lapply(pool_aboveNA_mainD, function(x) {
        return(Dpool.list[[x]])
      }))
      D_mainP_imp <- unlist(lapply(pool_mainD_imp, function(x) {
        return(Dpool.list[[x]])
      })) ## Distance in main pool will be imputed

      ## Extract any pools that will be imputed
      temp.table <- scHiC.table
      temp.table <- temp.table[((temp.table$region2 - temp.table$region1) / res)
      %in% D_mainP_imp, ]
      ## Impute
      temp.table_imp <- pools_impute(
        scHiC.table = temp.table, n.imputation = n.imputation,
        outlier.rm = outlier.rm, pool.style = pool.style
      )

      ## Now, impute D_set1 to its mean of that distance
      new_table_list <- list()
      for (j in seq_along(pool_aboveNA_NOTmainD)) {
        message(sprintf("pooled band %s, ", pool_aboveNA_NOTmainD[j]))
        data <- do.call(rbind, lapply(Dpool.list[[
          pool_aboveNA_NOTmainD[j]
        ]], function(x) {
          data <- predictorMatrixNP_sc_D(scHiC.table, distance = x)
          return(data)
        }))
        data <- data[order(data$region1), ] # sort data
        data_input <- data[, -c(2, 3, 4)]
        extreme.value_pos <- data_input$IF %in% graphics::boxplot(data_input$IF,
          plot = FALSE
        )$out # remove any extreme value
        extreme.value <- data_input$IF[data_input$IF %in%
          graphics::boxplot(data_input$IF,
            plot = FALSE
          )$out]
        data_input$IF[extreme.value_pos] <- NA
        data_input$IF[is.na(data_input$IF)] <- round(mean(data_input$IF,
          na.rm = TRUE
        ))
        agg_new_if2 <- data_input$IF
        agg_new_if2[extreme.value_pos] <- extreme.value

        # Create new scHicTable
        new_data <- data.frame(
          Single_cell = data$Single_cell, region1 = data$region1,
          region2 = data$region2, IF = agg_new_if2
        )
        ## Transform new_data into wide format and Re-format to scHicTable
        update_new_table <- tidyr::pivot_wider(
          new_data,
          names_from = Single_cell,
          values_from = IF
        )

        update_new_table <- data.frame(
          region1 = update_new_table$region1,
          region2 = update_new_table$region2,
          cell = rep(unique(scHiC.table$cell), nrow(update_new_table)),
          chr = rep(unique(scHiC.table$chr), nrow(update_new_table)),
          update_new_table[, -c(1, 2)]
        )

        # Rename columns to have "IF_" prefix for each imputed column
        colnames(update_new_table)[5:ncol(update_new_table)] <- sprintf(
          "IF_%d",
          seq_len(ncol(update_new_table) - 4)
        )

        # Append to list
        new_table_list[[j]] <- update_new_table
      }

      if (max(D) %in% D_aboveNA_NOTmainD) {
        # Check if any pool contain backward D (NA replacement).
        # If so, then need to remove it/them
        ## Where does the last D locate in last pool
        lastD.poolPosition <- which(Dpool.list[[length.Dpool.list]] == max(D))
        lastPool <- Dpool.list[[length.Dpool.list]]
        # check if last D stay at the end of last pool
        if (lastD.poolPosition < length(lastPool)) {
          # identify backward D
          backward_D <- lastPool[(lastD.poolPosition + 1):length(lastPool)]
          # identify backward D locate in new_table_list of last pool
          new_table_lastPool <- new_table_list[[length((new_table_list))]]
          backward_D_position <- which(abs((new_table_lastPool$region1 -
            new_table_lastPool$region2) / res) %in% backward_D)
          new_table_list[[length(pool_aboveNA_NOTmainD)]] <-
            new_table_lastPool[-backward_D_position, ]
        }
      }

      new_table_D1 <- do.call(rbind, new_table_list)

      ## Merge for final table
      new_table <- rbind(temp.table_imp, new_table_D1)
    }
  }

  return(new_table)
}

################## Checking collinear ################## 
find.collinear <- function(x, threshold = 0.999, ...) {
  nvar <- ncol(x)
  x <- data.matrix(x)
  r <- !is.na(x)
  nr <- apply(r, 2, sum, na.rm = TRUE)
  ord <- order(nr, decreasing = TRUE)
  xo <- x[, ord, drop = FALSE] ## SvB 24mar2011
  varnames <- dimnames(xo)[[2]]
  z <- suppressWarnings(cor(xo, use = "pairwise.complete.obs"))
  hit <- outer(seq_len(nvar), seq_len(nvar), "<") & (abs(z) >= threshold)
  out <- apply(hit, 2, any, na.rm = TRUE)
  return(varnames[out])
}





##################### Pooling + MICE RF imputation  ##################### 

################# progressive pooling ########################
.all_progressive_pooling <- function(vector_distance){
  ## input: 
  ##    + vector_distance = vector of distance D (can include D=0)
  ## output: all poolings list
  
  #exclude distance 0
  vector_distance <- vector_distance[-1] 
  #length of distance vector
  l <- length(vector_distance)  
  # total of all possible pooling groups
  all.pooling.group <- ceiling((sqrt(8 * (l - 1) + 1) - 1) / 2 )
  # extract number members of all pooling groups
  all.group.n <- 1:all.pooling.group
  # Initialize pooled_distance as an empty vector
  pooled_distance <- c()
  # Use lapply to create the poolings_list
  poolings_list <- lapply(all.group.n, function(i) {
    D_pooling_n <- vector_distance[i]
    # Exclude elements in pooled_distance
    poolings_elements <- vector_distance[!(vector_distance %in% pooled_distance)][1:D_pooling_n]
    # Update pooled_distance with the new elements
    pooled_distance <<- c(poolings_elements, pooled_distance)
    return(poolings_elements)
  })
  
  # if any pool D contain NA, replace these NA with backward D
  checkNA.pool <- sapply(poolings_list, function(x) any(is.na(x)))
  if ( T %in% checkNA.pool ){
    which.pool.NA <- which(checkNA.pool == T)
    pool.contain.NA <- poolings_list[[ which.pool.NA ]]
    n_NA <- length( which(is.na(pool.contain.NA)) )
    backforward_Ds <- seq(min(pool.contain.NA, na.rm = T) - n_NA, min(pool.contain.NA, na.rm = T) - 1)
    pool.contain.NA[is.na(pool.contain.NA)] <- backforward_Ds
    poolings_list[[ which.pool.NA ]] <- pool.contain.NA
  }
  
  # Add back Distance 0
  poolings_list <- c(0, poolings_list)
  return(poolings_list)
}


################# Fibonacci pooling ########################
.fibonacciR <- function(x){
  if (x < 2){
    out = x
  } else {
    out = fibonacciR(x - 1) + fibonacciR(x - 2)
  }
  return(out)
}

.all_fib_pooling <- function(vector_distance){
  ## input: 
  ##    + vector_distance = vector of distance D (can include D=0)
  ## output: all poolings list
  
  #length of distance vector
  l <- length(vector_distance)  
  allD_pooled <- NULL; poolings_list <- list()
  for (i in 1:l){ #set a list of trial iteration -> set large
    # extract number of D will be pool for each iteration
    nElement_apool <- fibonacciR(i)
    # identify which D will be pool in iteration i
    left.D.i <- vector_distance[!vector_distance %in% allD_pooled]
    whichD_apool <- left.D.i[1 : nElement_apool]
    # all pooled D have been pooled
    allD_pooled <- c(allD_pooled, whichD_apool)
    # store pool element for each group (i) in list
    poolings_list[[i]] <- whichD_apool
    # break if whichD_apool contain last D
    if (max(vector_distance) %in% whichD_apool){
      break
    }
  }
  
  # if any pool D contain NA, replace these NA with backward D
  checkNA.pool <- sapply(poolings_list, function(x) any(is.na(x)))
  if ( T %in% checkNA.pool ){
    which.pool.NA <- which(checkNA.pool == T)
    pool.contain.NA <- poolings_list[[ which.pool.NA ]]
    n_NA <- length( which(is.na(pool.contain.NA)) )
    backforward_Ds <- seq(min(pool.contain.NA, na.rm = T) - n_NA, min(pool.contain.NA, na.rm = T) - 1)
    pool.contain.NA[is.na(pool.contain.NA)] <- backforward_Ds
    poolings_list[[ which.pool.NA ]] <- pool.contain.NA
  }
  
  #poolings_list[[1]] <- 1
  return(poolings_list)
  
}



################# progressive RF imputation process ################# 

############ Function to design data matrix 
predictorMatrixNP_sc_D <- function(scHiC.table, distance){
  ## input: schic.table and assigned distance
  ## output: table of single cell name, region1, region2, bin, and if
  
  ## Add D - distance column into scHiC table
  res  = min( abs( diff(unique(scHiC.table$region1)) ) )
  D = abs(scHiC.table$region1 - scHiC.table$region2)/res
  scHiC.table.up <- cbind(D,scHiC.table)
  
  ## Extract IF vector at corresponding distance 
  entire.if = scHiC.table[,-c(1:4)]
  row.index.D = which(D == distance)
  IF =  as.vector( unlist(c(entire.if[row.index.D,])) )
  IF[IF == 0] <- NA
  
  ## Extract name of level 1: single cells for each bin
  matrix_position = ncol(scHiC.table.up) - 5
  single_cell <- rep(c(1:matrix_position), each= length(row.index.D))
  
  ## Extract name of level 2: bins
  region1 = region2 = bin = NULL
  for(i in 1:length(row.index.D)){
    rg1 <- scHiC.table$region1[row.index.D[i]]
    rg2 <- scHiC.table$region2[row.index.D[i]]
    b <- paste0('(', rg1, ',', rg2, ')')
    # get vector name of bin
    region1 = c(region1, rg1)
    region2 = c(region2, rg2)
    bin = c(bin,b)
  }
  
  ## Combining to data.frame
  report <- data.frame( Single_cell = single_cell, region1, region2, bin, IF = IF)
  
  return(report)
}



######## Function to set up RF imputation 
mice.rf_impute <- function(data_input, n.imputation = 5, maxit = 1, outlier.rm = TRUE, seed = 123){
  require(rstatix)
  require(mice)
  require(lattice)
  library(miceadds)
  # Temporally remove outlier
  if(outlier.rm == TRUE){
    outlier.threshold <- quantile(data_input$IF, 0.85, na.rm = T) + 3 * IQR(data_input$IF, na.rm = T) ## extreme value [??? The outlier.threshold with a 0.85 quantile seems to be a pretty important assumption. Perhaps describe this and pass this as a parameter further upstream?]
    outlier.threshold <- max(2,outlier.threshold)
    outlier <-  unique(na.omit(data_input$IF)[na.omit(data_input$IF) > outlier.threshold] )
    outlier.pos <- which(data_input$IF %in% na.omit(outlier))
    outlier.value <- data_input$IF[outlier.pos]
    data_input$IF[outlier.pos] <- NA #replace outlier with NA for now
  }
  
  ## Classify variable class
  data_input$Single_cell <- as.numeric(data_input$Single_cell)
  #data_input$Cell_type= as.numeric(as.factor(data_input$Cell_type))
  
  ## Classify method and predictor matrix
  # get initial default imputation setting
  set.seed(123)
  ini <- suppressWarnings(mice(data_input, maxit = 0))
  #set up method
  meth= ini$meth
  meth['IF'] <- "rf"
  #set up predictor matrix
  pred = ini$pred
  #pred[,'Single_cell'] <- 1 #random intercept
  
  ## Imputation
  #simulate for n time of multiple imputations, with 5 iteration
  imp <- suppressWarnings(mice(data_input, method = meth, predictorMatrix = pred, print = FALSE,maxit = maxit, seed = seed, m=n.imputation) )
  #returns the long format of all multiple imputation
  imp_data = complete(imp, action = 'long', include = F)
  #if the vector has all if>1, aggregate mean of all imputed complete data
  agg_new_if2 = round(aggregate(imp_data[,4] , by = list(imp_data$.id),FUN= mean))
  agg_new_if2 = agg_new_if2$x
  
  # Add back outlier, if they were removed in previous steps
  if(outlier.rm == TRUE){
    agg_new_if2[outlier.pos] <- outlier.value
  }
  return(agg_new_if2)
}



#### Impute process for each pool band
pools_impute <-  function(scHiC.table, n.imputation = 5, outlier.rm = TRUE,
                          pool.style = 'progressive', maxit = 1){
  # input: scHiC.table, n.imputation, and, option for outlier remover, option for impute at main closer distance
  # output: schic table (all single cell in columns) with imputed if
  
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_all.distance = length(D)
  
  ## List of all Distance Pool
  if(pool.style == 'progressive'){
    Dpool.list.full <- .all_progressive_pooling(vector_distance = c(min(D):max(D)))
  } else if ( pool.style == 'fibonancci'){
    Dpool.list.full <- all_fib_pooling(vector_distance = c(min(D):max(D)) )
  }
  
  ## only pool list from the input in correct original order
  Dpool.list.pos <- which(sapply(Dpool.list.full, function(x) all(x %in% D)))
  Dpool.list <- Dpool.list.full[Dpool.list.pos]
  length.Dpool.list <- length(Dpool.list)
  
  new_table <- NULL
  new_table_list <- rep(list(NULL), length.Dpool.list)
  
  
  ## Impute all pools by RF
  process_pool <- function(i) {
    cat(paste0(' pooled band ', Dpool.list.pos[i] ,', ') )
    
    require(tidyr)
    ## Pooling data
    data.pool <- do.call(rbind, lapply(Dpool.list[[i]], function(x) { 
      data <- predictorMatrixNP_sc_D(scHiC.table, distance = x)
      return(data)
    }))
    
    ## Impute by pools
    data <- unique(data.pool)
    data <- data[order(data$region1), ] #sort data 
    data_input <- data[,-c(2,3,4)]
    rm.na.data.input <- na.omit(data_input)
    
    if (all(rm.na.data.input$IF == rm.na.data.input$IF[1]) || 
        anyDuplicated(rm.na.data.input$Single_cell) == 0  || 
        length(find.collinear(data_input)) > 0) {
      agg_new_if2 <- data_input$IF
      agg_new_if2[is.na(agg_new_if2)] <- 0
      agg_new_if2[agg_new_if2 == 0] <- round(mean(agg_new_if2, na.rm = TRUE))
    } else if (nrow(rm.na.data.input) == 0) {
      agg_new_if2 <- data_input$IF
      agg_new_if2[is.na(agg_new_if2)] <- rep(0, nrow(data_input))
    } else { 
      agg_new_if2 <- mice.rf_impute(data_input = data_input, n.imputation = n.imputation, outlier.rm = outlier.rm, maxit = maxit) 
    }
    
    # Create new scHicTable
    new_data <- data.frame(Single_cell = data$Single_cell , region1 = data$region1,
                           region2 = data$region2, IF = agg_new_if2) 
    new_data <- unique(new_data)
    
    ## Transform new_data into wide format and Re-format to scHicTable
    update_new_table <- pivot_wider(new_data, 
                                    names_from = c(Single_cell), 
                                    values_from = IF)
    update_new_table <- data.frame(region1 = update_new_table$region1,
                                   region2 = update_new_table$region2,
                                   cell = rep(unique(scHiC.table$cell), nrow(update_new_table)),
                                   chr  = rep(unique(scHiC.table$chr), nrow(update_new_table)), 
                                   update_new_table[,-c(1:2)])
    colnames(update_new_table)[5:ncol(update_new_table)] <- paste0("IF_", 1:(ncol(update_new_table)-4))
    
    return(update_new_table)
  }
  
  ## run lapply for faster result
  new_table_list <- lapply(1:length.Dpool.list, process_pool)
  
  
  # Check if any pool contains backward D (NA replacement). If so, then need to remove it/them
  ## Where does the last D locate in last pool
  lastD.poolPosition <- which(Dpool.list[[length.Dpool.list]] == max(D))
  lastPool <- Dpool.list[[length.Dpool.list]]
  # check if last D stays at the end of the last pool
  if (lastD.poolPosition < length(lastPool)) {
    # identify backward D
    backward_D <- lastPool[(lastD.poolPosition + 1):length(lastPool)]
    # identify backward D location in new_table_list of last pool
    new_table_lastPool <- new_table_list[[length.Dpool.list]]
    backward_D_position <- which(abs((new_table_lastPool$region1 - new_table_lastPool$region2) / res) %in% backward_D) 
    new_table_list[[length.Dpool.list]] <- new_table_lastPool[-backward_D_position,]
    
  }
  ## Merge all elements of new_table_list into 1 data frame
  new_table <- do.call(rbind, new_table_list)
  
  return(new_table)
}







#####################  RF without Pooling imputation  ##################### 
## Function to design data matrix
predictorMatrix_sc_D <- function(scHiC.table, distance){
  ## Add D - distance column into scHiC table
  res  = min( abs( diff(unique(scHiC.table$region1)) ) )
  D = abs(scHiC.table$region1 - scHiC.table$region2)/res
  scHiC.table.up <- cbind(D,scHiC.table)
  
  ## Extract IF vector at corresponding distance 
  entire.if = scHiC.table[,-c(1:4)]
  row.index.D = which(D == distance)
  IF =  as.vector( unlist(c(entire.if[row.index.D,])) )
  IF[IF == 0] <- NA
  
  ## Extract name of level 1: single cells for each bin
  matrix_position = ncol(scHiC.table.up) - 5
  single_cell <- rep(c(1:matrix_position), each= length(row.index.D))
  
  ## Extract name of level 2: bins
  region1 = region2 = NULL
  for(i in 1:length(row.index.D)){
    rg1 <- scHiC.table$region1[row.index.D[i]]
    rg2 <- scHiC.table$region2[row.index.D[i]]
    # get vector name of bin
    region1 = c(region1, rg1)
    region2 = c(region2, rg2)
    
  }
  
  ## Combining to data.frame
  report <- data.frame( Single_cell = single_cell, region1, region2, IF = IF)
  
  return(report)
}



#### Implement full range RF
RF_impute.outrm.schic <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE, maxit = 1){
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_distance = length(D)
  
  l = n_distance-1
  new_table = NULL
  # only impute on distance data have single cell>1
  for( i in 1:n_distance){
    D_pos = 0:l
    cat(paste0('band ', D_pos[i] ,', ') )
    data = predictorMatrix_sc_D(scHiC.table, distance = D_pos[i])
    data =data[order(data$region1),]
    data_input = data[,-c(2,3)]
    
    # Temporally remove outlier
    if(outlier.rm == TRUE){
      outlier <- unique(boxplot.stats(data_input$IF)$out)
      outlier.pos <- which(data_input$IF %in% outlier)
      outlier.value <- data_input$IF[outlier.pos]
      data_input$IF[outlier.pos] <- NA #replace outlier with NA for now
    }
    
    require(mice)
    require(lattice)
    library(miceadds)
    
    ## Classify variable class
    data_input$Single_cell <- as.character(data_input$Single_cell)
    #data_input$Cell_type= as.numeric(as.factor(data_input$Cell_type))
    
    ## Classify method and predictor matrix
    rm.na.data = na.omit(data)
    rm.na.data.input = na.omit(data_input)
    if( all(rm.na.data.input$IF == rm.na.data.input$IF[1]) || anyDuplicated(rm.na.data.input$Single_cell) == 0 || all(data.frame(table(rm.na.data$bin))$Freq == 1) || length(find.collinear(data_input)) > 0 ){ # if all vector only have IF = 1 or if there is only 1 value per bin_single cell or if there is collinearity happening in IF
      agg_new_if2 <- data_input$IF
      agg_new_if2[is.na(agg_new_if2)] <- round(mean(rm.na.data.input$IF, na.rm = T))
    } else if( nrow(rm.na.data.input) ==0 ){ #if the distance is entirely missing, we input 0 for now
      agg_new_if2 <- data_input$IF
      agg_new_if2[is.na(agg_new_if2)] <- rep(1, nrow(data_input))
      
    } else {
      # get initial default imputation setting
      set.seed(123)
      
      ini <- suppressWarnings(mice(data_input, maxit = 0))
      #set up method
      meth= ini$meth
      meth['IF'] <- "rf" 
      #set up predictor matrix
      pred = ini$pred
      #pred[,'Single_cell'] <- -2 #random intercept
      
      ## Imputation
      #simulate for n time of multiple imputations, with 5 iteration
      imp <- suppressWarnings(mice(data_input, method = meth, predictorMatrix = pred, print = FALSE,maxit = maxit, set.seed(123), m=n_imputation) )
      #returns the long format of all multiple imputation
      imp_data = complete(imp, action = 'long', include = F) 
      #if the vector has all if>1, aggregate mean of all imputed complete data
      agg_new_if2 = round(aggregate(imp_data[,4] , by = list(imp_data$.id),FUN= mean))
      agg_new_if2 = agg_new_if2$x}
    
    # Add back outlier, if they were removed in previous steps
    if(outlier.rm == TRUE){
      agg_new_if2[outlier.pos] <- outlier.value
    }
    
    # Create new scHicTable
    new_data <- data.frame(Single_cell = data$Single_cell, region1 = data$region1,
                           region2 = data$region2, IF = agg_new_if2)
    ## Transform new_data into wide format and Re-format to new scHicTable
    library(tidyr)
    update_new_table <- pivot_wider(new_data, names_from = Single_cell, values_from = IF)
    update_new_table <- data.frame(region1 = update_new_table$region1,
                                   region2 = update_new_table$region2,
                                   cell = rep(unique(scHiC.table$cell), nrow(update_new_table)),
                                   chr  = rep(unique(scHiC.table$chr), nrow(update_new_table)), 
                                   update_new_table[,-c(1:2)])
    colnames(update_new_table)[5:ncol(update_new_table)] <- paste0("IF_", 1:(ncol(update_new_table)-4))
    
    new_table <- rbind(new_table, update_new_table)
    
    
  }
  return(new_table)
}  





###### Entire of Random Forest Process
RF_process <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE, maxit = 1,
                                   main_Distances = 1:10000000, missPerc.threshold = 95){
  
  
  ############################# Identify important Pools and Parameter 
  
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_all.distance = length(D)
  n_distance = length(D)
  
  ## Identify main Distances
  maxD = max(main_Distances)
  maxD_res = maxD/res
  main_D_range = 0:maxD_res
  
  
  
  ############################# Imputation process 
  new_table = NULL
  na_perc_all = NULL
  for (i in D){
    
    ## Test na percent
    #print(i)
    data.pool <- predictorMatrixNP_sc_D(scHiC.table, distance = i)
    data_input = data.pool$IF
    na_perc = ( sum(is.na(data_input)) / length(data_input) ) *100
    na_perc_all = c(na_perc_all, na_perc)
  }
  
  
  ### If all Distance have NA < 95%
  if(length(which(na_perc_all > missPerc.threshold)) == 0){
    new_table = RF_impute.outrm.schic(scHiC.table = scHiC.table, n_imputation = n_imputation,
                                                      outlier.rm = outlier.rm, maxit = maxit)
    
    
  } else { ### if any pool have NA > 95%
    
    ## Check which pool have above 95% -> impute mean only
    D_outside.mainD = D[!(D %in% main_D_range)]
    D_above.t = D[na_perc_all > missPerc.threshold]
    D_outside.mainD_obove.t = D_outside.mainD[D_outside.mainD %in% D_above.t]
    
    ## Impute
    require(tidyr)
    l = n_distance-1
    new_table = NULL
    # only impute on distance data have single cell>1
    for( i in 1:n_distance){
      D_pos = 0:l
      cat(paste0('band ', D_pos[i] ,', ') )
      data = predictorMatrix_sc_D(scHiC.table, distance = D_pos[i])
      data =data[order(data$region1),]
      data_input = data[,-c(2,3)]
      
      # Temporally remove outlier
      if(outlier.rm == TRUE){
        outlier <- unique(boxplot.stats(data_input$IF)$out)
        outlier.pos <- which(data_input$IF %in% outlier)
        outlier.value <- data_input$IF[outlier.pos]
        data_input$IF[outlier.pos] <- NA #replace outlier with NA for now
      }
      
      require(mice)
      require(lattice)
      library(miceadds)
      
      ## Classify variable class
      data_input$Single_cell <- as.character(data_input$Single_cell)
      #data_input$Cell_type= as.numeric(as.factor(data_input$Cell_type))
      
      ## Classify method and predictor matrix
      rm.na.data = na.omit(data)
      rm.na.data.input = na.omit(data_input)
      if( D_pos[i] %in% D_outside.mainD_obove.t ||
          all(rm.na.data.input$IF == rm.na.data.input$IF[1]) ||
          anyDuplicated(rm.na.data.input$Single_cell) == 0 ||
          all(data.frame(table(rm.na.data$bin))$Freq == 1) || 
          length(find.collinear(data_input)) > 0 ){ # if all vector only have IF = 1 or if there is only 1 value per bin_single cell or if there is collinearity happening in IF
        if( nrow(rm.na.data.input) ==0 ){ #if the distance is entirely missing, we input 0 for now
          agg_new_if2 <- data_input$IF
          agg_new_if2[is.na(agg_new_if2)] <- rep(1, nrow(data_input))
        } else {
        agg_new_if2 <- data_input$IF
        agg_new_if2[is.na(agg_new_if2)] <- round(mean(rm.na.data.input$IF, na.rm = T))
        }
        
      } else {
        # get initial default imputation setting
        set.seed(123)
        
        ini <- suppressWarnings(mice(data_input, maxit = 0))
        #set up method
        meth= ini$meth
        meth['IF'] <- "rf" 
        #set up predictor matrix
        pred = ini$pred
        #pred[,'Single_cell'] <- -2 #random intercept
        
        ## Imputation
        #simulate for n time of multiple imputations, with 5 iteration
        imp <- suppressWarnings(mice(data_input, method = meth, predictorMatrix = pred, print = FALSE,maxit = maxit, set.seed(123), m=n_imputation) )
        #returns the long format of all multiple imputation
        imp_data = complete(imp, action = 'long', include = F) 
        #if the vector has all if>1, aggregate mean of all imputed complete data
        agg_new_if2 = round(aggregate(imp_data[,4] , by = list(imp_data$.id),FUN= mean))
        agg_new_if2 = agg_new_if2$x}
      
      # Add back outlier, if they were removed in previous steps
      if(outlier.rm == TRUE){
        agg_new_if2[outlier.pos] <- outlier.value
      }
      
      # Create new scHicTable
      new_data <- data.frame(Single_cell = data$Single_cell, region1 = data$region1,
                             region2 = data$region2, IF = agg_new_if2)
      ## Transform new_data into wide format and Re-format to new scHicTable
      library(tidyr)
      update_new_table <- pivot_wider(new_data, names_from = Single_cell, values_from = IF)
      update_new_table <- data.frame(region1 = update_new_table$region1,
                                     region2 = update_new_table$region2,
                                     cell = rep(unique(scHiC.table$cell), nrow(update_new_table)),
                                     chr  = rep(unique(scHiC.table$chr), nrow(update_new_table)), 
                                     update_new_table[,-c(1:2)])
      colnames(update_new_table)[5:ncol(update_new_table)] <- paste0("IF_", 1:(ncol(update_new_table)-4))
      
      new_table <- rbind(new_table, update_new_table)
    }
    
  }
  
  return(new_table)
}









#'Random Forest Imputation with Pooling options for scHi-C Data
#'
#' This function performs imputation of single-cell Hi-C (scHi-C) interaction frequencies (IF) using Random Forest imputation methods with option of distance-based pooling strategies.
#'
#' @param scHiC.table A data frame containing interaction frequencies across single cells, created by the `scHiC_table` function. 
#'                    The first four columns should represent 'cell', 'chr', 'region1', and 'region2', 
#'                    followed by columns representing interaction frequencies ('IF') for individual cells.
#' @param n.imputation An integer specifying the number of imputations to be performed. Default is 5.
#' @param maxit An integer specifying the number of iterations for the internal refinement process within a single imputation cycle. Default is 1.
#' @param outlier.rm A logical value indicating whether to remove outliers during the imputation process. Default is TRUE.
#' @param main.Distances A vector of integers or 'full' representing the scHiC data in main distance range to focus the imputation on, in bp units (e.g., 1:1,000,000). Genomic distances (in bp) is the number of base pairs between two regions in the genome (e.g., loci or bins).
#' Default is from 1 to 10,000,000. 
#' @param pool.style A string specifying the pooling technique to use. Options are 'none', 'progressive' or 'fibonacci'. If 'none' is specified, the method apply imputation with individual genomic distance (consider a band include 1 Distance)
#' Default is 'progressive'.
#' @param missPerc.theshold An integer specifying the missing value percentage threshold in each pool band. 
#' 
#' @return A table in the format of an scHiC table (same structure as the output of the `scHiC_table` function) with imputed interaction frequencies (IF) across all single cells. The output table is formatted with regions and single cells 
#' in wide format, with one column per single cell containing imputed IF values.
#'
#' @details 
#' The function first identifies important pools based on the given scHi-C genomic distance effect by 
#' pooling distance data according to the chosen method. For progressive pooling, pools of distances are consecutively combined 
#' to form larger sets, while Fibonacci pooling uses a Fibonacci sequence to combine distances. If `none` pooling style is selected, the band contain individual 1 genomic distance.
#' During the imputation process, the function imputes all missing values (NAs) within each pool within the main distance range. For distances outside this main focus range, if any pool contains more than `missPerc.theshold` missing values, it triggers an alternative imputation 
#' method, filling in missing values based on the mean for distances.
#' 
#' @references
#' 
#' Doove, L.L., van Buuren, S., Dusseldorp, E. (2014), Recursive partitioning for missing data imputation in the presence of interaction Effects. Computational Statistics & Data Analysis, 72, 92-104.
#' 
#' Shah, A.D., Bartlett, J.W., Carpenter, J., Nicholas, O., Hemingway, H. (2014), Comparison of random forest and parametric imputation models for imputing missing data using MICE: A CALIBER study. American Journal of Epidemiology, doi:10.1093/aje/kwt312.
#' 
#' Van Buuren, S. (2018). Flexible Imputation of Missing Data. Second Edition. Chapman & Hall/CRC. Boca Raton, FL.
#' 
#' @examples
#' # Load MG data folder example
#' load_example_MGfolder()
#' #Create scHicCompare table to be used in scHicCompare
#' IF_table <- scHiC_table(file.path = "MGs_example", cell.type = 'MG', position.dataset =  1:50, type='txt', select.chromosome = 'chr22')
#' # Example usage of Pooling_RF_impute
#' imputed_table <- scHiCcompare_impute(IF_table, n.imputation = 5, outlier.rm = TRUE, 
#'                                   main.Distances = 1:10000000, pool.style = 'progressive')
#'
#' @export




##################### Entire imputation process ##################### 

scHiCcompare_impute <-  function(scHiC.table, n.imputation = 5,  maxit = 1, outlier.rm = TRUE, 
                               main.Distances = 1:10000000, pool.style = 'progressive', missPerc.threshold = 95
                               ){

  ###################################################################################################
  ############################# Identify important Pools and Parameter ############################# 
  
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region2 - scHiC.table$region1)/res )
  n_all.distance = length(D)
  
  ## Identify main Distances
  if(is.infinite(max(main.Distances)) ){
    main_D_range = 1:max(D)
  } else {
    maxD = max(main.Distances)
    maxD_res = maxD/res
    main_D_range = 1:maxD_res
  }
  
  
  ###################################################################################################
  ############################# Imputation process ##################################################
  
  ## If Pooling styles is applied
  if(pool.style == 'none'){     # If None pooling is selected  
    new_table = RF_process(scHiC.table = scHiC.table, n_imputation = n.imputation,
                                       outlier.rm = outlier.rm , main_Distances = main.Distances,
                                       missPerc.threshold = missPerc.threshold)
    
    
  } else { # If Pooling is selected  
    
    ## List of all Distance Pool
    if(pool.style == 'progressive'){
      Dpool.list <- .all_progressive_pooling(vector_distance = D)
    } else if ( pool.style == 'fibonancci'){
      Dpool.list <- .all_fib_pooling(vector_distance = D)
    }
    length.Dpool.list <- length(Dpool.list)
  
    new_table = NULL; new_table_list = rep(list(NULL), (length.Dpool.list) )
    na_perc_all = NULL
    for (i in 1 : length.Dpool.list){
      ## Test na percent
      #print(i)
      data.pool <- do.call(rbind, lapply(Dpool.list[[i]], function(x) { 
        data <- predictorMatrixNP_sc_D(scHiC.table, distance = x)
        return(data)
      }))
      data_input = data.pool$IF
      na_perc = ( sum(is.na(data_input)) / length(data_input) ) *100
      na_perc_all = c(na_perc_all, na_perc)
      ## Test if any element in pools D within main Distances range
      which.in.mainD = Dpool.list[[i]] %in% main_D_range
      mainD.range = any(which.in.mainD)
    }
    
    
    ### If all Distance have NA < `missPerc.theshold`
    if(length(which(na_perc_all> missPerc.threshold)) == 0){
      new_table =  pools_impute(scHiC.table = scHiC.table, n.imputation = n.imputation,
                                outlier.rm = outlier.rm, pool.style = pool.style)
      
      
    } else { ### if any pool have NA > `missPerc.theshold`
      ## Check which pool have above `missPerc.theshold` -> impute mean only
      which_pool_aboveNA = which(na_perc_all>missPerc.threshold)
      pool_aboveNA_Dlist = lapply(which_pool_aboveNA, function(x) return(Dpool.list[[x]])) # list of detail D in which_pool_aboveNA
      pool_aboveNA_D = unlist(pool_aboveNA_Dlist)
      
      ## check which pool above NA level and in main distance
      pool_aboveNA_mainD = which_pool_aboveNA[  unlist(lapply(which_pool_aboveNA, function(x){ any(main_D_range %in% Dpool.list[[x]])}) )  ]
      pool_aboveNA_NOTmainD = which_pool_aboveNA[!which_pool_aboveNA %in% pool_aboveNA_mainD]
      pool_mainD_imp = c(1 : length.Dpool.list)[ !c(1 : length.Dpool.list) %in% pool_aboveNA_NOTmainD]
      D_aboveNA_NOTmainD =  unlist( lapply(pool_aboveNA_NOTmainD, function(x) return(Dpool.list[[x]])) )# which D above NA level and NOT mainD
      D_aboveNA_mainD =  unlist( lapply(pool_aboveNA_mainD, function(x) return(Dpool.list[[x]])) )
      D_mainP_imp = unlist( lapply(pool_mainD_imp, function(x) return(Dpool.list[[x]])) ) ## Distance in main pool will be imputed
      
      ## Extract any pools that will be imputed
      temp.table = scHiC.table 
      temp.table = temp.table[ ((temp.table$region2 - temp.table$region1)/res) %in% D_mainP_imp, ] 
      ## Impute
      library(tidyr)
      temp.table_imp = pools_impute(scHiC.table = temp.table, n.imputation = n.imputation,
                                    outlier.rm = outlier.rm, pool.style = pool.style)
      
      
      ## Now, impute D_set1 to its mean of that distance
      new_table_list = list()
      for (j in 1:length(pool_aboveNA_NOTmainD)){
        cat(paste(' ,pool bans', pool_aboveNA_NOTmainD[j]))
        data = do.call(rbind, lapply(Dpool.list[[pool_aboveNA_NOTmainD[j]]], function(x) { 
          data <- predictorMatrixNP_sc_D(scHiC.table, distance = x)
          return(data)
        }))
        data = data[order(data$region1), ] #sort data
        data_input = data[,-c(2,3,4)]
        extreme.value_pos = data_input$IF %in% boxplot(data_input$IF, plot = F)$out # remove any extreme value
        extreme.value = data_input$IF[data_input$IF %in% boxplot(data_input$IF, plot = F)$out]
        data_input$IF[extreme.value_pos] = NA
        data_input$IF[is.na(data_input$IF)] <- round(mean(data_input$IF, na.rm = T))
        agg_new_if2 = data_input$IF
        agg_new_if2[extreme.value_pos] = extreme.value
        
        # Create new scHicTable
        new_data <- data.frame(Single_cell = data$Single_cell, region1 = data$region1,
                               region2 = data$region2, IF = agg_new_if2)
        ## Transform new_data into wide format and Re-format to scHicTable
        library(tidyr)
        update_new_table <- pivot_wider(new_data, names_from = Single_cell, values_from = IF)
        update_new_table <- data.frame(region1 = update_new_table$region1,
                                       region2 = update_new_table$region2,
                                       cell = rep(unique(scHiC.table$cell),
                                                  nrow(update_new_table)),
                                       chr  = rep(unique(scHiC.table$chr),
                                                  nrow(update_new_table)),
                                       update_new_table[,-c(1:2)])
        colnames(update_new_table)[5:ncol(update_new_table)] <- paste0("IF_",
                                                                       1:(ncol(update_new_table)-4))
        new_table_list[[j]] = update_new_table
      }
      
      if(max(D) %in% D_aboveNA_NOTmainD){
        # Check if any pool contain backward D (NA replacement). If so, then need to remove it/them
        ## Where does the last D locate in last pool
        lastD.poolPosition <- which(Dpool.list[[length.Dpool.list]] ==  max(D))
        lastPool <- Dpool.list[[length.Dpool.list]]
        # check if last D stay at the end of last pool
        if (lastD.poolPosition < length(lastPool)){
          # identify backward D
          backward_D <- lastPool[ (lastD.poolPosition + 1) : length(lastPool)]
          # identify backward D location in new_table_list of last pool
          new_table_lastPool <- new_table_list[[length((new_table_list))]]
          backward_D_position <- which(abs((new_table_lastPool$region1 - new_table_lastPool$region2) / res) %in% backward_D)
          new_table_list[[length(pool_aboveNA_NOTmainD)]] <- new_table_lastPool[-backward_D_position,]
        }
      }
      
      new_table_D1 <- do.call(rbind, new_table_list)
      
      
      ## Merge for final table
      new_table = rbind(temp.table_imp, new_table_D1)
    }
  }
  
  return(new_table)
}



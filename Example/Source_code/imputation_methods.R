

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


mice_impute <- function(data_input, n_imputation = 5, outlier.rm = TRUE){
  require(rstatix)
  require(mice)
  require(lattice)
  library(miceadds)
  # Temporally remove outlier
  if(outlier.rm == TRUE){
    outlier.threshold <- quantile(data_input$IF, 0.85, na.rm = T) + 3 * IQR(data_input$IF, na.rm = T) ## extreme value
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
  meth['IF'] <- "2l.pmm" 
  #set up predictor matrix
  pred = ini$pred
  pred[,'Single_cell'] <- -2 #random intercept
  
  ## Imputation
  #simulate for n time of multiple imputations, with 5 iteration
  imp <- suppressWarnings(mice(data_input, method = meth, predictorMatrix = pred, print = FALSE,maxit = 1, set.seed(123), m=n_imputation) )
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






# MICE4 + outlier removed ----


## mice: check collinear
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



## Function for MICE imputation with outlier.rm for all Distance data of each cell type
MICE_impute.outrm.schic <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE){
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_distance = length(D)
  
  l = n_distance-1
  new_table = NULL
  # only impute on distance data have single cell>1
  for( i in 1:n_distance){
    D_pos = 0:l
    data = predictorMatrix_sc_D(scHiC.table, distance = D_pos[i])
    data =data[order(data$region1),]
    data_input = data[,-c(2,3)]
    print(i)
    
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
    data_input$Single_cell <- as.numeric(data_input$Single_cell)
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
      meth['IF'] <- "2l.pmm" 
      #set up predictor matrix
      pred = ini$pred
      pred[,'Single_cell'] <- -2 #random intercept
      
      ## Imputation
      #simulate for n time of multiple imputations, with 5 iteration
      imp <- suppressWarnings(mice(data_input, method = meth, predictorMatrix = pred, print = FALSE,maxit = 1, set.seed(123), m=n_imputation) )
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






# MICE4 + mainD + outlier removed ----


MICE_impute_mainD <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE, imputed.main.distance = TRUE){
  # input: scHiC.table, n_imputation, neighbor.d - how far is the number of neighbors that imputed Distance is centered, option for outlier remover, option for impute at main closer distance
  # output: schic table (all single cell in columns) with imputed if
  
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_all.distance = length(D)
  
  l = n_all.distance-1 #index of Distance
  new_table = NULL
  
  ## check if option imputation at main Distance is selected
  if (imputed.main.distance == TRUE) {
    n_distance <- max(1000000 / res, 10) * 2
  } else {
    n_distance <- n_all.distance
  }
  
  # Imputing process
  
  
  ############# Imputed main Distance if the option is selected 
  for (i in 1 : n_distance){
    D_pos = 0:l
    print(i)
    
    ## Impute distance by its self 
    data_i = predictorMatrixNP_sc_D(scHiC.table, distance = D_pos[i])
    data = data_i
    data =data[order(data$region1), ] #sort data 
    data_input = data[,-c(2,3,4)]
    agg_new_if2 = mice_impute(data_input, n_imputation = 5, outlier.rm = outlier.rm) #imputed distance vector
    
    # Create new scHicTable
    new_data <- data.frame(Single_cell = data$Single_cell, region1 = data$region1,
                           region2 = data$region2, IF = agg_new_if2) #all table that have been used to impute
    
    ## Transform new_data into wide format and Re-format to scHicTable
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
  
  
  
  ############# Imputed remaining Distance if the option is selected 
  if (imputed.main.distance == TRUE) {
    remaining_n_distance <- (n_distance+1) : length(D)
    for (i in remaining_n_distance){
      print(i)
      data_i = predictorMatrixNP_sc_D(scHiC.table, distance = D_pos[i])
      data = data_i
      data = data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      # Temporally remove outlier
      if(outlier.rm == TRUE){
        outlier.threshold <- quantile(data_input$IF, 0.85, na.rm = T) + 3 * 
          IQR(data_input$IF, na.rm = T) ## extreme value
        outlier.threshold <- max(2,outlier.threshold)
        outlier <-  unique(na.omit(data_input$IF)[na.omit(data_input$IF) > outlier.threshold] )
        outlier.pos <- which(data_input$IF %in% na.omit(outlier))
        outlier.value <- data_input$IF[outlier.pos]
        data_input$IF[outlier.pos] <- NA #replace outlier with NA for now
      }
      data_input$IF[is.na(data_input$I)] <- round(mean(data_input$IF, na.rm = T))
      agg_new_if2 = data_input$IF
      
      # Create new scHicTable
      new_data <- data.frame(Single_cell = data$Single_cell, region1 = data$region1,
                             region2 = data$region2, IF = agg_new_if2) #all table that have been used to impute
      
      ## Transform new_data into wide format and Re-format to scHicTable
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











# MICE4 + neighbor pooling ----

############# Function for MICE imputation process 

MICE_neighbor_impute.outrm <-  function(scHiC.table, n_imputation = 5, neighbor.d = 3, outlier.rm = TRUE, imputed.main.distance = TRUE){
  # input: scHiC.table, n_imputation, neighbor.d - how far is the number, option for outlier remover, option for impute at main closer distance
  # output: schic table (all single cell in columns) with imputed if
  
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_distance = length(D)
  
  l = n_distance-1 #index of Distance
  new_table = NULL
  
  # Imputing process
  for (i in 1 : n_distance){
    D_pos = 0:l
    print(i)
    
    ## Impute distance by its self and neighbor vectors
    if (i %in% (3 : (n_distance-1))){ 
      ### if D > 1, and NOT last distance
      data_i = predictorMatrixNP_sc_D(scHiC.table, distance = D_pos[i])
      # Neighbor Distances
      d1 = min(i-2, neighbor.d) # number of previous data index used to imputed given distance
      data1 =  do.call(rbind, lapply(c(1:d1), function(x){ #previous distance 
        data <- predictorMatrixNP_sc_D(new_table, distance = D_pos[i-x])
        return(data)
      }))
      d3 = min(n_distance-i, neighbor.d) # number of following data index used to imputed given distance
      data3 = do.call(rbind, lapply(c(1:d3), function(x){ #following distance 
        data <- predictorMatrixNP_sc_D(scHiC.table, distance = D_pos[i+x])
        data[is.na(data$IF), 'IF'] <- round(mean(data$IF, na.rm=T))
        return(data)
      }))
      # Combine given D and its neighbors
      data = rbind(data1, data_i, data3)
      data =data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      agg_new_if2 = mice_impute(data_input, n_imputation = 5, outlier.rm = outlier.rm) #imputed distance vector
    } else if (i == 1 || i == 2) { 
      
      
      ### if D = 0 or D = 1
      data_i = predictorMatrixNP_sc_D(scHiC.table, distance = D_pos[i])
      data = data_i
      data =data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      agg_new_if2 = mice_impute(data_input, n_imputation = 5, outlier.rm = outlier.rm) #imputed distance vector
    } else if (i == n_distance){
      
      
      ### if D = last distance
      d1.l = neighbor.d
      data1.l = do.call(rbind, lapply(c(1:d1.l), function(x){ #previous distance 
        data <- predictorMatrixNP_sc_D(new_table, distance = D_pos[i-x])
        return(data)
      }))
      data_i = predictorMatrixNP_sc_D(scHiC.table, distance = D_pos[i])
      data = rbind(data1.l, data_i)
      data =data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      agg_new_if2 = mice_impute(data_input, n_imputation = 5, outlier.rm = outlier.rm) #imputed distance vector
    }  
    
    
    # Create new scHicTable
    all_data <- data.frame(Single_cell = data$Single_cell, region1 = data$region1,
                           region2 = data$region2, IF = agg_new_if2) #all table that have been used to impute
    new_data <- merge(data_i[ , names(data_i) %in% c('Single_cell', 'region1', 'region2')], all_data)
    
    ## Transform new_data into wide format and Re-format to scHicTable
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





# MICE4 + progressive pooling ----

all_progressive_pooling <- function(vector_distance){
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




MICE_progressivePooling_impute.outrm <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE){
  # input: scHiC.table, n_imputation, and, option for outlier remover, option for impute at main closer distance
  # output: schic table (all single cell in columns) with imputed if
  
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_all.distance = length(D)
  
  ## List of all Distance Pool
  Dpool.list <- all_progressive_pooling(vector_distance = D)
  length.Dpool.list <- length(Dpool.list)
  
  new_table = NULL; new_table_list = rep(list(NULL), (length.Dpool.list) )
  # Imputing process
  for (i in 1 : (length.Dpool.list)){
    print(i)
    
    ## Impute by pools
    data.pool = do.call(rbind, lapply(Dpool.list[[i]], function(x){ #previous distance 
      data <- predictorMatrixNP_sc_D(scHiC.table, distance = x)
      return(data)
    }))
    data = data.pool
    data = data[order(data$region1), ] #sort data 
    data_input = data[,-c(2,3,4)]
    agg_new_if2 = mice_impute(data_input, n_imputation = 5, outlier.rm = outlier.rm) 
    
    # Create new scHicTable
    new_data <- data.frame(Single_cell = data$Single_cell, region1 = data$region1,
                           region2 = data$region2, IF = agg_new_if2) 
    
    ## Transform new_data into wide format and Re-format to scHicTable
    library(tidyr)
    update_new_table <- pivot_wider(new_data, names_from = Single_cell, values_from = IF)
    update_new_table <- data.frame(region1 = update_new_table$region1,
                                   region2 = update_new_table$region2,
                                   cell = rep(unique(scHiC.table$cell), nrow(update_new_table)),
                                   chr  = rep(unique(scHiC.table$chr), nrow(update_new_table)), 
                                   update_new_table[,-c(1:2)])
    colnames(update_new_table)[5:ncol(update_new_table)] <- paste0("IF_", 1:(ncol(update_new_table)-4))
    
    new_table_list[[i]] <- update_new_table
  }
  
  # Check if any pool contain backward D (NA replacement). If so, then need to remove it/them
  ## Where does the last D locate in last pool
  lastD.poolPosition <- which(Dpool.list[[length.Dpool.list]] ==  max(D))
  lastPool <- Dpool.list[[length.Dpool.list]]
  # check if last D stay at the end of last pool
  if (lastD.poolPosition < length(lastPool)){
    # identify backward D
    backward_D <- lastPool[ (lastD.poolPosition + 1) : length(lastPool)]
    # identify backward D location in new_table_list of last pool
    new_table_lastPool <- new_table_list[[length.Dpool.list]]
    backward_D_position <- which(abs((new_table_lastPool$region1 - new_table_lastPool$region2) / res) %in% backward_D) 
    new_table_list[[length.Dpool.list]] <- new_table_lastPool[-backward_D_position,]
  }
  
  ## Merge all element of new_table_list to 1 data frame
  new_table <- do.call(rbind, new_table_list)
  
  return(new_table)
}  




# MICE4 +  fibonacci pooling ----

fibonacciR <- function(x){
  if (x < 2){
    out = x
  } else {
    out = fibonacciR(x - 1) + fibonacciR(x - 2)
  }
  return(out)
}

all_fib_pooling <- function(vector_distance){
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
  
  poolings_list[[1]] <- 1
  return(poolings_list)
  
}





MICE_fibonacciPooling_impute.outrm <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE){
  # input: scHiC.table, n_imputation, option for impute at main closer distance
  # output: schic table (all single cell in columns) with imputed if
  
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_all.distance = length(D)
  
  ## List of all Distance Pool
  Dpool.list <- all_fib_pooling(vector_distance = D)
  length.Dpool.list <- length(Dpool.list)
  
  new_table = NULL; new_table_list = rep(list(NULL), (length.Dpool.list) )
  # Imputing process
  for (i in 1 : (length.Dpool.list)){
    print(i)
    
    ## Impute by pools
    data.pool = do.call(rbind, lapply(Dpool.list[[i]], function(x){ #previous distance 
      data <- predictorMatrixNP_sc_D(scHiC.table, distance = x)
      return(data)
    }))
    data = data.pool
    data = data[order(data$region1), ] #sort data 
    data_input = data[,-c(2,3,4)]
    agg_new_if2 = mice_impute(data_input, n_imputation = 5, outlier.rm = outlier.rm) 
    
    # Create new scHicTable
    new_data <- data.frame(Single_cell = data$Single_cell, region1 = data$region1,
                           region2 = data$region2, IF = agg_new_if2) 
    
    ## Transform new_data into wide format and Re-format to scHicTable
    library(tidyr)
    update_new_table <- pivot_wider(new_data, names_from = Single_cell, values_from = IF)
    update_new_table <- data.frame(region1 = update_new_table$region1,
                                   region2 = update_new_table$region2,
                                   cell = rep(unique(scHiC.table$cell), nrow(update_new_table)),
                                   chr  = rep(unique(scHiC.table$chr), nrow(update_new_table)), 
                                   update_new_table[,-c(1:2)])
    colnames(update_new_table)[5:ncol(update_new_table)] <- paste0("IF_", 1:(ncol(update_new_table)-4))
    
    new_table_list[[i]] <- update_new_table
  }
  
  # Check if any pool contain backward D (NA replacement). If so, then need to remove it/them
  ## Where does the last D locate in last pool
  lastD.poolPosition <- which(Dpool.list[[length.Dpool.list]] ==  max(D))
  lastPool <- Dpool.list[[length.Dpool.list]]
  # check if last D stay at the end of last pool
  if (lastD.poolPosition < length(lastPool)){
    # identify backward D
    backward_D <- lastPool[ (lastD.poolPosition + 1) : length(lastPool)]
    # identify backward D location in new_table_list of last pool
    new_table_lastPool <- new_table_list[[length.Dpool.list]]
    backward_D_position <- which(abs((new_table_lastPool$region1 - new_table_lastPool$region2) / res) %in% backward_D) 
    new_table_list[[length.Dpool.list]] <- new_table_lastPool[-backward_D_position,]
  }
  
  ## Merge all element of new_table_list to 1 data frame
  new_table <- do.call(rbind, new_table_list)
  
  return(new_table)
}  



# Mean Distance imputation ----

meanD_impute <-  function(scHiC.table){
  # input: scHiC.table, n_imputation, option for impute at main closer distance
  # output: schic table (all single cell in columns) with imputed if
  
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D = unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_all.distance = length(D)
  D_pos = 1: n_all.distance
  new_table = NULL; new_table_list = rep(list(NULL), (n_all.distance) )
  # Imputing process
  for (i in 1 : n_all.distance){
    print(i)
    
    ## Impute by mean of D
    data_i = predictorMatrixNP_sc_D(scHiC.table, distance = D[i])
    data = data_i
    data = data[order(data$region1), ] #sort data 
    data_input = data[,-c(2,3,4)]
    data_input$IF[is.na(data_input$IF)] <- round(mean(data_input$IF, na.rm = T))
    agg_new_if2 = data_input$IF
    
    # Create new scHicTable
    new_data <- data.frame(Single_cell = data$Single_cell, region1 = data$region1,
                           region2 = data$region2, IF = agg_new_if2) 
    
    ## Transform new_data into wide format and Re-format to scHicTable
    library(tidyr)
    update_new_table <- pivot_wider(new_data, names_from = Single_cell, values_from = IF)
    update_new_table <- data.frame(region1 = update_new_table$region1,
                                   region2 = update_new_table$region2,
                                   cell = rep(unique(scHiC.table$cell), nrow(update_new_table)),
                                   chr  = rep(unique(scHiC.table$chr), nrow(update_new_table)), 
                                   update_new_table[,-c(1:2)])
    colnames(update_new_table)[5:ncol(update_new_table)] <- paste0("IF_", 1:(ncol(update_new_table)-4))
    
    new_table_list[[i]] <- update_new_table
  }
  
  new_table <- do.call(rbind, new_table_list)
  return(new_table)
}  







# MICE5 RF + outlier rm ----


## mice: check collinear
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



## Function for MICE imputation with outlier.rm for all Distance data of each cell type
MICE5_impute.outrm.schic <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE){
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_distance = length(D)
  
  l = n_distance-1
  new_table = NULL
  # only impute on distance data have single cell>1
  for( i in 1:n_distance){
    D_pos = 0:l
    data = predictorMatrix_sc_D(scHiC.table, distance = D_pos[i])
    data =data[order(data$region1),]
    data_input = data[,-c(2,3)]
    print(i)
    
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
      imp <- suppressWarnings(mice(data_input, method = meth, predictorMatrix = pred, print = FALSE,maxit = 1, set.seed(123), m=n_imputation) )
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










# MICE5 RF + neighbor pooling ----

mice.rf_impute <- function(data_input, n_imputation = 5, outlier.rm = TRUE){
  require(rstatix)
  require(mice)
  require(lattice)
  library(miceadds)
  # Temporally remove outlier
  if(outlier.rm == TRUE){
    outlier.threshold <- quantile(data_input$IF, 0.85, na.rm = T) + 3 * IQR(data_input$IF, na.rm = T) ## extreme value
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
  imp <- suppressWarnings(mice(data_input, method = meth, predictorMatrix = pred, print = FALSE,maxit = 1, set.seed(123), m=n_imputation) )
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




##################### Function for MICE imputation process 

MICE5_neighbor_impute.outrm <-  function(scHiC.table, n_imputation = 5, neighbor.d = 3, outlier.rm = TRUE){
  # input: scHiC.table, n_imputation, neighbor.d - how far is the number, option for outlier remover, option for impute at main closer distance
  # output: schic table (all single cell in columns) with imputed if
  
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_distance = length(D)
  
  l = n_distance-1 #index of Distance
  new_table = NULL
  
  # Imputing process
  for (i in 1 : n_distance){
    D_pos = 0:l
    print(i)
    
    ## Impute distance by its self and neighbor vectors
    if (i %in% (3 : (n_distance-1))){ 
      ### if D > 1, and NOT last distance
      data_i = predictorMatrixNP_sc_D(scHiC.table, distance = D_pos[i])
      # Neighbor Distances
      d1 = min(i-2, neighbor.d) # number of previous data index used to imputed given distance
      data1 =  do.call(rbind, lapply(c(1:d1), function(x){ #previous distance 
        data <- predictorMatrixNP_sc_D(new_table, distance = D_pos[i-x])
        return(data)
      }))
      d3 = min(n_distance-i, neighbor.d) # number of following data index used to imputed given distance
      data3 = do.call(rbind, lapply(c(1:d3), function(x){ #following distance 
        data <- predictorMatrixNP_sc_D(scHiC.table, distance = D_pos[i+x])
        data[is.na(data$IF), 'IF'] <- round(mean(data$IF, na.rm=T))
        return(data)
      }))
      # Combine given D and its neighbors
      data = rbind(data1, data_i, data3)
      data =data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      agg_new_if2 = mice.rf_impute(data_input, n_imputation = 5, outlier.rm = outlier.rm) #imputed distance vector
    } else if (i == 1 || i == 2) { 
      
      
      ### if D = 0 or D = 1
      data_i = predictorMatrixNP_sc_D(scHiC.table, distance = D_pos[i])
      data = data_i
      data =data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      agg_new_if2 = mice.rf_impute(data_input, n_imputation = 5, outlier.rm = outlier.rm) #imputed distance vector
    } else if (i == n_distance){
      
      
      ### if D = last distance
      d1.l = neighbor.d
      data1.l = do.call(rbind, lapply(c(1:d1.l), function(x){ #previous distance 
        data <- predictorMatrixNP_sc_D(new_table, distance = D_pos[i-x])
        return(data)
      }))
      data_i = predictorMatrixNP_sc_D(scHiC.table, distance = D_pos[i])
      data = rbind(data1.l, data_i)
      data =data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      agg_new_if2 = mice.rf_impute(data_input, n_imputation = 5, outlier.rm = outlier.rm) #imputed distance vector
    }  
    
    
    # Create new scHicTable
    all_data <- data.frame(Single_cell = data$Single_cell, region1 = data$region1,
                           region2 = data$region2, IF = agg_new_if2) #all table that have been used to impute
    new_data <- merge(data_i[ , names(data_i) %in% c('Single_cell', 'region1', 'region2')], all_data)
    
    ## Transform new_data into wide format and Re-format to scHicTable
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




# MICE5 RF + progressive pooling ########################
all_progressive_pooling <- function(vector_distance){
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



##################### Function for MICE imputation process 


MICE5_progressivePooling_impute.outrm <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE, parallel = FALSE){
  # input: scHiC.table, n_imputation, and, option for outlier remover, option for impute at main closer distance
  # output: schic table (all single cell in columns) with imputed if
  
  ## distance of scHiC table
  res <- min(abs(diff(unique(scHiC.table$region1))))
  D <- unique(abs(scHiC.table$region1 - scHiC.table$region2) / res)
  n_all.distance <- length(D)
  
  ## List of all Distance Pool
  Dpool.list <- all_fib_pooling(vector_distance = D)
  length.Dpool.list <- length(Dpool.list)
  
  new_table <- NULL
  new_table_list <- rep(list(NULL), length.Dpool.list)
  
  process_pool <- function(i) {
    print(i)
    
    ## Pooling data
    data.pool <- do.call(rbind, lapply(Dpool.list[[i]], function(x) { 
      data <- predictorMatrixNP_sc_D(scHiC.table, distance = x)
      return(data)
    }))
    
    ## Impute by pools
    data <- data.pool
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
      agg_new_if2 <- mice.rf_impute(data_input = data_input, n_imputation = n_imputation, outlier.rm = outlier.rm) 
    }
    
    # Create new scHicTable
    new_data <- data.frame(Single_cell = data$Single_cell, region1 = data$region1,
                           region2 = data$region2, IF = agg_new_if2) 
    
    ## Transform new_data into wide format and Re-format to scHicTable
    update_new_table <- pivot_wider(new_data, names_from = Single_cell, values_from = IF)
    update_new_table <- data.frame(region1 = update_new_table$region1,
                                   region2 = update_new_table$region2,
                                   cell = rep(unique(scHiC.table$cell), nrow(update_new_table)),
                                   chr  = rep(unique(scHiC.table$chr), nrow(update_new_table)), 
                                   update_new_table[,-c(1:2)])
    colnames(update_new_table)[5:ncol(update_new_table)] <- paste0("IF_", 1:(ncol(update_new_table)-4))
    
    return(update_new_table)
  }

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
    
    ## Merge all elements of new_table_list into 1 data frame
    new_table <- do.call(rbind, new_table_list)
    
    
  } else {
    new_table_list <- bplapply(1:length.Dpool.list, process_pool)
    
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
  } 

  
  return(new_table)
}





# MICE5 RF + fibonacci pooling ########################
fibonacciR <- function(x){
  if (x < 2){
    out = x
  } else {
    out = fibonacciR(x - 1) + fibonacciR(x - 2)
  }
  return(out)
}

all_fib_pooling <- function(vector_distance){
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
  
  poolings_list[[1]] <- 1
  return(poolings_list)
  
}




################ Function for MICE imputation process 


MICE5_fibonacciPooling_impute.outrm <- function(scHiC.table, n_imputation = 5, outlier.rm = TRUE, parallel = FALSE){
  # input: scHiC.table, n_imputation, and, option for outlier remover, option for impute at main closer distance
  # output: schic table (all single cell in columns) with imputed if
  
  ## distance of scHiC table
  res <- min(abs(diff(unique(scHiC.table$region1))))
  D <- unique(abs(scHiC.table$region1 - scHiC.table$region2) / res)
  n_all.distance <- length(D)
  
  ## List of all Distance Pool
  Dpool.list <- all_fib_pooling(vector_distance = D)
  length.Dpool.list <- length(Dpool.list)
  
  new_table <- NULL
  new_table_list <- rep(list(NULL), length.Dpool.list)
  
  process_pool <- function(i) {
    print(i)
    
    ## Pooling data
    data.pool <- do.call(rbind, lapply(Dpool.list[[i]], function(x) { 
      data <- predictorMatrixNP_sc_D(scHiC.table, distance = x)
      return(data)
    }))
    
    ## Impute by pools
    data <- data.pool
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
      agg_new_if2 <- mice.rf_impute(data_input = data_input, n_imputation = n_imputation, outlier.rm = outlier.rm) 
    }
    
    # Create new scHicTable
    new_data <- data.frame(Single_cell = data$Single_cell, region1 = data$region1,
                           region2 = data$region2, IF = agg_new_if2) 
    
    ## Transform new_data into wide format and Re-format to scHicTable
    update_new_table <- pivot_wider(new_data, names_from = Single_cell, values_from = IF)
    update_new_table <- data.frame(region1 = update_new_table$region1,
                                   region2 = update_new_table$region2,
                                   cell = rep(unique(scHiC.table$cell), nrow(update_new_table)),
                                   chr  = rep(unique(scHiC.table$chr), nrow(update_new_table)), 
                                   update_new_table[,-c(1:2)])
    colnames(update_new_table)[5:ncol(update_new_table)] <- paste0("IF_", 1:(ncol(update_new_table)-4))
    
    return(update_new_table)
  }
  
  if(parallel == FALSE){
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
    
    
  } else {
    new_table_list <- bplapply(1:length.Dpool.list, process_pool)
    
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
  } 
  
  
  return(new_table)
}









# MICE5 RF progressive + 95% NA rule #####

progressivePooling_impute.outrm <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE, main_Distances = 1:10000000){
  # input: scHiC.table, n_imputation, and, option for outlier remover, option for impute at main closer distance
  # output: schic table (all single cell in columns) with imputed if
  
  
  ###################################################################################################
  ############################# Identify important Pools and Parameter ############################# 
  
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_all.distance = length(D)
  
  ## Identify main Distances
  maxD = max(main_Distances)
  maxD_res = maxD/res
  main_D_range = 1:maxD_res
  
  
  ## List of all Distance Pool
  Dpool.list <- all_progressive_pooling(vector_distance = D)
  length.Dpool.list <- length(Dpool.list)
  
  
  ###################################################################################################
  ############################# Imputation process ##################################################
  new_table = NULL; new_table_list = rep(list(NULL), (length.Dpool.list) )
  for (i in 1 : length.Dpool.list){
    print(i)
    
    ## Test na percent
    data.pool <- do.call(rbind, lapply(Dpool.list[[i]], function(x) { 
      data <- predictorMatrixNP_sc_D(scHiC.table, distance = x)
      return(data)
    }))
    data_input = data.pool$IF
    na_perc = ( sum(is.na(data_input)) / length(data_input) ) *100
    
    ## Test if any element in pools D within main Distances range
    which.in.mainD = Dpool.list[[i]] %in% main_D_range
    mainD.range = any(which.in.mainD)
    
    
    
    ############################# Impute main Distance  
    ## Separate into 3 cases
    if(mainD.range == TRUE){
      
      data = data.pool
      data = data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      rm.na.data.input = na.omit(data_input)
      
      if( all(rm.na.data.input$IF == rm.na.data.input$IF[1]) ||
          anyDuplicated(rm.na.data.input$Single_cell) == 0  ||
          length(find.collinear(data_input)) > 0 ){ # if all vector only have constant IF or if there is only 1 value per bin_single cell or if there is colinearity happening in IF
        agg_new_if2 <- data_input$IF
        agg_new_if2[is.na(agg_new_if2)] <- 0
        agg_new_if2[agg_new_if2 == 0] <- round(mean(agg_new_if2, na.rm = T))
      } else if( nrow(rm.na.data.input) ==0 ){ #if the distance is entirely missing, we input 0 for now
        agg_new_if2 <- data_input$IF
        agg_new_if2[is.na(agg_new_if2)] <- rep(0, nrow(data_input))
        
      } else { 
        agg_new_if2 = mice.rf_impute(data_input=data_input, n_imputation = n_imputation,
                                     outlier.rm = outlier.rm) 
      }
      
      
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
      
      new_table_list[[i]] <- update_new_table
      
      
    } else if(mainD.range == FALSE & na_perc < 90){
      ############################# Impute Remmaining Distance that NA percent < 90%  
      data = data.pool
      data = data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      rm.na.data.input = na.omit(data_input)
      
      if( all(rm.na.data.input$IF == rm.na.data.input$IF[1]) ||
          anyDuplicated(rm.na.data.input$Single_cell) == 0  ||
          length(find.collinear(data_input)) > 0 ){ # if all vector only have constant IF or if there is only 1 value per bin_single cell or if there is colinearity happening in IF
        agg_new_if2 <- data_input$IF
        agg_new_if2[is.na(agg_new_if2)] <- 0
        agg_new_if2[agg_new_if2 == 0] <- round(mean(agg_new_if2, na.rm = T))
      } else if( nrow(rm.na.data.input) ==0 ){ #if the distance is entirely missing, we input 0 for now
        agg_new_if2 <- data_input$IF
        agg_new_if2[is.na(agg_new_if2)] <- rep(0, nrow(data_input))
        
      } else { 
        agg_new_if2 = mice.rf_impute(data_input=data_input, n_imputation = n_imputation,
                                     outlier.rm = outlier.rm) 
      }
      
      
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
      
      new_table_list[[i]] <- update_new_table
      
      
      
    } else if(mainD.range == FALSE & na_perc > 90){
      ############################# Impute Remmaining Distance that NA percent > 90%  
      data = data.pool
      data = data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      data_input$IF[is.na(data_input$IF)] <- round(mean(data_input$IF, na.rm = T))
      agg_new_if2 = data_input$IF
      
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
      
      new_table_list[[i]] <- update_new_table
      
    }
    
  }
  
  # Check if any pool contain backward D (NA replacement). If so, then need to remove it/them
  ## Where does the last D locate in last pool
  lastD.poolPosition <- which(Dpool.list[[length.Dpool.list]] ==  max(D))
  lastPool <- Dpool.list[[length.Dpool.list]]
  # check if last D stay at the end of last pool
  if (lastD.poolPosition < length(lastPool)){
    # identify backward D
    backward_D <- lastPool[ (lastD.poolPosition + 1) : length(lastPool)]
    # identify backward D location in new_table_list of last pool
    new_table_lastPool <- new_table_list[[length.Dpool.list]]
    backward_D_position <- which(abs((new_table_lastPool$region1 - new_table_lastPool$region2) / res) %in% backward_D) 
    new_table_list[[length.Dpool.list]] <- new_table_lastPool[-backward_D_position,]
  }
  
  new_table <- do.call(rbind, new_table_list)
  
  return(new_table)
}






# MICE5 RF fibonacci + 95% NA rule #####

fibPooling_impute.outrm <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE, main_Distances = 1:10000000){
  # input: scHiC.table, n_imputation, and, option for outlier remover, option for impute at main closer distance
  # output: schic table (all single cell in columns) with imputed if
  
  
  ###################################################################################################
  ############################# Identify important Pools and Parameter ############################# 
  
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_all.distance = length(D)
  
  ## Identify main Distances
  maxD = max(main_Distances)
  maxD_res = maxD/res
  main_D_range = 1:maxD_res
  
  
  ## List of all Distance Pool
  Dpool.list <- all_fib_pooling(vector_distance = D)
  length.Dpool.list <- length(Dpool.list)
  
  
  ###################################################################################################
  ############################# Imputation process ##################################################
  new_table = NULL; new_table_list = rep(list(NULL), (length.Dpool.list) )
  for (i in 1 : length.Dpool.list){
    print(i)
    
    ## Test na percent
    data.pool <- do.call(rbind, lapply(Dpool.list[[i]], function(x) { 
      data <- predictorMatrixNP_sc_D(scHiC.table, distance = x)
      return(data)
    }))
    data_input = data.pool$IF
    na_perc = ( sum(is.na(data_input)) / length(data_input) ) *100
    
    ## Test if any element in pools D within main Distances range
    which.in.mainD = Dpool.list[[i]] %in% main_D_range
    mainD.range = any(which.in.mainD)
    
    
    
    ############################# Impute main Distance  
    ## Separate into 3 cases
    if(mainD.range == TRUE){
      
      data = data.pool
      data = data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      rm.na.data.input = na.omit(data_input)
      
      if( all(rm.na.data.input$IF == rm.na.data.input$IF[1]) ||
          anyDuplicated(rm.na.data.input$Single_cell) == 0  ||
          length(find.collinear(data_input)) > 0 ){ # if all vector only have constant IF or if there is only 1 value per bin_single cell or if there is colinearity happening in IF
        agg_new_if2 <- data_input$IF
        agg_new_if2[is.na(agg_new_if2)] <- 0
        agg_new_if2[agg_new_if2 == 0] <- round(mean(agg_new_if2, na.rm = T))
      } else if( nrow(rm.na.data.input) ==0 ){ #if the distance is entirely missing, we input 0 for now
        agg_new_if2 <- data_input$IF
        agg_new_if2[is.na(agg_new_if2)] <- rep(0, nrow(data_input))
        
      } else { 
        agg_new_if2 = mice.rf_impute(data_input=data_input, n_imputation = n_imputation,
                                     outlier.rm = outlier.rm) 
      }
      
      
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
      
      new_table_list[[i]] <- update_new_table
      
      
    } else if(mainD.range == FALSE & na_perc < 90){
      ############################# Impute Remmaining Distance that NA percent < 90%  
      data = data.pool
      data = data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      rm.na.data.input = na.omit(data_input)
      
      if( all(rm.na.data.input$IF == rm.na.data.input$IF[1]) ||
          anyDuplicated(rm.na.data.input$Single_cell) == 0  ||
          length(find.collinear(data_input)) > 0 ){ # if all vector only have constant IF or if there is only 1 value per bin_single cell or if there is colinearity happening in IF
        agg_new_if2 <- data_input$IF
        agg_new_if2[is.na(agg_new_if2)] <- 0
        agg_new_if2[agg_new_if2 == 0] <- round(mean(agg_new_if2, na.rm = T))
      } else if( nrow(rm.na.data.input) ==0 ){ #if the distance is entirely missing, we input 0 for now
        agg_new_if2 <- data_input$IF
        agg_new_if2[is.na(agg_new_if2)] <- rep(0, nrow(data_input))
        
      } else { 
        agg_new_if2 = mice.rf_impute(data_input=data_input, n_imputation = n_imputation,
                                     outlier.rm = outlier.rm) 
      }
      
      
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
      
      new_table_list[[i]] <- update_new_table
      
      
      
    } else if(mainD.range == FALSE & na_perc > 90){
      ############################# Impute Remmaining Distance that NA percent > 90%  
      data = data.pool
      data = data[order(data$region1), ] #sort data 
      data_input = data[,-c(2,3,4)]
      data_input$IF[is.na(data_input$IF)] <- round(mean(data_input$IF, na.rm = T))
      agg_new_if2 = data_input$IF
      
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
      
      new_table_list[[i]] <- update_new_table
      
    }
    
  }
  
  # Check if any pool contain backward D (NA replacement). If so, then need to remove it/them
  ## Where does the last D locate in last pool
  lastD.poolPosition <- which(Dpool.list[[length.Dpool.list]] ==  max(D))
  lastPool <- Dpool.list[[length.Dpool.list]]
  # check if last D stay at the end of last pool
  if (lastD.poolPosition < length(lastPool)){
    # identify backward D
    backward_D <- lastPool[ (lastD.poolPosition + 1) : length(lastPool)]
    # identify backward D location in new_table_list of last pool
    new_table_lastPool <- new_table_list[[length.Dpool.list]]
    backward_D_position <- which(abs((new_table_lastPool$region1 - new_table_lastPool$region2) / res) %in% backward_D) 
    new_table_list[[length.Dpool.list]] <- new_table_lastPool[-backward_D_position,]
  }
  
  new_table <- do.call(rbind, new_table_list1)
  
  return(new_table)
}


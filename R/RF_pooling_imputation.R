################## Checking collinear ################## 
.find.collinear <- function(x, threshold = 0.999, ...) {
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
pools_impute <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE,
                          pool_style = 'progressive'){
  # input: scHiC.table, n_imputation, and, option for outlier remover, option for impute at main closer distance
  # output: schic table (all single cell in columns) with imputed if
  
  ## distance of scHiC table
  res = min( abs( diff(unique(scHiC.table$region1)) ) )
  D =unique( abs(scHiC.table$region1 - scHiC.table$region2)/res )
  n_all.distance = length(D)
  
  ## List of all Distance Pool
  if(pool_style == 'progressive'){
    Dpool.list <- all_progressive_pooling(vector_distance = D)
  } else if ( pool_style == 'fibonancci'){
    Dpool.list <- all_fib_pooling(vector_distance = D)
  }
 
  length.Dpool.list <- length(Dpool.list)
  
  new_table <- NULL
  new_table_list <- rep(list(NULL), length.Dpool.list)
  
  
  ## Impute all pools by RF
  process_pool <- function(i) {
    print(i)
    
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
      agg_new_if2 <- mice.rf_impute(data_input = data_input, n_imputation = n_imputation, outlier.rm = outlier.rm) 
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














##################### Entire Pooling RF imputation process ##################### 

Pooling_RF_impute <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE, 
                               main_Distances = 1:10000000, pool_style = 'progressive'
                               ){
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
  if(pool_style == 'progressive'){
    Dpool.list <- all_progressive_pooling(vector_distance = D)
  } else if ( pool_style == 'fibonancci'){
    Dpool.list <- all_fib_pooling(vector_distance = D)
  }
  length.Dpool.list <- length(Dpool.list)
  
  
  ###################################################################################################
  ############################# Imputation process ##################################################
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
  
  
  ### If all Distance have NA < 95%
  if(length(which(na_perc_all>95)) == 0){
    new_table =  pools_impute(scHiC.table = scHiC.table, n_imputation = n_imputation,
                                                      outlier.rm = outlier.rm, pool_style = pool_style)
    
    
  } else { ### if any pool have NA > 95%
    ## Check which pool have above 95% -> impute mean only
    which_pool_aboveNA = which(na_perc_all>95)
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
    temp.table_imp = pools_impute(scHiC.table = temp.table, n_imputation = n_imputation,
                                                           outlier.rm = outlier.rm, pool_style = pool_style)
    
    
    ## Now, impute D_set1 to its mean of that distance
    new_table_list = list()
    for (j in 1:length(pool_aboveNA_NOTmainD)){
      print(pool_aboveNA_NOTmainD[j])
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
    
    new_table_D1 <- do.call(rbind, new_table_list)
    
    
    ## Merge for final table
    new_table = rbind(temp.table_imp, new_table_D1)
    
    
  }
  
  return(new_table)
}



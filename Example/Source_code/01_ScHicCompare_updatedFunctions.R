
## Create scHiC.table -----
scHiC.table <- function(cell_type, n_sc , select_chromosome){
  #input - single cell data has column: chr1,start1,chr2,start2,if
  ### n_sc - number of single cell
  #output: table with columns of IF from given single cells: cell, chr, region1, region2, if_i
  regions =NULL
  options(scipen = 999)
  for(i in 1:n_sc){
    #input dataset
    dataset_name <- paste0(cell_type,"_", i) #dataset name
    data <- get(dataset_name) #get dataset have name above
    # Select chromosome
    data <- data[data[,1] == select_chromosome,] # by using information in col1
    # Transform dataset into sparse
    data <- data[,c(2,4,5)]
    names(data) <- c('region1', 'region2','IF')
    # Remove rows with 0 values in any column
    data <- data[!data[,1]==0,]
    data <- data[!data[,2]==0,]
    # Extract single cells regions
    region1sc <- unique(c(data$region1, data$region2))
    regions <- unique(c(regions,region1sc))
  }
  
  # Extract start and end regions
  start.region = min(regions, na.rm = T); end.region = max(regions, na.rm = T); bin = min(abs(diff(regions)), na.rm = T)
  regions.seq = seq(start.region, end.region, by = bin)
  
  # coordination of pair of bin
  grid1 <- data.frame(X1 = 1:length(regions.seq), X2=1:length(regions.seq)) # diagonal line
  grid2 <- data.frame(t(combn(1:length(regions.seq),2))) # off diagnoal line value
  grid <- rbind(grid1, grid2)
  region1 <- regions.seq[grid[,1]]
  region2 <- regions.seq[grid[,2]]
  cordination <- cbind(region1, region2)
  
  
  # Preallocate memory for the table
  table <- data.frame(
    region1 = cordination[,'region1'],
    region2 = cordination[,'region2']
  )
  
  for (i in 1:n_sc) {
    # Input dataset
    dataset_name <- paste0(cell_type, "_", i)
    data <- get(dataset_name)
    
    # Filter rows based on chromosome
    data <- data[data[, 1] == select_chromosome, c(2, 4, 5)]
    names(data)[c(1, 2)] <- c('region1', 'region2')
    names(data)[3] <- paste0('IF_', i)
    
    # Remove rows with 0 values in any column
    data <- data[rowSums(data[, c(1, 2)] == 0) == 0, ]
    
    # Merge data into table
    table <- dplyr::full_join(table, data, by = c('region1', 'region2'), suffix = c('', paste0('_IF_', i)))
  }
  
  # If needed, replace NA values with 0
  table[is.na(table)] <- 0
  
  up.table <- data.frame(
    cell = rep(cell_type, nrow(table)),
    chr = rep(select_chromosome, nrow(table)),
    table
  )
  
  return(up.table)
}



## Create PseudoBulk Sparse ----
pseudo_bulkHic_sparse <- function(scHiC.table){ #output: pseudobulk sparse data
  # sparse matrix of the pseudo bulk
  if_scs = scHiC.table[,-c(1,2,3,4)]
  IF = rowSums(if_scs)
  bulk.table <- cbind(scHiC.table[,c(1,2,3,4)],IF )
  
  #remove regions with IF_bulk = 0
  bulk.table <- bulk.table[!bulk.table$IF==0,]
  
  #remove 'chr' column
  bulk.table <- bulk.table[,!names(bulk.table) %in% c('cell','chr')]
  return(bulk.table)
}


## Create Distance vector from scHic table
### report for bin distance
bin_distance.up <- function(search_distance, scHiC.table){
  resolution <-  min( abs( diff(unique(scHiC.table$region1)) ) )
  table <- cbind(scHiC.table,
                 D = abs((scHiC.table$region1 - scHiC.table$region2)/resolution) )
  distance.table = table[table$D == search_distance,]
  distance.table = distance.table[, !(names(distance.table) %in% c('region1', 'region2', 'cell', 'chr')),]
  cell_list <- c(distance.table)
  vector <- as.vector(unlist(cell_list))
  return(vector)
}




## Control Change ----

## Function to make change DCI on individual single cell
sc_controlDCI.schicTable <- function(groundtruth_control.report, cell_type1.scHiC_table,
                          cell_type2.scHiC_table, select_chromosome,
                          FC=3, sign_interaction='up'){ #output: scHiC_table with DCI changes
  ## extract changes DCI position from ground truth report
  rt.report = groundtruth_control.report
  changes.pos = rt.report[rt.report[['truth']]==1,]
  region1 = changes.pos[['start1']]; region2 = changes.pos[['start2']] #extract which pair of regions will make change
  
  if(sign_interaction=='up') {
    for(j in 1:length(region1)){
      data = cell_type2.scHiC_table
      DCI.pos = which(data$region1 == region1[j] & data$region2 == region2[j])
      data[DCI.pos,-c(1:4)] <- data[DCI.pos,-c(1:4)]* FC %>% as.integer()
    }
    cc.cells_list = list(cell_type1.scHiC_table, data)
    
    
    
  } else if(sign_interaction=='down'){
    for(j in 1:length(region1)){
      data = cell_type1.scHiC_table
      DCI.pos = which(data$region1 == region1[j] & data$region2 == region2[j])
      data[DCI.pos,-c(1:4)] <- data[DCI.pos,-c(1:4)]* FC %>% as.integer()
    }
    cc.cells_list = list(data, cell_type2.scHiC_table)
    
    
    
    
  } else if(sign_interaction=='up_down'){
    up_region1 = changes.pos$'start1'[changes.pos$sign_interaction=='up']
    up_region2 = changes.pos$'start2'[changes.pos$sign_interaction=='up']
    down_region1 = changes.pos$'start1'[changes.pos$sign_interaction=='down']
    down_region2 = changes.pos$'start2'[changes.pos$sign_interaction=='down']
    
    ## Multiply FC to cell_type2 for up-regulated interaction
    dataup = cell_type2.scHiC_table  
    for(j in 1:length(up_region1)){
      DCI.up.pos = which(dataup$region1 == up_region1[j] & dataup$region2 == up_region2[j])
      dataup[DCI.up.pos,-c(1:4)] <- dataup[DCI.up.pos,-c(1:4)]* FC %>% as.integer()
    }
    ## Multiply FC to cell_type1 for down-regulated interaction
    datadown =  cell_type1.scHiC_table  
    for(j in 1:length(down_region1)){
      DCI.down.pos = which(datadown$region1 == down_region1[j] &
                             datadown$region2 == down_region2[j])
      datadown[DCI.down.pos,-c(1:4)] <- datadown[DCI.down.pos,-c(1:4)]* FC %>% as.integer()
    }
    cc.cells_list = list(datadown, dataup)
  } ## for else if { 
  
  return(cc.cells_list)
} ## for function(){







# MICE imputation with outlier removed ----


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


### 2l.pmm ----
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






### 2l.pmm + neighbor pooling----

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









############ Function to apply MICE impute with outlier removed option 

mice_impute <- function(data_input, n_imputation = n_imputation, outlier.rm = outlier.rm){
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




### Rf ----
## Function for MICE imputation with outlier.rm for all Distance data of each cell type
MICE_impute.rf.outrm.schic <-  function(scHiC.table, n_imputation = 5, outlier.rm = TRUE){
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
    require('CALIBERrfimpute')
    
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
      meth['IF'] <- "rf" 
      #set up predictor matrix
      pred = ini$pred
      pred[,'Single_cell'] <- 0
      
      ## Imputation
      #simulate for n time of multiple imputations, with 5 iteration
      imp <- suppressWarnings(mice(data_input, method = meth, predictorMatrix = pred, print = FALSE,maxit = 1, set.seed(123), m=n_imputation,  ntree = 3) )
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












# Non-parametric Tests ----

# scHiC.table.group1 = ODCs1.table.1mb; scHiC.table.group2 = ODCs2.table.1mb; ith.cell = 5
x.vector <- function(scHiC.table.group1, scHiC.table.group2, ith.cell){
  ## input - scHiC.table.group1
  ##       - scHiC.table.group2
  ##       - ith.cell = i position (row th in scHiC.table) of cell
  ## output - matrix of 2 IF vector (gr1, gr2)
  
  # IF vector of cell ith in Group 1
  x1.coordinations <- scHiC.table.group1[ith.cell, c('region1', 'region2')]
  x1.vector <- scHiC.table.group1[ith.cell, -c(1:4)]
  # IF vector of cell ith in Group 2
  x2.coordinations.pos <- which(scHiC.table.group2$region1 == x1.coordinations$region1[1] & 
                                  scHiC.table.group2$region2 == x1.coordinations$region2[1]  )
  x2.vector <- scHiC.table.group2[x2.coordinations.pos, -c(1:4)]
  #Append vectors from 2 cell types
  x.vector <- append( as.vector(t(x1.vector)), as.vector(t(x2.vector)) )
  # Design matrix
  ref.group <- c(rep('Group1',ncol(x1.vector)), rep('Group2',ncol(x2.vector)) )
  design.matrix <- data.frame(ref.group, x.vector)
  design.matrix[ ,2] <- as.numeric(design.matrix[ ,2])
  # report
  report <- data.frame(chr = unique(scHiC.table.group1$chr), start1 = x1.coordinations$region1,
                       start2 = x1.coordinations$region2, design.matrix)
  return(report)
}


## Wilcoxon Rank test
WR.scHiC.table <- function(scHiC.table.group1, scHiC.table.group2){
  n_cell = nrow(scHiC.table.group1)
  report = NULL
  for ( i in 1:n_cell){
    # Extract design matrix
    m <- x.vector(scHiC.table.group1 = scHiC.table.group1, scHiC.table.group2 = scHiC.table.group2, ith.cell = i)
    # KS Test
    test <- wilcox.test(x.vector ~ ref.group, data = m)
    # Report
    report.cell <- data.frame(unique(m[, -c(4,5)]), wr_p.value = test$p.value, wr_p.adj = NA)
    report <- rbind(report, report.cell)
  }
  report$wr_p.adj <- p.adjust(report$wr_p.value,method = 'fdr' )
  return(report)
}


## KS test
KS.scHiC.table <- function(scHiC.table.group1, scHiC.table.group2){
  n_cell = nrow(scHiC.table.group1)
  report = NULL
  for ( i in 1:n_cell){
    # Extract design matrix
    m <- x.vector(scHiC.table.group1 = scHiC.table.group1, scHiC.table.group2 = scHiC.table.group2, ith.cell = i)
    # KS Test
    test <- ks.test(m$x.vector[m$ref.group == 'Group1'], m$x.vector[m$ref.group == 'Group2'])
    # Report
    report.cell <- data.frame(unique(m[, -c(4,5)]), ks_p.value = test$p.value, ks_p.adj = NA)
    report <- rbind(report, report.cell)
  }
  report$ks_p.adj <- p.adjust(report$ks_p.value,method = 'fdr' )
  return(report)
}


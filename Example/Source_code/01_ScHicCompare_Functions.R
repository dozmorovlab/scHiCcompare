##################### scHicCompare_Functions ####################


## Extract data of single cell -----

### Funtion import data Lee(2019)
### Create function to Extract name of dataset in the file path ----

sc_txt <- function(file_path, cell, position_dataset, type='txt'){
  #list name of selected dataset in the file path
  data_names <- list.files(path = file_path  ,full.names=TRUE,recursive=TRUE) [position_dataset] 
  
  if (type == 'txt'){
    # To write all dataset - txt file.
    datasets <- list()
    for(i in 1: length(position_dataset)){
      dataset <- read.delim(data_names[i])
      # Assign a name to each dataset
      dataset_name <- paste0(cell,"_", i)
      ## assign each dataset into global environment
      assign(dataset_name, dataset, envir = .GlobalEnv) 
    }
  } else if (type == 'cool'){
    # To write all dataset - cool file.
    library(HiCcompare)
    datasets <- list()
    for(i in 1:length(position_dataset)){
      dataset <- cooler2bedpe(data_names[i])
      # Assign a name to each dataset
      dataset_name <- paste0(cell,"_", i)
      ## assign each dataset into global environment
      assign(dataset_name, dataset, envir = .GlobalEnv) 
    }
  }
  return(datasets)
}


sc_list <- function(file_path, cell, position_dataset, type='txt'){
  #list name of selected dataset in the file path
  data_names <- list.files(path = file_path  ,full.names=TRUE,recursive=TRUE) [position_dataset] 
  
  if (type == 'txt'){
    # To write all dataset - txt file.
    datasets <- list()
    datasets <- lapply(data_names, read.delim)
  } else if (type == 'cool'){
    # To write all dataset - cool file.
    library(HiCcompare)
    datasets <- list()
    datasets <- lapply(data_names, cooler2bedpe)
  }
  return(datasets)
}


### Create function to Extract name of dataset in the datalist ----
sc.list_txt <- function(cell, datalist, position_dataset){
  for(i in 1: length(position_dataset) ){
    dataset = datalist[[position_dataset[i]]]
    # Assign a name to each dataset
    dataset_name <- paste0(cell,"_", i)
    ## assign each dataset into global environment
    assign(dataset_name, dataset, envir = .GlobalEnv) 
  }
}




## Transform dataset into sparse and full format ----

### Sparse - Function transform data into sparse format ----
#### Create function to choose specific chromosome and transform into sparse 
chr_sparse <- function(cell_type, n , select_chromosome, type='txt' ){
  if (type == 'txt'){
    options(scipen = 999)
    for(i in 1:n){
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
      # save data with new name
      dataset_name <- paste0(cell_type,"_", i, "_",select_chromosome)
      assign(dataset_name, data, envir = .GlobalEnv) ## assign each dataset into global environment
    }
  } else if (type == 'cool'){
    options(scipen = 999)
    for(i in 1:n){
      #input dataset
      dataset_name <- paste0(cell_type,"_", i) #dataset name
      data <- get(dataset_name) #get dataset have name above
      # Select chromosome
      data <- data$cis[[select_chromosome]] # by using information in col1
      # Transform dataset into sparse
      data <- data[,c('start1', 'start2','IF')]
      names(data) <- c('region1', 'region2','IF')
      # Remove region1 and region2 with 0 values in any column
      data <- data[!data[,1]==0,]
      data <- data[!data[,2]==0,]
      # save data with new name
      dataset_name <- paste0(cell_type,"_", i, "_",select_chromosome)
      assign(dataset_name, data, envir = .GlobalEnv) ## assign each dataset into global environment
    }
  }
}


### Full - Function transform data into full format ----
### Create function to transform into sparse format

full_format <- function(cell_type, n , select_chromosome){
  for(i in 1:n){
    #input dataset
    dataset_name <-paste0(cell_type,"_", i, "_",select_chromosome) #dataset name
    data <- get(dataset_name) #get dataset have name above
    # Transform dataset into full
    library(HiCcompare)
    data_full <-  sparse2full(data) 
    # save data with new name
    datafull_name <- paste0(cell_type,"_", i, "_",select_chromosome,'.full' )
    assign(datafull_name, data_full, envir = .GlobalEnv) ## assign each dataset into global environment
  }
}




## Pseudo bulk Hi-C matrix ----

### Pseudo Bulk HiC - Funtion transform to pseudo bulk HiC ----
### Create function to transform into speudo bulk HiC

pseudo_bulkHiC <- function(cell_type, n , select_chromosome, impute = 'full'){
  data_list = list()
  for(i in 1:n){
    #input dataset
    dataset_name <-paste0(cell_type,"_", i, "_",select_chromosome,'.', impute) #dataset name
    data <- get(dataset_name) #get dataset have name above
    # create a list of all dataset
     data_list <- append(data_list, list(data))
  }
  bulkHiC <- Reduce('+', data_list)
  return(bulkHiC)
}


### Check Sparse odd matrix ----
### Get all possible regions for the resolution
number.region <- function( n, cell_type, chromosome){
  report=NULL
  for (i in 1:n) {
  #input dataset
    dataset_name <- paste0(cell_type,'_', i, '_',chromosome) #dataset name
    data <- get(dataset_name) #get dataset have name above
  #check number of region-- 
    report <- c(as.vector(data[,1]), as.vector(data[,2]), report)
  }
  report <- unique(report)
  return(report)
} #output: all possible regions from all single cell


### Check which ODC matrix have different dimension, since ODC have some single cell with missing regions
matrix_odd_dim <- function(reference, n, cell_type, chromosome){
  # Get the dimensions of the first matrix
  # reference_dim <- dim(ODC_1_chr22)
  # Iterate over the datasets and identify matrices with different dimensions
  differing_matrices <- list()
  for (i in 1:n) {
   #input dataset
    dataset_name <-paste0(cell_type,'_', i, '_',chromosome) #dataset name
    data <- get(dataset_name) #get dataset have name above
  #check region
    reference_region <- reference
    data_region <- unique(c(data[,1], data[,2]))
   if (!identical(reference_region, data_region)) {
      differing_matrices <- c(differing_matrices, i)
    }
  }
  return(differing_matrices)
} #output = order of the differing matrix


correct_odd_sparse <- function(cell_type, select_chromosome, n_sc ){
  ## Input: cell_type; chromosome, etc 'chr22'; n_sc = number of single cells
  ## Output: corrected sparse single cell data (directly input to environment)
  
  #Check any sparse full matrix that has different dimension than the rest
  ## reference dataset 
  ref_number_region <-  number.region(n=n_sc, cell_type = cell_type, chromosome = select_chromosome)
  differing_matrices = matrix_odd_dim(reference =  ref_number_region, n =n_sc, 
                                      cell_type = cell_type, chromosome =  select_chromosome)
  ## Correct adds matrix
  if (length(differing_matrices) != 0){
    for(i in 1: length(differing_matrices)){
      data.sparse_name <- paste0(cell_type, '_',differing_matrices[[i]], '_', select_chromosome)
      data.sparse <- get(data.sparse_name)
      reference_region <- ref_number_region
      data_region <- unique(c(data.sparse[,1], data.sparse[,2]))
      missing_region <-  reference_region[!reference_region %in% data_region]
      # add 0 into missing region in the sparse matrix
      col1 <- c(data.sparse[,1], missing_region) 
      col2<- c(data.sparse[,2], missing_region) 
      col3 <- c(data.sparse[,3],rep(0,length(missing_region)))
      new.data <- data.frame(region1 = col1, region2 = col2, IF=col3)
      # save sparse data with new added regions
      new.data_name <- paste0(cell_type,'_',differing_matrices[[i]], '_', select_chromosome)
      assign(new.data_name, new.data, envir = .GlobalEnv) 
    }
  }
}




## Number of interaction in the transformed pseudo matrix ----

### Function identify \# IF vs \# single cell matrix ----
IF.vs.sc <- function(matrix_name, n){
  pseudo <- get(matrix_name) #get matrix
  n_if <- length(which(pseudo!=0))
  report <- data.frame(matrix_name,
                       n_sc = n,
                       n_IF = n_if)
  return(report)
}


### Function identify sum of IF vs # single cell matrix ----
sumIF.vs.sc <- function(matrix_name, n){
  pseudo <- get(matrix_name) #get matrix
  sum_if <- sum(pseudo)
  report <- data.frame(matrix_name,
                       n_sc = n,
                       sum_IF = sum_if)
  return(report)
}



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



## Find optimal A ----

### Randomize_IFs - function to add noise to IFs of one matrix to create a similar 2 matrix----
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


### function to find best A
best_A <- function (hic.table, SD = 2, numChanges = 35, FC = 3, alpha = 0.05, 
                    Plot = FALSE) 
{
  if (is(hic.table, "list")) {
    stop("Enter a single hic.table object, not a list of hic.tables")
  }
  new.table <- randomize_IFs(hic.table, SD)
  new.table <- new.table[abs(M) < 2, ]
  sample_space <- 1:nrow(new.table)
  tmp_A <- (new.table$IF1 + new.table$IF2)/2
  low_A <- which(tmp_A < quantile(tmp_A, 0.1))
  changes <- sample(sample_space[-low_A], numChanges)
  meanIF <- ((new.table[changes, ]$IF1 + new.table[changes, 
  ]$IF2)/2) %>% round() %>% as.integer()
  suppressWarnings(new.table[changes, `:=`(IF1, meanIF)])
  suppressWarnings(new.table[changes, `:=`(IF2, meanIF)])
  midpoint <- floor(numChanges/2)
  newIF1 <- new.table[changes[1:midpoint], ]$IF1 * FC %>% as.integer()
  newIF2 <- new.table[changes[(midpoint + 1):numChanges], ]$IF2 * 
    FC %>% as.integer()
  new.table[changes[1:midpoint], `:=`(IF1, newIF1)]
  new.table[changes[(midpoint + 1):numChanges], `:=`(IF2, newIF2)]
  new.table = new.table[, `:=`(M, log2(IF2/IF1))]
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
    TP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 
                   1)
    FP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 
                   0)
    FN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 
                   1)
    TN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 
                   0)
  }
  MCC <- ((TP * TN) - (FP * FN))/(sqrt((TP + FP)) * sqrt((TP + 
                                                            FN)) * sqrt((TN + FP)) * sqrt((TN + FN)))
  FPR <- FP/(FP + TP)
  FNR <- FN/(FN + TN)
  TPR <- TP/(TP + FP)
  
  #Best A in term of MCC, TPR, FPR
  MCC.A = which(MCC == max(MCC))
  TPR.A = which(TPR == max(TPR))
  FPR.A = which(FPR == min(FPR))
  intersect.MCC_TPR.A <- intersect(MCC.A,TPR.A)
  best_A <- intersect(intersect.MCC_TPR.A, FPR.A)[1]
  return(best_A)
}




## Vector of each interacting bin on all single cells matrices ----

### function to extract values at an assigned cell in all single cell matrices
vector_1bin  <- function(cell_type, matrix_position , select_chromosome,
                         row, column){
  options(scipen = 999)
  bin_value = rep(NA, length(matrix_position))
  for(i in 1:length(matrix_position)){
    #input matrix
    sc_matrix_name <- paste0(cell_type,"_", i, "_",select_chromosome,'.full' )
    data <- get(sc_matrix_name) #get dataset have name above
    # extract value of a cell
    bin_value[i] = data[row, column]
  }
  return(bin_value)
}

### function to extract values at all cells in all single cell matrices
vector_bins <- function(cell_type, matrix_position , select_chromosome){
  #input matrix
  sc_matrix_name <- paste0(cell_type,"_", matrix_position[1], "_",select_chromosome,'.full' )
  data <- get(sc_matrix_name) #get dataset have name above
  
  # coordination of bins in the triangle of spseudo bulk hi-c
  grid1 <- data.frame(X1 = 1:nrow(data), X2=1:nrow(data)) # diagonal line
  grid2 <- data.frame(t(combn(1:nrow(data),2))) # off diagnoal line value
  grid <- rbind(grid1, grid2)
  
  grid.row <- grid[,1]
  grid.col <- grid[,2]
  # all vector 
  for(i in 1:length(grid.row)){
    vector = vector_1bin(cell_type, matrix_position , select_chromosome,
                         grid.row[i], grid.col[i])
    vector = data.frame(vector)
    ## assign each vector into global
    dataset_name <- paste0(cell_type,"_(",grid.row[i] ,
                           ",",grid.col[i], ')_',select_chromosome)
    assign(dataset_name, vector, envir = .GlobalEnv) 
  }
}




# Vector of bins on each distance (D) in half of speudo bulk matrices ----

bin_distance <- function(cell_type, matrix_position , select_chromosome){ #report for bins information and distance
  # get reference matrix for dimension
  sc_matrix_name <- paste0(cell_type,"_",matrix_position[1] ,
                           "_",select_chromosome,'.full' )
  data <- get(sc_matrix_name) #get dataset have name above
  
  # coordination of bins in the triangle of spseudo bulk hi-c
  grid1 <- data.frame(X1 = 1:nrow(data), X2=1:nrow(data)) # diagonal line
  grid2 <- data.frame(t(combn(1:nrow(data),2))) # off diagnoal line value
  grid <- rbind(grid1, grid2)
  
  # set up the report format
  n_cell = nrow(grid)
  report = data.frame(bin = rep(NA, n_cell), region1 = rep(NA, n_cell), 
                      region2 =rep(NA, n_cell), D = rep(NA, n_cell) )
  # report
  for(i in 1:n_cell){
    
    grid.row <- grid[,1]
    grid.col <- grid[,2]
    # get vector
    dataset_name <- paste0(cell_type,"_(",grid.row[i] ,
                           ",",grid.col[i], ')_','chr22')
    vector <- get(dataset_name); vector <- as.numeric(vector$vector)
    # derive report
    report[i,'bin'] = dataset_name
    report[i,'region1'] =  rownames(ODC_1_chr22.full)[grid.row[i]]
    report[i,'region2'] = colnames(ODC_1_chr22.full)[grid.col[i]] 
    report[i,'D'] = abs(grid.row[i] - grid.col[i])
  }
  return(report)
}



vector_distance <- function(cell_type, matrix_position , select_chromosome){
  D.report <- bin_distance(cell_type, matrix_position , select_chromosome)
  D <- unique(D.report$D)
  for(i in 1:length(D)){
    #data at each cell names correspnding to each distance
    data_names <- D.report$bin[D.report$D == D[i]] 
    vector=NULL
    for(j in 1: length(data_names)){
      data <- get(data_names[j])
      data = as.numeric(data$vector)
      vector = c(vector, data)
    }
    ## assign each vector into global
    vector <- data.frame(vector)
    dataset_name <- paste0(cell_type,"_D",D[i],'_',select_chromosome)
    assign(dataset_name, vector, envir = .GlobalEnv)
  }
}

#vector_distance(cell_type='ODC', matrix_position=1:100 , select_chromosome='chr22')




## Control Change ----

## Function to Create hic.table include changes in IF and columns recording the changes 
groundtruth_control<- function(hic_table, DCI_prop = 0.1, D_interval=0:10,seed =123,
                               p='p.adj', sign_DCI='up_down'){
  ## hic.table = table result after run FPC test
  ## DCI_prop = proportion of true differential chromatin interaction (DCI)
  ## D_interval = Chosen Distance interval to make changes DCI (non-zero value)
  ## seed = random seed for random pick Distance changes postion
  ##p-value = the true/false positive result from hic.table's test based on which before or after FDR ajusted p-value (p.value/p.adjusted)
  ## sign_interaction = for sign of (IF2-IF1), IF2>IF1 -> 'up'; IF2<IF1 -> 'down'; mix of IF2>IF1 and IF2<IF1 -> 'up_down'
  
  
  ## True positive region position based on p-value from FPC test
  set.seed(seed)
  new.table = hic_table
  true_region.pos = which(new.table[[p]] >0.05)
  
  ## Extract bins position within D_interval
  region_Dinterval.pos = which(new.table$D %in% D_interval)
  
  ## Extract bin position with M close from 0, sd=1 -> -1<if<1
  region_M0.pos = which(new.table$adj.M>(-1) & new.table$adj.M< 1)
  
  ## Extract high average expression bin position
  tmp_A <- (new.table$adj.IF1 + new.table$adj.IF2)/2
  high_A <- which(tmp_A > quantile(tmp_A, 0.1)) #remove sample space where A<10th percentile
  
  ## Selected all bin positions
  select_reg.pos <- intersect(intersect(true_region.pos, high_A), region_Dinterval.pos)
  new.table$select_DCI_loci <- 0
  new.table$select_DCI_loci[select_reg.pos] <- 1
  
  ## Select Number of sample to assign changes 
  changeSample <- round(nrow(new.table)*DCI_prop) #number of DCI
  # library(purrr)
  # distribute_D <- rdunif(changeSample,b=D_interval[length(D_interval)],a=D_interval[1]) #selected frequency DCI changes happen at each distance
  
  ## Assign position of change sample (assume changes equally distributed among Distance)
  #distribute_D.report = data.frame(table(distribute_D))
  changes.pos <- NULL
  for(i in 1:length(D_interval)){
    ## bins positions at each distance
    bins.D.pos <- which(new.table$D == D_interval[i] & new.table$select_DCI_loci ==1)
    ## same proportion of changed DCI at each Distance
    changeSample_D <- round(length(bins.D.pos)*DCI_prop)
    changes.D.pos <- sample(bins.D.pos,changeSample_D)
    changes.pos <- c(changes.pos, changes.D.pos)
  }

  
  ## Create ground truth report
  truth <- rep(0, nrow(new.table)) 
  truth[changes.pos] <- 1
  new.table$truth <- truth
  
  if (sign_DCI=='up'){
    ## Create sign of regulated interaction (up vs down) truth report
    sign_interaction <- rep(0, nrow(new.table)) 
    sign_interaction[changes.pos] <- 'up'
    new.table$sign_interaction <- sign_interaction
    ## reconstruct new table
    new.table <- new.table[,-c('mc','A','Z','p.value','p.adj')]
    
    
  } else if (sign_DCI=='down'){
    ## Create sign of regulated interaction (up vs down) truth report
    sign_interaction <- rep(0, nrow(new.table)) 
    sign_interaction[changes.pos] <- 'down'
    new.table$sign_interaction <- sign_interaction
    ## reconstruct new table
    new.table <- new.table[,-c('mc','A','Z','p.value','p.adj')]
    
    
  } else if (sign_DCI=='up_down'){
    ## Create sign of regulated interaction (up vs down) truth report
    sign_interaction <- rep(0, nrow(new.table)) 
    changes.pos2 <- sample(changes.pos) #shuffle elements of changes position
    changeSample <- length(changes.pos2) #total DCI length
    midpoint <- floor(changeSample/2)
    sign_interaction[changes.pos2[1:midpoint]] <- 'up'
    sign_interaction[changes.pos2[(midpoint + 1):changeSample]] <- 'down'
    new.table$sign_interaction <- sign_interaction
    ## reconstruct new table
    new.table <- new.table[,-c('mc','A','Z','p.value','p.adj')]
  }
  
  new.table <- new.table[!new.table$D == 0, ]
  
  return(new.table)
}



## Function to make change DCI on individual single cell
## Function to make change DCI on individual single cell
sc_controlDCI <- function(groundtruth_control.report, cell_type1, cell_type2,
                          select_chromosome, sc_position =1:100, FC=3, sign_interaction='up'){
  ## extract changes DCI position from ground truth report
  rt.report = groundtruth_control.report
  changes.pos = rt.report[rt.report[['truth']]==1,]
  region1 = changes.pos[['start1']]; region2 = changes.pos[['start2']] #extract which pair of regions will make change
  
  if(sign_interaction=='up') {
    for(i in 1:length(sc_position)){
      ## extract single cell(sc) full matrix at the sc_position
      datafull_name <- paste0(cell_type2,"_", i, "_",select_chromosome,'.full' )
      data <- get(datafull_name)
      for(j in 1:length(region1)){
        rg1 = which(rownames(data) == region1[j])
        rg2 = which(rownames(data) == region2[j])
        if(rg1==rg2){ # rg1 = rg2 -> diagnal line in the full matrix
          data[rg1,rg2] <- data[rg1,rg2] * FC %>% as.integer()
        } else { # other pair of region, matrix is symmetric
          data[rg1,rg2] <- data[rg1,rg2] * FC %>% as.integer()
          data[rg2,rg1] <- data[rg2,rg1] * FC %>% as.integer()
        }
      }
      ## Assign new simulation data matrix
      new_name <- paste0(cell_type2,".simCC_", i, "_",select_chromosome,'.full' )
      assign(new_name, data , envir = .GlobalEnv)
    }
    
    
  } else if(sign_interaction=='down'){
    for(i in 1:length(sc_position)){
      ## extract single cell(sc) full matrix at the sc_position
      datafull_name <- paste0(cell_type1,"_", i, "_",select_chromosome,'.full' )
      data <- get(datafull_name)
      for(j in 1:length(region1)){
        rg1 = which(rownames(data) == region1[j])
        rg2 = which(rownames(data) == region2[j])
        if(rg1==rg2){ # rg1 = rg2 -> diagnal line in the full matrix
          data[rg1,rg2] <- data[rg1,rg2] * FC %>% as.integer()
        } else { # other pair of region, matrix is symmetric
          data[rg1,rg2] <- data[rg1,rg2] * FC %>% as.integer()
          data[rg2,rg1] <- data[rg2,rg1] * FC %>% as.integer()
        }
      }
      ## Assign new simulation data matrix
      new_name <- paste0(cell_type1,".simCC_", i, "_",select_chromosome,'.full' )
      assign(new_name, data, envir = .GlobalEnv)
    }
    
    
  } else if(sign_interaction=='up_down'){
    up_region1 = changes.pos$'start1'[changes.pos$sign_interaction=='up']
    up_region2 = changes.pos$'start2'[changes.pos$sign_interaction=='up']
    down_region1 = changes.pos$'start1'[changes.pos$sign_interaction=='down']
    down_region2 = changes.pos$'start2'[changes.pos$sign_interaction=='down']
    
    ## Multiply FC to cell_type2 for up-regulated interaction
    for(i in 1:length(sc_position)){
      ## extract single cell(sc) full matrix at the sc_position
      dataup_name <- paste0(cell_type2,"_", i, "_",select_chromosome,'.full' )
      dataup <- get(dataup_name)
      for(j in 1:length(up_region1)){
        up.rg1 = which(rownames(dataup) == up_region1[j])
        up.rg2 = which(rownames(dataup) == up_region2[j])
        if(up.rg1==up.rg2){ # rg1 = rg2 -> diagnal line in the full matrix
          dataup[up.rg1, up.rg2] <- dataup[up.rg1, up.rg2] * FC %>% as.integer()
        } else { # other pair of region, matrix is symmetric
          dataup[up.rg1, up.rg2] <- dataup[up.rg1, up.rg2] * FC %>% as.integer()
          dataup[up.rg2,up.rg1] <- dataup[up.rg2,up.rg1] * FC %>% as.integer()
        }
      }
      ## Assign new simulation data matrix
      new_name.up <- paste0(cell_type2,".simCC_", i, "_",select_chromosome,'.full' )
      assign(new_name.up, dataup, envir = .GlobalEnv)
    } 
    
    ## Multiply FC to cell_type1 for down-regulated interaction
    for(i in 1:length(sc_position)){
      ## extract single cell(sc) full matrix at the sc_position
      datadown_name <- paste0(cell_type1,"_", i, "_",select_chromosome,'.full' )
      datadown <- get(datadown_name)
      for(j in 1:length(down_region1)){
        down.rg1 = which(rownames(datadown) == down_region1[j])
        down.rg2 = which(rownames(datadown) == down_region2[j])
        if(down.rg1==down.rg2){ # rg1 = rg2 -> diagnal line in the full matrix
          datadown[down.rg1, down.rg2] <- datadown[down.rg1,down.rg2] * FC %>% as.integer()
        } else { # other pair of region, matrix is symmetric
          datadown[down.rg1, down.rg2] <- datadown[down.rg1, down.rg2] * FC %>% as.integer()
          datadown[down.rg2, down.rg1] <- datadown[down.rg2, down.rg1] * FC %>% as.integer()
        }
      }
      ## Assign new simulation data matrix
      new_name.down <- paste0(cell_type1,".simCC_", i, "_",select_chromosome,'.full' )
      assign(new_name.down, datadown, envir = .GlobalEnv)
    } ## for (i in 1: ___){}
    
  } ## for else if { 
} ## for function(){



### Function to Calculate measurement metric
metric_controlDCI <-  function(hic_table, groundtruth, p='p.value'){
  tmp.table<- merge(hic_table, groundtruth, by = c('start1','start2'))
  ## MCC, FPR, FNR, TPR - double check the formula
  TP <- sum(tmp.table[[p]] < 0.05 & tmp.table$truth == 1)
  FP <- sum(tmp.table[[p]] < 0.05 & tmp.table$truth == 0)
  FN <- sum(tmp.table[[p]] >= 0.05 & tmp.table$truth == 1)
  TN <- sum(tmp.table[[p]] >= 0.05 & tmp.table$truth == 0)
  
  MCC <- ((TP * TN) - (FP * FN))/(sqrt((TP + FP)) * sqrt((TP +  FN)) * sqrt((TN + FP)) * sqrt((TN + FN)))
  Sensitivity <- TP/(TP + FN) # (recall) true positive predictions among all actual positive 
  Specificity <- TN/(TN + FP) #true negative predictions among all actual negative
  Precision <- TP/(TP+FP) #true positive predictions among all positive predictions
  F1 <- (2* Precision * Sensitivity)/ (Precision + Sensitivity)
  return(data.frame(MCC,Sensitivity, Specificity,Precision, F1))
}



## Other tests compares bin vector between groups ----

### Design IF vector for 2 group
x.vector <- function(cell_type1, cell_type2, select_chromosome,n_Distance,
                     scale = FALSE, round = FALSE){
  # coordination of bins in the triangle of spseudo bulk hi-c
  grid1 <- data.frame(X1 = 1:n_Distance, X2=1:n_Distance) # diagonal line
  grid2 <- data.frame(t(combn(1:n_Distance,2))) # off diagnoal line value
  grid <- rbind(grid1, grid2)
  
  grid.row <- grid[,1]
  grid.col <- grid[,2]
  #Extract bins of cell type 1
  for(i in 1:length(grid.row)){
    x1.name <- paste0(cell_type1,"_(",grid.row[i] ,
                      ",",grid.col[i], ')_',select_chromosome)
    x1.vector <- get(x1.name)
    #Extract bins of cell type 2
    x2.name <- paste0(cell_type2,"_(",grid.row[i] ,
                      ",",grid.col[i], ')_',select_chromosome)
    x2.vector <- get(x2.name)
    if(scale == TRUE){
      x1 = as.numeric(x1.vector$vector); x2 = as.numeric(x2.vector$vector)
      x1.vector <- round((x1 - mean(x1))/sd(x1))
      x2.vector <- round((x2 - mean(x2))/sd(x2))
    } else if(round == TRUE){
      x1 = as.numeric(x1.vector$vector); x2 = as.numeric(x2.vector$vector)
      x1.vector <- round(x1)
      x2.vector <- round(x2)
    } else {
      x1.vector <- as.numeric(x1.vector$vector)
      x2.vector <- as.numeric(x2.vector$vector)
    }
    #Append vectors from 2 cell types
    x.vector <- append(x1.vector, x2.vector)
    x.name <- paste0(cell_type1, '_',cell_type2,"_(",grid.row[i] ,
                     ",",grid.col[i], ')_',select_chromosome)
    assign(x.name, x.vector, envir = .GlobalEnv) 
  }
}


## Function Negative Binomial for each bin vector between 2 groups
nb.sc <- function(cell_type1, cell_type2, select_chromosome,n_Distance,
                  start_region, resolution){
  library(MASS)
  ## extract bins vector at each cell matrix
  x.vector(cell_type1, cell_type2, select_chromosome, n_Distance,  scale = FALSE, round = TRUE)
  
  ## ref data to extract length of x.vector
  ref.name <- paste0(cell_type1, '_',cell_type2,"_(1,1)_",select_chromosome)
  ref.vector <- get(ref.name)
  n <- length(ref.vector)
  ref.group <- c(rep(cell_type1, n/2), rep(cell_type2,n/2))
  
  ## report table
  ##### number of vector with pair of regions
  grid1 <- data.frame(X1 = 1:n_Distance, X2=1:n_Distance) # diagonal line
  grid2 <- data.frame(t(combn(1:n_Distance,2))) # off diagnoal line value
  grid <- rbind(grid1, grid2)
  grid.row <- grid[,1]
  grid.col <- grid[,2]
  l <- nrow(grid) ## number of x,vectors
  ##### create data frame for report table
  report <- data.frame(chr = rep(select_chromosome,l), start1 = rep(NA, l),
                       start2 = rep(NA,l), D = rep(NA,l), nb_p.value= rep(NA, l),
                       nb_p.adj= rep(NA, l))
  
  ## run the test an extract p.value to table
  for(i in 1:length(grid.row)){
    x.name <- paste0(cell_type1, '_',cell_type2,"_(",grid.row[i],
                     ",",grid.col[i], ')_',select_chromosome)
    x.vector <- get(x.name)
    ## run negative binomial test
    if (var(x.vector) == 0) { 
      report$nb_p.value[i] <- 1
    } else {
      nb.model <- glm.nb(x.vector ~ ref.group)
      p.value <- coef(summary(nb.model))[2,4] #extract p.value
      report$nb_p.value[i] <- p.value
    }
    ## assign result into report table
    report$start1[i] <- start_region + (grid.row[i] -1)*resolution
    report$start2[i] <- start_region+ (grid.col[i] -1)*resolution
    report$D[i] <- abs(grid.row[i] - grid.col[i])
  }
  report$nb_p.adj = p.adjust(report$nb_p.value,method = 'fdr' )
  return(report)
}



## Function t-test for each bin vector between 2 groups
ttest.sc <- function(cell_type1, cell_type2, select_chromosome,n_Distance,
                     start_region, resolution){
  
  ## extract bins vector at each cell matrix
  x.vector(cell_type1, cell_type2, select_chromosome, n_Distance)
  
  ## ref data to extract length of x.vector
  ref.name <- paste0(cell_type1, '_',cell_type2,"_(1,1)_",select_chromosome)
  ref.vector <- get(ref.name)
  n <- length(ref.vector)
  ref.group <- c(rep(cell_type1, n/2), rep(cell_type2,n/2))
  
  ## report table
  ##### number of vector with pair of regions
  grid1 <- data.frame(X1 = 1:n_Distance, X2=1:n_Distance) # diagonal line
  grid2 <- data.frame(t(combn(1:n_Distance,2))) # off diagnoal line value
  grid <- rbind(grid1, grid2)
  grid.row <- grid[,1]
  grid.col <- grid[,2]
  l <- nrow(grid) ## number of x,vectors
  ##### create data frame for report table
  report <- data.frame(chr = rep(select_chromosome,l), start1 = rep(NA, l),
                       start2 = rep(NA,l), D = rep(NA,l), tt_p.value= rep(NA, l),
                       tt_p.adj= rep(NA, l))
  
  ## run the test an extract p.value to table
  for(i in 1:length(grid.row)){
    x.name <- paste0(cell_type1, '_',cell_type2,"_(",grid.row[i],
                     ",",grid.col[i], ')_',select_chromosome)
    x.vector <- get(x.name)
    ## run t.test
    if (var(x.vector) == 0 | length(unique(x.vector[1:(n/2)])) == 1 | length(unique(x.vector[(n/2+1):n])) == 1 ) { 
      report$tt_p.value[i] <- 1
    } else {
      ttest <- t.test(x.vector ~ ref.group)
      p.value <- ttest$p.value #extract p.value
      report$tt_p.value[i] <- p.value
    }
    ## assign result into report table
    report$start1[i] <- start_region + (grid.row[i] -1)*resolution
    report$start2[i] <- start_region+ (grid.col[i] -1)*resolution
    report$D[i] <- abs(grid.row[i] - grid.col[i])
  }
  report$tt_p.adj = p.adjust(report$tt_p.value,method = 'fdr' )
  return(report)
}







## Function distance standardized t-test for each bin vector between 2 groups
D.standardize <- function(cell_type, select_chromosome,n_Distance){
  # Input: cell_type = name of cell type; select_chromosome = ex. 'chr22'; n_Distance = # of total Distance
  # Output: standardized vector of bin with old name (ex. ODC_(1,1)_chr22)
  
  
  ## Get distance data
  D_mean = D_sd = rep(NA, n_Distance)
  for(d in 1:n_Distance){
    D = 0:(n_Distance-1)
    dataset_name <- paste0(cell_type,"_D",D[d],'_',select_chromosome)
    D_data = get(dataset_name)
  ## Get mean and sd at each Distance
    D_mean[d] = mean(D_data$vector)
    D_sd[d] = sd(D_data$vector)
  }
  
  ## Standardized Cell Bin vector
  ##### number of vector with pair of regions
  grid1 <- data.frame(X1 = 1:n_Distance, X2=1:n_Distance) # diagonal line
  grid2 <- data.frame(t(combn(1:n_Distance,2))) # off diagnoal line value
  grid <- rbind(grid1, grid2)
  grid.row <- grid[,1]
  grid.col <- grid[,2]
  #### Standardize bin vector
  for(g in 1:length(grid.row)){
    x.name <- paste0(cell_type,"_(",grid.row[g],
                     ",",grid.col[g], ')_',select_chromosome)
    x.vector <- get(x.name)
    x.d <- abs(grid.row[g] - grid.col[g]) #get distance index
    x.sd.vector <- (x.vector - D_mean[x.d+1])/D_sd[x.d+1] #z-score tranform
   #### assign each vector into global
    vector <- data.frame(x.sd.vector)
    assign(x.name, vector, envir = .GlobalEnv)  
  }
}


D.standard.ttest.sc <- function(cell_type1, cell_type2, select_chromosome,n_Distance,
                     start_region, resolution){
  
  # Distance standardize bin vector
  D.standardize(cell_type = cell_type1, select_chromosome = select_chromosome ,n_Distance = n_Distance)
  D.standardize(cell_type = cell_type2, select_chromosome = select_chromosome ,n_Distance = n_Distance)
  ## extract bins vector at each cell matrix
  x.vector(cell_type1, cell_type2, select_chromosome, n_Distance)
  
  ## ref data to extract length of x.vector
  ref.name <- paste0(cell_type1, '_',cell_type2,"_(1,1)_",select_chromosome)
  ref.vector <- get(ref.name)
  n <- length(ref.vector)
  ref.group <- c(rep(cell_type1, n/2), rep(cell_type2,n/2))
  
  ## report table
  ##### number of vector with pair of regions
  grid1 <- data.frame(X1 = 1:n_Distance, X2=1:n_Distance) # diagonal line
  grid2 <- data.frame(t(combn(1:n_Distance,2))) # off diagnoal line value
  grid <- rbind(grid1, grid2)
  grid.row <- grid[,1]
  grid.col <- grid[,2]
  l <- nrow(grid) ## number of x,vectors
  ##### create data frame for report table
  report <- data.frame(chr = rep(select_chromosome,l), start1 = rep(NA, l),
                       start2 = rep(NA,l), D = rep(NA,l), tt_p.value= rep(NA, l),
                       tt_p.adj= rep(NA, l))
  
  ## run the test an extract p.value to table
  for(i in 1:length(grid.row)){
    x.name <- paste0(cell_type1, '_',cell_type2,"_(",grid.row[i],
                     ",",grid.col[i], ')_',select_chromosome)
    x.vector <- get(x.name)
    ## run t.test
    if (var(x.vector) == 0 | length(unique(x.vector[1:(n/2)])) == 1 | length(unique(x.vector[(n/2+1):n])) == 1 ) { 
      report$tt_p.value[i] <- 1
    } else {
      ttest <- t.test(x.vector ~ ref.group)
      p.value <- ttest$p.value #extract p.value
      report$tt_p.value[i] <- p.value
    }
    ## assign result into report table
    report$start1[i] <- start_region + (grid.row[i] -1)*resolution
    report$start2[i] <- start_region+ (grid.col[i] -1)*resolution
    report$D[i] <- abs(grid.row[i] - grid.col[i])
  }
  report$tt_p.adj = p.adjust(report$tt_p.value,method = 'fdr' )
  return(report)
}




## Function Wilcoxon Rank test for each bin vector between 2 groups
wr.sc <- function(cell_type1, cell_type2, select_chromosome,n_Distance,
                     start_region, resolution){
  
  ## extract bins vector at each cell matrix
  x.vector(cell_type1, cell_type2, select_chromosome, n_Distance)
  
  ## ref data to extract length of x.vector
  ref.name <- paste0(cell_type1, '_',cell_type2,"_(1,1)_",select_chromosome)
  ref.vector <- get(ref.name)
  n <- length(ref.vector)
  ref.group <- c(rep(cell_type1, n/2), rep(cell_type2,n/2))
  
  ## report table
  ##### number of vector with pair of regions
  grid1 <- data.frame(X1 = 1:n_Distance, X2=1:n_Distance) # diagonal line
  grid2 <- data.frame(t(combn(1:n_Distance,2))) # off diagnoal line value
  grid <- rbind(grid1, grid2)
  grid.row <- grid[,1]
  grid.col <- grid[,2]
  l <- nrow(grid) ## number of x,vectors
  ##### create data frame for report table
  report <- data.frame(chr = rep(select_chromosome,l), start1 = rep(NA, l),
                       start2 = rep(NA,l), D = rep(NA,l), wr_p.value= rep(NA, l),
                       wr_p.adj= rep(NA, l))
  
  ## run the test an extract p.value to table
  for(i in 1:length(grid.row)){
    x.name <- paste0(cell_type1, '_',cell_type2,"_(",grid.row[i],
                     ",",grid.col[i], ')_',select_chromosome)
    x.vector <- get(x.name)
    ## run wr.test
    if (var(x.vector) == 0 | length(unique(x.vector[1:(n/2)])) == 1 | length(unique(x.vector[(n/2+1):n])) == 1 ) { 
      report$wr_p.value[i] <- 1
    } else {
      wrest <- wilcox.test(x.vector ~ ref.group)
      p.value <- wrest$p.value #extract p.value
      report$wr_p.value[i] <- p.value
    }
    ## assign result into report table
    report$start1[i] <- start_region + (grid.row[i] -1)*resolution
    report$start2[i] <- start_region+ (grid.col[i] -1)*resolution
    report$D[i] <- abs(grid.row[i] - grid.col[i])
  }
  report$wr_p.adj = p.adjust(report$wr_p.value,method = 'fdr' )
  return(report)
}





## Function KS test for each bin vector between 2 groups
ks.sc <- function(cell_type1, cell_type2, select_chromosome,n_Distance,
                  start_region, resolution){
  
  ## extract bins vector at each cell matrix
  x.vector(cell_type1, cell_type2, select_chromosome, n_Distance)
  
  ## ref data to extract length of x.vector
  ref.name <- paste0(cell_type1, '_',cell_type2,"_(1,1)_",select_chromosome)
  ref.vector <- get(ref.name)
  n <- length(ref.vector)
  ref.group <- c(rep(cell_type1, n/2), rep(cell_type2,n/2))
  
  ## report table
  ##### number of vector with pair of regions
  grid1 <- data.frame(X1 = 1:n_Distance, X2=1:n_Distance) # diagonal line
  grid2 <- data.frame(t(combn(1:n_Distance,2))) # off diagnoal line value
  grid <- rbind(grid1, grid2)
  grid.row <- grid[,1]
  grid.col <- grid[,2]
  l <- nrow(grid) ## number of x,vectors
  ##### create data frame for report table
  report <- data.frame(chr = rep(select_chromosome,l), start1 = rep(NA, l),
                       start2 = rep(NA,l), D = rep(NA,l), ks_p.value= rep(NA, l),
                       ks_p.adj= rep(NA, l))
  
  ## run the test an extract p.value to table
  for(i in 1:length(grid.row)){
    x.name <- paste0(cell_type1, '_',cell_type2,"_(",grid.row[i],
                     ",",grid.col[i], ')_',select_chromosome)
    x.vector <- get(x.name)
    ## run ks.test
    if (var(x.vector) == 0 | length(unique(x.vector[1:(n/2)])) == 1 | length(unique(x.vector[(n/2+1):n])) == 1 ) { 
      report$ks_p.value[i] <- 1
    } else {
      ksest <- ks.test(x.vector[ref.group == cell_type1], x.vector[ref.group == cell_type2])
      p.value <- ksest$p.value #extract p.value
      report$ks_p.value[i] <- p.value
    }
    ## assign result into report table
    report$start1[i] <- start_region + (grid.row[i] -1)*resolution
    report$start2[i] <- start_region+ (grid.col[i] -1)*resolution
    report$D[i] <- abs(grid.row[i] - grid.col[i])
  }
  report$ks_p.adj = p.adjust(report$ks_p.value,method = 'fdr' )
  return(report)
}




## Mice Imputation ----

### MICE 1 - Algorithm with predictor matrix for each Distance imputed both group test ----

## Function to design data matrix

impute_sc_D <- function(cell_type1, cell_type2, matrix_position, select_chromosome, distance){
  D.report1 = bin_distance(cell_type = cell_type1, matrix_position , select_chromosome)
  D.report2 = bin_distance(cell_type = cell_type2, matrix_position , select_chromosome)
  
  ## Extract name of level 2: bins corresponding to pair of interacting regions
  #### Cell type 1
  data_names1 <- D.report1$bin[D.report1$D == distance] #find all bins vector same in corresponding distance
  bins1_name <- rep(data_names1, each = length(matrix_position))
  #### Cell type 2
  data_names2 <- D.report2$bin[D.report2$D == distance] #find all bins vector same in corresponding distance
  bins2_name <- rep(data_names2, each = length(matrix_position))
  bins_name <- c(bins1_name, bins2_name)
  
  
  ## Extract name of level 1: single cells for each in above level 
  single_cell1 <- rep(matrix_position, length(data_names1))
  single_cell2 <- rep(matrix_position, length(data_names2))
  single_cell <- c(single_cell1,single_cell2 )
  
  ## Extract IF vector at corresponding distance 
  vector_distance(cell_type1, matrix_position, select_chromosome) #### Cell type 1
  vector_name1 <- paste0(cell_type1,"_D",distance,'_',select_chromosome)
  IF1 <- get(vector_name1)
  vector_distance(cell_type2, matrix_position, select_chromosome) #### Cell type 2
  vector_name2 <- paste0(cell_type2,"_D",distance,'_',select_chromosome)
  IF2 <- get(vector_name2)
  IF <- c(IF1$vector, IF2$vector)
  IF[IF == 0] <- NA
  
  ## Extract name of cell_type
  type1 = rep(cell_type1, length(IF1$vector))
  type2 = rep(cell_type2, length(IF2$vector))
  type = c(type1, type2)
  
  ## Combining to data.frame
  matrix <- data.frame(Bins = bins_name, Single_cell = single_cell, IF = IF, Cell_type = type)
  
  return(matrix)
}



## Function for MICE imputation for all Distance data of each cell type
MICE_impute1 <-  function(cell_type1, cell_type2, matrix_position, select_chromosome,
                          n_imputation = 5, n_dis){
  #n_dis = number of distance
  l = n_dis-1
  for( i in 1:n_dis){
    D_pos = 0:l
    data = impute_sc_D(cell_type1 = cell_type1, cell_type2 = cell_type2,
                       matrix_position = 1:100, select_chromosome = select_chromosome,
                       distance =  D_pos[i])
    data_input = data[,-1]
    require(mice)
    require(lattice)
    library(miceadds)
    
    ## Classify variable class
    data_input$Single_cell <- as.numeric(data_input$Single_cell)
    data_input$Cell_type= as.numeric(as.factor(data_input$Cell_type))
    
    ## Classify method and predictor matrix
    # get initial default imputation setting
    set.seed(123)
    ini <- mice(data_input, maxit = 0)
    #set up method
    meth= ini$meth
    meth['IF'] <- "2l.pmm" 
    #set up predictor matrix
    pred = ini$pred
    pred[,'Single_cell'] <- -2 #random intercept
    pred[,'Cell_type'] <- c(0,0,0) #fix effect
    pred[,'IF'] <- c(1,0,1) #level 2 imputation
    
    ## Imputation
    #simulate for n time of multiple imputations, with 5 iteration
    imp <- mice(data_input, method = meth, predictorMatrix = pred, print = FALSE,maxit = 1,
                set.seed(123), m=n_imputation) 
    #returns the long format of all multiple imputation
    imp_data = complete(imp, action = 'long', include = F) 
    
    #if the vector has all if=1, the NA will be 1 too
    if(anyNA(imp_data) == T) {
      imp_data$IF[is.na(imp_data$IF)] <- 1
      agg_new_if = imp_data$IF
    } else{
      #if the vector has all if>1, aggregate mean of all imputed complete data
      agg_new_if = round(aggregate(imp_data[,4] , by = list(imp_data$.id),FUN= mean))
      agg_new_if = agg_new_if$x}
    
    
    ## Assign each vector into global
    pos1 <- which(data$Cell_type == cell_type1) 
    new_data1 <- data.frame(Bins = data$Bins[pos1], Single_cell = data$Single_cell[pos1],
                            IF = agg_new_if[pos1], Cell_type = data$Cell_type[pos1])
    dataset_name1 <- paste0(cell_type1,"_impD",D_pos[i],'_',select_chromosome)
    assign(dataset_name1, new_data1 , envir = .GlobalEnv)
    
    pos2 <- which(data$Cell_type == cell_type2) 
    new_data2 <- data.frame(Bins = data$Bins[pos2], Single_cell = data$Single_cell[pos2],
                            IF = agg_new_if[pos2], Cell_type = data$Cell_type[pos2])
    dataset_name2 <- paste0(cell_type2,"_impD",D_pos[i],'_',select_chromosome)
    assign(dataset_name2, new_data2 , envir = .GlobalEnv)
  }
}  




## Function to Merge these imputed vector into sparse matrices
imp_scKth_maxtrix1 <- function(cell_type , select_chromosome, n_distance,
                              Kth_matrix, start_region, end_region, bin,
                              n_sc){
  l = n_distance-1
  region1 = region2 = IF = NULL
  for(i in 1:length(1:n_distance)){
    # get imputed distance vector
    D = 0:l
    dis_matrix_name <-  paste0(cell_type,"_impD",D[i],'_',select_chromosome)
    dis_data <- get(dis_matrix_name)
    # extract distance in i th single cell
    ## extract IF element in simulated distance vector corresponding to position of Kth single cell matrix
    D_length <- n_distance - D[i]
    IF_position <- seq(Kth_matrix, D_length*n_sc, by = n_sc)
    if_ = dis_data[IF_position,3]
    If_pos = dis_data[IF_position,1]
    ## region1 and region2 at each distance D's IF vector above
    reg1 = seq(start_region, end_region - D[i]*bin, by = bin)
    reg2 = reg1 + D[i]*bin
    
    ## create frame for result
    l_matrix = (n_distance/2) * (1 + n_distance) #length for sc sparse length
    ## Input result to sc_matrix
    region1 = c(region1, reg1)
    region2 = c(region2, reg2)
    IF = c(IF, if_)
  }
  
  ## assign sc sparse format into global
  sc_matrix = data.frame(region1, region2, IF)
  dataset_name <- paste0(cell_type,"_imp_", Kth_matrix, "_",select_chromosome)
  assign(dataset_name,sc_matrix, envir = .GlobalEnv)
  return(sc_matrix)
}  



### MICE 2 - Algorithm with predictor matrix for total bins imputed individual group test with Distane predictor ----

## Function to design data matrix

predictorMatrix_sc <- function(cell_type, matrix_position, select_chromosome){
  ## Report bins location and distance
  D.report = bin_distance(cell_type = 'ODC', matrix_position , select_chromosome)
  
  ## Extract IF vector for each bin (compiling from `matrix_position` sc)
  IF.vector <- NULL
  for(i in 1:nrow(D.report)){
    name <- D.report$bin[i]
    IF <- get(name)
    IF.vector <- c(IF.vector, IF$vector)
  }
  IF.vector[IF.vector == 0] = NA
  
  ## Build the designed matrix
  bin <- rep(D.report$bin,  each= length(matrix_position)) #bin
  sc <- rep(matrix_position,  nrow(D.report)) #repeated single cell
  d <-  rep(D.report$D, each= length(matrix_position)) #distance 
  matrix <- data.frame(Bin = bin, Single_cell = sc, D = d, IF = IF.vector)
  return(matrix)
}




## Function for MICE imputation for all Distance data of each cell type
MICE_impute2 <-  function(cell_type, matrix_position, select_chromosome, n_imputation = 5){
  data = predictorMatrix_sc(cell_type = cell_type, matrix_position = matrix_position,
                            select_chromosome = select_chromosome)
  data_input = data[,-1]
  require(mice)
  require(lattice)
  library(miceadds)
  
  ## Classify variable class
  data_input$Single_cell <- as.numeric(data_input$Single_cell)
  data_input$D= as.numeric(data_input$D)
  
  ## Classify method and predictor matrix
  # get initial default imputation setting
  set.seed(101)
  ini <- mice(data_input, maxit = 0)
  #set up method
  meth= ini$meth
  meth['IF'] <- "2l.pmm" 
  #set up predictor matrix
  pred = ini$pred
  pred[,'Single_cell'] <- -2 #random intercept
  pred[,'D'] <- c(0,0,1) #fix effect
  pred[,'IF'] <- c(0,1,0) #level 2 imputation
  
  ## Imputation
  #simulate for n time of multiple imputations, with 5 iteration
  imp <- mice(data_input, method = meth, predictorMatrix = pred, print = FALSE,maxit = 1,
              set.seed(123), m=n_imputation) 
  #returns the long format of all multiple imputation
  imp_data = complete(imp, action = 'long', include = F) 
  #iaggregate mean of all imputed complete data
  agg_new_if = round(aggregate(imp_data[,'IF'] , by = list(imp_data$.id),FUN= mean))
  agg_new_if = agg_new_if$x
  
  
  ## Assign each vector into global
  new_data <- data.frame(Bin = data$Bin, Single_cell = data$Single_cell,
                         D = data$D, IF = agg_new_if)
  return(new_data)
}  



#### Func to seperate imputed Distance vector 
MICE_D_impute <-  function(cell_type, matrix_position, select_chromosome, 
                           n_imputation = 5, n_distance){
  #impute data in all single cells and distances
  imp_data =  MICE_impute2(cell_type, matrix_position, select_chromosome, n_imputation = 5)
  
  l = n_distance-1 
  region1 = region2 = IF = NULL
  for(i in 1:length(1:n_distance)){
    # get imputed distance vector
    D = 0:(n_distance-1)
    dis_data <-  imp_data[imp_data$D == D[i],] 
    # assign dataset for each imputed distance
    name_dis_data <- paste0(cell_type,"_impD",D[i],'_',select_chromosome)
    assign(name_dis_data, dis_data, envir = .GlobalEnv)
  }
}



## Function to Merge these imputed vector into sparse matrices
imp_scKth_maxtrix2 <- function(cell_type , select_chromosome, n_distance,
                               Kth_matrix, start_region, end_region, bin,
                               n_sc){
  l = n_distance-1
  region1 = region2 = IF = NULL
  for(i in 1:length(1:n_distance)){
    # get imputed distance vector
    D = 0:l
    dis_matrix_name <-  paste0(cell_type,"_impD",D[i],'_',select_chromosome)
    dis_data <- get(dis_matrix_name)
    # extract distance in i th single cell
    ## extract IF element in simulated distance vector corresponding to position of Kth single cell matrix
    D_length <- n_distance - D[i]
    IF_position <- seq(Kth_matrix, D_length*n_sc, by = n_sc)
    if_ = dis_data[IF_position,]$IF
    If_pos = dis_data[IF_position,1]
    ## region1 and region2 at each distance D's IF vector above
    reg1 = seq(start_region, end_region - D[i]*bin, by = bin)
    reg2 = reg1 + D[i]*bin
    
    ## create frame for result
    l_matrix = (n_distance/2) * (1 + n_distance) #length for sc sparse length
    ## Input result to sc_matrix
    region1 = c(region1, reg1)
    region2 = c(region2, reg2)
    IF = c(IF, if_)
  }
  
  ## assign sc sparse format into global
  sc_matrix = data.frame(region1, region2, IF)
  dataset_name <- paste0(cell_type,"_imp_", Kth_matrix, "_",select_chromosome)
  assign(dataset_name,sc_matrix, envir = .GlobalEnv)
  return(sc_matrix)
}  




### MICE 3 - Algorithm with predictor matrix for each Distance imputed individual group test with no predictor ----

## Function to design data matrix

predictorMatrix_sc_D <- function(cell_type, matrix_position, select_chromosome, distance){
  D.report = bin_distance(cell_type = cell_type, matrix_position , select_chromosome)
  
  ## Extract name of level 1: single cells for each bin
  single_cell <- rep(matrix_position, nrow(D.report[D.report$D==distance,]) )
  
  ## Extract name of level 2: bins
  bin <- rep(D.report$bin[D.report$D==distance],  each= length(matrix_position))
  
  ## Extract IF vector at corresponding distance 
  vector_name <- paste0(cell_type,"_D",distance,'_',select_chromosome)
  IF <- get(vector_name)[,1]
  IF[IF == 0] <- NA
  
  ## Extract name of cell_type
  type = rep(cell_type, length(IF))
  
  ## Combining to data.frame
  matrix <- data.frame(Bins = bin, Single_cell = single_cell, IF = IF, Cell_type= type)
  
  return(matrix)
}


## Function for MICE imputation for all Distance data of each cell type
MICE_impute3 <-  function(cell_type,  matrix_position, select_chromosome,
                          n_imputation = 5, n_distance){
  #n_dis = number of distance
  l = n_distance-1
  for( i in 1:n_distance){
    D_pos = 0:l
    data = predictorMatrix_sc_D(cell_type, matrix_position, select_chromosome, distance = D_pos[i])
    data_input = data[,-1]
    require(mice)
    require(lattice)
    library(miceadds)
    
    ## Classify variable class
    data_input$Single_cell <- as.numeric(data_input$Single_cell)
    data_input$Cell_type= as.numeric(as.factor(data_input$Cell_type))
    
    ## Classify method and predictor matrix
    # get initial default imputation setting
    set.seed(123)
    ini <- mice(data_input, maxit = 0)
    #set up method
    meth= ini$meth
    meth['IF'] <- "2l.pmm" 
    #set up predictor matrix
    pred = ini$pred
    pred[,'Single_cell'] <- -2 #random intercept
    pred[,'Cell_type'] <- 0 #no effect
    pred[,'IF'] <- c(1,0,1) #level 2 imputation
    
    ## Imputation
    #simulate for n time of multiple imputations, with 5 iteration
    imp <- mice(data_input, method = meth, predictorMatrix = pred, print = FALSE,maxit = 1, set.seed(123), m=n_imputation) 
    #returns the long format of all multiple imputation
    imp_data = complete(imp, action = 'long', include = F) 
    
    #if the vector has all if=1, the NA will be 1 too
    if(anyNA(imp_data) == T) {
      imp_data$IF[is.na(imp_data$IF)] <- 1
      agg_new_if2 = imp_data$IF
    } else{
      #if the vector has all if>1, aggregate mean of all imputed complete data
      agg_new_if2 = round(aggregate(imp_data[,4] , by = list(imp_data$.id),FUN= mean))
      agg_new_if2 = agg_new_if2$x}
    
    
    ## Assign each vector into global
    new_data <- data.frame(Bins = data$Bins, Single_cell = data$Single_cell,
                           IF = agg_new_if2, Cell_type = data$Cell_type)
    dataset_name <- paste0(cell_type,"_impD",D_pos[i],'_',select_chromosome)
    assign(dataset_name, new_data , envir = .GlobalEnv)
  }
}  


### Function to Merge these imputed vector into sparse matrices
imp_scKth_maxtrix3 <- function(cell_type , select_chromosome, n_distance,
                               Kth_matrix, start_region, end_region, bin,
                               n_sc){
  l = n_distance-1
  region1 = region2 = IF = NULL
  for(i in 1:length(1:n_distance)){
    # get imputed distance vector
    D = 0:l
    dis_matrix_name <-  paste0(cell_type,"_impD",D[i],'_',select_chromosome)
    dis_data <- get(dis_matrix_name)
    # extract distance in i th single cell
    ## extract IF element in simulated distance vector corresponding to position of Kth single cell matrix
    D_length <- n_distance - D[i]
    IF_position <- seq(Kth_matrix, D_length*n_sc, by = n_sc)
    if_ = dis_data[IF_position,3]
    If_pos = dis_data[IF_position,1]
    ## region1 and region2 at each distance D's IF vector above
    reg1 = seq(start_region, end_region - D[i]*bin, by = bin)
    reg2 = reg1 + D[i]*bin
    
    ## create frame for result
    l_matrix = (n_distance/2) * (1 + n_distance) #length for sc sparse length
    ## Input result to sc_matrix
    region1 = c(region1, reg1)
    region2 = c(region2, reg2)
    IF = c(IF, if_)
  }
  
  ## assign sc sparse format into global
  sc_matrix = data.frame(region1, region2, IF)
  dataset_name <- paste0(cell_type,"_imp_", Kth_matrix, "_",select_chromosome)
  assign(dataset_name,sc_matrix, envir = .GlobalEnv)
  return(sc_matrix)
}  



### MICE 4 - Algorithm with predictor matrix for each Distance imputed individual group test with no predictor and remove outlier ----

## Function for MICE imputation with outlier.rm for all Distance data of each cell type
MICE_impute.outrm <-  function(cell_type,  matrix_position, select_chromosome,
                               n_imputation = 5, n_distance, outlier.rm = TRUE){
  l = n_distance-1
  
  # only imput on distance data have single cell>1
  for( i in 1:n_distance){
    D_pos = 0:l
    data = predictorMatrix_sc_D(cell_type, matrix_position, select_chromosome, distance = D_pos[i])
    data_input = data[,-c(1,4)]
    
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
      agg_new_if2[is.na(agg_new_if2)] <- round(mean(rm.na.data.input$IF))
    } else {
      # get initial default imputation setting
      set.seed(123)
      ini <- mice(data_input, maxit = 0)
      #set up method
      meth= ini$meth
      meth['IF'] <- "2l.pmm" 
      #set up predictor matrix
      pred = ini$pred
      pred[,'Single_cell'] <- -2 #random intercept
      
      ## Imputation
      #simulate for n time of multiple imputations, with 5 iteration
      imp <- mice(data_input, method = meth, predictorMatrix = pred, print = FALSE,maxit = 1, set.seed(123), m=n_imputation) 
      #returns the long format of all multiple imputation
      imp_data = complete(imp, action = 'long', include = F) 
      #if the vector has all if>1, aggregate mean of all imputed complete data
      agg_new_if2 = round(aggregate(imp_data[,4] , by = list(imp_data$.id),FUN= mean))
      agg_new_if2 = agg_new_if2$x}
    
    # Add back outlier, if they were removed in previous steps
    if(outlier.rm == TRUE){
      agg_new_if2[outlier.pos] <- outlier.value
    }
    
    
    ## Assign each vector into global
    new_data <- data.frame(Bins = data$Bins, Single_cell = data$Single_cell,
                           IF = agg_new_if2, Cell_type = data$Cell_type)
    dataset_name <- paste0(cell_type,"_impD",D_pos[i],'_',select_chromosome)
    assign(dataset_name, new_data , envir = .GlobalEnv)
  }
}  


## Function to Merge these imputed vector into sparse matrices
imp_scKth_maxtrix <- function(cell_type , select_chromosome, n_distance,
                              Kth_matrix, start_region, end_region, bin,
                              n_sc){
  l = n_distance-1
  region1 = region2 = IF = NULL
  for(i in 1:length(1:n_distance)){
    # get imputed distance vector
    D = 0:l
    dis_matrix_name <-  paste0(cell_type,"_impD",D[i],'_',select_chromosome)
    dis_data <- get(dis_matrix_name)
    # extract distance in i th single cell
    ## extract IF element in simulated distance vector corresponding to position of Kth single cell matrix
    D_length <- n_distance - D[i]
    IF_position <- seq(Kth_matrix, D_length*n_sc, by = n_sc)
    if_ = dis_data[IF_position,3]
    If_pos = dis_data[IF_position,1]
    ## region1 and region2 at each distance D's IF vector above
    reg1 = seq(start_region, end_region - D[i]*bin, by = bin)
    reg2 = reg1 + D[i]*bin
    
    ## create frame for result
    l_matrix = (n_distance/2) * (1 + n_distance) #length for sc sparse length
    ## Input result to sc_matrix
    region1 = c(region1, reg1)
    region2 = c(region2, reg2)
    IF = c(IF, if_)
  }
  
  ## assign sc sparse format into global
  sc_matrix = data.frame(region1, region2, IF)
  dataset_name <- paste0(cell_type,"_imp_", Kth_matrix, "_",select_chromosome)
  assign(dataset_name,sc_matrix, envir = .GlobalEnv)
  return(sc_matrix)
}  





## scHicNorm ----
scHicNorm <- function(dat_feature, dat_HiC, method = "NBH") {
  # Remove rows with zero values in density, GC, or map
  rows_del_feature <- which(dat_feature$density == 0 | dat_feature$GC == 0 | dat_feature$map == 0)
  rows_del_HiC <- which(apply(dat_HiC, 1, sum) == 0)
  rows_del <- c(rows_del_feature, rows_del_HiC)
  rows_del <- rows_del[!duplicated(rows_del)]
  dat_feature <- dat_feature[-rows_del, ]
  dat_HiC <- dat_HiC[-rows_del, -rows_del]
  
  if (length(rows_del) > 0) {
    cat("The following rows are removed: ", rows_del, ".\n")
  }
  
  mat_HiC <- as.matrix(dat_HiC)
  vec_HiC <- mat_HiC[upper.tri(mat_HiC, diag = FALSE)]
  
  # Get the feature matrix and Z-score
  mat_density <- as.matrix(log(dat_feature[, 3] %o% dat_feature[, 3]))
  mat_GC <- as.matrix(log(dat_feature[, 4] %o% dat_feature[, 4]))
  mat_map <- as.matrix(log(dat_feature[, 5] %o% dat_feature[, 5]))
  mat_density <- (mat_density - mean(c(mat_density))) / sd(c(mat_density))
  mat_GC <- (mat_GC - mean(c(mat_GC))) / sd(c(mat_GC))
  
  vec_density <- mat_density[upper.tri(mat_density, diag = FALSE)]
  vec_GC <- mat_GC[upper.tri(mat_GC, diag = FALSE)]
  vec_map <- mat_map[upper.tri(mat_map, diag = FALSE)]
  
  fit_model <- NULL
  
  if (method == "Poisson") {
    cat("Using Poisson. \n")
    fit_model <- glm(vec_HiC ~ vec_density + vec_GC + offset(vec_map), family = "poisson")
  } else if (method == "NB") {
    cat("Using Negative Binomial. \n")
    fit_model <- glm.nb(vec_HiC ~ vec_density + vec_GC + offset(vec_map))
  } else if (method == "ZIP") {
    cat("Using Zero-inflated Poisson. \n")
    fit_model <- zeroinfl(vec_HiC ~ vec_density + vec_GC + offset(vec_map) | vec_density + vec_GC + vec_map, dist = "poisson")
  } else if (method == "ZINB") {
    cat("Using Zero-inflated Negative Binomial. \n")
    fit_model <- zeroinfl(vec_HiC ~ vec_density + vec_GC + offset(vec_map) | vec_density + vec_GC + vec_map, dist = "negbin")
  } else if (method == "PH") {
    cat("Using Poisson Hurdle. \n")
    fit_model <- hurdle(vec_HiC ~ vec_density + vec_GC + offset(vec_map) | vec_density + vec_GC + vec_map)
  } else if (method == "NBH") {
    cat("Using Negative Ninomial Hurdle. \n")
    fit_model <- hurdle(vec_HiC ~ vec_density + vec_GC + offset(vec_map) | vec_density + vec_GC + vec_map, dist = "negbin")
  } else {
    stop("Please select a valid method.\n", call. = FALSE)
  }
  
  # Output normalized matrix
  fitted_vals <- fit_model$fitted.values
  fitted_vals_mat <- matrix(0, nrow(mat_HiC), ncol(mat_HiC))
  fitted_vals_mat[upper.tri(fitted_vals_mat, diag = FALSE)] <- fitted_vals
  fitted_vals_mat <- t(fitted_vals_mat)
  fitted_vals_mat[upper.tri(fitted_vals_mat, diag = FALSE)] <- fitted_vals
  res <- round(mat_HiC / fitted_vals_mat, 4)
  diag(res) <- 0
  # write.table(res, file = output_file, row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)
  return(res)
}

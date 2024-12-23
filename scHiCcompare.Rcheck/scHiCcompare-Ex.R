pkgname <- "scHiCcompare"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('scHiCcompare')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("plot_HiCmatrix_heatmap")
### * plot_HiCmatrix_heatmap

flush(stderr()); flush(stdout())

### Name: plot_HiCmatrix_heatmap
### Title: Plot Hi-C Interaction Matrix Heatmap
### Aliases: plot_HiCmatrix_heatmap

### ** Examples

data("ODC.bandnorm_chr20_1")

# Call the function with this sparse matrix to generate the heatmap
plot_HiCmatrix_heatmap(
  scHiC.sparse = ODC.bandnorm_chr20_1,
  zlim = c(0, 7), # Log scale color limits
  color_low = "white", # Color for low values
  color_high = "red", # Color for high values
  main = "Single-Cell Hi-C Heatmap", # Title of the plot
  figure_name = "Example Heatmap" # Subtitle for the plot
)



cleanEx()
nameEx("plot_imputed_distance_diagnostic")
### * plot_imputed_distance_diagnostic

flush(stderr()); flush(stdout())

### Name: plot_imputed_distance_diagnostic
### Title: Plot Imputed Distance Diagnostic
### Aliases: plot_imputed_distance_diagnostic

### ** Examples

data("scHiC.table_MG_chr22")
## Impute data above
scHiC.table_MG_imp = scHiCcompare_impute(scHiC.table = scHiC.table_MG_chr22)

# Call the function with this sparse matrix to generate the plot
plot_imputed_distance_diagnostic(raw_sc_data = scHiC.table_MG_chr22,
 imp_sc_data = scHiC.table_MG_imp, D = 1)




cleanEx()
nameEx("pseudo_bulkHic")
### * pseudo_bulkHic

flush(stderr()); flush(stdout())

### Name: pseudo_bulkHic
### Title: Generate Pseudo-bulk Hi-C Data
### Aliases: pseudo_bulkHic

### ** Examples

data("scHiC.table_MG_chr22")
data("scHiC.table_ODC_chr22")
pseudo_bulk_MG = pseudo_bulkHic(scHiC.table_MG_chr22)
pseudo_bulk_ODC = pseudo_bulkHic(scHiC.table_ODC_chr22)
head(pseudo_bulk_MG)
head(pseudo_bulk_ODC)




cleanEx()
nameEx("read_files")
### * read_files

flush(stderr()); flush(stdout())

### Name: read_files
### Title: Read Files for scHi-C Dataset
### Aliases: read_files

### ** Examples

# Load MG data folder example
MGs_example <- system.file("MGs_example", package = "scHiCcompare")
datasets <- read_files(file.path = MGs_example, position.dataset =  c(1, 2, 3, 4, 5),
 txt.sparse.heads.position = c(1, 2, 3, 4, 5)
 )



cleanEx()
nameEx("scHiC.table_MG_chr22")
### * scHiC.table_MG_chr22

flush(stderr()); flush(stdout())

### Name: scHiC.table_MG_chr22
### Title: scHi-C table from microglia (MG) cell type - chromosome 22 at 1
###   MB resolution
### Aliases: scHiC.table_MG_chr22
### Keywords: datasets

### ** Examples

data("scHiC.table_MG_chr22")
head(scHiC.table_MG_chr22)



cleanEx()
nameEx("scHiC_bulk_compare")
### * scHiC_bulk_compare

flush(stderr()); flush(stdout())

### Name: scHiC_bulk_compare
### Title: Compare Bulk Hi-C Data
### Aliases: scHiC_bulk_compare

### ** Examples

# Load data folder example to the current working directory
ODCs_example <- system.file("ODCs_example", package = "scHiCcompare")
MGs_example <- system.file("MGs_example", package = "scHiCcompare")
# Input single-cell Hi-C in sparse format (.txt) from a path
scHiC.table_ODC <- scHiC_table(
  file.path = ODCs_example,
  cell.type = "ODC",
  select.chromosome = "chr20"
)
scHiC.table_MG <- scHiC_table(
  file.path = MGs_example,
  cell.type = "MG",
  select.chromosome = "chr20"
)
# Bulk matrix in sparse format
bulk.sparse.1 <- na.omit(pseudo_bulkHic(scHiC.table = scHiC.table_ODC, out = "sparse"))
bulk.sparse.2 <- na.omit(pseudo_bulkHic(scHiC.table = scHiC.table_MG, out = "sparse"))
# Create the `hic.table` object
library(HiCcompare)
bulk.hic.table <- create.hic.table(bulk.sparse.1, bulk.sparse.2, chr = "chr20", scale = FALSE)
# Jointly normalize data for a single chromosome
hic.table_normalize <- hic_loess(bulk.hic.table, Plot = TRUE, Plot.smooth = FALSE)
# Example usage of the BulkHiC_compare function
result <- scHiC_bulk_compare(hic.table_normalize, D.interval = c(1:10), fprControl.logfc = 0.8)




cleanEx()
nameEx("scHiC_table")
### * scHiC_table

flush(stderr()); flush(stdout())

### Name: scHiC_table
### Title: Create scHiC Interaction Frequency Table
### Aliases: scHiC_table

### ** Examples

# Load MG data folder example
MGs_example <- system.file("MGs_example", package = "scHiCcompare")

# Create scHiC table to be used in scHiCcompare
IF_table <- scHiC_table(file.path = MGs_example, cell.type ='MG', select.chromosome = "chr20")




cleanEx()
nameEx("scHiCcompare")
### * scHiCcompare

flush(stderr()); flush(stdout())

### Name: scHiCcompare
### Title: ScHiCcompare: Differential Analysis of Single-Cell Hi-C Data
### Aliases: scHiCcompare

### ** Examples

## Load example data for ODC and MG file paths
ODCs_example <- system.file("ODCs_example", package = "scHiCcompare")
MGs_example <- system.file("MGs_example", package = "scHiCcompare")

## Run scHiCcompare on example data
result <- scHiCcompare(
  file.path.1 = MGs_example,
  file.path.2 = ODCs_example,
  select.chromosome = "chr20"
)
print(result)




cleanEx()
nameEx("scHiCcompare_impute")
### * scHiCcompare_impute

flush(stderr()); flush(stdout())

### Name: scHiCcompare_impute
### Title: Random Forest Imputation with Pooling options for scHi-C Data
### Aliases: scHiCcompare_impute

### ** Examples

# Load MG data folder example
MGs_example <- system.file("MGs_example", package = "scHiCcompare")
# Create scHicCompare table to be used in scHicCompare
IF_table <- scHiC_table(file.path = MGs_example, cell.type = "MG", select.chromosome = "chr20")
# Example usage of Pooling_RF_impute
library(tidyr)
imputed_table <- scHiCcompare_impute(IF_table,
  n.imputation = 5, outlier.rm = TRUE,
  main.Distances = 1:10000000, pool.style = "progressive"
)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

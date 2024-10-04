setwd("/Users/mynguyen/Documents/GitHub/scHiCcompare")
# Install required packages if you don't have them already
install.packages("usethis")
install.packages("devtools")
install.packages("roxygen2")

library(roxygen2)


# Create a package
usethis::create_package("/Users/mynguyen/Documents/GitHub/scHiCcompare/scHiCcompare")

# Add importance package
usethis::use_package("dplyr") # Adding dplyr to Imports
usethis::use_package("dplyr", type = "Suggests") # Adding dplyr to Suggests
usethis::use_package("tidyr")
usethis::use_package("HiCcompare")
usethis::use_package("mice")
usethis::use_package("miceadds")
usethis::use_package("rstatix")
usethis::use_package("lattice")

# Add function
usethis::use_r("scHiC_table")
usethis::use_r("PseudoBulk_sparse")
usethis::use_r("RF_pooling_imputation")
#usethis::use_r("scHiC_Bulk_loess")
usethis::use_r("ScHiC_bulk_compare")
usethis::use_r("Load_example_MGfolder")
usethis::use_r("Load_example_ODCfolder")

# Add doument for function
usethis::use_roxygen_md()
roxygen2::roxygenise()


# Build Vignete
devtools::build_vignettes()






# Add example data

# Loop to save each dataset in the ODC data/ folder
ODC_names = list.files('/Users/mynguyen/Documents/GitHub/scHiCcompare/Example/Data/ODC', full.names=TRUE, recursive=TRUE)

for (i in 1:length(ODC_names)) {
  dataset <- read.delim(ODC_names[i])
  
  # Assign the dataset to a variable with the desired name
  dataset_name <- paste0("ODC_", i)
  assign(dataset_name, dataset)
  
  # Save each dataset to the 'data/' directory, ensuring the correct name is saved
  save(list = dataset_name, file = paste0("data/", dataset_name, ".rda"))
}


# Loop to save each dataset in the MG data/ folder
MG_names = list.files('/Users/mynguyen/Documents/GitHub/scHiCcompare/Example/Data/MG', full.names=TRUE, recursive=TRUE)

for (i in 1:length(MG_names)) {
  dataset <- read.delim(MG_names[i])
  
  # Assign the dataset to a variable with the desired name
  dataset_name <- paste0("MG_", i)
  assign(dataset_name, dataset)
  
  # Save each dataset to the 'data/' directory, ensuring the correct name is saved
  save(list = dataset_name, file = paste0("data/", dataset_name, ".rda"))
}



# Save scHiC table of MG and ODC
scHiC.table_ODC = scHiC_table(file.path = "Example/Data/ODC", cell.type = 'ODC', 
                             position.dataset =  1:50, type='txt', select.chromosome = 'chr22')
scHiC.table_MG = scHiC_table(file.path = "Example/Data/MG", cell.type = 'MG', 
                              position.dataset =  1:50, type='txt', select.chromosome = 'chr22')
save(scHiC.table_ODC, file = paste0("data/scHiC.table_ODC.rda"))
save(scHiC.table_MG, file = paste0("data/scHiC.table_MG.rda"))




### Test functions
# Set up testthat
usethis::use_testthat()



# Build the package
library(devtools)
build()


# Replace "yourPackageName_0.0.0.9000.tar.gz" with your actual package filename
install.packages("/Users/mynguyen/Documents/GitHub/scHiCcompare_0.0.0.9000.tar.gz", repos = NULL, type = "source")


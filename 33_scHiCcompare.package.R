# Install required packages if you don't have them already
install.packages("usethis")
install.packages("devtools")

# Create a package
usethis::create_package("/Users/mynguyen/Documents/GitHub/scHiCcompare/scHiCcompare")

# Add importance package
usethis::use_package("dplyr") # Adding dplyr to Imports
usethis::use_package("dplyr", type = "Suggests") # Adding dplyr to Suggests
usethis::use_package("HiCcompare")

# Add function
usethis::use_r("scHiC_table")
usethis::use_r("PseudoBulk_sparse")
usethis::use_r("RF_pooling_imputation")

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
usethis::use_r("BulkHiC_loess")
usethis::use_r("BulkHiC_compare")


# Add doument for function
usethis::use_roxygen_md()
roxygen2::roxygenise()


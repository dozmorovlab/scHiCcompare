% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scHiCcompare_impute.R
\name{scHiCcompare_impute}
\alias{scHiCcompare_impute}
\title{Random Forest Imputation with Pooling options for scHi-C Data}
\usage{
scHiCcompare_impute(
  scHiC.table,
  n.imputation = 5,
  maxit = 1,
  outlier.rm = TRUE,
  main.Distances = seq(1e+07),
  pool.style = "progressive",
  missPerc.threshold = 95
)
}
\arguments{
\item{scHiC.table}{A data frame containing interaction frequencies across
single cells, created by the `scHiC_table` function. The first four columns
should represent 'cell', 'chr', 'region1', and 'region2', followed by
columns representing interaction frequencies ('IF') for individual cells.}

\item{n.imputation}{An integer specifying the number of imputations to be
performed. The default is 5.}

\item{maxit}{An integer specifying the number of iterations for the internal
refinement process within a single imputation cycle. The default is 1.}

\item{outlier.rm}{A logical value indicate whether to remove outliers during
the imputation process. The default is TRUE.}

\item{main.Distances}{A vector of integers or 'full' representing the scHiC
data in the main distance range to focus the imputation on, in bp units
(e.g., 1:1,000,000). Genomic distances (in bp) are the number of base pairs
between two regions in the genome (e.g., loci or bins).
The default is from 1 to 10,000,000.}

\item{pool.style}{A string specifying the pooling technique to use. Options
are 'none', 'progressive', or 'Fibonacci'. Default is 'progressive'.}

\item{missPerc.threshold}{An integer specifying the missing value percentage
threshold in each pool band.}
}
\value{
A table in the format of a scHiC table (same structure as the output
 of the `scHiC_table` function) with imputed interaction frequencies (IF)
 across all single cells. The output table is formatted with regions and
 single cells in wide format, with one column per single cell containing
 imputed IF values.
}
\description{
This function performs imputation of single-cell Hi-C (scHi-C) interaction
 frequencies (IF) using Random Forest imputation methods with different
 options for distance-based pooling strategies.
}
\details{
The function first identifies important pools based on the given scHi-C
 genomic distance effect by pooling distance data according to the chosen
 method. For progressive pooling, pools of distances are consecutively
 combined to form larger sets, while Fibonacci pooling uses a Fibonacci
 sequence to combine distances. If the pooling style `none` is selected, the
 band contains individual 1 genomic distance. During the imputation process,
 the function imputes all missing values (NAs) within each pool within the
 main distance range. For distances outside this main focus range, if any
 pool contains more than `missPerc.threshold` missing values, it triggers an
 alternative imputation method, filling in missing values based on the mean
 for distances.
}
\examples{
# Load MG data folder example
MGs_example <- system.file("extdata/MGs_example", package = "scHiCcompare")
# Create scHicCompare table to be used in scHicCompare
IF_table <- scHiC_table(
  file.path = MGs_example, cell.type = "MG",
  select.chromosome = "chr20"
)
# Example usage of Pooling_RF_impute
library(tidyr)
imputed_table <- scHiCcompare_impute(IF_table,
  n.imputation = 5, outlier.rm = TRUE,
  main.Distances = 1:10000000, pool.style = "progressive"
)

}
\references{
Doove, L.L., van Buuren, S., Dusseldorp, E. (2014), Recursive partitioning
for missing data imputation in the presence of interaction Effects.
Computational Statistics & Data Analysis, 72, 92-104.

Shah, A.D., Bartlett, J.W., Carpenter, J., Nicholas, O., Hemingway, H.
(2014), Comparison of random forest and parametric imputation models for
imputing missing data using MICE: A CALIBER study. American Journal of
Epidemiology, doi:10.1093/aje/kwt312.

Van Buuren, S. (2018). Flexible Imputation of Missing Data. Second Edition.
Chapman & Hall/CRC. Boca Raton, FL.
}

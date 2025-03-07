test_that("scHiCcompare runs successfully with valid inputs", {
  # Mock data for scHiCcompare inputs
  # You can replace these with file paths containing test data
  file.path.1 <- system.file("extdata/ODCs_example", package = "scHiCcompare")
  file.path.2 <- system.file("extdata/MGs_example", package = "scHiCcompare")

  # Mock chromosome selection (e.g., "chr22")
  select.chromosome <- "chr20"

  # Run scHiCcompare function with mock parameters
  result <- scHiCcompare(
    file.path.1 = file.path.1,
    file.path.2 = file.path.2,
    select.chromosome = select.chromosome,
    imputation = "RF",
    normalization = "LOESS",
    differential.detect = "MD.cluster",
    main.Distances = 1:1000000,
    pool.style = "progressive",
    n.imputation = 2, # Reduced imputations for testing speed
    maxit = 1,
    outlier.rm = TRUE,
    missPerc.threshold = 95,
    A.min = NULL,
    fprControl.logfc = 0.8,
    alpha = 0.05,
    Plot = FALSE,
    Plot.normalize = FALSE,
    save.output.path = NULL
  )

  # Check if the result is a list with expected components
  expect_true("Differential_Analysis" %in% names(result))
  expect_true("Intermediate" %in% names(result))
  expect_true("Imputation" %in% names(result$Intermediate))
})

library(testthat)
library(SpatialExperiment)

# Test parameters
input_file <- "resources_test/xenium/xenium_tiny.h5mu"
output_file <- tempfile(fileext = ".rds")

test_that("H5MU to SpatialExperiment conversion works", {
  # Check that input file exists
  expect_true(file.exists(input_file))

  # Set parameters for the conversion
  par <- list(
    input = input_file,
    output = output_file,
    modality = "rna",
    var_gene_names = NULL,
    obs_cell_id = "cell_id",
    obsm_spatial_coordinates = "spatial",
    obsm_reduced_dims = NULL
  )

  # Source and run the script
  source("script.R")

  # Check that output file was created
  expect_true(file.exists(output_file))

  # Load the SpatialExperiment object
  spe <- readRDS(output_file)

  # Test the object structure
  expect_s4_class(spe, "SpatialExperiment")
  expect_gt(nrow(spe), 0) # Should have genes
  expect_gt(ncol(spe), 0) # Should have cells

  # Test that assays are present
  expect_gt(length(assays(spe)), 0)
  expect_true("X" %in% names(assays(spe)))

  # Test spatial coordinates
  expect_true(!is.null(spatialCoords(spe)))
  expect_equal(ncol(spatialCoords(spe)), 2) # Should have x, y coordinates
  expect_equal(nrow(spatialCoords(spe)), ncol(spe)) # One coord per cell

  # Test metadata
  expect_gt(ncol(colData(spe)), 0) # Should have cell metadata
  expect_gt(ncol(rowData(spe)), 0) # Should have gene metadata

  # Clean up
  unlink(output_file)
})

cat("All tests passed!\n")

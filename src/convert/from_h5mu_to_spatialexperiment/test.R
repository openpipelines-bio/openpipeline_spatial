library(testthat)
library(SpatialExperiment)
library(hdf5r)
library(Matrix)
library(reticulate)

## VIASH START
meta <- list(
  resources_dir = "resources_test",
  executable = "./target/executable/convert/from_h5mu_to_spatialexperiment/from_h5mu_to_spatialexperiment",
  config = "./src/convert/from_h5mu_to_spatialexperiment/config.vsh.yaml"
)
## VIASH END

mu <- reticulate::import("mudata")
ad <- reticulate::import("anndata")

# Helper function to create mock H5MU test data
create_mock_h5mu <- function(path) {

  n_obs <- 5
  n_var_mod1 <- 4
  n_var_mod2 <- 3

  # ============== MOD1 MODALITY ==============

  mod1_X_data <- matrix(c(
    1, 2, 3, 0,
    4, 5, 6, 2,
    0, 1, 2, 3,
    2, 0, 1, 4,
    1, 3, 0, 2
  ), nrow = n_obs, ncol = n_var_mod1, byrow = TRUE)


  # Create obs dataframe
  mod1_obs <- data.frame(
    Obs1 = c("A", "B", "A", "C", "B"),
    Obs2 = c(0.9, 0.8, 0.95, 0.7, 0.85),
    Obs3 = c(FALSE, FALSE, TRUE, FALSE, FALSE),
    row.names = paste0("cell_", 1:n_obs),
    stringsAsFactors = FALSE
  )
  # Create var dataframe
  mod1_var <- data.frame(
    Feat1 = c("A", "B", "C", "D"),
    Feat2 = c(TRUE, FALSE, TRUE, FALSE),
    Feat3 = c(1.6, 2.2, 1.2, 1.8),
    row.names = paste0("gene_", 1:n_var_mod1),
    stringsAsFactors = FALSE
  )

  # Create layers
  mod1_layers <- list(
    counts = mod1_X_data * 2
  )

  # Create obsm
  obsm_1 <- matrix(c(
    100.5, 200.3,
    150.2, 180.7,
    120.8, 220.1,
    180.4, 160.9,
    200.1, 190.5
  ), nrow = n_obs, ncol = 2, byrow = TRUE)

  obsm_2 <- matrix(c(
    -1.2, 0.8, 0.3,
    1.1, -0.5, -0.2,
    0.3, 1.2, 0.7,
    -0.8, -0.3, 1.1,
    0.9, 0.2, -0.9
  ), nrow = n_obs, ncol = 3, byrow = TRUE)

  mod1_obsm <- list(
    Obsm1 = obsm_1,
    Obsm2 = obsm_2
  )

  # Create uns (unstructured metadata)
  mod1_uns <- list(
    experiment_info = "metadata"
  )

  # Create AnnData object for mod1 using AnnDataR
  ad_mod1 <- ad$AnnData(
    X = mod1_X_data,
    obs = mod1_obs,
    var = mod1_var,
    layers = mod1_layers,
    obsm = mod1_obsm,
    uns = mod1_uns
  )

  # ============== MOD2 MODALITY ==============

  # Create expression matrix
  mod2_X_data <- matrix(c(
    10, 20, 15,
    25, 30, 18,
    12, 22, 20,
    18, 25, 12,
    20, 28, 16
  ), nrow = n_obs, ncol = n_var_mod2, byrow = TRUE)

  # Create obs dataframe
  mod2_obs <- data.frame(
    Obs = c("C", "D", "C", "E", "D"),
    row.names = paste0("cell_", 1:n_obs),
    stringsAsFactors = FALSE
  )

  # Create var dataframe
  mod2_var <- data.frame(
    Feat = c("d", "e", "g"),
    row.names = paste0("protein_", 1:n_var_mod2),
    stringsAsFactors = FALSE
  )

  # Create AnnData object for mod2
  ad_mod2 <- ad$AnnData(
    X = mod2_X_data,
    obs = mod2_obs,
    var = mod2_var
  )

  # ============== CREATE MUDATA ==============

  # Create MuData object using reticulate
  mdata <- mu$MuData(list(
    mod1 = ad_mod1,
    mod2 = ad_mod2
  ))

  # Write Mudata to path
  mdata$write_h5mu(path)
  path
}

# Main test

cat("> Creating mock H5MU file\n")

# Create mock H5MU file
test_h5mu <- tempfile(fileext = ".h5mu")
create_mock_h5mu(test_h5mu)

# Output file
out_rds <- tempfile(fileext = ".rds")

# Run conversion
cat("> Running conversion\n")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", test_h5mu,
    "--modality", "mod1",
    "--output", out_rds,
    "--obsm_spatial_coordinates", "Obsm1"
  )
)

cat("> Checking execution status\n")
expect_equal(out$status, 0)
expect_true(file.exists(out_rds))

cat("> Reading output file\n")
spe <- readRDS(file = out_rds)
expect_s4_class(spe, "SpatialExperiment")

cat("> Opening input file for comparison\n")
mod1 <- mu$read_h5ad(test_h5mu, mod = "mod1")

cat("> Testing dimensions\n")
dim_spe <- dim(spe)
dim_h5mu <- dim(mod1$X)

expect_equal(dim_spe[1], dim_h5mu[2])
expect_equal(dim_spe[2], dim_h5mu[1])
expect_equal(nrow(spe), 4)
expect_equal(ncol(spe), 5)

cat("> Testing colData (obs) transfer and data types\n")
col_data <- colData(spe)
coldata_cols <- colnames(col_data)
obs_cols <- colnames(mod1$obs)
expect_true(all(obs_cols %in% coldata_cols))

# Test data types in colData
expect_true(is.factor(col_data$Obs1))
expect_true(is.numeric(col_data$Obs2))
expect_true(is.logical(col_data$Obs3))

cat("> Testing rowData (var) transfer and data types\n")
row_data <- rowData(spe)
row_names <- colnames(row_data)
var_cols <- colnames(mod1$var)
expect_true(all(var_cols %in% row_names))

# Test data types in rowData
expect_true(is.character(row_data$Feat1))
expect_true(is.logical(row_data$Feat2))
expect_true(is.numeric(row_data$Feat3))

cat("> Testing spatialCoords\n")
spatial_coords <- spatialCoords(spe)
expect_false(is.null(spatial_coords))
expect_equal(ncol(spatial_coords), 2)
expect_equal(nrow(spatial_coords), ncol(spe))
expect_identical(colnames(spatial_coords), c("x", "y"))

# Test spatial coordinate data types and values
expect_true(is.numeric(spatial_coords[, "x"]))
expect_true(is.numeric(spatial_coords[, "y"]))

# Compare with original spatial coordinates
original_spatial <- mod1$obsm[["Obsm1"]]
expect_equal(as.numeric(original_spatial), as.numeric(spatial_coords))

cat("> Testing assay data\n")
counts_matrix <- assays(spe)[["counts"]]
expect_true(is(counts_matrix, "Matrix") || is.matrix(counts_matrix))
expect_true(all(counts_matrix >= 0))
expect_equal(dim(counts_matrix), c(4, 5))

cat("> Testing reducedDims\n")
# PCA should not be in reducedDims since we only specified spatial
red_dims <- reducedDims(spe)
expect_false(is.null(red_dims))
expect_equal(names(red_dims), c("Obsm2"))
expect_equal(dim(red_dims$Obsm2), c(5, 3))
expect_true(is.numeric(red_dims$Obsm2))

# Compare with original spatial coordinates
original_dimred <- mod1$obsm[["Obsm2"]]
expect_equal(as.numeric(red_dims$Obsm2), as.numeric(original_dimred))

# Clean up
unlink(c(test_h5mu, out_rds))

cat("All tests completed!\n")

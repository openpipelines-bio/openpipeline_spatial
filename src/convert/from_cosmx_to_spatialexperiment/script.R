library(SpatialExperimentIO)

### VIASH START
par <- list(
  input = "resources_test/cosmx/test2.zip",
  add_tx_path = TRUE,
  add_polygon_path = FALSE,
  add_fov_positions = TRUE,
  alternative_experiment_features = c(
    "NegPrb", "Negative", "SystemControl", "FalseCode"
  ),
  output = "spe_cosmx_test.rds"
)
meta <- list(
  resources_dir = "src/utils/"
)
### VIASH END

source(paste0(meta$resources_dir, "/unzip_archived_folder.R"))

cat("Reading input data...")
if (tools::file_ext(par$input) == "zip") {
  expected_file_patterns <- c(
    "*.csv",
    "*.parquet"
  )
  tmp_dir <- extract_selected_files(
    par$input,
    members = expected_file_patterns
  )
  cosmx_output_bundle <- file.path(
    tmp_dir,
    tools::file_path_sans_ext(basename(par$input))
  )
} else {
  cosmx_output_bundle <- par$input
}

cat("Setting parameters...")
if (par$add_polygon_path == FALSE && par$add_tx_path == FALSE) {
  add_parquet_paths <- FALSE
} else {
  add_parquet_paths <- TRUE
}

cat("Converting to SpatialExperiment...")
spe <- readCosmxSXE(
  dirName = cosmx_output_bundle,
  returnType = "SPE",
  countMatPattern = "exprMat_file.csv",
  metaDataPattern = "metadata_file.csv",
  coordNames = c("CenterX_global_px", "CenterY_global_px"),
  addFovPos = par$add_fov_positions,
  fovPosPattern = "fov_positions_file.csv",
  addParquetPaths = add_parquet_paths,
  addPolygon = par$add_polygon_path,
  addTx = par$add_tx_path,
  altExps = par$alternative_experiment_features
)

cat("Saving output...")
saveRDS(spe, file = par$output)

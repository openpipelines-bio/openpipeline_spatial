library(SpatialExperimentIO)

### VIASH START
par <- list(
  input = "resources_test/cosmx/Lung5_Rep2_tiny",
  add_tx_path = TRUE,
  add_polygon_path = FALSE,
  add_fov_positions = TRUE,
  alternative_experiment_features = c(
    "NegPrb", "Negative", "SystemControl", "FalseCode"
  ),
  output = "spe_cosmx_test.rds"
)
### VIASH END

if (par$add_polygon_path == FALSE && par$add_tx_path == FALSE) {
  add_parquet_paths <- FALSE
} else {
  add_parquet_paths <- TRUE
}

spe <- readCosmxSXE(
  dirName = par$input,
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

saveRDS(spe, file = par$output)

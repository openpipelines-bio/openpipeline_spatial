library(SpatialExperimentIO)

### VIASH START
par <- list(
  input = "resources_test/cosmx/Lung5_Rep2_tiny",
  add_parquet_paths = TRUE,
  add_fov_positions = TRUE,
  alternative_experiment_features = c("NegPrb", "Negative", "SystemControl", "FalseCode"),
  output = "spe_cosmx_test.rds"
)
### VIASH END


spe <- readCosmxSXE(
  dirName = par$input,
  returnType = "SPE",
  countMatPattern = "exprMat_file.csv",
  metaDataPattern = "metadata_file.csv",
  coordNames = c("CenterX_global_px", "CenterY_global_px"),
  addFovPos = par$add_fov_positions,
  fovPosPattern = "fov_positions_file.csv",
  addParquetPaths = par$add_parquet_paths,
  altExps = par$alternative_experiment_features
)

saveRDS(spe, file = par$output)
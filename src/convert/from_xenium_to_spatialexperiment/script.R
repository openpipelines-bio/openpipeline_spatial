library(SpatialExperimentIO)

### VIASH START
par <- list(
  input = "resources_test/xenium/xenium_tiny",
  add_experiment_xenium = TRUE,
  add_parquet_paths = TRUE,
  alternative_experiment_features = c("NegControlProbe", "UnassignedCodeword", "NegControlCodeword", "antisense", "BLANK"),
  output = "spe_test.rds"
)
### VIASH END


spe <- readXeniumSXE(
  dirName = par$input,
  returnType = "SPE",
  countMatPattern = "cell_feature_matrix.h5",
  metaDataPattern = "cells.parquet",
  coordNames = c("x_centroid", "y_centroid"),
  addExperimentXenium = par$add_experiment_xenium,
  addParquetPaths = par$add_parquet_paths,
  altExps = par$alternative_experiment_features
)

saveRDS(spe, file = par$output)
library(SpatialExperimentIO)

### VIASH START
par <- list(
  input = "resources_test/xenium/xenium_tiny",
  add_experiment_xenium = TRUE,
  add_parquet_paths = TRUE,
  alternative_experiment_features = c(
    "NegControlProbe", "UnassignedCodeword",
    "NegControlCodeword", "antisense", "BLANK"
  ),
  output = "spe_test.rds"
)
meta <- list(
  resources_dir = "src/utils/"
)
### VIASH END

source(paste0(meta$resources_dir, "/unzip_archived_folder.R"))

xenium_output_bundle <- par$input
if (grepl("\\.zip$", xenium_output_bundle)) {
  expected_file_patterns <- c(
    "cell_feature_matrix.h5",
    "*.parquet",
    "experiment.xenium"
  )

  xenium_output_bundle <- extract_selected_files(
    xenium_output_bundle,
    members = expected_file_patterns
  )
}

spe <- readXeniumSXE(
  dirName = xenium_output_bundle,
  returnType = "SPE",
  countMatPattern = "cell_feature_matrix.h5",
  metaDataPattern = "cells.parquet",
  coordNames = c("x_centroid", "y_centroid"),
  addExperimentXenium = par$add_experiment_xenium,
  addParquetPaths = par$add_parquet_paths,
  altExps = par$alternative_experiment_features
)

saveRDS(spe, file = par$output)

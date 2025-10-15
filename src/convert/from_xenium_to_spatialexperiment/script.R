library(SpatialExperimentIO)

### VIASH START
par <- list(
  input = "output-XETG00150__0031015__slidearray0085__20241023__195946",
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

cat("Reading input data...")
if (tools::file_ext(par$input) == "zip") {
  required_file_patterns <- c(
    "**/cell_feature_matrix.h5",
    "**/*.parquet",
    "**/experiment.xenium"
  )
  tmp_dir <- extract_selected_files(
    par$input,
    members = required_file_patterns
  )
  xenium_output_bundle <- file.path(
    tmp_dir,
    tools::file_path_sans_ext(basename(par$input))
  )
} else {
  xenium_output_bundle <- par$input
}

cat("Converting to SpatialExperiment")
spe <- readXeniumSXE(
  dirName = xenium_output_bundle,
  returnType = "SPE",
  countMatPattern = "cell_feature_matrix.h5",
  metaDataPattern = "cells.parquet",
  coordNames = c("x_centroid", "y_centroid"),
  addExperimentXenium = par$add_experiment_xenium,
  addParquetPaths = par$add_parquet_paths,
  altExps = par$alternative_experiment_features,
  addCellBound = TRUE,
  addNucBound = TRUE
)

cat("Saving output...")
saveRDS(spe, file = par$output)

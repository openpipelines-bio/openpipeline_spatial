library(testthat, warn.conflicts = FALSE)
library(SpatialExperimentIO)
library(SpatialExperiment)

## VIASH START
meta <- list(
  executable = "./from_xenium_to_spatialexperiment",
  resources_dir = "resources_test/xenium",
  name = "from_xenium_to_spatial_experiment"
)
## VIASH END

cat("> Checking simple execution\n")

spe <- paste0(
  meta[["resources_dir"]],
  "/xenium_tiny"
)
out_rds <- "output.rds"

cat("> Running ", meta[["name"]], "\n", sep = "")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", spe,
    "--output", out_rds
  )
)

cat("> Checking whether output file exists\n")
expect_equal(out$status, 0)
expect_true(file.exists(out_rds))

cat("> Reading output file\n")
obj <- readRDS(file = out_rds)

cat("> Checking whether Seurat object is in the right format\n")
# Object type
expect_is(obj, "SpatialExperiment")
# Assay structure
expect_equal(names(slot(obj, "assays")), "counts")
# Spatial coordinates
expect_equal(spatialCoordsNames(obj), c("x_centroid", "y_centroid"))
# Alternative experiments
expect_equal(
  altExpNames(obj),
  c("NegControlProbe", "UnassignedCodeword", "NegControlCodeword")
)
# Metadata components
metadata_components <- c(
  "experiment.xenium", "transcripts", "cell_boundaries", "nucleus_boundaries"
)
expect_named(
  metadata(obj),
  metadata_components,
  ignore.order = TRUE
)
# Parquet paths
parquet_components <- c("transcripts", "cell_boundaries", "nucleus_boundaries")
for (component in parquet_components) {
  expect_true(grepl("\\.parquet$", metadata(obj)[[component]]))
}
# Dimensions
input <- readXeniumSXE(
  dirName = spe,
  returnType = "SPE"
)
dim_rds <- dim(obj)
dim_input <- dim(input)

expect_equal(dim_rds, dim_input)



cat("> Checking execution with compressed input\n")

spe <- paste0(
  meta[["resources_dir"]],
  "/xenium_tiny"
)
out_rds <- "output.rds"

create_folder_archive <- function(folder_path, archive = "temp_dir.zip") {
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(folder_path)
  system2("zip", c("-r", archive, "."))
  paste0(folder_path, "/", archive)
}

zipped_spe <- create_folder_archive(spe)

cat("> Running ", meta[["name"]], "\n", sep = "")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", zipped_spe,
    "--output", out_rds
  )
)

cat("> Checking whether output file exists\n")
expect_equal(out$status, 0)
expect_true(file.exists(out_rds))

cat("> Reading output file\n")
obj <- readRDS(file = out_rds)

cat("> Checking whether Seurat object is in the right format\n")
# Object type
expect_is(obj, "SpatialExperiment")
# Assay structure
expect_equal(names(slot(obj, "assays")), "counts")
# Spatial coordinates
expect_equal(spatialCoordsNames(obj), c("x_centroid", "y_centroid"))
# Alternative experiments
expect_equal(
  altExpNames(obj),
  c("NegControlProbe", "UnassignedCodeword", "NegControlCodeword")
)
# Metadata components
metadata_components <- c(
  "experiment.xenium", "transcripts", "cell_boundaries", "nucleus_boundaries"
)
expect_named(
  metadata(obj),
  metadata_components,
  ignore.order = TRUE
)
# Parquet paths
parquet_components <- c("transcripts", "cell_boundaries", "nucleus_boundaries")
for (component in parquet_components) {
  expect_true(grepl("\\.parquet$", metadata(obj)[[component]]))
}
# Dimensions
input <- readXeniumSXE(
  dirName = spe,
  returnType = "SPE"
)
dim_rds <- dim(obj)
dim_input <- dim(input)

expect_equal(dim_rds, dim_input)


cat("> Checking parameter functionality\n")

out_rds_ext <- "output_ext.rds"

cat("> Running ", meta[["name"]], "\n", sep = "")
out_ext <- processx::run(
  meta[["executable"]],
  c(
    "--input", spe,
    "--add_experiment_xenium", FALSE,
    "--add_parquet_paths", FALSE,
    "--alternative_experiment_features", c("NegControlProbe"),
    "--output", out_rds_ext
  )
)

cat("> Checking whether output file exists\n")
expect_equal(out_ext$status, 0)
expect_true(file.exists(out_rds_ext))

cat("> Reading output file\n")
obj_ext <- readRDS(file = out_rds_ext)

cat("> Checking whether Seurat object is in the right format\n")
# Object type
expect_is(obj_ext, "SpatialExperiment")
# Assay structure
expect_equal(names(slot(obj_ext, "assays")), "counts")
# Spatial coordinates
expect_equal(spatialCoordsNames(obj_ext), c("x_centroid", "y_centroid"))
# Alternative experiments
expect_equal(altExpNames(obj_ext), c("NegControlProbe"))
# Metadata components
expect_true(length(metadata(obj_ext)) == 0)

dim_rds_ext <- dim(obj_ext)
expect_true(identical(dim_rds_ext[2], dim_input[2]))
expect_false(identical(dim_rds_ext[1], dim_input[1]))

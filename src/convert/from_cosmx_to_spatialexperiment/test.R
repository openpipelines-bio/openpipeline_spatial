library(testthat, warn.conflicts = FALSE)
library(SpatialExperimentIO)
library(SpatialExperiment)

## VIASH START
meta <- list(
  executable = "./from_cosmx_to_spatialexperiment",
  resources_dir = "resources_test/cosmx/",
  name = "from_cosmx_to_spatialexperiment"
)
## VIASH END

cat("> Checking simple execution\n")

spe <- paste0(
  meta[["resources_dir"]],
  "/Lung5_Rep2_tiny"
)
out_rds <- "output.rds"

cat("> Running ", meta[["name"]], "\n", sep = "")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", spe,
    "--add_tx_path", TRUE,
    "--add_polygon_path", FALSE,
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
expect_equal(
  spatialCoordsNames(obj),
  c("CenterX_global_px", "CenterY_global_px")
)
# Alternative experiments
expect_equal(altExpNames(obj), c("NegPrb"))
# Metadata components
expect_named(
  metadata(obj),
  c("fov_positions", "transcripts"),
  ignore.order = TRUE
)
# Parquet paths
expect_true(grepl("\\.parquet$", metadata(obj)[["transcripts"]]))
# Dimensions
input <- readCosmxSXE(
  dirName = spe,
  addParquetPaths = FALSE,
  returnType = "SPE"
)

dim_rds <- dim(obj)
dim_input <- dim(input)

expect_equal(dim_rds, dim_input)


cat("> Checking execution with compressed input\n")

spe <- paste0(meta[["resources_dir"]], "/Lung5_Rep2_tiny")
out_rds <- "output.rds"

create_folder_archive <- function(
    folder_path,
    archive = "Lung5_Rep2_tiny.zip") {
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(meta$resources_dir)
  system2("zip", c("-r", archive, "Lung5_Rep2_tiny"))
  paste0(meta$resources_dir, "/", archive)
}

zipped_spe <- create_folder_archive(spe)

cat("> Running ", meta[["name"]], "\n", sep = "")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", zipped_spe,
    "--add_tx_path", TRUE,
    "--add_polygon_path", FALSE,
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
expect_equal(
  spatialCoordsNames(obj),
  c("CenterX_global_px", "CenterY_global_px")
)
# Alternative experiments
expect_equal(altExpNames(obj), c("NegPrb"))
# Metadata components
expect_named(
  metadata(obj),
  c("fov_positions", "transcripts"),
  ignore.order = TRUE
)
# Parquet paths
expect_true(grepl("\\.parquet$", metadata(obj)[["transcripts"]]))
# Dimensions
input <- readCosmxSXE(
  dirName = spe,
  addParquetPaths = FALSE,
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
    "--add_fov_positions", FALSE,
    "--add_tx_path", FALSE,
    "--add_polygon_path", FALSE,
    "--alternative_experiment_features", c("Negative"),
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
expect_equal(
  spatialCoordsNames(obj_ext),
  c("CenterX_global_px", "CenterY_global_px")
)
# Alternative experiments
expect_length(altExpNames(obj_ext), 0)
# Metadata components
expect_length(metadata(obj_ext), 0)

dim_rds_ext <- dim(obj_ext)
expect_true(identical(dim_rds_ext[2], dim_input[2]))
expect_false(identical(dim_rds_ext[1], dim_input[1]))

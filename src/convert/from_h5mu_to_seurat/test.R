library(testthat, warn.conflicts = FALSE)
library(hdf5r)
library(Seurat)

## VIASH START
meta <- list(
  executable = "target/executable/convert/from_h5ad_to_spatial_seurat",
  resources_dir = "resources_test",
  name = "from_h5ad_to_spatial_seurat"
)
## VIASH END


# ---- No FOV ----------------------------------------------------------
cat("> Test conversion without adding FOV\n")

in_h5mu <- paste0(
  meta[["resources_dir"]],
  "/xenium_tiny.h5mu"
)
out_rds <- "output.rds"

cat("> Running ", meta[["name"]], "\n", sep = "")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", in_h5mu,
    "--output", out_rds,
    "--modality", "rna",
    "--assay", "Xenium"
  )
)

cat("> Checking whether output file exists\n")
expect_equal(out$status, 0)
expect_true(file.exists(out_rds))

cat("> Reading output file\n")
obj <- readRDS(file = out_rds)
adata <- H5File$new(in_h5mu, mode = "r")[["/mod/rna/X"]]

cat("> Checking whether Seurat object is in the right format\n")
expect_equal(Assays(obj), "Xenium")
expect_true(all(Layers(obj) == c("counts")))

dim_rds <- dim(obj)
dim_ad <- adata$attr_open("shape")$read()

expect_equal(dim_rds[1], dim_ad[2])
expect_equal(dim_rds[2], dim_ad[1])

expect_false("fov" %in% names(obj))

# # ---- Xenium ----------------------------------------------------------
cat("> Test conversion Xenium\n")

in_h5mu <- paste0(
  meta[["resources_dir"]],
  "/xenium_tiny.h5mu"
)
out_rds <- "output.rds"

cat("> Running ", meta[["name"]], "\n", sep = "")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", in_h5mu,
    "--output", out_rds,
    "--modality", "rna",
    "--assay", "Xenium",
    "--obsm_centroid_coordinates", "spatial"
  )
)

cat("> Checking whether output file exists\n")
expect_equal(out$status, 0)
expect_true(file.exists(out_rds))

cat("> Reading output file\n")
obj <- readRDS(file = out_rds)
adata <- H5File$new(in_h5mu, mode = "r")[["/mod/rna/X"]]

cat("> Checking whether Seurat object is in the right format\n")
expect_equal(Assays(obj), "Xenium")
expect_true(all(Layers(obj) == c("counts")))

dim_rds <- dim(obj)
dim_ad <- adata$attr_open("shape")$read()

expect_equal(dim_rds[1], dim_ad[2])
expect_equal(dim_rds[2], dim_ad[1])

cat("> Checking FOV object\n")
expect_true("fov" %in% names(obj))
expect_true("fov" %in% Images(obj))

fov <- obj[["fov"]]
expect_equal(fov@assay, "Xenium")
expect_equal(fov@key, "Xenium_")

centroids <- fov@boundaries$centroids
expect_equal(nrow(centroids@coords), dim_rds[2])

centroid_coords <- centroids@coords
expect_true(is.numeric(centroid_coords[, 1]))
expect_true(is.numeric(centroid_coords[, 2]))
expect_false(any(is.na(centroid_coords)))

# # ---- Xenium with args-------------------------------------------------
cat("> Test conversion Xenium with centroid arguments\n")

in_h5mu <- paste0(
  meta[["resources_dir"]],
  "/xenium_tiny.h5mu"
)
out_rds <- "output.rds"

cat("> Running ", meta[["name"]], "\n", sep = "")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", in_h5mu,
    "--output", out_rds,
    "--modality", "rna",
    "--assay", "Xenium",
    "--obsm_centroid_coordinates", "spatial",
    "--centroid_nsides", "8",
    "--centroid_radius", "3",
    "--centroid_theta", "0.1"
  )
)

cat("> Checking whether output file exists\n")
expect_equal(out$status, 0)
expect_true(file.exists(out_rds))

cat("> Reading output file\n")
obj <- readRDS(file = out_rds)
adata <- H5File$new(in_h5mu, mode = "r")[["/mod/rna/X"]]

cat("> Checking FOV object\n")
fov <- obj[["fov"]]
centroids <- fov@boundaries$centroids
expect_equal(centroids@nsides, 8)
expect_equal(centroids@radius, 3)
expect_equal(centroids@theta, 0.1)


# ---- CosMx ----------------------------------------------------------

cat("> Test conversion CosMx\n")

in_h5mu <- paste0(
  meta[["resources_dir"]],
  "/Lung5_Rep2_tiny.h5mu"
)
out_rds <- "output.rds"

cat("> Running ", meta[["name"]], "\n", sep = "")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", in_h5mu,
    "--output", out_rds,
    "--modality", "rna",
    "--assay", "CosMx",
    "--obsm_centroid_coordinates", "spatial"
  )
)

cat("> Checking whether output file exists\n")
expect_equal(out$status, 0)
expect_true(file.exists(out_rds))

cat("> Reading output file\n")
obj <- readRDS(file = out_rds)
adata <- H5File$new(in_h5mu, mode = "r")[["/mod/rna/X"]]

cat("> Checking whether Seurat object is in the right format\n")
expect_equal(Assays(obj), "CosMx")
expect_true(all(Layers(obj) == c("counts")))

dim_rds <- dim(obj)
dim_ad <- adata$attr_open("shape")$read()

expect_equal(dim_rds[1], dim_ad[2])
expect_equal(dim_rds[2], dim_ad[1])

cat("> Checking FOV object\n")
expect_true("fov" %in% names(obj))
expect_true("fov" %in% Images(obj))

fov <- obj[["fov"]]
expect_equal(fov@assay, "CosMx")
expect_equal(fov@key, "CosMx_")

centroids <- fov@boundaries$centroids
expect_equal(nrow(centroids@coords), dim_rds[2])

centroid_coords <- centroids@coords
expect_true(is.numeric(centroid_coords[, 1]))
expect_true(is.numeric(centroid_coords[, 2]))
expect_false(any(is.na(centroid_coords)))

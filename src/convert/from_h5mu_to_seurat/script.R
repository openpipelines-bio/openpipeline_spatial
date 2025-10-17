library(anndataR)
library(hdf5r)
library(Seurat)

### VIASH START
par <- list(
  input = "resources_test/xenium/xenium_tiny_processed.h5mu",
  output = "test.rds",
  obsm_centroid_coordinates = "spatial",
  assay = "RNA",
  centroid_nsides = 8,
  centroid_radius = 3,
  centroid_theta = 0.1,
  modality = "rna"
)
### VIASH END


h5mu_to_h5ad <- function(h5mu_path, modality_name) {
  tmp_path <- tempfile(fileext = ".h5ad")
  mod_location <- paste("mod", modality_name, sep = "/")
  h5src <- hdf5r::H5File$new(h5mu_path, "r")
  h5dest <- hdf5r::H5File$new(tmp_path, "w")
  # Copy over the child objects and the child attributes from root
  children <- hdf5r::list.objects(h5src,
    path = mod_location,
    full.names = FALSE, recursive = FALSE
  )
  for (child in children) {
    h5dest$obj_copy_from(
      h5src, paste(mod_location, child, sep = "/"),
      paste0("/", child)
    )
  }
  # Also copy the root attributes
  root_attrs <- hdf5r::h5attr_names(x = h5src)
  for (attr in root_attrs) {
    h5a <- h5src$attr_open(attr_name = attr)
    robj <- h5a$read()
    h5dest$create_attr_by_name(
      attr_name = attr,
      obj_name = ".",
      robj = robj,
      space = h5a$get_space(),
      dtype = h5a$get_type()
    )
  }
  h5src$close()
  h5dest$close()

  tmp_path
}

# Read in H5AD
h5ad_path <- h5mu_to_h5ad(par$input, par$modality)

# Convert to Seurat
seurat_obj <- read_h5ad(
  h5ad_path,
  mode = "r",
  as = "Seurat",
  assay_name = par$assay
)

# Create Centroids object
if (!is.null(par$obsm_centroid_coordinates)) {
  reductions <- seurat_obj@reductions[[par$obsm_centroid_coordinates]]
  spatial_coords <- as.data.frame(reductions@cell.embeddings)
  colnames(spatial_coords) <- c("x_coord", "y_coord")

  if (is.null(par$centroid_nsides)) {
    par$centroid_nsides <- Inf
  }

  if (is.null(par$centroid_theta)) {
    par$centroid_theta <- 0L
  }

  centroids <- CreateCentroids(
    coords = spatial_coords,
    nsides = par$centroid_nsides,
    radius = par$centroid_radius,
    theta = par$centroid_theta
  )

  # Create FOV object
  fov <- CreateFOV(coords = centroids, assay = par$assay)
  seurat_obj[["fov"]] <- fov
  seurat_obj@reductions[[par$obsm_centroid_coordinates]] <- NULL
}


saveRDS(seurat_obj, file = par$output)

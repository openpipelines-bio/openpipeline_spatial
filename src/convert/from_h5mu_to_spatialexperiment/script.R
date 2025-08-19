library(SpatialExperiment)
library(SingleCellExperiment)
library(anndataR)
library(hdf5r)

## VIASH START
par <- list(
  input = "resources_test/xenium/xenium_tiny.h5mu",
  output = "xenium_test.rds",
  modality = "rna",
  obsm_spatial_coordinates = "spatial"
)
## VIASH END


h5mu_to_h5ad <- function(h5mu_path, modality_name) {
  tmp_path <- tempfile(fileext = ".h5ad")
  mod_location <- paste("mod", modality_name, sep = "/")
  h5src <- hdf5r::H5File$new(h5mu_path, "r")
  h5dest <- hdf5r::H5File$new(tmp_path, "w")
  # Copy over the child objects and the child attributes from root
  # Root cannot be copied directly because it always exists and
  # copying does not allow overwriting.
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

read_spatial_coordinates <- function(sce, spatial_coordinates_name) {
  if (par$obsm_spatial_coordinates %in% names(reducedDims(sce))) {
    spatial_coords <- reducedDims(sce)[[par$obsm_spatial_coordinates]]
    if (ncol(spatial_coords) != 2) {
      stop(
        "Spatial coordinates must have 2 columns, but found ",
        ncol(spatial_coords), " columns"
      )
    }
    # Set proper column names for spatial coordinates
    colnames(spatial_coords) <- c("x", "y")
  } else {
    warning(
      "Spatial coordinates '", par$obsm_spatial_coordinates,
      "' not found in reducedDims. Available dimensions: ",
      paste(names(reducedDims(sce)), collapse = ", ")
    )
    spatial_coords <- NULL
  }
  spatial_coords
}

main <- function() {
  # Convert to AnnData
  cat("Converting H5MU file to H5AD...\n")
  h5file <- h5mu_to_h5ad(par$input, par$modality)

  # Convert to SpatialExperiment
  cat("Converting to SingleCellExperiment...\n")
  sce <- read_h5ad(h5file, as = "SingleCellExperiment")

  # Extract spatial coordinates if specified
  if (
    !is.null(par$obsm_spatial_coordinates) &&
      length(par$obsm_spatial_coordinates) > 0
  ) {
    cat("Reading in spatial coordinates...\n")
    spatial_coords <- read_spatial_coordinates(
      sce, par$obsm_spatial_coordinates
    )
    reducedDims(sce)[[par$obsm_spatial_coordinates]] <- NULL
  } else {
    spatial_coords <- NULL
  }

  # Converting SingleCellExperiment to SpatialExperiment
  cat("Converting to SpatialExperiment...\n")
  spe <- as(sce, "SpatialExperiment")
  spatialCoords(spe) <- spatial_coords

  # Saving SpatialExperiment object
  cat("Saving SpatialExperiment object to:", par$output, "\n")
  saveRDS(spe, file = par$output, compress = FALSE)
}

main()

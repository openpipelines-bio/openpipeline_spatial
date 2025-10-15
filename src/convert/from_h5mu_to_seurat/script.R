library(anndataR)
library(hdf5r)
# library(Seurat)

### VIASH START
par <- list(
  input = "resources_test_sc/10x_5k_anticmv/5k_human_antiCMV_T_TBNK_connect.h5mu",
  output = "resources_test_sc/10x_5k_anticmv/5k_human_antiCMV_T_TBNK_connect.rds",
  assay = c("RNA", "ADT", "TCR"),
  modality = c("rna", "prot", "vdj_t")
)
### VIASH END

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

has_zero_dimension <- function(adata) {
  adata_dims <- adata$shape()
  any(adata_dims == 0)
}

map_obs_to_metadata <- function(adata, seurat_obj, mod, assay) {
  cells <- colnames(seurat_obj)
  mod_meta <- adata$obs
  obs_names <- rownames(mod_meta)
  colnames(mod_meta) <- paste0(assay, "_", make.names(colnames(mod_meta)))

  if (!all(obs_names %in% cells)) {
    stop(paste0(
      "Not all cells in the adata modality", mod,
      "are present in the Seurat object."
    )
    )
  }
  mod_meta <- mod_meta[match(cells, obs_names), , drop = FALSE]
  rownames(mod_meta) <- cells

  seurat_obj <- AddMetaData(seurat_obj, mod_meta)
}

map_modality_to_assay <- function(adata, seurat_obj, mod, assay) {


}

# Initialize seurat object
seurat_obj <- NULL
modalities_to_metadata <- list()

# Check that modalities and assays have the same length
if (length(par$modality) != length(par$assay)) {
  stop("The number of modalities should match the number of assays.")
}

# Loop through modalities and assays
for (i in seq_along(par$modality)) {
  cat("Processing modality:", par$modality[i], "as assay:", par$assay[i], "\n")

  # Read the specific modality from h5mu file
  adata_path <- h5mu_to_h5ad(par$input, par$modality[i])
  adata <- read_h5ad(adata_path, mode = "r+")

  # Check dimensions
  # Modalities with dimension 0 will be added later
  if (has_zero_dimension(adata)) {
    if (adata$shape()[1] == 0) {
      cat("Skipping modality", par$modality[i], "- has zero observations\n")
    } else {
      cat(
        "Modality",
        par$modality[i],
        "has zero features - moving to metadata\n"
      )
      modalities_to_metadata[[par$modality[i]]] <- par$assay[i]
    }
    next
  }

  # Convert to Seurat assay
  temp_seurat <- read_h5ad(
    adata_path,
    mode = "r",
    as = "Seurat",
    assay_name = par$assay[i]
  )

  # If this is the first modality, initialize the main Seurat object
  if (is.null(seurat_obj)) {
    seurat_obj <- temp_seurat
  } else {
    # Add as additional assay to existing Seurat object
    seurat_obj[[par$assay[i]]] <- temp_seurat[[par$assay[i]]]
  }
}

if (is.null(seurat_obj)) {
  stop("No valid modalities found to create a Seurat object. At least one modality must have non-zero dimensions.")
}

for (mod in names(modalities_to_metadata)) {
  assay <- modalities_to_metadata[[mod]]
  cat("Moving modality", mod, "to metadata\n")
  seurat_obj <- map_obs_to_metadata(adata, seurat_obj, mod, assay)
}

saveRDS(seurat_obj, file = par$output)

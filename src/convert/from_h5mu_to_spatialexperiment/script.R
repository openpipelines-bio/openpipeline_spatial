library(SpatialExperiment)
library(hdf5r)
library(Matrix)
library(S4Vectors)

### VIASH START
par <- list(
  input = "/Users/dorienroosen/code/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
  output = "xenium_tiny_qc.h5mu",
  modality = "rna",
  var_gene_names = NULL,
  obs_cell_id = NULL,
  obsm_spatial_coordinates = NULL,
  obsm_reduced_dims = c("X_pca", "X_umap")
)
### VIASH END


read_matrix_from_h5 <- function(h5file, matrix_path, obs_names, feat_names) {
  if (!h5file$exists(matrix_path)) {
    stop("Matrix not found at ", matrix_path)
  }

  cat("Reading matrix from:", matrix_path)
  counts <- h5file[[matrix_path]]

  # Check matrix encoding and dimensions
  if (counts$attr_exists("encoding-type")) {
    encoding <- counts$attr_open("encoding-type")$read()

    if (!encoding %in% c("csr_matrix", "csc_matrix")) {
      stop("Unsupported sparse matrix encoding: ", encoding)
    }

    cat("Reading sparse", encoding, "...")
    data <- counts[["data"]]$read()
    indices <- counts[["indices"]]$read()
    indptr <- counts[["indptr"]]$read()
    shape <- counts$attr_open("shape")$read()

    if (encoding == "csr_matrix") {
      # SPARSE MATRIX (CSR format)
      row_indices <- rep(seq_len(length(indptr) - 1), diff(indptr))
      counts_matrix <- sparseMatrix(
        i = row_indices,
        j = indices + 1,  # Convert to 1-based indexing
        x = data,
        dims = shape
      )
    } else if (encoding == "csc_matrix") {
      # SPARSE MATRIX (CSC format)
      counts_matrix <- sparseMatrix(
        j = indices + 1,
        p = indptr,
        x = data,
        dims = shape,
        index1 = FALSE
      )
    } 
  } else {
    # DENSE MATRIX - read directly
    cat("Reading dense matrix...")
    dense_matrix <- counts$read()
    counts_matrix <- as(dense_matrix, "sparseMatrix")
  }

  # Transpose and assign rownames/colnames
  counts_matrix <- t(counts_matrix)
  colnames(counts_matrix) <- obs_names
  rownames(counts_matrix) <- feat_names
  counts_matrix
}

read_h5_layers <- function(h5file, mod_location, obs_names, feat_names) {
  assay_list <- list()
  
  # Read X matrix (main assay)
  x_path <- paste0(mod_location, "/X")
  if (h5file$exists(x_path)) {
    cat("Reading X matrix...")
    assay_list[["X"]] <- read_matrix_from_h5(h5file, x_path, obs_names, feat_names)
  }
  
  # Read all layers
  layers_path <- paste0(mod_location, "/layers")
  if (h5file$exists(layers_path)) {
    layers_group <- h5file[[layers_path]]
    layer_names <- names(layers_group)
    
    for (layer_name in layer_names) {
      layer_path <- paste0(layers_path, "/", layer_name)
      cat("Reading layer:", layer_name)
      assay_list[[layer_name]] <- read_matrix_from_h5(h5file, layer_path, obs_names, feat_names)
    }
  }
  
  if (length(assay_list) == 0) {
    stop("No matrices found. Expected at least X or layers.")
  }
  
  return(assay_list)
}

read_h5_field <- function(h5file, base_path, field_name = NULL, index_type = NULL) {
  if (!is.null(field_name) && field_name != "") {
    field_path <- paste0(base_path, "/", field_name)
    if (!h5file$exists(field_path)) {
      stop("Field path not found: ", field_path)
    }
    return(h5file[[field_path]]$read())
  } else {
    # Read index values from the index field
    group <- h5file[[base_path]]
    
    # First, get the index field name from the _index attribute
    index_attr <- "_index"
    
    # Check if attribute exists using h5attr_names
    attr_names <- hdf5r::h5attr_names(group)
    if (!index_attr %in% attr_names) {
      stop("Index attribute '", index_attr, "' not found in group: ", base_path)
    }
    
    # Read the index field 
    h5a <- group$attr_open(attr_name = index_attr)
    index_field_name <- h5a$read()
    h5a$close()
    field_path <- paste0(base_path, "/", index_field_name)
    if (!h5file$exists(field_path)) {
      stop("Index field '", index_field_name, "' not found at path: ", field_path)
    }
    
    return(h5file[[field_path]]$read())
  }
}

read_dataframe_from_h5 <- function(h5file, group_path, index_names) {
  if (!h5file$exists(group_path)) return(DataFrame())
  
  # Get all child objects in the group
  children <- hdf5r::list.objects(h5file, 
    path = group_path,
    full.names = FALSE, 
    recursive = FALSE
  )
  
  # Filter out metadata columns
  
  df <- DataFrame(row.names = index_names)
  
  for (col in children) {
    tryCatch({
      col_path <- paste(group_path, col, sep = "/")
      col_obj <- h5file[[col_path]]
      
      # Check if it's a dataset and try to read
      if (inherits(col_obj, "H5D")) {
        col_data <- col_obj$read()
        
        # Handle NaN/missing values
        if (is.numeric(col_data) && any(is.nan(col_data))) {
          col_data[is.nan(col_data)] <- NA
        }
        
        df[[col]] <- col_data
      } else if (inherits(col_obj, "H5Group")) {
        # Handle categorical data (groups with codes + categories)
        if (col_obj$exists("codes") && col_obj$exists("categories")) {
          codes <- col_obj[["codes"]]$read()
          categories <- col_obj[["categories"]]$read()
          col_data <- categories[codes + 1]  # Convert 0-based to 1-based indexing
          df[[col]] <- col_data
        } else {
          warning("Unknown group structure for column '", col, "' - skipping")
        }
      } else {
        # Try direct read for other types
        col_data <- col_obj$read()
        df[[col]] <- col_data
      }
    }, error = function(e) {
      # More detailed error message for debugging
      warning("Could not read column '", col, "' at path '", col_path, "': ", e$message)
    })
  }
  
  df
}

read_spatial_coords <- function(h5file, mod_location, spatial_coord_path, obs_names) {
  if (!h5file$exists(spatial_coord_path)) {
    stop("Spatial coordinates path not found: ", spatial_coord_path)
  }
  
  coords <- t(h5file[[spatial_coord_path]]$read())

  # Ensure it's a 2D matrix
    if (!is.matrix(coords) || ncol(coords) != 2) {
      stop("Spatial coordinates must be a matrix with exactly 2 columns. Found: ", 
          if(is.matrix(coords)) ncol(coords) else "not a matrix")
    }

    # Set column names as x, y coordinates
    colnames(coords) <- c("x", "y")

    if (nrow(coords) != length(obs_names)) {
      stop("Spatial coordinates dimension mismatch: ", nrow(coords), " coords vs ", length(obs_names), " observations")
    }

    rownames(coords) <- obs_names
    coords
}

read_reduced_dims <- function(h5file, mod_location, obs_names, reduced_dim_names = NULL) {
  reduced_dims <- list()
  
  # If no reduced dims specified, return empty list
  if (is.null(reduced_dim_names) || length(reduced_dim_names) == 0) {
    cat("No reduced dimensions specified")
    return(reduced_dims)
  }
  
  for (obsm_name in reduced_dim_names) {
    obsm_field_path <- paste0(mod_location, "/obsm/", obsm_name)
    
    if (!h5file$exists(obsm_field_path)) {
      warning("Reduced dimension '", obsm_name, "' not found at path: ", obsm_field_path)
      next
    }
    
    tryCatch({
      cat("Reading reduced dimension:", obsm_name)
      obsm_data <- h5file[[obsm_field_path]]$read()
      
      # Ensure it's a matrix and transpose if needed (cells should be rows)
      if (is.matrix(obsm_data)) {
        # If rows don't match obs_names, try transposing
        if (nrow(obsm_data) != length(obs_names) && ncol(obsm_data) == length(obs_names)) {
          obsm_data <- t(obsm_data)
        }
        
        if (nrow(obsm_data) == length(obs_names)) {
          rownames(obsm_data) <- obs_names
          reduced_dims[[obsm_name]] <- obsm_data
        } else {
          warning("Dimension mismatch for reduced dimension '", obsm_name, "': ", 
                  nrow(obsm_data), " rows vs ", length(obs_names), " cells")
        }
      } else {
        warning("Reduced dimension '", obsm_name, "' is not a matrix")
      }
    }, error = function(e) {
      warning("Could not read reduced dimension '", obsm_name, "': ", e$message)
    })
  }
  
  return(reduced_dims)
}

main <- function() {

  cat("Converting H5MU to SpatialExperiment...")

  # Open H5 file and verify modality
  h5file <- hdf5r::H5File$new(par$input, "r")
  mod_location <- paste("mod", par$modality, sep = "/")
  
  tryCatch({
    if (!h5file$exists(mod_location)) {
      stop("Could not find modality '", par$modality)
    }

    # Get observation_and_feature_name
    cat("Reading observation and feature names...")
    var_names <- read_h5_field(h5file, paste0(mod_location, "/var"), par$var_gene_names)
    obs_names <- read_h5_field(h5file, paste0(mod_location, "/obs"), par$obs_cell_id)

    # Read all matrices (X and layers)
    cat("Reading matrices...")
    assay_list <- read_h5_layers(h5file, mod_location, obs_names, var_names)

    # Read metadata
    cat("Reading obs metadata...")
    obs_data <- read_dataframe_from_h5(h5file, paste0(mod_location, "/obs"), obs_names)

    cat("Reading var metadata...")  
    var_data <- read_dataframe_from_h5(h5file, paste0(mod_location, "/var"), var_names)

    # Read spatial coordinates (optional)
    if (!is.null(par$obsm_spatial_coordinates) && par$obsm_spatial_coordinates != "") {
      cat("Reading spatial coordinates...")
      spatial_coord_path <- paste0(mod_location, "/obsm/", par$obsm_spatial_coordinates)
      spatial_coords <- read_spatial_coords(h5file, mod_location, spatial_coord_path, obs_names)
    } else {
      spatial_coords <- NULL
    }

    # Read obsm dimred fields (optional)
    cat("Reading reduced dimensions...")
    if (!is.null(par$obsm_reduced_dims) && length(par$obsm_reduced_dims) > 0) {
      reduced_dims <- read_reduced_dims(h5file, mod_location, obs_names, par$obsm_reduced_dims)
    } else {
      reduced_dims <- NULL
    }


    # Create SpatialExperiment object
    cat("Creating SpatialExperiment object...")
    spe <- SpatialExperiment(
      assays = assay_list,
      spatialCoords = spatial_coords,
      reducedDims = reduced_dims,
      colData = obs_data,
      rowData = var_data
    )

    # Save with compression for large files
    cat("Saving SpatialExperiment object...")
    saveRDS(spe, file = par$output, compress = TRUE)
    
  }, finally = {
    h5file$close()
  })
}

# Execute main function
main()

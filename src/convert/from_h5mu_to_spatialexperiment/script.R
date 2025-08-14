library(SpatialExperiment)
library(hdf5r)
library(Matrix)
library(S4Vectors)

### VIASH START
par <- list(
  input = "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
  output = "resources_test/annotation_test_data/TS_Blood_filtered.rds",
  assay = "RNA",
  modality = "rna",
  layer = NULL,
  output_compression = "gzip"
)
### VIASH END


# STEP 2: Matrix reading with layer support
read_counts_matrix <- function(h5file, mod_location, layer_name = NULL) {
  # Determine matrix location: either X or in layers
  if (is.null(layer_name) || layer_name == "") {
    layer <- paste0(mod_location, "/X")
  } else {
    layer <- paste0(mod_location, "/layers/", layer_name)
  }

  if (!h5file$exists(layer)) {
    stop("Counts matrix not found at ", layer)
  }

  cat("Reading matrix from:", layer, "\n")
  counts <- h5file[[layer]]

  # Check matrix encoding and dimensions
  if (counts$attr_exists("encoding-type")) {
    encoding <- counts$attr_open("encoding-type")$read()

    if (encoding %in% c("csr_matrix", "csc_matrix")) {
      # SPARSE MATRIX (CSR or CSC format)
      cat("Reading sparse", encoding, "...\n")
      data <- counts[["data"]]$read()
      indices <- counts[["indices"]]$read()
      indptr <- counts[["indptr"]]$read()
      shape <- counts$attr_open("shape")$read()

      # Create sparse matrix (R uses 1-based indexing)
      # CSR uses 'i' parameter, CSC uses 'j' parameter
      if (encoding == "csr_matrix") {
        counts_matrix <- sparseMatrix(
          i = indices + 1,
          p = indptr,
          x = data,
          dims = shape,
          index1 = FALSE
        )
      } else {  # csc_matrix
        counts_matrix <- sparseMatrix(
          j = indices + 1,
          p = indptr,
          x = data,
          dims = shape,
          index1 = FALSE
        )
      }

    } else {
      stop("Unsupported sparse matrix encoding: ", encoding)
    }
  } else {
    # DENSE MATRIX - read directly
    cat("Reading dense matrix...\n")
    dense_matrix <- counts$read()
    counts_matrix <- as(dense_matrix, "sparseMatrix")
  }

  return(counts_matrix)
}

# STEP 3: Efficient metadata reading with error handling
read_metadata_safely <- function(h5file, path, default_names = NULL) {
  if (!h5file$exists(path)) {
    if (!is.null(default_names)) {
      return(default_names)
    }
    return(NULL)
  }

  tryCatch({
    h5file[[path]]$read()
  }, error = function(e) {
    warning("Failed to read metadata from ", path, ": ", e$message)
    return(default_names)
  })
}

read_dataframe_from_h5 <- function(h5file, group_path, index_names) {
  df <- DataFrame(row.names = index_names)
  
  if (!h5file$exists(group_path)) {
    return(df)
  }
  
  group <- h5file[[group_path]]
  column_names <- names(group)
  
  for (col_name in column_names) {
    if (col_name == "_index" || col_name == "__categories") {
      next  # Skip index and category metadata
    }
    
    tryCatch({
      col_data <- group[[col_name]]$read()
      
      # Handle different data types
      if (is.character(col_data) || is.factor(col_data)) {
        df[[col_name]] <- col_data
      } else if (is.numeric(col_data)) {
        df[[col_name]] <- col_data
      } else if (is.logical(col_data)) {
        df[[col_name]] <- col_data
      } else {
        # Try to convert to character for unknown types
        df[[col_name]] <- as.character(col_data)
      }
    }, error = function(e) {
      warning("Could not read column '", col_name, "': ", e$message)
    })
  }
  
  return(df)
}

# STEP 4: Spatial coordinates extraction
read_spatial_coords <- function(h5file, mod_location, cell_names) {
  # Try different common spatial coordinate locations
  spatial_paths <- c(
    paste0(mod_location, "/obsm/spatial"),
    paste0(mod_location, "/obsm/X_spatial"),
    paste0(mod_location, "/obsm/coordinates")
  )
  
  for (spatial_path in spatial_paths) {
    if (h5file$exists(spatial_path)) {
      tryCatch({
        coords <- h5file[[spatial_path]]$read()
        
        # Ensure it's a matrix with 2+ columns
        if (is.matrix(coords) && ncol(coords) >= 2) {
          # Use first two columns as x, y coordinates
          spatial_coords <- coords[, 1:2, drop = FALSE]
          colnames(spatial_coords) <- c("x", "y")
          
          if (nrow(spatial_coords) == length(cell_names)) {
            rownames(spatial_coords) <- cell_names
            cat("Found spatial coordinates at:", spatial_path, "\n")
            return(spatial_coords)
          }
        }
      }, error = function(e) {
        warning("Could not read spatial coordinates from ", spatial_path)
      })
    }
  }
  
  cat("No spatial coordinates found\n")
  return(NULL)
}

# STEP 5: Main conversion function
h5mu_to_spe <- function(h5mu_path, modality_name, assay_name = "counts", 
                       layer_name = NULL) {
  
  cat("Converting H5MU to SpatialExperiment...\n")
  cat("Input file:", h5mu_path, "\n")
  cat("Modality:", modality_name, "\n")
  if (!is.null(layer_name)) {
    cat("Layer:", layer_name, "\n")
  }
  
  # Open H5 file
  h5file <- H5File$new(h5mu_path, "r")
  mod_location <- paste("mod", modality_name, sep = "/")
  
  tryCatch({
    # Verify modality exists
    if (!h5file$exists(mod_location)) {
      available_mods <- h5file[["mod"]]$ls()$name
      stop("Modality '", modality_name, "' not found. Available modalities: ", 
           paste(available_mods, collapse = ", "))
    }
    
    cat("Reading counts matrix...\n")
    # Read counts matrix with layer support
    counts_matrix <- read_counts_matrix(h5file, mod_location, layer_name)
    
    # Transpose if necessary (H5MU typically stores as cells x genes)
    if (ncol(counts_matrix) > nrow(counts_matrix)) {
      cat("Transposing matrix (genes x cells)...\n")
      counts_matrix <- t(counts_matrix)
    }
    
    cat("Matrix dimensions:", nrow(counts_matrix), "genes x", ncol(counts_matrix), "cells\n")
    
    # Read feature and cell names
    cat("Reading feature and cell names...\n")
    gene_names <- read_metadata_safely(
      h5file, 
      paste0(mod_location, "/var/_index"),
      paste0("Gene_", seq_len(nrow(counts_matrix)))
    )
    
    cell_names <- read_metadata_safely(
      h5file,
      paste0(mod_location, "/obs/_index"), 
      paste0("Cell_", seq_len(ncol(counts_matrix)))
    )
    
    # Set matrix names
    rownames(counts_matrix) <- gene_names
    colnames(counts_matrix) <- cell_names
    
    # Read metadata
    cat("Reading cell metadata...\n")
    col_data <- read_dataframe_from_h5(h5file, paste0(mod_location, "/obs"), cell_names)
    
    cat("Reading gene metadata...\n")  
    row_data <- read_dataframe_from_h5(h5file, paste0(mod_location, "/var"), gene_names)
    
    # Read spatial coordinates
    cat("Reading spatial coordinates...\n")
    spatial_coords <- read_spatial_coords(h5file, mod_location, cell_names)
    
    # Create assay list
    assay_list <- list()
    assay_list[[assay_name]] <- counts_matrix
    
    # Create SpatialExperiment object
    cat("Creating SpatialExperiment object...\n")
    spe <- SpatialExperiment(
      assays = assay_list,
      spatialCoords = spatial_coords,
      colData = col_data,
      rowData = row_data
    )
    
    cat("Conversion completed successfully!\n")
    return(spe)
    
  }, finally = {
    h5file$close()
  })
}

# STEP 6: Main execution
main <- function() {
  # Validate required inputs
  if (is.null(par$modality) || par$modality == "") {
    stop("'modality' argument must be specified for H5MU files")
  }
  
  # Convert H5MU to SpatialExperiment
  spe <- h5mu_to_spe(
    h5mu_path = par$input,
    modality_name = par$modality, 
    assay_name = par$assay,
    layer_name = par$layer
  )
  
  # Save with compression for large files
  cat("Saving SpatialExperiment object...\n")
  if (!is.null(par$compression_level) && par$compression_level > 0) {
    saveRDS(spe, file = par$output, compress = "gzip")
  } else {
    saveRDS(spe, file = par$output, compress = TRUE)
  }
  
  cat("Output saved to:", par$output, "\n")
  cat("File size:", file.size(par$output), "bytes\n")
}

# Execute main function
main()


# ### VIASH START
# par <- list(
#   input = "resources_test/xenium/xenium_tiny",
#   add_experiment_xenium = TRUE,
#   add_parquet_paths = TRUE,
#   alternative_experiment_features = c("NegControlProbe", "UnassignedCodeword", "NegControlCodeword", "antisense", "BLANK"),
#   output = "spe_test.rds"
# )
# ### VIASH END

# obj <- readRDS(file = par$output)


# spe <- readXeniumSXE(
#   dirName = par$input,
#   returnType = "SPE",
#   countMatPattern = "cell_feature_matrix.h5",
#   metaDataPattern = "cells.parquet",
#   coordNames = c("x_centroid", "y_centroid"),
#   addExperimentXenium = par$add_experiment_xenium,
#   addParquetPaths = par$add_parquet_paths,
#   altExps = par$alternative_experiment_features
# )

# saveRDS(spe, file = par$output)

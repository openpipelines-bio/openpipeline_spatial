# openpipeline_spatial 0.2.0

## NEW FUNCTIONALITY

* `neighbors/spatial_neighborhood_graph`: Calculate the spatial neighborhood graph (PR #29).

* `convert/from_spaceranger_to_h5mu`: Added converter component for convert Spaceranger output to H5MU files (PR #33).

* `workflows/ingestion/spaceranger_mapping`: Added a workflow to ingest Visium data using Spaceranger and convert the count matrix to an H5MU file (PR #33).

## MINOR CHANGES

* Add `scope` to component and workflow configurations (PR #22).

* Bump version of spatialdata-io to 0.3.0 and spatialdata to 0.5.0. Pin version of pyarrow to 18.0.0 for compatibility (PR #24).

* `convert/from_xenium_to_spatialexperiment`: Add arrow with zstd codec support to handle I/O of zstd-compressed Xenium parquet files (PR #30).

* `mapping/spaceranger_count`: Allow providing individual FASTQ files instead of directories (PR #32).

* Bump anndata to 0.12.7 and mudata to 0.3.2 (PR #34).

* Bump spatialdata to 0.6.1 and spatialdata-io to 0.5.1 (PR #34).

## BUG FIXES

* `convert/from_cosmx_to_h5mu`: Fixed an issue where parent directories of the cosmx output bundle were duplicated when reading in data (PR #25).

* `mapping/spaceranger_count`: Fixed issue with long temporary folder paths causing write failures (PR #31).

# openpipeline_spatial 0.1.1

## MINOR CHANGES

* Add a README (PR #21).

## NEW FUNCTIONALITY

* `convert`: Updated multiple components to accept spatial output bundles in .zip format (for CosMx, Xenium and Aviti) as input (PR #19, PR #20).

* `convert/from_cosmx_to_h5mu`: Updated component to handle CosMx output bundles generated with AtoMx SIP versions < v1.3.2 (PR #25).

# openpipeline_spatial 0.1.0

## NEW FUNCTIONALITY

* `filter/subset_cosmx`: Added a component to subset COSMX data (PR #3, PR #9).

* `convert/from_cosmx_to_h5mu`: Added converter component for COSMX data (PR #3, PR #9).

* `mapping/spaceranger_count`: Added a spaceranger count component (PR #2).

* `convert/from_spatialdata_to_h5mu`, `convert/from_xenium_to_spatialdata`, `convert/from_xenium_to_h5mu`: Added converter components for xenium data (PR #1, #10).

* `convert/from_xenium_to_spatialexperiment`, `convert/from_cosmx_to_spatialexperiment`: Added converter components for Xenium or CosMx data to SpatialExperiment objects (PR #9).

* `convert/from_cells2stats_to_h5mu`: Added a component to convert data resulting from Aviti Teton sequencers processed by Cells2Stats into an H5MU file (PR #15).

* `workflows/qc/qc`: Added a pipeline for calculating qc metrics of spatial omics samples (PR #5).

* `workflows/multiomics/spatial_process_samples`: Added a pipeline to pre-process multiple spatial omics samples (PR #7).

* `convert/from_h5mu_to_spatialexperiment`: Added converter component for H5MU data to SpatialExperiment objects (PR #15).

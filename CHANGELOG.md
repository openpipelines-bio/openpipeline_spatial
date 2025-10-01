# openpipeline_spatial x.x.x

## MINOR CHANGES

* Add a README (PR #21).

## NEW FUNCTIONALITY

* `convert`: Updated multiple components to accept spatial output bundles in .zip format (for CosMx, Xenium and Aviti) as input (PR #19).

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

viash ns build -q spatial_sample_processing --setup cb

nextflow run . \
-resume \
-profile docker \
-c src/labels_ci.config \
-main-script target/nextflow/workflows/multiomics/spatial_sample_processing/main.nf \
--input resources_test/xenium/xenium_tiny.h5mu \
--output xenium_tiny_processed.h5mu \
--publish_dir resources_test/xenium/
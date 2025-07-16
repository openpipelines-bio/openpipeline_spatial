viash ns build -q spatial_process_samples --setup cb

nextflow run . \
-resume \
-profile docker \
-c src/workflows/utils/labels_ci.config \
-main-script target/nextflow/workflows/multiomics/spatial_process_samples/main.nf \
--input resources_test/xenium/xenium_tiny.h5mu \
--output xenium_tiny_processed.h5mu \
--publish_dir resources_test/xenium/
viash ns build -q spatial_qc_report --setup b

nextflow run . \
-resume \
-profile docker \
-c src/labels_ci.config \
-main-script target/nextflow/workflows/spatial_qc/spatial_qc_report/main.nf \
--input resources_test/xenium/xenium_tiny.h5mu \
--output_processed_h5mu test.h5mu \
--output_qc_report test.html \
--var_name_mitochondrial_genes mitochondrial \
--var_name_ribosomal_genes ribosomal \
--publish_dir test
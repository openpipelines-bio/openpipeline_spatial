viash ns build -q spatial_qc/spatial_qc --setup b

nextflow run . \
-resume \
-profile docker \
-c src/labels_ci.config \
-main-script target/nextflow/workflows/spatial_qc/spatial_qc/main.nf \
--input resources_test/xenium/xenium_tiny.h5mu \
--output resources_test/xenium/xenium_tiny_qc.h5mu \
--var_name_mitochondrial_genes mitochondrial \
--var_name_ribosomal_genes ribosomal

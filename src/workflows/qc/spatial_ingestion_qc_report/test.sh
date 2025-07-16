aws s3 sync s3://openpipelines-bio/openpipeline_spatial/resources_test/xenium resources_test/xenium

viash ns build -q spatial_ingestion_qc_report --setup cb

nextflow run . \
-resume \
-profile docker \
-c src/workflows/utils/labels_ci.config \
-main-script target/nextflow/workflows/qc/spatial_ingestion_qc_report/main.nf \
--input resources_test/xenium/xenium_tiny.h5mu \
--ingestion_method xenium \
--output_processed_h5mu test.h5mu \
--output_qc_report test.html \
--var_name_mitochondrial_genes mitochondrial \
--var_name_ribosomal_genes ribosomal \
--publish_dir test
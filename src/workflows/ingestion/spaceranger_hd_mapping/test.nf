nextflow.enable.dsl=2

include { spaceranger_hd_mapping } from params.rootDir + "/target/nextflow/workflows/ingestion/spaceranger_hd_mapping/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  // Human Visium HD FFPE dataset from nf-core/test-datasets (spatialvi). Its reads
  // are curated to match a tiny GRCh38 reference, so Space Ranger produces
  // non-empty counts and its HD pipeline (including cell segmentation) runs to
  // completion -- exercising the full spaceranger_count -> converter chain.
  hd = resources_test.resolve("visium_hd_ffpe")
  prefix = "Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2"

  output_ch = Channel.fromList([
      [
        id: "foo",
        input: hd.resolve("${prefix}_fastqs"),
        gex_reference: hd.resolve("GRCh38"),
        probe_set: hd.resolve("${prefix}_probe_set.csv"),
        cytaimage: hd.resolve("${prefix}_cytaimage.tif"),
        image: hd.resolve("${prefix}_image.btf"),
        slide: "H1-84QJZFR",
        area: "D1",
        create_bam: false,
        // Human sample -> the bundled Pan-Human Azimuth model runs locally (no
        // 10x Cloud token), so cell-type annotation is exercised end-to-end.
        cell_annotation_model: "human_pca_v1_beta",
        output_type: "filtered",
      ]
    ])
    | map{ state -> [state.id, state] }
    | spaceranger_hd_mapping
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, out]"
      assert output[1] instanceof Map : "Output should be a Map."
      assert output[1].containsKey("output_spatialdata") : "Output should contain output_spatialdata."
      assert file(output[1].output_spatialdata).isDirectory() : "output_spatialdata should be a SpatialData Zarr store."
      "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
}

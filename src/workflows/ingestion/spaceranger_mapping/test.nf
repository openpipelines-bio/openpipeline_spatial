nextflow.enable.dsl=2

include { spaceranger_mapping } from params.rootDir + "/target/nextflow/workflows/ingestion/spaceranger_mapping/main.nf"
include { cellranger_mapping_test } from params.rootDir + "/target/_test/nextflow/test_workflows/ingestion/cellranger_mapping_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [  
        id: "foo",
        input: resources_test.resolve("visium/subsampled"),
        gex_reference: resources_test.resolve("GRCh38"),
        probe_set: resources_test.resolve("visium/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv"),
        image: resources_test.resolve("visium/subsampled/Visium_FFPE_Human_Ovarian_Cancer_image.jpg"),
        unkown_slide: "visium-1",
        slide: "V10L13-020",
        area: "D1",
        output_type: "filtered",
      ]
    ])
    | map{ state -> [state.id, state] }
    | spaceranger_mapping
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, out]"
      assert output[1] instanceof Map : "Output should be a Map."
      "Output: $output"
    }

    | cellranger_mapping_test.run(
      fromState: ["input": "output_h5mu"]
    )

    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
}

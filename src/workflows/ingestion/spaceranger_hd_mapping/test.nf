nextflow.enable.dsl=2

include { spaceranger_hd_mapping } from params.rootDir + "/target/nextflow/workflows/ingestion/spaceranger_hd_mapping/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  // Feed a pre-computed Space Ranger HD output (real 10x segmented output) so the
  // workflow's conversion step is exercised without rerunning Space Ranger, which
  // needs a full reference + real data that are too large for CI. The live
  // Space Ranger alignment path is covered separately (spaceranger_count tests).
  output_ch = Channel.fromList([
      [
        id: "foo",
        spaceranger_output: resources_test.resolve("visium_hd/Visium_HD_Mouse_Brain_tiny_spaceranger"),
        dataset_id: "Visium_HD_Mouse_Brain",
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

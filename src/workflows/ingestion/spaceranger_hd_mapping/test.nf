nextflow.enable.dsl=2

include { spaceranger_hd_mapping } from params.rootDir + "/target/nextflow/workflows/ingestion/spaceranger_hd_mapping/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "foo",
        input: resources_test.resolve("visium_hd/Visium_HD_Mouse_Brain_tiny"),
        gex_reference: resources_test.resolve("mm10"),
        probe_set: resources_test.resolve("visium_hd/probe_set.csv"),
        cytaimage: resources_test.resolve("visium_hd/Visium_HD_Mouse_Brain_cytassist_tiny.tiff"),
        image: resources_test.resolve("visium_hd/Visium_HD_Mouse_Brain_image_tiny.jpg"),
        create_bam: "false",
        // Skip secondary analysis (per-bin clustering metrics come out null on
        // this tiny fixture and break the HD web-summary builder) and cell
        // annotation (requires a 10x Cloud token).
        nosecondary: true,
        disable_cell_annotation: true,
        output_type: "filtered",
      ]
    ])
    | map{ state -> [state.id, state] }
    | spaceranger_hd_mapping
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, out]"
      assert output[1] instanceof Map : "Output should be a Map."
      assert output[1].containsKey("output_spatialdata") : "Output should contain output_spatialdata."
      "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
}

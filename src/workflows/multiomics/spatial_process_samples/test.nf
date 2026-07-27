nextflow.enable.dsl=2
targetDir = params.rootDir + "/target/nextflow"

include { spatial_process_samples } from targetDir + "/workflows/multiomics/spatial_process_samples/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
    [
      id: "xenium",
      input: resources_test.resolve("xenium/xenium_tiny.h5mu"),
      publish_dir: "foo/",
      output: "test.h5mu",
    ]
  ])
  | map{ state -> [state.id, state] }
  | spatial_process_samples
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith("test.h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    assert output_list.size() == 1 : "output channel should contain one event"
    assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
  }
  
}
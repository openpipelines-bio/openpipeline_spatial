nextflow.enable.dsl=2
targetDir = params.rootDir + "/target/nextflow"

include { spatial_leiden_domains } from targetDir + "/workflows/domains/spatial_leiden_domains/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
    [
      id: "xenium",
      input: resources_test.resolve("xenium/xenium_tiny.qc.all_neighbors.pca.h5mu"),
      output: "output.h5mu",
      device_type: "cpu",
      resolution: [0.5, 1.0],
    ]
  ])
  | map { state -> [state.id, state] }
  | spatial_leiden_domains
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, state]"

    def id = output[0]
    assert id == "xenium"

    def state = output[1]
    assert state instanceof Map : "State should be a map. Found: ${state}"
    assert state.containsKey("output") : "Output should contain key 'output'."
    assert state.output.isFile() : "'output' should be a file."
    assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    assert output_list.size() == 1 : "output channel should contain one event"
  }
}

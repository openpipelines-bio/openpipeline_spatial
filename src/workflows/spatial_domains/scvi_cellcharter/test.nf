nextflow.enable.dsl=2
targetDir = params.rootDir + "/target/nextflow"

include { scvi_cellcharter } from targetDir + "/workflows/spatial_domains/scvi_cellcharter/main.nf"
include { scvi_cellcharter_test } from params.rootDir + "/target/_test/nextflow/test_workflows/spatial_domains/scvi_cellcharter_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
    [
      // two samples, generic coordinates with pruned Delaunay graph, automatic selection of the number of clusters
      id: "cosmx",
      input: resources_test.resolve("cosmx/Lung5_tiny_processed.h5mu"),
      // CosMx 'spatial' coordinates are relative to each FOV, 'spatial_fov' holds the global coordinates
      input_obsm_spatial_coords: "spatial_fov",
      scvi_max_epochs: 2,
      n_clusters_range: [2, 4],
      max_runs: 2,
      cluster_max_epochs: 10,
      expect_stability: true,
    ],
    [
      // grid coordinates, highly variable genes and a fixed number of clusters
      id: "visium",
      input: resources_test.resolve("visium/visium_tiny_processed.h5mu"),
      coord_type: "grid",
      var_input: "filter_with_hvg",
      scvi_max_epochs: 2,
      n_clusters: 3,
      cluster_max_epochs: 10,
      expect_stability: false,
      expected_n_clusters: 3,
    ],
  ])
  | map { state -> [state.id, state] }
  | scvi_cellcharter.run(
    toState: { id, output, state -> output + [
      n_dimensions_latent_space: 10,
      expect_stability: state.expect_stability,
      expected_n_clusters: state.expected_n_clusters,
    ]}
  )
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, state]"

    def id = output[0]
    assert id in ["cosmx", "visium"]

    def state = output[1]
    assert state instanceof Map : "State should be a map. Found: ${state}"
    assert state.containsKey("output") : "Output should contain key 'output'."
    assert state.output.isFile() : "'output' should be a file."
    assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"
    assert state.containsKey("output_scvi_model") : "Output should contain key 'output_scvi_model'."
    assert state.output_scvi_model.isDirectory() : "'output_scvi_model' should be a directory."

    "Output: $output"
  }
  | scvi_cellcharter_test.run(
    fromState: [
      "input": "output",
      "n_dimensions_latent_space": "n_dimensions_latent_space",
      "expected_n_clusters": "expected_n_clusters",
      "expect_stability": "expect_stability",
    ]
  )
  | toSortedList()
  | map { output_list ->
    assert output_list.size() == 2 : "output channel should contain two events"
  }
}

workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

    | map { id, state ->
      [id, state + [workflow_output: state.output]]
    }

    // Compute the spatial neighborhood graph from the spatial coordinates
    | spatial_neighborhood_graph.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_obsm_spatial_coords": "input_obsm_spatial_coords",
        "coord_type": "coord_type",
        "n_spatial_neighbors": "n_spatial_neighbors",
        "delaunay": "delaunay",
        "output_compression": "output_compression",
        "output": "workflow_output",
      ],
      toState: ["input": "output"]
    )

    // Fuse the expression and spatial neighborhood graphs
    | join_graphs.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "alpha": "alpha",
        "output_compression": "output_compression",
        "output": "workflow_output",
      ],
      args: [
        "input_obsp_expression_graph": "connectivities",
        "input_obsp_spatial_graph": "spatial_connectivities",
        "output_obsp_graph": "spatial_expression_connectivities",
      ],
      toState: ["input": "output"]
    )

    // Leiden clustering on the fused graph
    | leiden.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "resolution": "resolution",
        "n_iterations": "leiden_n_iterations",
        "seed": "leiden_seed",
        "obsm_name": "obsm_output",
        "output_compression": "output_compression",
        "output": "workflow_output",
      ],
      args: [
        "obsp_connectivities": "spatial_expression_connectivities",
      ],
      toState: ["output": "output"]
    )

    | setState(["output": "output"])

  emit:
    output_ch
}

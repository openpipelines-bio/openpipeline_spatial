workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

    | map { id, state ->
      [id, state + [workflow_output: state.output]]
    }

    // Compute the expression k-NN graph from the PCA embedding
    | expression_neighbors.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "input_obsm_pca",
        "device_type": "device_type",
        "num_neighbors": "expression_num_neighbors",
        "metric": "expression_metric",
        "random_state": "expression_random_state",
        "output_compression": "output_compression",
        "output": "workflow_output",
      ],
      toState: ["input": "output"]
    )

    // Compute the spatial neighborhood graph from the spatial coordinates
    | spatial_neighborhood_graph.run(
      fromState: { id, state -> [
        "input": state.input,
        "modality": state.modality,
        "input_obsm_spatial_coords": state.input_obsm_spatial_coords,
        "coord_type": state.coord_type ?: (state.technology in ["visium", "visium_hd"] ? "grid" : "generic"),
        "n_spatial_neighbors": state.n_spatial_neighbors,
        "delaunay": state.delaunay,
        "output_compression": state.output_compression,
        "output": state.workflow_output,
      ]},
      toState: ["input": "output"]
    )

    // Technology-specific spatial statistics
    | xenium_spatial_statistics.run(
      runIf: { id, state -> state.technology == "xenium" },
      fromState: [
        "input": "input",
        "modality": "modality",
        "output": "workflow_output",
      ],
      toState: ["input": "output"]
    )

    | visium_spatial_statistics.run(
      runIf: { id, state -> state.technology in ["visium", "visium_hd"] },
      fromState: { id, state -> [
        "input": state.input,
        "modality": state.modality,
        "tissue_edge_max_neighbors": state.technology == "visium_hd" ? 8 : 6,
        "output": state.workflow_output,
      ]},
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
        "device_type": "device_type",
        "resolution": "resolution",
        "n_iterations": "leiden_n_iterations",
        "random_state": "leiden_random_state",
        "obsm_name": "obsm_output",
        "output_compression": "output_compression",
        "output": "workflow_output",
      ],
      args: [
        "obsp_connectivities": "spatial_expression_connectivities",
      ],
      toState: ["input": "output"]
    )

    // Spatially variable gene detection via spatial autocorrelation
    | spatial_autocorr.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "device_type": "device_type",
        "mode": "svg_mode",
        "n_perms": "svg_n_perms",
        "output_compression": "output_compression",
        "output": "workflow_output",
      ],
      args: [
        "obsp_connectivities": "spatial_connectivities",
      ],
      toState: ["output": "output"]
    )

    | setState(["output": "output"])

  emit:
    output_ch
}

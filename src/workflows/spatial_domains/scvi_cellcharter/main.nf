workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

    | map { id, state ->
      // Validate the clustering arguments before training the scVI model
      if (!workflow.stubRun) {
        assert state.n_clusters != null || state.n_clusters_range :
          "At least one of --n_clusters or --n_clusters_range must be provided."
        assert !state.n_clusters_range || state.n_clusters_range.size() == 2 :
          "--n_clusters_range must contain exactly two values: the minimum and maximum number of clusters."
      }
      [id, state + [workflow_output: state.output]]
    }

    // Reduce the dimensionality of the gene expression with scVI
    | scvi.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "layer",
        "obs_batch": "obs_batch",
        "var_input": "var_input",
        "n_dimensions_latent_space": "n_dimensions_latent_space",
        "max_epochs": "scvi_max_epochs",
        "early_stopping": "early_stopping",
        "obsm_output": "obsm_scvi",
        "output_model": "output_scvi_model",
      ],
      toState: [
        "input": "output",
        "output_scvi_model": "output_model",
      ]
    )

    // Compute the spatial neighborhood graph per library, pruning long links
    | spatial_neighborhood_graph.run(
      fromState: { id, state -> [
        "input": state.input,
        "modality": state.modality,
        "input_obsm_spatial_coords": state.input_obsm_spatial_coords,
        "input_obs_library_key": state.input_obs_library_key,
        "coord_type": state.coord_type,
        "n_spatial_neighbors": state.n_spatial_neighbors,
        "delaunay": state.delaunay,
        // pruning long links is only supported for generic coordinates
        "percentile": state.coord_type == "generic" ? state.percentile : null,
      ]},
      toState: ["input": "output"]
    )

    // Aggregate the scVI embedding over the spatial neighborhood of every cell
    | aggregate_neighbors.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_obsm_features": "obsm_scvi",
        "input_obs_sample_key": "input_obs_library_key",
        "n_layers": "n_layers",
        "aggregations": "aggregations",
        "output_obsm_aggregated": "obsm_cellcharter",
      ],
      args: [
        "input_obsp_spatial_connectivities": "spatial_connectivities",
      ],
      toState: ["input": "output"]
    )

    // Cluster the aggregated features into spatial domains
    | cluster.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_obsm_features": "obsm_cellcharter",
        "n_clusters": "n_clusters",
        "n_clusters_range": "n_clusters_range",
        "max_runs": "max_runs",
        "stability_convergence_tolerance": "stability_convergence_tolerance",
        "covariance_type": "covariance_type",
        "max_epochs": "cluster_max_epochs",
        "batch_size": "cluster_batch_size",
        "accelerator": "accelerator",
        "random_state": "cluster_random_state",
        "output_obs_cluster": "obs_cluster",
        "output_uns_stability": "uns_stability",
        "output_compression": "output_compression",
        "output": "workflow_output",
      ],
      toState: ["output": "output"]
    )

    | setState(["output", "output_scvi_model"])

  emit:
    output_ch
}

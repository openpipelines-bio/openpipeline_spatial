workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

    | map { id, state ->
      [id, state + [
        workflow_output: state.output
      ]]
    }
    // Compute the spatial neighborhood graph per sample
    | spatial_neighborhood_graph.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_obsm_spatial_coords": "input_obsm_spatial_coords",
        "input_obs_library_key": "input_obs_batch_key",
        "coord_type": "coord_type",
        "n_spatial_neighbors": "n_spatial_neighbors",
        "delaunay": "delaunay"
      ],
      toState: ["input": "output"]
    )

    | nichecompass.run(
      fromState: {id, state -> [
        "input": state.input,
        "input_gp_mask": state.input_gp_mask,
        "input_obs_covariates": state.input_obs_covariates,
        "var_input": state.var_input,
        "modality": state.modality,
        "layer": state.layer,
        "min_genes_per_gp": state.min_genes_per_gp,
        "min_source_genes_per_gp": state.min_source_genes_per_gp,
        "min_target_genes_per_gp": state.min_target_genes_per_gp,
        "max_genes_per_gp": state.max_genes_per_gp,
        "max_source_genes_per_gp": state.max_source_genes_per_gp,
        "max_target_genes_per_gp": state.max_target_genes_per_gp,
        "filter_genes_not_in_masks": state.filter_genes_not_in_masks,
        "covariate_edges": state.covariate_edges,
        "gene_expr_recon_dist": state.gene_expr_recon_dist,
        "log_variational": state.log_variational,
        "node_label_method": state.node_label_method,
        "active_gp_thresh_ratio": state.active_gp_thresh_ratio,
        "active_gp_type": state.active_gp_type,
        "n_addon_gp": state.n_addon_gp,
        "cat_covariates_embeds_nums": state.cat_covariates_embeds_nums,
        "random_state": state.random_state,
        "n_epochs": state.n_epochs,
        "n_epochs_all_gps": state.n_epochs_all_gps,
        "n_epochs_no_edge_recon": state.n_epochs_no_edge_recon,
        "n_epochs_no_cat_covariates_contrastive": state.n_epochs_no_cat_covariates_contrastive,
        "lr": state.lr,
        "weight_decay": state.weight_decay,
        "edge_val_ratio": state.edge_val_ratio,
        "node_val_ratio": state.node_val_ratio,
        "edge_batch_size": state.edge_batch_size,
        "node_batch_size": state.node_batch_size,
        "n_sampled_neighbors": state.n_sampled_neighbors,
        "output_obsm_embedding": state.output_obsm_embedding,
        "output_model": state.output_model
      ]},
      args: [
        "input_obsp_spatial_connectivities": "spatial_connectivities"
      ],
      toState: [
        "input": "output",
        "output_model": "output_model"
      ]
    )

    | neighbors_leiden_umap.run(
      fromState: { id, state -> [
        "input": state.input,
        "modality": state.modality,
        "obsm_input": state.output_obsm_embedding,
        "output": state.workflow_output,
        "uns_neighbors": state.uns_neighbors,
        "obsp_neighbor_distances": state.obsp_neighbor_distances,
        "obsp_neighbor_connectivities": state.obsp_neighbor_connectivities,
        "leiden_resolution": state.leiden_resolution,
        "obs_cluster": state.obs_cluster,
        "obsm_umap": state.obsm_umap,
      ]},
      toState: ["output": "output"]
    )

    | setState(["output": "output", "output_model": "output_model"])

  emit:
    output_ch
}

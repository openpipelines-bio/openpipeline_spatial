workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

    | map { id, state ->
      [id, state + [
        workflow_output: state.output,
        _meta: [join_id: id]
      ]]
    }
    // If requested, add the id of the events (samples) to a column in .obs. 
    // Also allows to make .obs_names (the .obs index) unique, by prefixing the values with an unique id per .h5mu file.
    // The latter is usefull to avoid duplicate observations during concatenation.
    | add_id.run(
      filter: {id, state -> state.add_id_to_obs },
      fromState: {id, state -> 
        def newState = [
          "input": state.input,
          "input_id": id,
          "make_observation_keys_unique": state.add_id_make_observation_keys_unique,
          "obs_output": state.add_id_obs_output,
          "add_id_to_obs": state.add_id_to_obs
        ]
        newState
      },
      toState: {id, output, state ->
        def keysToRemove = ["add_id_to_obs", "add_id_obs_output", "add_id_make_observation_keys_unique"]
        def newState = state.findAll{it.key !in keysToRemove}
        newState + ["input": output.output]
      }
    )

    | spatial_neighborhood_graph.run(
      fromState: {id, state -> [
        "input": state.input,
        "modality": state.modality,
        "layer": state.layer,
        "input_obsm_spatial_coords": state.input_obsm_spatial_coords,
        "coord_type": state.coord_type,
        "n_spatial_neighbors": state.n_spatial_neighbors,
        "delaunay": state.delaunay
      ]},
      toState: {id, output, state -> 
          def keysToRemove = ["input_obsm_spatial_coords", "coord_type", "n_spatial_neighbors", "delaunay"]
          def newState = state.findAll{it.key !in keysToRemove}
          newState + ["input": output.output]
        }
    )

    | joinStates { ids, states ->
      def newId = "merged"
      // gather keys with unique values across states that should be combined
      def new_state_non_unique_values = [
        input: states.collect{it.input},
        input_id: ids,
        _meta: [join_id: ids[0]]
      ]
      // gather keys from different states
      def all_state_keys = states.inject([].toSet()){ current_keys, state ->
          def new_keys = current_keys + state.keySet()
          return new_keys
      }.minus(["output", "id", "input", "_meta"])
      // Create the new state from the keys, values should be the same across samples
      def new_state = all_state_keys.inject([:]){ old_state, argument_name ->
          argument_values = states.collect{it.get(argument_name)}.unique()
          assert argument_values.size() == 1, "Arguments should be the same across samples. Argument name: $argument_name, \
                                                argument value: $argument_values"
          // take the unique value from the set (there is only one)
          def argument_value
          argument_values.each { argument_value = it }
          def current_state = old_state + [(argument_name): argument_value]
          return current_state
      }
      def data_state = new_state_non_unique_values + new_state
      [ newId, data_state ]
    }

    | obsp_block_concatenation.run(
      fromState: { id, state -> [
        "input": state.input,
        "modality": state.modality,
        "input_id": state.input_id
      ]},
      toState: {id, output, state -> 
          def keysToRemove = ["input_id"]
          def newState = state.findAll{it.key !in keysToRemove}
          newState + ["input": output.output]
        }
    )

    | nichecompass.run(
      fromState: {id, state -> [
        "input": state.input,
        "input_gp_mask": state.input_gp_mask,
        "input_obs_covariates": state.input_obs_covariates,
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
        "gene_expr_recon_distribution": state.gene_expr_recon_dist,
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
        "n_epochs_no_cat_covariates_contrastive_loss": state.n_epochs_no_cat_covariates_contrastive_loss,
        "lr": state.lr,
        "weight_decay": state.weight_decay,
        "edge_val_ratio": state.edge_val_ratio,
        "node_val_ratio": state.node_val_ratio,
        "edge_batch_size": state.edge_batch_size,
        "node_batch_size": state.node_batch_size,
        "n_sampled_neighbors": state.n_sampled_neighbors,
        "output_model": state.output_model
      ]},
      args: [
        "input_obsm_spatial_connectivities": "spatial_connectivities"
      ],
      toState: [
        "input": "output"
      ]
    )

    | leiden.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_name": "obs_cluster",
        "resolution": "leiden_resolution"
      ],
      args: [
        "obsp_connectivities": "spatial_connectivities"
      ],
      toState: ["input": "output"]
    )

    | move_obsm_to_obs.run(
      runIf: {id, state -> state.leiden_resolution},
      fromState: [
        "input": "input",
        "obsm_key": "obs_cluster",
        "modality": "modality",
      ],
      toState: ["input": "output"]
    )
    | umap.run(
      runIf: {id, state -> !state.obsm_umap?.trim()?.isEmpty()},
      fromState: [
          "input": "input",
          "output": "workflow_output",
          "uns_neighbors": "uns_neighbors",
          "obsm_output": "obsm_umap",
          "modality": "modality",
        ],
      args: ["output_compression": "gzip"],
      toState: ["output": "output"]
    )
    | setState(["output": "output", "output_model": "output_model"])
    | view()

  emit:
    output_ch
}
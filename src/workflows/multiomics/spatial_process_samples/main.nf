workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      | map { id, state ->
        def newEvent = [state.id, state + ["_meta": ["join_id": id]]]
        newEvent
      }
      | spatial_sample_processing.run(
        fromState: { id, state -> [
          "id": id,
          "input": state.input,
          "rna_layer": state.rna_layer,
          "prot_layer": state.prot_layer,
          "add_id_to_obs": state.add_id_to_obs,
          "add_id_obs_output": state.add_id_obs_output,
          "add_id_make_observation_keys_unique": state.add_id_make_observation_keys_unique,
          "rna_min_counts": state.rna_min_counts,
          "rna_max_counts": state.rna_max_counts,
          "rna_min_genes_per_cell": state.rna_min_genes_per_cell,
          "rna_max_genes_per_cell": state.rna_max_genes_per_cell,
          "rna_min_cells_per_gene": state.rna_min_cells_per_gene,
          "rna_min_fraction_mito": state.rna_min_fraction_mito,
          "rna_max_fraction_mito": state.rna_max_fraction_mito,
          "rna_min_fraction_ribo": state.rna_min_fraction_ribo,
          "rna_max_fraction_ribo": state.rna_max_fraction_ribo,
          "prot_min_counts": state.prot_min_counts,
          "prot_max_counts": state.prot_max_counts,
          "prot_min_proteins_per_cell": state.prot_min_proteins_per_cell,
          "prot_max_proteins_per_cell": state.prot_max_proteins_per_cell,
          "prot_min_cells_per_protein": state.prot_min_cells_per_protein,
          "highly_variable_features_var_output": state.highly_variable_features_var_output,
          "highly_variable_features_obs_batch_key": state.highly_variable_features_obs_batch_key,
          "var_gene_names": state.var_gene_names,
          "var_name_mitochondrial_genes": state.var_name_mitochondrial_genes,
          "obs_name_mitochondrial_fraction": state.obs_name_mitochondrial_fraction,
          "mitochondrial_gene_regex": state.mitochondrial_gene_regex,
          "var_name_ribosomal_genes": state.var_name_ribosomal_genes,
          "obs_name_ribosomal_fraction": state.obs_name_ribosomal_fraction,
          "ribosomal_gene_regex": state.ribosomal_gene_regex,
          "var_qc_metrics": state.var_qc_metrics,
          "top_n_vars": state.top_n_vars,
          "pca_overwrite": state.pca_overwrite,
          "clr_axis": state.clr_axis,
          "rna_enable_scaling": state.rna_enable_scaling,
          "rna_scaling_output_layer": state.rna_scaling_output_layer,
          "rna_scaling_pca_obsm_output": state.rna_scaling_pca_obsm_output,
          "rna_scaling_pca_loadings_varm_output": state.rna_scaling_pca_loadings_varm_output,
          "rna_scaling_pca_variance_uns_output": state.rna_scaling_pca_variance_uns_output,
          "rna_scaling_umap_obsm_output": state.rna_scaling_umap_obsm_output,
          "rna_scaling_max_value": state.rna_scaling_max_value,
          "rna_scaling_zero_center": state.rna_scaling_zero_center,
        ]},
        args: [
          "skip_scrublet_filtering": "true",
        ],
        toState: [
          "output": "output"
        ]
      )

      | setState( 
        [
          "_meta": "_meta",
          "output": "output"
        ]
      )

      // | map {combined_id, state ->
      //   def resultState = [
      //     "output": state.output,
      //     // The join ID is the same across all samples from the same run
      //     "_meta": ["join_id": state._meta.join_id]
      //   ]
      //   return [combined_id, resultState]
      // }

  emit:
    output_ch
}

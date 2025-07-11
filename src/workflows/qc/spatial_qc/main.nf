workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      | spatial_qc_workflow.run(
        fromState: { id, state -> [
          "id": id,
          "input": state.input,
          "modality": state.modality,
          "layer": state.layer,
          "var_gene_names": state.var_gene_names,
          "var_name_mitochondrial_genes": state.var_name_mitochondrial_genes,
          "obs_name_mitochondrial_fraction": state.obs_name_mitochondrial_fraction,
          "mitochondrial_gene_regex": state.mitochondrial_gene_regex,
          "var_name_ribosomal_genes": state.var_name_ribosomal_genes,
          "obs_name_ribosomal_fraction": state.obs_name_ribosomal_fraction,
          "ribosomal_gene_regex": state.ribosomal_gene_regex,
          "var_qc_metrics": state.var_qc_metrics,
          "top_n_vars": state.top_n_vars,
          "output_obs_num_nonzero_vars": state.output_obs_num_nonzero_vars,
          "output_obs_total_counts_vars": state.output_obs_total_counts_vars,
          "output_var_num_nonzero_obs": state.output_var_num_nonzero_obs,
          "output_var_total_counts_obs": state.output_var_total_counts_obs,
          "output_var_obs_mean": state.output_var_obs_mean,
          "output_var_pct_dropout": state.output_var_pct_dropout
        ]},
        toState: [
          "output": "output"
        ]
      )

      | setState(["output"])
    
  emit:
    output_ch
}
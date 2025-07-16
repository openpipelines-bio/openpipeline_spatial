workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
    // store join id
      | map { id, state ->
        [id, state + [_meta: [join_id: id]]]
      }
      | spatial_qc_report.run(
        fromState: { id, state -> [
          "id": id,
          "input": state.input,
          "sample_metadata": state.sample_metadata,
          "max_samples_per_report": state.max_samples_per_report,
          "var_gene_names": state.var_gene_names,
          "obs_metadata": state.obs_metadata,
          "var_name_mitochondrial_genes": state.var_name_mitochondrial_genes,
          "var_name_ribosomal_genes": state.var_name_ribosomal_genes,

        ]},
        args: [
            "ingestion_method": "xenium",
            "run_cellbender": false
        ],
        toState: {id, output, state -> [
          "output_processed_h5mu": output.output_processed_h5mu,
          "output_qc_report": output.output_qc_report,
          "_meta": state._meta
        ]}
      )

      | setState(["output_processed_h5mu", "output_qc_report"])
      | niceView()
    
  emit:
    output_ch
}
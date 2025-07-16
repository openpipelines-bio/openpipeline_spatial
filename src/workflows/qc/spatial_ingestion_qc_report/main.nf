workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      | map { id, state ->
        def new_state = [
          state.id, 
          state + ["_meta": ["join_id": id]]
        ]
        new_state
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
          "output_processed_h5mu": state.output_processed_h5mu,
          "output_qc_report": state.output_qc_report

        ]},
        args: [
            "ingestion_method": "xenium",
            "run_cellbender": false
        ],
        toState: {id, output, state -> [
          "output_processed_h5mu": output.output_processed_h5mu,
          "output_qc_report": output.output_qc_report
        ]}
      )

      | setState( 
        [
          "_meta": "_meta",
          "output_processed_h5mu": "output_processed_h5mu",
          "output_qc_report": "output_qc_report"
        ]
      )
    
  emit:
    output_ch
}
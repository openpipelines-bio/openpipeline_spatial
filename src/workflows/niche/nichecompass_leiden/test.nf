nextflow.enable.dsl=2

include { nichecompass_leiden } from params.rootDir + "/target/nextflow/workflows/niche/nichecompass_leiden/main.nf"
include { nichecompass_leiden_test } from params.rootDir + "/target/_test/nextflow/test_workflows/niche/nichecompass_leiden_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = 
    Channel.fromList([
      [
        id: "test",
        input: resources_test.resolve("cosmx/Lung5_tiny_processed.h5mu"),
        input_gp_mask: resources_test.resolve("niche/prior_knowledge_gp_mask.json"),
        n_epochs: 1,
        n_epochs_all_gps: 0,
        n_epochs_no_edge_recon: 0,
        n_epochs_no_cat_covariates_contrastive_loss: 0,
        output_model: "simple_execution_test_model"
      ]
    ])
    | map { state -> [state.id, state] }
    | nichecompass_leiden.run(
      toState: { id, output, state -> output + [og_input: state.input] }
    )

    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id == "test"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"
      
      // check model_output
      assert state.containsKey("output_model") : "Output should contain key 'output_model'."
      assert state.output_model.isDirectory() : "'output_model' should be a directory."
      assert state.output_model.toString().endsWith("_model") : "Model output directory should end with '_model'. Found: ${state.output_model}"
      
      "Output: $output"
    }
    | nichecompass_leiden_test.run(
      fromState: [
        "input": "output"
      ]
    )
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 events"
    }
}
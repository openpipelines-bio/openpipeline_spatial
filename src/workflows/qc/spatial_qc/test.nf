nextflow.enable.dsl=2

include { spatial_qc } from params.rootDir + "/target/nextflow/workflows/qc/spatial_qc/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = 
    Channel.fromList([
      [
        id: "xenium_test",
        input: resources_test.resolve("xenium/xenium_tiny.h5mu"),
        var_name_mitochondrial_genes: "mitochondrial",
        var_name_ribosomal_genes: "ribosomal",
      ]
    ])
    | map { state -> [state.id, state] }
    | spatial_qc.run(
      toState: { id, output, state -> output + [og_input: state.input] }
    )

    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id.endsWith("_test")

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

    }
}
workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      | spatial_sample_processing.run(
        fromState: { id, state -> [
          "id": id,
          "input": state.input
        ]},
        toState: [
          "output": "output"
        ]
      )

      | setState(["output"])

  emit:
    output_ch
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | spaceranger_count.run(
      fromState: [
        "input": "input",
        "gex_reference": "gex_reference",
        "probe_set": "probe_set",
        "cytaimage": "cytaimage",
        "image": "image",
        "slide": "slide",
        "area": "area",
        "unkown_slide": "unkown_slide",
        "slidefile": "slidefile",
        "override_id": "override_id",
        "darkimage": "darkimage",
        "colorizedimage": "colorizedimage",
        "dapi_index": "dapi_index",
        "image_scale": "image_scale",
        "reorient_images": "reorient_images",
        "create_bam": "create_bam",
        "nosecondary": "nosecondary",
        "r1_length": "r1_length",
        "r2_length": "r2_length",
        "filter_probes": "filter_probes",
        "custom_bin_size": "custom_bin_size",
        "output": "output_raw",
      ],
      toState: [
        "input": "output",
        "output_raw": "output"
      ]
    )
    // split output dir into map
    | spaceranger_count_split.run(
      fromState: {id, state -> 
        def stateMapping = [
          "input": state.input,
        ]
        outputType = state.output_type == "raw" ? "raw_h5" : "filtered_h5"
        stateMapping += [outputType: "\$id.\$key.${outputType}.h5"]
        stateMapping += ["metrics_summary": "\$id.\$key.metrics_summary.csv"]
        return stateMapping
      },
      toState: {id, output, state -> 
        def outputFile = state.output_type == "raw" ? output.raw_h5 : output.filtered_h5
        def newState = state + [ "input": outputFile ] 
        return newState
      }
    )
    // convert to h5mu
    | from_10xh5_to_h5mu.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "output_compression": "gzip",
          "output": state.output_h5mu,
          "uns_metrics": state.uns_metrics,
          "input_metrics_summary": state.metrics_summary
        ]
      },
      toState: ["output_h5mu": "output"]
    )
    | setState(["output_raw", "output_h5mu"])

  emit:
  output_ch
}
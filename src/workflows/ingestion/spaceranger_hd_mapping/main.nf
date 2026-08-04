workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map { id, state ->
      assert state.input != null || state.spaceranger_output != null :
        "ERROR [$id]: provide either --input (FASTQs to align) or --spaceranger_output (a pre-computed Space Ranger output)."
      [id, state]
    }
    // Align reads with Space Ranger, unless a pre-computed output is provided.
    | spaceranger_count.run(
      runIf: { id, state -> state.spaceranger_output == null },
      fromState: { id, state -> [
        "input": state.input,
        "gex_reference": state.gex_reference,
        "probe_set": state.probe_set,
        "cytaimage": state.cytaimage,
        "image": state.image,
        "slide": state.slide,
        "area": state.area,
        "unknown_slide": state.unknown_slide,
        "create_bam": state.create_bam,
        "nosecondary": state.nosecondary,
        "custom_bin_size": state.custom_bin_size,
        "nucleus_segmentation": state.nucleus_segmentation,
        "disable_cell_annotation": state.disable_cell_annotation,
        "output": state.output_raw,
      ]},
      toState: [
        "output_raw": "output"
      ]
    )
    // Convert the segmented cells to SpatialData.
    | from_spaceranger_hd_to_spatialdata.run(
      fromState: { id, state -> [
        "input": state.spaceranger_output ?: state.output_raw,
        "output_type": state.output_type,
        "output_layer": state.output_layer,
        // Space Ranger writes a bare "feature_slice.h5"; the converter infers
        // dataset_id from a "<id>_feature_slice.h5" name, so default it to the
        // sample id to name the SpatialData elements.
        "dataset_id": state.dataset_id ?: id,
        "output": state.output_spatialdata,
      ]},
      toState: [
        "output_spatialdata": "output"
      ]
    )
    // When Space Ranger was skipped, report the provided output as output_raw.
    | map { id, state ->
      [id, state + ["output_raw": state.spaceranger_output ?: state.output_raw]]
    }
    | setState(["output_raw", "output_spatialdata"])

  emit:
  output_ch
}

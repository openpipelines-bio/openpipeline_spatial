workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | spaceranger_count.run(
      fromState: { id, state -> [
        "input": state.input,
        "gex_reference": state.gex_reference,
        "probe_set": state.probe_set,
        "cytaimage": state.cytaimage,
        "image": state.image,
        "slide": state.slide,
        "area": state.area,
        "unknown_slide": state.unknown_slide,
        "loupe_alignment": state.loupe_alignment,
        "create_bam": state.create_bam,
        "nosecondary": state.nosecondary,
        "custom_bin_size": state.custom_bin_size,
        "nucleus_segmentation": state.nucleus_segmentation,
        "cell_annotation_model": state.cell_annotation_model,
        "tenx_cloud_token_path": state.tenx_cloud_token_path,
        "disable_cell_annotation": state.disable_cell_annotation,
        "output": state.output_raw,
      ]},
      toState: [
        "output_raw": "output"
      ]
    )
    // convert the segmented cells to SpatialData
    | from_spaceranger_hd_to_spatialdata.run(
      fromState: { id, state -> [
        "input": state.output_raw,
        "mode": state.mode,
        "bin_size": state.bin_size,
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
    | setState(["output_raw", "output_spatialdata"])

  emit:
  output_ch
}

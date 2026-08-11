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
        "slidefile": state.slidefile,
        "override_id": state.override_id,
        "darkimage": state.darkimage,
        "colorizedimage": state.colorizedimage,
        "dapi_index": state.dapi_index,
        "image_scale": state.image_scale,
        "reorient_images": state.reorient_images,
        "loupe_alignment": state.loupe_alignment,
        "create_bam": state.create_bam,
        "nosecondary": state.nosecondary,
        "r1_length": state.r1_length,
        "r2_length": state.r2_length,
        "filter_probes": state.filter_probes,
        "custom_bin_size": state.custom_bin_size,
        "include_introns": state.include_introns,
        "nucleus_segmentation": state.nucleus_segmentation,
        "cell_annotation_model": state.cell_annotation_model,
        "tenx_cloud_token_path": state.tenx_cloud_token_path,
        "disable_cell_annotation": state.disable_cell_annotation,
        "custom_segmentation_file": state.custom_segmentation_file,
        "nucleus_expansion_distance_micron": state.nucleus_expansion_distance_micron,
        "max_nucleus_diameter_px": state.max_nucleus_diameter_px,
        "umi_registration": state.umi_registration,
        "umi_to_image_offset": state.umi_to_image_offset,
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

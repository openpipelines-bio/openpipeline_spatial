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
        "create_bam": state.create_bam,
        "nosecondary": state.nosecondary,
        "r1_length": state.r1_length,
        "r2_length": state.r2_length,
        "filter_probes": state.filter_probes,
        "custom_bin_size": state.custom_bin_size,
        "include_introns": state.include_introns,
        "nucleus_segmentation": state.nucleus_segmentation,
        "custom_segmentation_file": state.custom_segmentation_file,
        "nucleus_expansion_distance_micron": state.nucleus_expansion_distance_micron,
        "max_nucleus_diameter_px": state.max_nucleus_diameter_px,
        "umi_registration": state.umi_registration,
        "umi_to_image_offset": state.umi_to_image_offset,
        "cell_annotation_model": state.cell_annotation_model,
        "tenx_cloud_token_path": state.tenx_cloud_token_path,
        "disable_cell_annotation": state.disable_cell_annotation,
        "output": state.output_raw,
      ]},
      toState: [
        "input": "output",
        "output_raw": "output"
      ]
    )
    // convert to h5mu
    | from_spaceranger_to_h5mu.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "output_compression": state.output_compression,
          "output": state.output_h5mu,
          "uns_metrics": state.uns_metrics,
          "uns_probe_set": state.uns_probe_set,
          "obsm_coordinates": state.obsm_coordinates,
          "output_type": state.output_type,
          "output_compression": state.output_compression,
        ]
      },
      toState: ["output_h5mu": "output"]
    )
    | setState(["output_raw", "output_h5mu"])

  emit:
  output_ch
}
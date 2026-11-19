#!/bin/bash

set -eo pipefail

## VIASH START
par_xenium_bundle='resources_test/xenium_multichannel/xenium_multicellseg_tiny_cropped'
par_units='pixels'
par_expansion_distance=5
par_nuclei='resources_test/xenium_multichannel/xenium_multicellseg_tiny_cropped/cells.zarr.zip'
par_output='xeniumranger_import_segmentation_test'
## VIASH END

# Resolve paths before changing directory so relative inputs keep working
par_xenium_bundle=$(realpath "$par_xenium_bundle")
par_output=$(realpath -m "$par_output")
[[ -n "${par_nuclei:-}" ]] && par_nuclei=$(realpath "$par_nuclei")
[[ -n "${par_cells:-}" ]] && par_cells=$(realpath "$par_cells")
[[ -n "${par_coordinate_transform:-}" ]] && par_coordinate_transform=$(realpath "$par_coordinate_transform")
[[ -n "${par_viz_polygons:-}" ]] && par_viz_polygons=$(realpath "$par_viz_polygons")
[[ -n "${par_transcript_assignment:-}" ]] && par_transcript_assignment=$(realpath "$par_transcript_assignment")

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

cd "$tmpdir"

temp_id="xeniumranger_import_segmentation_run"

# Disable anonymized telemetry collection
export TENX_DISABLE_TELEMETRY=1

xeniumranger import-segmentation \
  --id="$temp_id" \
  --xenium-bundle="$par_xenium_bundle" \
  --disable-ui=true \
  ${par_units:+--units="$par_units"} \
  ${par_expansion_distance:+--expansion-distance="$par_expansion_distance"} \
  ${par_nuclei:+--nuclei="$par_nuclei"} \
  ${par_cells:+--cells="$par_cells"} \
  ${par_coordinate_transform:+--coordinate-transform="$par_coordinate_transform"} \
  ${par_viz_polygons:+--viz-polygons="$par_viz_polygons"} \
  ${par_transcript_assignment:+--transcript-assignment="$par_transcript_assignment"} \
  ${meta_cpus:+--localcores="$meta_cpus"} \
  ${meta_memory_gb:+--localmem=$(($meta_memory_gb-2))}

mkdir -p "$par_output"
mv -f "$temp_id"/outs/* "$par_output"/

#!/bin/bash

set -eo pipefail

## VIASH START
par_xenium_bundle='resources_test/xenium/xenium_tiny'
par_id='xenium_tiny_import_segmentation'
par_units='pixels'
par_expansion_distance=5
par_nuclei='resources_test/xenium/xenium_tiny/cells.zarr.zip'
par_output='xeniumranger_import_segmentation_test'
## VIASH END

par_xenium_bundle=`realpath $par_xenium_bundle`
par_output=`realpath $par_output`

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

cd "$tmpdir"

xeniumranger import-segmentation \
  --id="$par_id" \
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
mv -f "$par_id"/outs/* "$par_output"/
rm -rf "$par_id"/outs

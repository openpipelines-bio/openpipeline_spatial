#!/bin/bash

set -eo pipefail

## VIASH START
par_xenium_bundle='resources_test/xenium/xenium_tiny'
par_id='xenium_tiny_resegment'
par_boundary_stain='ATP1A1/CD45/E-Cadherin'
par_interior_stain='18S'
par_expansion_distance=5
par_dapi_filter=15
par_output='xeniumranger_resegment_test'
## VIASH END

unset_if_false=(
   par_resegment_nuclei
   par_segment_large_cells
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

par_xenium_bundle=`realpath $par_xenium_bundle`
par_output=`realpath $par_output`

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

cd "$tmpdir"

xeniumranger resegment \
  --id="$par_id" \
  --xenium-bundle="$par_xenium_bundle" \
  --disable-ui=true \
  ${par_boundary_stain:+--boundary-stain="$par_boundary_stain"} \
  ${par_interior_stain:+--interior-stain="$par_interior_stain"} \
  ${par_segment_large_cells:+--segment-large-cells} \
  ${par_expansion_distance:+--expansion-distance="$par_expansion_distance"} \
  ${par_dapi_filter:+--dapi-filter="$par_dapi_filter"} \
  ${par_resegment_nuclei:+--resegment-nuclei="$par_resegment_nuclei"} \
  ${meta_cpus:+--localcores="$meta_cpus"} \
  ${meta_memory_gb:+--localmem=$(($meta_memory_gb-2))}

mkdir -p "$par_output"
mv -f "$par_id"/outs/* "$par_output"/
rm -rf "$par_id"/outs

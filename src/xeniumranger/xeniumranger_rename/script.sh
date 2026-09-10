#!/bin/bash

set -eo pipefail

## VIASH START
par_xenium_bundle='resources_test/xenium/xenium_tiny'
par_id='xenium_tiny_rename'
par_region_name=''
par_cassette_name=''
par_output='xeniumranger_rename_test'
## VIASH END

# Make sure paths are absolute, since we cd into a tempdir before running xeniumranger
par_xenium_bundle=`realpath $par_xenium_bundle`
par_output=`realpath $par_output`

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# cd into tempdir
cd "$tmpdir"

xeniumranger rename \
  --id="$par_id" \
  --xenium-bundle="$par_xenium_bundle" \
  ${par_region_name:+--region-name="$par_region_name"} \
  ${par_cassette_name:+--cassette-name="$par_cassette_name"} \
  --disable-ui=true \
  ${meta_cpus:+--localcores="$meta_cpus"} \
  ${meta_memory_gb:+--localmem=$(($meta_memory_gb-2))}

mkdir -p "$par_output"
mv -f "$par_id"/outs/* "$par_output"/
rm -rf "$par_id"/outs


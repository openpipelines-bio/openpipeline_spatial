#!/bin/bash

set -eo pipefail

## VIASH START
par_xenium_bundle="$par_xenium_bundle"
par_region_name=''
par_cassette_name=''
par_output="xeniumranger_rename_test"
## VIASH END

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# Resolve paths before changing directory so relative inputs keep working
par_xenium_bundle=$(realpath "$par_xenium_bundle")
par_output=$(realpath -m "$par_output")

cd "$tmpdir"

temp_id="xeniumranger_rename_run"

# Disable anonymized telemetry collection
export TENX_DISABLE_TELEMETRY=1

xeniumranger rename \
  --id="$temp_id" \
  --xenium-bundle="$par_xenium_bundle" \
  ${par_region_name:+--region-name="$par_region_name"} \
  ${par_cassette_name:+--cassette-name="$par_cassette_name"} \
  --disable-ui=true \
  ${meta_cpus:+--localcores="$meta_cpus"} \
  ${meta_memory_gb:+--localmem=$(($meta_memory_gb-2))}

mkdir -p "$par_output"
mv -f "$temp_id"/outs/* "$par_output"/

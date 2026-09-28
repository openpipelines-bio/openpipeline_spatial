#!/bin/bash

set -eo pipefail

## VIASH START
par_xenium_bundle='resources_test/xenium/xenium_tiny'
par_panel='resources_test/xenium/xenium_tiny/gene_panel.json'
par_output='xeniumranger_relabel_test'
## VIASH END

par_xenium_bundle=`realpath $par_xenium_bundle`
par_panel=`realpath $par_panel`
par_output=`realpath $par_output`

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

cd "$tmpdir"

temp_id="xeniumranger_relabel_run"

# Disable anonymized telemetry collection
export TENX_DISABLE_TELEMETRY=1

xeniumranger relabel \
  --id="$temp_id" \
  --xenium-bundle="$par_xenium_bundle" \
  --panel="$par_panel" \
  --disable-ui=true \
  ${meta_cpus:+--localcores="$meta_cpus"} \
  ${meta_memory_gb:+--localmem=$(($meta_memory_gb-2))}

mkdir -p "$par_output"
mv -f "$temp_id"/outs/* "$par_output"/

#!/bin/bash

set -eo pipefail

## VIASH START
par_xenium_bundle='resources_test/xenium/xenium_tiny'
par_id='xenium_tiny_relabel'
par_panel='resources_test/xenium/xenium_tiny/gene_panel.json'
par_output='xeniumranger_relabel_test'
## VIASH END 

# Make sure paths are absolute, since we cd into a tempdir before running xeniumranger
par_xenium_bundle=`realpath $par_xenium_bundle`
par_panel=`realpath $par_panel`
par_output=`realpath $par_output`

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# cd into tempdir
cd "$tmpdir"

xeniumranger relabel \
  --id="$par_id" \
  --xenium-bundle="$par_xenium_bundle" \
  --panel="$par_panel" \
  --disable-ui \
  ${meta_cpus:+--localcores="$meta_cpus"} \
  ${meta_memory_gb:+--localmem=$(($meta_memory_gb-2))}

mkdir -p "$par_output"
mv -f "$par_id"/outs/* "$par_output"/
rm -rf "$par_id"/outs


#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DIR="resources_test/niche"
ID="nichecompass"

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# Define file names
orthologue_file="human_mouse_gene_orthologs.csv"
enzymes_file="mouse_metabolite_enzymes.tsv"
sensors_file="mouse_metabolite_sensors.tsv"

if [ ! -d "$DIR" ]; then
    mkdir -p "$DIR"

    orthologue_url="https://raw.githubusercontent.com/Lotfollahi-lab/nichecompass/refs/tags/0.3.2/data/gene_annotations/human_mouse_gene_orthologs.csv"
    wget "$orthologue_url" -O "$DIR/$orthologue_file"

    enzymes_url="https://raw.githubusercontent.com/Lotfollahi-lab/nichecompass/refs/tags/0.3.2/data/gene_programs/metabolite_enzyme_sensor_gps/mouse_metabolite_enzymes.tsv"
    wget "$enzymes_url" -O "$DIR/$enzymes_file"

    sensors_url="https://raw.githubusercontent.com/Lotfollahi-lab/nichecompass/refs/tags/0.3.2/data/gene_programs/metabolite_enzyme_sensor_gps/mouse_metabolite_sensors.tsv"
    wget "$sensors_url" -O "$DIR/$sensors_file" 
fi

# Generate omnipath and collectri network files for testing
# These files are used to avoid API rate limits when running tests in parallel
omnipath_file="omnipath_lr_network.csv"
collectri_file="collectri_tf_network.csv"

gp_mask="prior_knowledge_gp_mask.json"
viash run src/nichecompass/gene_program_mask/config.vsh.yaml -- \
    --input_gene_orthologs_mapping_file "$DIR/$orthologue_file" \
    --input_metabolite_enzymes "$DIR/$enzymes_file" \
    --input_metabolite_sensors "$DIR/$sensors_file" \
    --output "${DIR}/${gp_mask}" \
    --output_omnipath_lr_network "${DIR}/${omnipath_file}" \
    --output_collectri_tf_network "${DIR}/${collectri_file}"

# Sync to S3
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/niche \
    --delete \
    --dryrun

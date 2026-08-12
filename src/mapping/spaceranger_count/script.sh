#!/bin/bash

set -eo pipefail

## VIASH START
par_input='resources_test/visium/Visium_FFPE_Human_Ovarian_Cancer_fastqs'
par_image='resources_test/visium/Visium_FFPE_Human_Ovarian_Cancer_image.jpg'
par_output='spaceranger_test'
par_gex_reference='resources_test/GRCh38'
par_probe_set='resources_test/visium/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv'
par_slide='V10L13-020'
par_area='D1'
par_create_bam='false'
## VIASH END

unset_if_false=(
    par_override_id
    par_nosecondary
    par_disable_cell_annotation
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Make sure paths are absolute, since we cd into a tempdir before running spaceranger
par_gex_reference=`realpath $par_gex_reference`
par_output=`realpath $par_output`
[[ -n "${par_probe_set:-}" ]] && par_probe_set=$(realpath "$par_probe_set")
[[ -n "${par_image:-}" ]] && par_image=$(realpath "$par_image")
[[ -n "${par_cytaimage:-}" ]] && par_cytaimage=$(realpath "$par_cytaimage")
[[ -n "${par_slidefile:-}" ]] && par_slidefile=$(realpath "$par_slidefile")
[[ -n "${par_loupe_alignment:-}" ]] && par_loupe_alignment=$(realpath "$par_loupe_alignment")
[[ -n "${par_darkimage:-}" ]] && par_darkimage=$(realpath "$par_darkimage")
[[ -n "${par_colorizedimage:-}" ]] && par_colorizedimage=$(realpath "$par_colorizedimage")
[[ -n "${par_custom_segmentation_file:-}" ]] && par_custom_segmentation_file=$(realpath "$par_custom_segmentation_file")
[[ -n "${par_tenx_cloud_token_path:-}" ]] && par_tenx_cloud_token_path=$(realpath "$par_tenx_cloud_token_path")

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# process inputs
# for every fastq file found, make a symlink into the tempdir
fastq_dir="$tmpdir/fastqs"
mkdir -p "$fastq_dir"
IFS=";"
for var in $par_input; do
  unset IFS
  abs_path=`realpath $var`
  if [ -d "$abs_path" ]; then
    find "$abs_path" -name *.fastq.gz -exec ln -s {} "$fastq_dir" \;
  else
    ln -s "$abs_path" "$fastq_dir"
  fi
done

# process reference: untar if it is a gzipped archive rather than a directory
# (`gzip -t` avoids depending on the `file` command, which the image lacks)
if [ -f "$par_gex_reference" ] && gzip -t "$par_gex_reference" 2>/dev/null; then
  echo "Untarring genome"
  reference_dir="$tmpdir/fastqs"
  mkdir -p "$reference_dir"
  tar -xvf "$par_gex_reference" -C "$reference_dir" --strip-components=1
  par_gex_reference="$reference_dir"
fi

# cd into tempdir
cd "$tmpdir"

temp_id="spaceranger_run"

# Disable anonymized telemetry collection (enabled by default since Space Ranger 4.0)
export TENX_DISABLE_TELEMETRY=1

spaceranger count \
  --id="$temp_id" \
  --disable-ui \
  --fastqs="$fastq_dir" \
  --transcriptome="$par_gex_reference" \
  ${par_probe_set:+--probe-set="$par_probe_set"} \
  ${par_cytaimage:+--cytaimage="$par_cytaimage"} \
  ${par_image:+--image="$par_image"} \
  ${par_slide:+--slide="$par_slide"} \
  ${par_area:+--area="$par_area"} \
  ${par_unknown_slide:+--unknown-slide="$par_unknown_slide"} \
  ${par_slidefile:+--slidefile="$par_slidefile"} \
  ${par_loupe_alignment:+--loupe-alignment="$par_loupe_alignment"} \
  ${par_override_id:+--override-id} \
  ${par_darkimage:+--darkimage="$par_darkimage"} \
  ${par_colorizedimage:+--colorizedimage="$par_colorizedimage"} \
  ${par_dapi_index:+--dapi-index="$par_dapi_index"} \
  ${par_image_scale:+--image-scale="$par_image_scale"} \
  ${par_reorient_images:+--reorient-images="$par_reorient_images"} \
  ${par_create_bam:+--create-bam="$par_create_bam"} \
  ${par_nosecondary:+--nosecondary} \
  ${par_r1_length:+--r1-length="$par_r1_length"} \
  ${par_r2_length:+--r2-length="$par_r2_length"} \
  ${par_filter_probes:+--filter-probes="$par_filter_probes"} \
  ${par_custom_bin_size:+--custom-bin-size="$par_custom_bin_size"} \
  ${par_include_introns:+--include-introns="$par_include_introns"} \
  ${par_nucleus_segmentation:+--nucleus-segmentation="$par_nucleus_segmentation"} \
  ${par_custom_segmentation_file:+--custom-segmentation-file="$par_custom_segmentation_file"} \
  ${par_nucleus_expansion_distance_micron:+--nucleus-expansion-distance-micron="$par_nucleus_expansion_distance_micron"} \
  ${par_max_nucleus_diameter_px:+--max-nucleus-diameter-px="$par_max_nucleus_diameter_px"} \
  ${par_umi_registration:+--umi-registration="$par_umi_registration"} \
  ${par_umi_to_image_offset:+--umi-to-image-offset="$par_umi_to_image_offset"} \
  ${par_cell_annotation_model:+--cell-annotation-model="$par_cell_annotation_model"} \
  ${par_tenx_cloud_token_path:+--tenx-cloud-token-path="$par_tenx_cloud_token_path"} \
  ${par_disable_cell_annotation:+--disable-cell-annotation} \
  ${par_project:+--project="$par_project"} \
  ${par_sample:+--sample="$par_sample"} \
  ${par_lanes:+--lanes="$par_lanes"} \
  ${meta_cpus:+--localcores="$meta_cpus"} \
  ${meta_memory_gb:+--localmem=$(($meta_memory_gb-2))}

mkdir -p "$par_output"
mv -f "$temp_id"/outs/* "$par_output"/
rm -rf "$temp_id"/outs
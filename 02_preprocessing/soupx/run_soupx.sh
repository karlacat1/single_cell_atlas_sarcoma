#!/usr/bin/env bash
#
# Run SoupX ambient RNA removal for a list of samples.
#
# Usage:
#   ./run_soupx.sh <sample_list> <data_dir> <output_dir> <soupx_script>
#
#   sample_list   Text file with one sample name per line.
#   data_dir      Directory containing the aligned samples, with one
#                 <sample_name>_aligned/outs directory per sample.
#   output_dir    Directory for the corrected count matrices. One
#                 subdirectory per sample is created.
#   soupx_script  Either soupx_auto.R (automatic contamination estimation,
#                 used for all entities except Ewing sarcoma) or
#                 soupx_fixed_cf.R (contamination fraction fixed at 0.4,
#                 used for Ewing sarcoma).
#
# Example:
#   ./run_soupx.sh sample_list.txt data/aligned /data/soupx soupx_auto.R

set -euo pipefail

SAMPLE_LIST="${1:?sample_list required}"
DATA_DIR="${2:?data_dir required}"
OUTPUT_DIR="${3:?output_dir required}"

SOUPX_SCRIPT="soupx.R"

# Set mode here:
SOUPX_MODE=(--auto)
# SOUPX_MODE=(--fixed 0.4)

while read -r sample_name; do
    [ -z "$sample_name" ] && continue
    echo "$sample_name"

    Rscript "$SOUPX_SCRIPT" \
        "${DATA_DIR}/${sample_name}_aligned/outs" \
        "$sample_name" \
        "${OUTPUT_DIR}/${sample_name}" \
        "${SOUPX_MODE[@]}"

done < "$SAMPLE_LIST"

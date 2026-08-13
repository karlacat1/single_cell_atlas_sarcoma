#!/usr/bin/env bash

#
# Run SoupX ambient RNA removal for a list of samples.
#
# Usage:
# ./run_soupx.sh <sample_list> <data_dir> <output_dir> <mode> [contamination_fraction]
#
# sample_list          Text file with one sample name per line.
# data_dir             Directory containing the aligned samples, with one
#                      <sample_name>_aligned/outs directory per sample.
# output_dir           Directory for the corrected count matrices. One
#                      subdirectory per sample is created.
# mode                 Either "auto" for automatic contamination estimation
#                      or "fixed" for a manually specified contamination fraction.
# contamination_fraction
#                      Contamination fraction used with "fixed" mode.
#
# Examples:
# ./run_soupx.sh sample_list.txt data/aligned /data/soupx auto
# ./run_soupx.sh sample_list.txt data/aligned /data/soupx fixed 0.4

set -euo pipefail

SAMPLE_LIST="${1:?sample_list required}"
DATA_DIR="${2:?data_dir required}"
OUTPUT_DIR="${3:?mode required}"
MODE="${4:?mode required (auto or fixed)}"

SOUPX_SCRIPT="soupx.R"

if [[ "$MODE" == "fixed" ]]; then
    CF="${5:?contamination fraction required for fixed mode}"
elif [[ "$MODE" != "auto" ]]; then
    echo "Error: mode must be 'auto' or 'fixed'."
    exit 1
fi

while read -r sample_name; do
    [ -z "$sample_name" ] && continue
    echo "$sample_name"

    if [[ "$MODE" == "auto" ]]; then
        Rscript "$SOUPX_SCRIPT" \
            "${DATA_DIR}/${sample_name}_aligned/outs" \
            "$sample_name" \
            "${OUTPUT_DIR}/${sample_name}" \
            --auto
    else
        Rscript "$SOUPX_SCRIPT" \
            "${DATA_DIR}/${sample_name}_aligned/outs" \
            "$sample_name" \
            "${OUTPUT_DIR}/${sample_name}" \
            --fixed "$CF"
    fi

done < "$SAMPLE_LIST"

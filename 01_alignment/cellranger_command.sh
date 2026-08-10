#!/usr/bin/env bash
#
# Align single-nucleus RNA-seq reads and generate count matrices with Cell Ranger.
#
# Runs `cellranger count` for one sample. Intronic reads are included
# (--include-introns=true) to capture unspliced transcripts, which make up a large
# fraction of single-nucleus libraries. Output is a standard Cell Ranger run
# directory containing the raw and filtered feature-barcode matrices used for all
# downstream analysis.
#
# Usage:
#   ./run_cellranger_count.sh <sample_id> <fastq_dir> <fastq_prefixes> [output_dir]
#
# Arguments:
#   sample_id       Identifier for this sample; used as the Cell Ranger run ID.
#   fastq_dir       Directory containing the demultiplexed FASTQ files.
#   fastq_prefixes  Comma-separated FASTQ prefixes for this sample, as passed to
#                   --sample. Libraries sequenced across several lanes or runs have
#                   more than one prefix (e.g. "21106a084_01,21106a084_02").
#   output_dir      Optional. Where the run directory is written. Default: current
#                   working directory.
#
# Example:
#   ./run_cellranger_count.sh STS_204_002 /data/fastq/STS_204_002 lib_01,lib_02 ./aligned
#
# Reference:
#   Set CELLRANGER_REF to the 10x Genomics human reference. This study used
#   refdata-gex-GRCh38-2020-A, available from
#   https://www.10xgenomics.com/support/software/cell-ranger/downloads
#
# Requirements:
#   Cell Ranger v6.1.1 on PATH; approximately 80 GB RAM and 10 cores per sample.

set -euo pipefail

SAMPLE_ID="${1:?sample_id required}"
FASTQ_DIR="${2:?fastq_dir required}"
FASTQ_PREFIXES="${3:?fastq_prefixes required}"
OUTPUT_DIR="${4:-$(pwd)}"

REF="${CELLRANGER_REF:?Set CELLRANGER_REF to the path of refdata-gex-GRCh38-2020-A}"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

cellranger count \
    --id="${SAMPLE_ID}_aligned" \
    --transcriptome="$REF" \
    --fastqs="$FASTQ_DIR" \
    --sample="$FASTQ_PREFIXES" \
    --include-introns=true 
#!/usr/bin/env python3
"""
Post-CNV biological annotation of sarcoma single-cell RNA-seq data.

This script performs the annotation step after CNV/SCEVAN-based
malignant/non-malignant classification.

The preceding CNV workflow must provide:

    adata.obs["cnv_annotation"]

with values:

    "malignant"
    "non-malignant"

Tumor-specific annotation
-------------------------
Each tumor type is associated with EXACTLY ONE malignant marker list.

Supported tumor types:

    - osteosarcoma
    - ewing
    - clear_cell_sarcoma
    - embryonal_rhabdomyosarcoma
    - synovial_sarcoma
    - epithelioid_sarcoma
    - undifferentiated_sarcoma
    - alveolar_soft_part_sarcoma
    - myoepithelial_carcinoma

For malignant cells, the marker panel corresponding to the selected
tumor type is scored and stored as a tumor-specific score.

Common non-malignant annotation
-------------------------------
The same marker panels are used for non-malignant cells across all
tumor types:

    - immune
    - myeloid
    - macrophages
    - monocytes
    - endothelial
    - pericytes
    - stromal

For non-malignant cells, the highest-scoring common marker set is used
for the final cell annotation.

Outputs
-------
The following columns are added to adata.obs:

    cnv_annotation
        Original CNV-based malignant/non-malignant classification.

    final_annotation
        Final biological annotation.

    annotation_score
        Score of the marker set used for the final annotation.

    annotation_marker_set
        Marker set used for the final annotation.

All individual Scanpy gene-set scores are retained in adata.obs.

Notes
-----
- CNV classification is NOT recalculated by this script.
- Marker scores are calculated with scanpy.tl.score_genes.
- Only marker genes present in the dataset are used.
- Missing marker genes are reported.
- The minimum score threshold is configurable.
- The supplied marker name "COLA1A" is retained exactly as provided in
  the original analysis. Verify this gene symbol before final repository
  deposition if necessary.

Example
-------
python annotate_sarcoma.py \
    --input-h5ad results/sample_cnv_classified.h5ad \
    --output-h5ad results/sample_final.h5ad \
    --output-annotations results/sample_final_annotations.txt \
    --tumor-type ewing
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


# =============================================================================
# Common non-malignant marker sets
# =============================================================================
#
# These marker sets are applied to all tumor types.
# =============================================================================

COMMON_NON_MALIGNANT_MARKERS = {
    "immune": [
        "PTPRC",
        "CD4",
        "IL7R",
    ],
    "myeloid": [
        "CD14",
        "CD163",
        "CD68",
        "CD80",
    ],
    "macrophages": [
        "CD68",
        "CD163",
        "CD14",
        "ITGAM",
    ],
    "monocytes": [
        "ITGAX",
        "CD80",
        "CD86",
    ],
    "endothelial": [
        "CD34",
        "VWF",
        "PECAM1",
        "KDR",
    ],
    "pericytes": [
        "DCN",
        "FBN1",
        "FBLN1",
        "THBS2",
    ],
    "stromal": [
        "PDGFRB",
        "RGS5",
        "CD248",
        "MCAM",
        "COLA1A",
    ],
}


# =============================================================================
# Tumor-specific malignant marker sets
# =============================================================================
#
# IMPORTANT:
# Each tumor type has ONE marker list.
# =============================================================================

MALIGNANT_MARKERS = {
    "osteosarcoma": [
        "RUNX2",
    ],
    "ewing": [
        "CD99",
        "CAV1",
        "CCND1",
        "HES1",
        "KDSR",
        "PAPPA",
    ],
    "clear_cell_sarcoma": [
        "MLANA",
        "S100B",
        "MITF",
        "CD99",
        "EWSR1",
        "ATF1",
    ],
    "embryonal_rhabdomyosarcoma": [
        "DES",
        "ACTA1",
        "PAX7",
    ],
    "synovial_sarcoma": [
        "SS18",
        "SSX1",
    ],
    "epithelioid_sarcoma": [
        "SMARCB1",
        "TFE3",
        "TFEB",
        "TP53",
    ],
    "undifferentiated_sarcoma": [
        "SMARCB1",
        "TFE3",
        "TFEB",
        "TP53",
    ],
    "alveolar_soft_part_sarcoma": [
        "ASPSCR1",
        "TFE3",
    ],
    "myoepithelial_carcinoma": [
        "EWSR1",
        "POU5F1",
    ],
}


# =============================================================================
# Input validation
# =============================================================================

def validate_input(adata: sc.AnnData) -> None:
    """
    Validate the input AnnData object.

    The input must contain a CNV-based binary malignant/non-malignant
    classification.
    """
    if "cnv_annotation" not in adata.obs.columns:
        raise ValueError(
            "The input AnnData object does not contain "
            "`adata.obs['cnv_annotation']`.\n"
            "Run the CNV classification script first."
        )

    valid_annotations = {
        "malignant",
        "non-malignant",
    }

    observed_annotations = set(
        adata.obs["cnv_annotation"]
        .dropna()
        .astype(str)
        .unique()
    )

    unexpected = observed_annotations - valid_annotations

    if unexpected:
        raise ValueError(
            "Unexpected values found in `adata.obs['cnv_annotation']`: "
            + ", ".join(sorted(unexpected))
            + ". Expected only 'malignant' and 'non-malignant'."
        )


# =============================================================================
# Marker preparation
# =============================================================================

def prepare_marker_sets(
    adata: sc.AnnData,
    marker_sets: dict[str, list[str]],
) -> dict[str, list[str]]:
    """
    Restrict marker sets to genes present in the dataset.

    Missing genes are reported. Marker sets with no genes present in the
    dataset are excluded from scoring.
    """
    available_genes = set(adata.var_names)
    usable_marker_sets: dict[str, list[str]] = {}

    for name, genes in marker_sets.items():
        present = [
            gene
            for gene in genes
            if gene in available_genes
        ]

        missing = [
            gene
            for gene in genes
            if gene not in available_genes
        ]

        if missing:
            print(
                f"Warning: marker set '{name}' is missing "
                f"{len(missing)} gene(s): "
                f"{', '.join(missing)}"
            )

        if not present:
            print(
                f"Warning: marker set '{name}' has no genes "
                "present in the dataset and will not be scored."
            )
            continue

        usable_marker_sets[name] = present

    return usable_marker_sets


# =============================================================================
# Marker scoring
# =============================================================================

def score_marker_set(
    adata: sc.AnnData,
    genes: list[str],
    score_name: str,
) -> str:
    """
    Calculate a Scanpy gene-set score.

    Returns
    -------
    str
        Name of the resulting score column.
    """
    print(
        f"Scoring '{score_name}' "
        f"using {len(genes)} gene(s): "
        f"{', '.join(genes)}"
    )

    sc.tl.score_genes(
        adata,
        gene_list=genes,
        score_name=score_name,
    )

    return score_name


def score_common_non_malignant_markers(
    adata: sc.AnnData,
    marker_sets: dict[str, list[str]],
) -> list[str]:
    """
    Score all common non-malignant marker sets.
    """
    score_columns = []

    for annotation, genes in marker_sets.items():
        score_name = f"{annotation}_score"

        score_marker_set(
            adata=adata,
            genes=genes,
            score_name=score_name,
        )

        score_columns.append(score_name)

    return score_columns


# =============================================================================
# Malignant annotation
# =============================================================================

def annotate_malignant_cells(
    adata: sc.AnnData,
    tumor_type: str,
    malignant_markers: list[str],
) -> None:
    """
    Annotate CNV-defined malignant cells using the single marker panel
    associated with the selected tumor type.

    The tumor-specific score is stored in:

        adata.obs["<tumor_type>_malignant_score"]

    Malignant cells are annotated as the selected tumor type.
    """
    malignant_mask = (
        adata.obs["cnv_annotation"].astype(str)
        == "malignant"
    )

    score_column = f"{tumor_type}_malignant_score"

    if not malignant_mask.any():
        print("No malignant cells detected.")
        return

    score_marker_set(
        adata=adata,
        genes=malignant_markers,
        score_name=score_column,
    )

    # The CNV classification establishes malignant status.
    # The tumor-specific marker score provides the biological annotation.
    adata.obs.loc[
        malignant_mask,
        "final_annotation",
    ] = tumor_type

    adata.obs.loc[
        malignant_mask,
        "annotation_score",
    ] = adata.obs.loc[
        malignant_mask,
        score_column,
    ]

    adata.obs.loc[
        malignant_mask,
        "annotation_marker_set",
    ] = tumor_type


# =============================================================================
# Non-malignant annotation
# =============================================================================

def annotate_non_malignant_cells(
    adata: sc.AnnData,
    marker_sets: dict[str, list[str]],
) -> None:
    """
    Annotate CNV-defined non-malignant cells using the common marker panels.

    The marker set with the highest score is assigned to each cell.

    Cells for which no marker set can be scored retain the generic
    "non-malignant" annotation.
    """
    non_malignant_mask = (
        adata.obs["cnv_annotation"].astype(str)
        == "non-malignant"
    )

    if not non_malignant_mask.any():
        print("No non-malignant cells detected.")
        return

    score_columns = score_common_non_malignant_markers(
        adata=adata,
        marker_sets=marker_sets,
    )

    if not score_columns:
        print(
            "No common non-malignant marker sets could be scored."
        )
        return

    scores = adata.obs.loc[
        non_malignant_mask,
        score_columns,
    ]

    best_score_column = scores.idxmax(axis=1)
    best_score = scores.max(axis=1)

    # Convert "<cell_type>_score" into the final cell-type label.
    final_labels = (
        best_score_column
        .str.replace(
            "_score",
            "",
            regex=False,
        )
    )

    adata.obs.loc[
        non_malignant_mask,
        "final_annotation",
    ] = final_labels.values

    adata.obs.loc[
        non_malignant_mask,
        "annotation_score",
    ] = best_score.values

    adata.obs.loc[
        non_malignant_mask,
        "annotation_marker_set",
    ] = final_labels.values


# =============================================================================
# Summary
# =============================================================================

def print_annotation_summary(
    adata: sc.AnnData,
) -> None:
    """
    Print the final annotation distribution.
    """
    counts = (
        adata.obs["final_annotation"]
        .value_counts()
    )

    total = counts.sum()

    print("\nFinal annotation summary")
    print("=" * 72)

    for annotation, count in counts.items():
        percentage = (
            count / total * 100
            if total > 0
            else 0.0
        )

        print(
            f"{annotation:40s}"
            f"{count:>8,} cells "
            f"({percentage:5.1f}%)"
        )

    print("=" * 72)


# =============================================================================
# Save results
# =============================================================================

def save_results(
    adata: sc.AnnData,
    output_h5ad: Path,
    output_annotations: Path,
) -> None:
    """
    Save the annotated AnnData object and final cell-level annotations.
    """
    output_h5ad.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_annotations.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Save the complete AnnData object, including all marker scores.
    adata.write(output_h5ad)

    # Save the main annotation table.
    annotation_table = adata.obs[
        [
            "cnv_annotation",
            "final_annotation",
            "annotation_score",
            "annotation_marker_set",
        ]
    ].copy()

    annotation_table.to_csv(
        output_annotations,
        sep="\t",
        index=True,
    )

    print("\nResults saved:")
    print(f"  AnnData:     {output_h5ad}")
    print(f"  Annotations: {output_annotations}")


# =============================================================================
# Main workflow
# =============================================================================

def run_annotation(
    input_h5ad: Path,
    output_h5ad: Path,
    output_annotations: Path,
    tumor_type: str,
) -> None:
    """
    Execute the post-CNV annotation workflow.
    """
    # -------------------------------------------------------------------------
    # 1. Read CNV-classified data
    # -------------------------------------------------------------------------
    print(f"Reading AnnData: {input_h5ad}")

    adata = sc.read_h5ad(input_h5ad)

    print(
        f"Dataset contains "
        f"{adata.n_obs:,} cells × "
        f"{adata.n_vars:,} genes."
    )

    # -------------------------------------------------------------------------
    # 2. Validate CNV classification
    # -------------------------------------------------------------------------
    validate_input(adata)

    # -------------------------------------------------------------------------
    # 3. Select the single malignant marker panel for this tumor type
    # -------------------------------------------------------------------------
    malignant_markers = MALIGNANT_MARKERS[tumor_type]

    print(f"\nTumor type: {tumor_type}")
    print(
        "Malignant marker panel: "
        + ", ".join(malignant_markers)
    )

    # -------------------------------------------------------------------------
    # 4. Prepare tumor-specific malignant markers
    # -------------------------------------------------------------------------
    usable_malignant_markers = prepare_marker_sets(
        adata,
        {tumor_type: malignant_markers},
    )

    if tumor_type not in usable_malignant_markers:
        raise ValueError(
            f"None of the malignant marker genes for '{tumor_type}' "
            "are present in the dataset."
        )

    malignant_markers = usable_malignant_markers[tumor_type]

    # -------------------------------------------------------------------------
    # 5. Prepare the common non-malignant marker panels
    # -------------------------------------------------------------------------
    usable_non_malignant_markers = prepare_marker_sets(
        adata,
        COMMON_NON_MALIGNANT_MARKERS,
    )

    # -------------------------------------------------------------------------
    # 6. Initialize the final annotation from the CNV classification
    # -------------------------------------------------------------------------
    adata.obs["final_annotation"] = (
        adata.obs["cnv_annotation"]
        .astype(str)
    )

    adata.obs["annotation_score"] = pd.NA
    adata.obs["annotation_marker_set"] = pd.NA

    # -------------------------------------------------------------------------
    # 7. Annotate malignant cells using the ONE tumor-specific marker list
    # -------------------------------------------------------------------------
    print("\nAnnotating malignant cells...")

    annotate_malignant_cells(
        adata=adata,
        tumor_type=tumor_type,
        malignant_markers=malignant_markers,
    )

    # -------------------------------------------------------------------------
    # 8. Annotate non-malignant cells using the common marker panels
    # -------------------------------------------------------------------------
    print("\nAnnotating non-malignant cells...")

    annotate_non_malignant_cells(
        adata=adata,
        marker_sets=usable_non_malignant_markers,
    )

    # -------------------------------------------------------------------------
    # 9. Print final annotation summary
    # -------------------------------------------------------------------------
    print_annotation_summary(adata)

    # -------------------------------------------------------------------------
    # 10. Save results
    # -------------------------------------------------------------------------
    save_results(
        adata=adata,
        output_h5ad=output_h5ad,
        output_annotations=output_annotations,
    )


# =============================================================================
# Command-line interface
# =============================================================================

def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Perform post-CNV biological annotation of "
            "sarcoma single-cell RNA-seq data."
        )
    )

    parser.add_argument(
        "--input-h5ad",
        type=Path,
        required=True,
        help=(
            "Input AnnData file produced by the CNV classification "
            "workflow. Must contain adata.obs['cnv_annotation']."
        ),
    )

    parser.add_argument(
        "--output-h5ad",
        type=Path,
        required=True,
        help="Output path for the annotated AnnData file.",
    )

    parser.add_argument(
        "--output-annotations",
        type=Path,
        required=True,
        help="Output path for the final cell-level annotation table.",
    )

    parser.add_argument(
        "--tumor-type",
        required=True,
        choices=[
            "osteosarcoma",
            "ewing",
            "clear_cell_sarcoma",
            "embryonal_rhabdomyosarcoma",
            "synovial_sarcoma",
            "epithelioid_sarcoma",
            "undifferentiated_sarcoma",
            "alveolar_soft_part_sarcoma",
            "myoepithelial_carcinoma",
        ],
        help=(
            "Tumor type. This selects exactly one "
            "tumor-specific malignant marker panel."
        ),
    )

    return parser.parse_args()


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    args = parse_arguments()

    run_annotation(
        input_h5ad=args.input_h5ad,
        output_h5ad=args.output_h5ad,
        output_annotations=args.output_annotations,
        tumor_type=args.tumor_type,
    )

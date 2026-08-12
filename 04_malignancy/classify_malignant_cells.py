#!/usr/bin/env python3
"""
Classify single-cell RNA-seq cells as malignant or non-malignant
using SCEVAN CNV-based predictions.

Overview
--------
This script integrates cell-level SCEVAN classifications with a processed
single-cell AnnData object and generates a cluster-aware malignant
annotation.

For each Leiden cluster:

1. The most frequent SCEVAN class is identified.
2. A cluster is classified as malignant ("malignant") when:
   - more than 40% of its cells are classified as "tumor", OR
   - the combined fraction of "tumor" and "filtered" cells exceeds 40%.
3. All cells in malignant clusters are labelled "malignant".
4. Cells in remaining clusters are labelled "non-malignant".

The resulting cell-level classification is stored in:
    adata.obs["cnv_annotation"]

The script also stores the intermediate SCEVAN result in:
    adata.obs["scevan_class"]

and the cluster-level SCEVAN assignment in:
    adata.obs["cluster_scevan_class"].

Input
-----
An AnnData file containing a Leiden clustering in:

    adata.obs["leiden"]

and a SCEVAN result table containing:

    - cell IDs as the index
    - a column named "class"

The SCEVAN "class" column is expected to contain categories such as
"tumor", "filtered", and non-tumor/non-filtered classes.

Output
------
The updated AnnData object and a tab-delimited cell-level annotation file.

Example
-------
python classify_malignant_by_cnv.py \
    --input-h5ad data/processed/sample.h5ad \
    --scevan-results results/sample/scevan_results.csv \
    --output-h5ad results/sample/sample_cnv_classified.h5ad \
    --output-annotations results/sample/sample_cnv_annotation.txt \
    --tumor-threshold 40

Notes
-----
The 40% threshold reproduces the criterion used in the original analysis.
No marker-based annotation is performed in this script.
"""


from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


def read_scevan_results(scevan_file: Path) -> pd.DataFrame:
    """
    Read and validate SCEVAN classification results.

    Parameters
    ----------
    scevan_file
        Path to the SCEVAN CSV result file.

    Returns
    -------
    pandas.DataFrame
        SCEVAN results indexed by cell ID.
    """
    if not scevan_file.exists():
        raise FileNotFoundError(
            f"SCEVAN results file does not exist: {scevan_file}"
        )

    scevan = pd.read_csv(scevan_file, index_col=0)

    if "class" not in scevan.columns:
        raise ValueError(
            "The SCEVAN results file must contain a column named 'class'."
        )

    # Missing SCEVAN values are treated as zero to reproduce the
    # original analysis workflow.
    scevan["class"] = scevan["class"].fillna(0).astype(str)

    return scevan


def add_scevan_predictions(
    adata: sc.AnnData,
    scevan: pd.DataFrame,
) -> int:
    """
    Add cell-level SCEVAN classifications to AnnData.

    Only cells present in both datasets are assigned a prediction.

    Parameters
    ----------
    adata
        Main AnnData object.
    scevan
        SCEVAN classification table.

    Returns
    -------
    int
        Number of cells shared between the two datasets.
    """
    common_cells = adata.obs.index.intersection(scevan.index)

    if len(common_cells) == 0:
        raise ValueError(
            "No cell IDs are shared between the AnnData object "
            "and the SCEVAN results."
        )

    adata.obs["scevan_class"] = pd.NA

    adata.obs.loc[common_cells, "scevan_class"] = (
        scevan.loc[common_cells, "class"].values
    )

    print(
        f"Matched SCEVAN classifications for "
        f"{len(common_cells):,} of {adata.n_obs:,} cells."
    )

    return len(common_cells)


def classify_clusters_by_cnv(
    adata: sc.AnnData,
    tumor_threshold: float = 40.0,
) -> dict[str, str]:
    """
    Assign a CNV-based malignant/non-malignant classification to each
    Leiden cluster.

    A cluster is classified as malignant when either:

    1. The percentage of SCEVAN "tumor" cells is greater than the
       specified threshold, or

    2. The combined percentage of "tumor" and "filtered" cells is
       greater than the threshold.

    Otherwise, the most frequent SCEVAN class is retained at the
    cluster level and subsequently classified as non-malignant.

    Parameters
    ----------
    adata
        AnnData object containing "leiden" and "scevan_class".
    tumor_threshold
        Percentage threshold for malignant cluster assignment.

    Returns
    -------
    dict
        Mapping from Leiden cluster IDs to cluster-level classification.
    """
    required_columns = {"leiden", "scevan_class"}
    missing_columns = required_columns - set(adata.obs.columns)

    if missing_columns:
        raise ValueError(
            "Missing required columns in adata.obs: "
            + ", ".join(sorted(missing_columns))
        )

    cluster_to_classification: dict[str, str] = {}

    print("\nCluster-level CNV classification:")
    print("-" * 72)

    for cluster, cells in adata.obs.groupby(
        "leiden",
        observed=True,
    ):
        # Remove cells without a SCEVAN prediction.
        scevan_classes = cells["scevan_class"].dropna()

        if scevan_classes.empty:
            cluster_to_classification[str(cluster)] = "non-malignant"

            print(
                f"Cluster {cluster}: no SCEVAN predictions "
                f"-> non-malignant"
            )
            continue

        counts = scevan_classes.value_counts()
        total_cells = counts.sum()

        tumor_fraction = (
            counts.get("tumor", 0) / total_cells * 100
        )

        filtered_fraction = (
            counts.get("filtered", 0) / total_cells * 100
        )

        tumor_plus_filtered = (
            tumor_fraction + filtered_fraction
        )

        # Preserve the original CNV-based thresholding strategy.
        is_malignant = (
            tumor_fraction > tumor_threshold
            or tumor_plus_filtered > tumor_threshold
        )

        classification = (
            "malignant" if is_malignant else "non-malignant"
        )

        cluster_to_classification[str(cluster)] = classification

        print(
            f"Cluster {cluster}: "
            f"n={total_cells:,}; "
            f"tumor={tumor_fraction:.1f}%; "
            f"filtered={filtered_fraction:.1f}%; "
            f"tumor+filtered={tumor_plus_filtered:.1f}% "
            f"-> {classification}"
        )

    print("-" * 72)

    return cluster_to_classification


def assign_cell_level_classification(
    adata: sc.AnnData,
    cluster_classification: dict[str, str],
) -> None:
    """
    Transfer cluster-level CNV classifications to individual cells.

    The resulting annotation is stored in:
        adata.obs["cnv_annotation"]

    Parameters
    ----------
    adata
        AnnData object.
    cluster_classification
        Mapping from Leiden cluster ID to malignant/non-malignant status.
    """
    adata.obs["cluster_scevan_class"] = (
        adata.obs["leiden"]
        .astype(str)
        .map(cluster_classification)
    )

    adata.obs["cnv_annotation"] = (
        adata.obs["cluster_scevan_class"]
        .fillna("non-malignant")
        .astype(str)
    )

    # Make the classification categorical for efficient storage and
    # consistent downstream analysis.
    adata.obs["cnv_annotation"] = pd.Categorical(
        adata.obs["cnv_annotation"],
        categories=["non-malignant", "malignant"],
    )


def print_summary(adata: sc.AnnData) -> None:
    """
    Print a summary of the final cell-level classification.
    """
    counts = adata.obs["cnv_annotation"].value_counts()

    total = counts.sum()

    print("\nFinal cell-level classification:")
    print("-" * 40)

    for label in ["malignant", "non-malignant"]:
        count = int(counts.get(label, 0))
        percentage = (
            count / total * 100
            if total > 0
            else 0.0
        )

        print(
            f"{label:15s}: "
            f"{count:>8,} cells "
            f"({percentage:5.1f}%)"
        )

    print("-" * 40)


def save_results(
    adata: sc.AnnData,
    output_h5ad: Path,
    output_annotations: Path,
) -> None:
    """
    Save the annotated AnnData object and cell-level classifications.
    """
    output_h5ad.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_annotations.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Save the complete AnnData object with the new CNV-based
    # classification stored in .obs.
    adata.write(output_h5ad)

    # Save a simple two-column annotation table:
    # cell_id <tab> malignant/non-malignant
    annotations = adata.obs[
        ["cnv_annotation"]
    ].copy()

    annotations.to_csv(
        output_annotations,
        sep="\t",
        header=True,
        index=True,
    )

    print("\nOutput files:")
    print(f"  AnnData:     {output_h5ad}")
    print(f"  Annotations: {output_annotations}")


def run_cnv_classification(
    input_h5ad: Path,
    scevan_results: Path,
    output_h5ad: Path,
    output_annotations: Path,
    tumor_threshold: float = 40.0,
) -> None:
    """
    Run the complete CNV-based malignant classification workflow.
    """
    # ------------------------------------------------------------------
    # 1. Read processed single-cell data
    # ------------------------------------------------------------------
    print(f"Reading AnnData: {input_h5ad}")

    adata = sc.read_h5ad(input_h5ad)

    print(
        f"Dataset contains "
        f"{adata.n_obs:,} cells × {adata.n_vars:,} genes."
    )

    if "leiden" not in adata.obs.columns:
        raise ValueError(
            "The AnnData object must contain "
            "`adata.obs['leiden']` for cluster-level classification."
        )

    # ------------------------------------------------------------------
    # 2. Read SCEVAN CNV predictions
    # ------------------------------------------------------------------
    print(f"\nReading SCEVAN results: {scevan_results}")

    scevan = read_scevan_results(scevan_results)

    print(
        f"SCEVAN results contain "
        f"{len(scevan):,} cells."
    )

    # ------------------------------------------------------------------
    # 3. Add SCEVAN predictions to AnnData
    # ------------------------------------------------------------------
    add_scevan_predictions(
        adata,
        scevan,
    )

    # ------------------------------------------------------------------
    # 4. Classify Leiden clusters using SCEVAN predictions
    # ------------------------------------------------------------------
    cluster_classification = classify_clusters_by_cnv(
        adata,
        tumor_threshold=tumor_threshold,
    )

    # ------------------------------------------------------------------
    # 5. Transfer cluster classification to individual cells
    # ------------------------------------------------------------------
    assign_cell_level_classification(
        adata,
        cluster_classification,
    )

    # ------------------------------------------------------------------
    # 6. Print final classification summary
    # ------------------------------------------------------------------
    print_summary(adata)

    # ------------------------------------------------------------------
    # 7. Save results
    # ------------------------------------------------------------------
    save_results(
        adata=adata,
        output_h5ad=output_h5ad,
        output_annotations=output_annotations,
    )


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Classify single-cell RNA-seq cells as malignant or "
            "non-malignant using SCEVAN CNV predictions."
        )
    )

    parser.add_argument(
        "--input-h5ad",
        type=Path,
        required=True,
        help="Processed AnnData file containing Leiden clusters.",
    )

    parser.add_argument(
        "--scevan-results",
        type=Path,
        required=True,
        help=(
            "SCEVAN classification CSV file. "
            "Cell IDs must be the index and a 'class' column must be present."
        ),
    )

    parser.add_argument(
        "--output-h5ad",
        type=Path,
        required=True,
        help="Path for the output annotated AnnData file.",
    )

    parser.add_argument(
        "--output-annotations",
        type=Path,
        required=True,
        help="Path for the output cell-level annotation table.",
    )

    parser.add_argument(
        "--tumor-threshold",
        type=float,
        default=40.0,
        help=(
            "Percentage threshold used to classify a cluster as malignant "
            "based on tumor or tumor+filtered SCEVAN classifications. "
            "Default: 40."
        ),
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    run_cnv_classification(
        input_h5ad=args.input_h5ad,
        scevan_results=args.scevan_results,
        output_h5ad=args.output_h5ad,
        output_annotations=args.output_annotations,
        tumor_threshold=args.tumor_threshold,
    )

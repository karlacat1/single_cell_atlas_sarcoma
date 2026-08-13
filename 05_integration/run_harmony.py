#!/usr/bin/env python3
"""
Harmony integration of multi-sample single-cell RNA-seq data.

Selects highly variable genes in a sample-aware manner, regresses out total
counts and mitochondrial fraction, scales, computes PCA, and corrects the PCA
representation with Harmony using sample and 10x chip as batch variables. The
neighborhood graph, UMAP and Leiden clustering are computed on the corrected
representation in adata.obsm["X_pca_harmony"].

Input: normalized, log-transformed .h5ad with obs columns sample, 10x_chip,
total_counts and pct_counts_mt.

Output: integrated .h5ad.

Usage:
  python harmony_integration.py --input-h5ad processed.h5ad \
      --output-h5ad harmony.h5ad [--n-hvgs 3000] [--n-pcs 50] \
      [--n-neighbors 15] [--min-dist 0.5] [--leiden-resolution 0.5]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
import scanpy.external as sce


# =============================================================================
# Default parameters
# =============================================================================

DEFAULT_N_HVGS = 3000
DEFAULT_N_PCS = 50
DEFAULT_N_NEIGHBORS = 15
DEFAULT_MIN_DIST = 0.5
DEFAULT_LEIDEN_RESOLUTION = 0.5
DEFAULT_MAX_SCALE = 10.0
DEFAULT_RANDOM_SEED = 42

DEFAULT_MIN_MEAN = 0.0125
DEFAULT_MAX_MEAN = 3.0
DEFAULT_MIN_DISP = 0.5
DEFAULT_SPAN = 1.0


# =============================================================================
# Input validation
# =============================================================================

def validate_input(
    adata: sc.AnnData,
    sample_key: str,
    chip_key: str,
    total_counts_key: str,
    pct_mt_key: str,
) -> None:
    """
    Validate that all required metadata fields are present.
    """
    required_columns = {
        sample_key,
        chip_key,
        total_counts_key,
        pct_mt_key,
    }

    missing_columns = (
        required_columns - set(adata.obs.columns)
    )

    if missing_columns:
        raise ValueError(
            "The input AnnData object is missing the following "
            "required obs columns: "
            + ", ".join(sorted(missing_columns))
        )

    if adata.n_obs == 0:
        raise ValueError(
            "The AnnData object contains no cells."
        )

    if adata.n_vars == 0:
        raise ValueError(
            "The AnnData object contains no genes."
        )


# =============================================================================
# Highly variable gene selection
# =============================================================================

def select_hvgs(
    adata: sc.AnnData,
    sample_key: str,
    n_hvgs: int,
    min_mean: float,
    max_mean: float,
    min_disp: float,
    span: float,
) -> sc.AnnData:
    """
    Identify batch-aware highly variable genes and subset the AnnData object.

    The original workflow stores the full expression matrix in `.raw`
    before subsetting to HVGs.
    """
    print(
        f"\nSelecting {n_hvgs:,} highly variable genes "
        f"using '{sample_key}' as the batch key..."
    )

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3_paper",
        batch_key=sample_key,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        n_top_genes=n_hvgs,
        span=span,
        check_values=False,
    )

    n_hvgs_found = int(
        adata.var["highly_variable"].sum()
    )

    print(
        f"Highly variable genes identified: "
        f"{n_hvgs_found:,}"
    )

    # Preserve the complete expression matrix before restricting the
    # working object to HVGs.
    adata.raw = adata

    adata = adata[
        :,
        adata.var["highly_variable"],
    ].copy()

    print(
        f"Working matrix after HVG selection: "
        f"{adata.n_obs:,} cells × "
        f"{adata.n_vars:,} genes."
    )

    return adata


# =============================================================================
# PCA preprocessing
# =============================================================================

def compute_pca(
    adata: sc.AnnData,
    total_counts_key: str,
    pct_mt_key: str,
    max_scale: float,
    random_seed: int,
) -> None:
    """
    Regress out technical covariates, scale the data, and compute PCA.
    """
    print(
        "\nRegressing out total counts and mitochondrial fraction..."
    )

    sc.pp.regress_out(
        adata,
        [
            total_counts_key,
            pct_mt_key,
        ],
    )

    print(
        f"Scaling expression values "
        f"(maximum absolute value = {max_scale})..."
    )

    sc.pp.scale(
        adata,
        max_value=max_scale,
    )

    print("Computing PCA...")

    sc.tl.pca(
        adata,
        svd_solver="arpack",
        random_state=random_seed,
    )

    print(
        f"PCA completed with "
        f"{adata.obsm['X_pca'].shape[1]} components."
    )


# =============================================================================
# Harmony integration
# =============================================================================

def run_harmony(
    adata: sc.AnnData,
    sample_key: str,
    chip_key: str,
) -> None:
    """
    Apply Harmony to the PCA representation.

    Harmony corrects the PCA embedding using sample and 10x-chip
    information. The corrected representation is stored by Scanpy as:

        adata.obsm["X_pca_harmony"]
    """
    print(
        "\nApplying Harmony integration using:"
    )
    print(
        f"  sample:   {sample_key}"
    )
    print(
        f"  10x chip: {chip_key}"
    )

    # Match the original workflow: use string-valued chip identifiers.
    adata.obs[chip_key] = (
        adata.obs[chip_key]
        .astype(str)
    )

    sce.pp.harmony_integrate(
        adata,
        [
            sample_key,
            chip_key,
        ],
    )

    if "X_pca_harmony" not in adata.obsm:
        raise RuntimeError(
            "Harmony integration completed without generating "
            "`adata.obsm['X_pca_harmony']`."
        )

    print(
        "Harmony integration completed."
    )


# =============================================================================
# Neighborhood graph, UMAP, and clustering
# =============================================================================

def compute_graph_and_clusters(
    adata: sc.AnnData,
    n_neighbors: int,
    n_pcs: int,
    min_dist: float,
    leiden_resolution: float,
    random_seed: int,
) -> None:
    """
    Compute neighbors, UMAP, and Leiden clustering from the Harmony embedding.
    """
    harmony_dimensions = adata.obsm[
        "X_pca_harmony"
    ].shape[1]

    if n_pcs > harmony_dimensions:
        raise ValueError(
            f"Requested {n_pcs} PCs, but the Harmony representation "
            f"contains only {harmony_dimensions} dimensions."
        )

    print(
        f"\nComputing neighbors using the first "
        f"{n_pcs} Harmony dimensions..."
    )

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep="X_pca_harmony",
        random_state=random_seed,
    )

    print("Computing UMAP...")

    sc.tl.umap(
        adata,
        random_state=random_seed,
        min_dist=min_dist,
    )

    print(
        f"Computing Leiden clusters "
        f"(resolution={leiden_resolution})..."
    )

    sc.tl.leiden(
        adata,
        random_state=random_seed,
        resolution=leiden_resolution,
    )

    print(
        f"Identified "
        f"{adata.obs['leiden'].nunique()} Leiden clusters."
    )


# =============================================================================
# Save results
# =============================================================================

def save_results(
    adata: sc.AnnData,
    output_h5ad: Path,
    output_annotations: Path | None,
) -> None:
    """
    Save the integrated AnnData object and optional cluster annotations.
    """
    output_h5ad.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        f"\nSaving integrated AnnData: "
        f"{output_h5ad}"
    )

    adata.write(
        output_h5ad,
    )

    if output_annotations is not None:
        output_annotations.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        cluster_annotations = adata.obs[
            ["sample", "10x_chip", "leiden"]
        ].copy()

        cluster_annotations.to_csv(
            output_annotations,
            sep="\t",
            index=True,
        )

        print(
            f"Cluster assignments saved to: "
            f"{output_annotations}"
        )

    print("Results saved successfully.")


# =============================================================================
# Main workflow
# =============================================================================

def run_harmony_integration(
    input_h5ad: Path,
    output_h5ad: Path,
    output_annotations: Path | None,
    sample_key: str,
    chip_key: str,
    total_counts_key: str,
    pct_mt_key: str,
    n_hvgs: int,
    n_pcs: int,
    n_neighbors: int,
    min_dist: float,
    leiden_resolution: float,
    max_scale: float,
    min_mean: float,
    max_mean: float,
    min_disp: float,
    span: float,
    random_seed: int,
) -> None:
    """
    Run the complete Harmony integration workflow.
    """
    # -------------------------------------------------------------------------
    # 1. Read input data
    # -------------------------------------------------------------------------
    print(
        f"Reading AnnData: {input_h5ad}"
    )

    adata = sc.read_h5ad(
        input_h5ad,
    )

    print(
        f"Loaded dataset: "
        f"{adata.n_obs:,} cells × "
        f"{adata.n_vars:,} genes."
    )

    # -------------------------------------------------------------------------
    # 2. Validate input
    # -------------------------------------------------------------------------
    validate_input(
        adata=adata,
        sample_key=sample_key,
        chip_key=chip_key,
        total_counts_key=total_counts_key,
        pct_mt_key=pct_mt_key,
    )

    # -------------------------------------------------------------------------
    # 3. Select batch-aware highly variable genes
    # -------------------------------------------------------------------------
    adata = select_hvgs(
        adata=adata,
        sample_key=sample_key,
        n_hvgs=n_hvgs,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        span=span,
    )

    # -------------------------------------------------------------------------
    # 4. Regress technical variables, scale, and compute PCA
    # -------------------------------------------------------------------------
    compute_pca(
        adata=adata,
        total_counts_key=total_counts_key,
        pct_mt_key=pct_mt_key,
        max_scale=max_scale,
        random_seed=random_seed,
    )

    # -------------------------------------------------------------------------
    # 5. Harmony integration
    # -------------------------------------------------------------------------
    run_harmony(
        adata=adata,
        sample_key=sample_key,
        chip_key=chip_key,
    )

    # -------------------------------------------------------------------------
    # 6. Neighbors, UMAP, and Leiden clustering
    # -------------------------------------------------------------------------
    compute_graph_and_clusters(
        adata=adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        min_dist=min_dist,
        leiden_resolution=leiden_resolution,
        random_seed=random_seed,
    )

    # -------------------------------------------------------------------------
    # 7. Save results
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
            "Perform Harmony integration and downstream graph-based "
            "clustering of multi-sample single-cell RNA-seq data."
        )
    )

    # -------------------------------------------------------------------------
    # Input/output
    # -------------------------------------------------------------------------

    parser.add_argument(
        "--input-h5ad",
        type=Path,
        required=True,
        help="Input processed AnnData file.",
    )

    parser.add_argument(
        "--output-h5ad",
        type=Path,
        required=True,
        help="Output AnnData file containing the Harmony integration.",
    )

    parser.add_argument(
        "--output-annotations",
        type=Path,
        default=None,
        help=(
            "Optional tab-delimited file containing sample, 10x-chip, "
            "and Leiden assignments."
        ),
    )

    # -------------------------------------------------------------------------
    # Metadata columns
    # -------------------------------------------------------------------------

    parser.add_argument(
        "--sample-key",
        default="sample",
        help="Sample/batch column. Default: sample.",
    )

    parser.add_argument(
        "--chip-key",
        default="10x_chip",
        help="10x sequencing chip column. Default: 10x_chip.",
    )

    parser.add_argument(
        "--total-counts-key",
        default="total_counts",
        help="Total-counts column. Default: total_counts.",
    )

    parser.add_argument(
        "--pct-mt-key",
        default="pct_counts_mt",
        help="Mitochondrial percentage column. Default: pct_counts_mt.",
    )

    # -------------------------------------------------------------------------
    # HVG parameters
    # -------------------------------------------------------------------------

    parser.add_argument(
        "--n-hvgs",
        type=int,
        default=DEFAULT_N_HVGS,
        help="Number of highly variable genes. Default: 2000.",
    )

    parser.add_argument(
        "--min-mean",
        type=float,
        default=DEFAULT_MIN_MEAN,
        help="Minimum mean expression for HVG selection.",
    )

    parser.add_argument(
        "--max-mean",
        type=float,
        default=DEFAULT_MAX_MEAN,
        help="Maximum mean expression for HVG selection.",
    )

    parser.add_argument(
        "--min-disp",
        type=float,
        default=DEFAULT_MIN_DISP,
        help="Minimum dispersion for HVG selection.",
    )

    parser.add_argument(
        "--span",
        type=float,
        default=DEFAULT_SPAN,
        help="Span parameter for HVG selection.",
    )

    # -------------------------------------------------------------------------
    # PCA / Harmony / clustering
    # -------------------------------------------------------------------------

    parser.add_argument(
        "--n-pcs",
        type=int,
        default=DEFAULT_N_PCS,
        help=(
            "Number of Harmony dimensions used for graph construction. "
            "Default: 50."
        ),
    )

    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=DEFAULT_N_NEIGHBORS,
        help="Number of neighbors. Default: 15.",
    )

    parser.add_argument(
        "--min-dist",
        type=float,
        default=DEFAULT_MIN_DIST,
        help="Minimum UMAP distance. Default: 0.5.",
    )

    parser.add_argument(
        "--leiden-resolution",
        type=float,
        default=DEFAULT_LEIDEN_RESOLUTION,
        help="Leiden resolution. Default: 0.5.",
    )

    parser.add_argument(
        "--max-scale",
        type=float,
        default=DEFAULT_MAX_SCALE,
        help="Maximum absolute scaled expression value. Default: 10.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_RANDOM_SEED,
        help="Random seed. Default: 42.",
    )

    return parser.parse_args()


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    args = parse_arguments()

    run_harmony_integration(
        input_h5ad=args.input_h5ad,
        output_h5ad=args.output_h5ad,
        output_annotations=args.output_annotations,
        sample_key=args.sample_key,
        chip_key=args.chip_key,
        total_counts_key=args.total_counts_key,
        pct_mt_key=args.pct_mt_key,
        n_hvgs=args.n_hvgs,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        leiden_resolution=args.leiden_resolution,
        max_scale=args.max_scale,
        min_mean=args.min_mean,
        max_mean=args.max_mean,
        min_disp=args.min_disp,
        span=args.span,
        random_seed=args.seed,
    )

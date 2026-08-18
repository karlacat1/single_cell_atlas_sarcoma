#!/usr/bin/env python3
"""
Harmony integration of multi-sample single-cell RNA-seq data.

This script performs sample-aware highly variable gene selection, prepares the
expression matrix for dimensionality reduction, computes PCA, and applies
Harmony to correct the PCA representation for specified batch variables.

The corrected embedding is stored in:

    adata.obsm["X_pca_harmony"]

Input
-----
A processed AnnData (.h5ad) file containing normalized, log-transformed
expression data and the metadata required for batch correction.

By default, the following columns are expected in `adata.obs`:

    sample
    10x_chip
    total_counts
    pct_counts_mt

Output
------
An AnnData (.h5ad) file containing the Harmony-corrected PCA representation.

Usage
-----
python harmony_integration.py \
    --input input.h5ad \
    --output output.h5ad

Metadata column names can be changed with the corresponding command-line
arguments if necessary.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import harmonypy as hm
import pandas as pd
import scanpy as sc
import scanpy.external as sce


# ---------------------------------------------------------------------------
# Default parameters
# ---------------------------------------------------------------------------

DEFAULT_N_HVGS = 3000
DEFAULT_N_PCS = 50
DEFAULT_MAX_SCALE = 10.0
DEFAULT_RANDOM_SEED = 42

DEFAULT_MIN_MEAN = 0.0125
DEFAULT_MAX_MEAN = 3.0
DEFAULT_MIN_DISP = 0.5
DEFAULT_SPAN = 1.0


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def validate_input(adata: sc.AnnData, sample_key: str, chip_key: str, total_counts_key: str, pct_mt_key: str) -> None:
    """Validate the input AnnData object and required metadata."""

    if adata.n_obs == 0:
        raise ValueError("The input AnnData object contains no cells.")

    if adata.n_vars == 0:
        raise ValueError("The input AnnData object contains no genes.")

    required_columns = {sample_key, chip_key, total_counts_key, pct_mt_key}
    missing_columns = required_columns.difference(adata.obs.columns)

    if missing_columns:
        raise ValueError(
            "The input AnnData object is missing the following required "
            f"metadata columns: {', '.join(sorted(missing_columns))}"
        )


# ---------------------------------------------------------------------------
# Highly variable gene selection
# ---------------------------------------------------------------------------

def select_highly_variable_genes(adata: sc.AnnData, sample_key: str, n_hvgs: int, min_mean: float, max_mean: float, min_disp: float, span: float) -> sc.AnnData:
    """Select highly variable genes using sample-aware HVG selection."""

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3_paper",
        batch_key=sample_key,
        n_top_genes=n_hvgs,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        span=span,
        check_values=False,
    )

    if "highly_variable" not in adata.var:
        raise RuntimeError(
            "Highly variable gene selection did not produce the expected "
            "`adata.var['highly_variable']` annotation."
        )

    n_selected = int(adata.var["highly_variable"].sum())

    if n_selected == 0:
        raise RuntimeError("No highly variable genes were identified.")

    adata.raw = adata

    return adata[:, adata.var["highly_variable"]].copy()


# ---------------------------------------------------------------------------
# PCA preparation
# ---------------------------------------------------------------------------

def compute_pca(adata: sc.AnnData, total_counts_key: str, pct_mt_key: str, n_pcs: int, max_scale: float, random_seed: int) -> None:
    """Regress technical covariates, scale expression values, and compute PCA."""

    sc.pp.regress_out(adata, [total_counts_key, pct_mt_key])
    sc.pp.scale(adata, max_value=max_scale)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", random_state=random_seed)


# ---------------------------------------------------------------------------
# Harmony integration
# ---------------------------------------------------------------------------

def apply_harmony(adata: sc.AnnData, batch_keys: list[str]) -> None:
    """Apply Harmony to the PCA representation."""

    for key in batch_keys:
        adata.obs[key] = adata.obs[key].astype(str)

    harmony_out = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        batch_keys,
    )

    adata.obsm["X_pca_harmony"] = harmony_out.Z_corr
    #sce.pp.harmony_integrate(adata, key=batch_keys)

    if "X_pca_harmony" not in adata.obsm:
        raise RuntimeError(
            "Harmony integration did not produce "
            "`adata.obsm['X_pca_harmony']`."
        )


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def run_harmony(input_path: Path, output_path: Path, sample_key: str, chip_key: str, total_counts_key: str, pct_mt_key: str, n_hvgs: int, n_pcs: int, max_scale: float, min_mean: float, max_mean: float, min_disp: float, span: float, random_seed: int) -> None:
    """Run preprocessing, PCA, and Harmony integration."""

    adata = sc.read_h5ad(input_path)

    # Ensure regression covariates are numeric.
    adata.obs[total_counts_key] = pd.to_numeric(adata.obs[total_counts_key].astype(str), errors="raise")
    adata.obs[pct_mt_key] = pd.to_numeric(adata.obs[pct_mt_key].astype(str), errors="raise")

    validate_input(
        adata,
        sample_key=sample_key,
        chip_key=chip_key,
        total_counts_key=total_counts_key,
        pct_mt_key=pct_mt_key,
    )

    adata = select_highly_variable_genes(
        adata,
        sample_key=sample_key,
        n_hvgs=n_hvgs,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        span=span,
    )

    compute_pca(
        adata,
        total_counts_key=total_counts_key,
        pct_mt_key=pct_mt_key,
        n_pcs=n_pcs,
        max_scale=max_scale,
        random_seed=random_seed,
    )

    apply_harmony(adata, batch_keys=[sample_key, chip_key])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(output_path)


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Apply Harmony batch correction to a PCA representation "
            "of multi-sample single-cell RNA-seq data."
        )
    )

    parser.add_argument("--input", type=Path, required=True, help="Input processed AnnData (.h5ad) file.")
    parser.add_argument("--output", type=Path, required=True, help="Output AnnData (.h5ad) file.")
    parser.add_argument("--sample-key", default="sample", help="Metadata column containing sample identifiers.")
    parser.add_argument("--chip-key", default="10x_chip", help="Metadata column containing sequencing batch/chip identifiers.")
    parser.add_argument("--total-counts-key", default="total_counts", help="Metadata column containing total counts per cell.")
    parser.add_argument("--pct-mt-key", default="pct_counts_mt", help="Metadata column containing mitochondrial fraction.")

    parser.add_argument(
        "--n-hvgs",
        type=int,
        default=DEFAULT_N_HVGS,
        help=f"Number of highly variable genes. Default: {DEFAULT_N_HVGS}.",
    )

    parser.add_argument(
        "--n-pcs",
        type=int,
        default=DEFAULT_N_PCS,
        help=f"Number of principal components. Default: {DEFAULT_N_PCS}.",
    )

    parser.add_argument(
        "--max-scale",
        type=float,
        default=DEFAULT_MAX_SCALE,
        help=f"Maximum absolute scaled value. Default: {DEFAULT_MAX_SCALE}.",
    )

    parser.add_argument("--min-mean", type=float, default=DEFAULT_MIN_MEAN, help="Minimum mean expression for HVG selection.")
    parser.add_argument("--max-mean", type=float, default=DEFAULT_MAX_MEAN, help="Maximum mean expression for HVG selection.")
    parser.add_argument("--min-disp", type=float, default=DEFAULT_MIN_DISP, help="Minimum dispersion for HVG selection.")
    parser.add_argument("--span", type=float, default=DEFAULT_SPAN, help="Span parameter for HVG selection.")
    parser.add_argument("--seed", type=int, default=DEFAULT_RANDOM_SEED, help=f"Random seed. Default: {DEFAULT_RANDOM_SEED}.")

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_arguments()

    run_harmony(
        input_path=args.input,
        output_path=args.output,
        sample_key=args.sample_key,
        chip_key=args.chip_key,
        total_counts_key=args.total_counts_key,
        pct_mt_key=args.pct_mt_key,
        n_hvgs=args.n_hvgs,
        n_pcs=args.n_pcs,
        max_scale=args.max_scale,
        min_mean=args.min_mean,
        max_mean=args.max_mean,
        min_disp=args.min_disp,
        span=args.span,
        random_seed=args.seed,
    )

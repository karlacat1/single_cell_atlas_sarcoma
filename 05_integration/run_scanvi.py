#!/usr/bin/env python3
"""
scVI and scANVI integration of the merged sarcoma dataset.

Normalizes and log-transforms the counts, selects highly variable genes in a
batch-aware manner, then trains scVI on the SoupX-corrected counts with sample
and 10x chip as categorical covariates and mitochondrial percentage and total
counts as continuous covariates. scANVI is initialized from the trained scVI
model and uses cell_types as the label field, with "unclear" as the unlabeled
category, to give the semi-supervised integrated representation.

Input: .h5ad with a soupX_counts layer and obs columns sample, 10x_chip,
pct_counts_mt, total_counts and cell_types.

Output: .h5ad with X_scVI and X_scANVI in .obsm, scvi_normalized in .layers,
annotation_scanvi in .obs, and the log-normalized matrix retained in .raw;
plus the trained scVI and scANVI models.

Usage:
  python integrate_scvi_scanvi.py --input-h5ad raw.h5ad \
      --output-h5ad integrated.h5ad \
      --scvi-model models/scvi --scanvi-model models/scanvi
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
import scvi


# =============================================================================
# Default analysis parameters
# =============================================================================

DEFAULT_TARGET_SUM = 1e4
DEFAULT_N_HVGS = 3000
DEFAULT_MIN_MEAN = 0.0125
DEFAULT_MAX_MEAN = 3.0
DEFAULT_MIN_DISP = 0.5
DEFAULT_SPAN = 1.0

DEFAULT_SAMPLE_KEY = "sample"
DEFAULT_CHIP_KEY = "10x_chip"
DEFAULT_LABELS_KEY = "cell_types"
DEFAULT_UNLABELED_CATEGORY = "unclear"

DEFAULT_PCT_MT_KEY = "pct_counts_mt"
DEFAULT_TOTAL_COUNTS_KEY = "total_counts"
DEFAULT_COUNTS_LAYER = "soupX_counts"


# =============================================================================
# Input validation
# =============================================================================

def validate_input(
    adata: ad.AnnData,
    sample_key: str,
    chip_key: str,
    labels_key: str,
    pct_mt_key: str,
    total_counts_key: str,
    counts_layer: str,
    unlabeled_category: str,
) -> None:
    """
    Validate that all fields required for SCVI/SCANVI are present.
    """
    required_obs_columns = {
        sample_key,
        chip_key,
        labels_key,
        pct_mt_key,
        total_counts_key,
    }

    missing_obs = (
        required_obs_columns
        - set(adata.obs.columns)
    )

    if missing_obs:
        raise ValueError(
            "The input AnnData object is missing the following "
            "required obs columns: "
            + ", ".join(sorted(missing_obs))
        )

    if counts_layer not in adata.layers:
        raise ValueError(
            f"The required count layer "
            f"'{counts_layer}' is not present in adata.layers."
        )

    labels = (
        adata.obs[labels_key]
        .astype(str)
    )

    if unlabeled_category not in set(labels.unique()):
        raise ValueError(
            f"The label column '{labels_key}' does not contain "
            f"the required unlabeled category "
            f"'{unlabeled_category}'."
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
# Preprocessing
# =============================================================================

def preprocess_data(
    adata: ad.AnnData,
    sample_key: str,
    target_sum: float,
    n_hvgs: int,
    min_mean: float,
    max_mean: float,
    min_disp: float,
    span: float,
) -> None:
    """
    Normalize, log-transform, and identify batch-aware HVGs.

    The full normalized/log-transformed matrix is retained in the AnnData
    object. HVG status is stored in adata.var rather than subsetting the
    dataset.
    """
    print("\nNormalizing expression data...")

    sc.pp.normalize_total(
        adata,
        target_sum=target_sum,
    )

    print("Log-transforming expression data...")

    sc.pp.log1p(adata)

    print(
        f"Selecting {n_hvgs:,} highly variable genes "
        f"using batch-aware selection with '{sample_key}'..."
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

    # Preserve the normalized/log-transformed expression matrix.
    adata.raw = adata

    n_hvgs_found = int(
        adata.var["highly_variable"].sum()
    )

    print(
        f"Highly variable genes identified: "
        f"{n_hvgs_found:,}"
    )


# =============================================================================
# SCVI setup
# =============================================================================

def setup_scvi(
    adata: ad.AnnData,
    counts_layer: str,
    sample_key: str,
    chip_key: str,
    pct_mt_key: str,
    total_counts_key: str,
) -> None:
    """
    Register the AnnData object for SCVI.

    The raw/SoupX count layer is used as the model input.
    Sample and 10x chip are treated as categorical covariates.
    Mitochondrial percentage and total counts are continuous covariates.
    """
    print("\nSetting up SCVI...")

    # Preserve the original workflow's explicit string conversion.
    adata.obs[chip_key] = (
        adata.obs[chip_key]
        .astype(str)
    )

    scvi.model.SCVI.setup_anndata(
        adata=adata,
        layer=counts_layer,
        categorical_covariate_keys=[
            sample_key,
            chip_key,
        ],
        continuous_covariate_keys=[
            pct_mt_key,
            total_counts_key,
        ],
    )

    print("SCVI AnnData setup completed.")


# =============================================================================
# SCVI model
# =============================================================================

def train_scvi(
    adata: ad.AnnData,
    max_epochs: int | None,
    accelerator: str,
    devices: int | str | None,
):
    """
    Initialize and train the SCVI model.
    """
    print("\nInitializing SCVI model...")

    model_scvi = scvi.model.SCVI(
        adata,
    )

    model_scvi.view_anndata_setup()

    print("\nTraining SCVI model...")

    train_kwargs = {}

    if max_epochs is not None:
        train_kwargs["max_epochs"] = max_epochs

    if accelerator != "auto":
        train_kwargs["accelerator"] = accelerator

    if devices is not None:
        train_kwargs["devices"] = devices

    model_scvi.train(
        **train_kwargs,
    )

    print("SCVI training completed.")

    return model_scvi


# =============================================================================
# SCVI outputs
# =============================================================================

def add_scvi_outputs(
    adata: ad.AnnData,
    model_scvi,
    library_size: float,
) -> None:
    """
    Add the SCVI latent representation and normalized expression.
    """
    print("\nComputing SCVI latent representation...")

    latent = (
        model_scvi
        .get_latent_representation()
    )

    adata.obsm["X_scVI"] = latent

    print(
        f"SCVI latent representation: "
        f"{latent.shape[0]:,} cells × "
        f"{latent.shape[1]:,} dimensions."
    )

    print(
        "Computing SCVI normalized expression..."
    )

    normalized_expression = (
        model_scvi
        .get_normalized_expression(
            library_size=library_size,
        )
    )

    # Ensure that the result is stored as a dense numeric matrix
    # compatible with AnnData layers.
    if isinstance(
        normalized_expression,
        pd.DataFrame,
    ):
        normalized_expression = (
            normalized_expression.loc[
                adata.obs_names,
                adata.var_names,
            ].to_numpy()
        )

    adata.layers["scvi_normalized"] = (
        normalized_expression
    )

    print(
        "SCVI normalized expression stored in "
        "`adata.layers['scvi_normalized']`."
    )


# =============================================================================
# SCANVI setup
# =============================================================================

def initialize_scanvi(
    model_scvi,
    labels_key: str,
    unlabeled_category: str,
):
    """
    Initialize SCANVI from a pretrained SCVI model.
    """
    print(
        "\nInitializing SCANVI from the trained SCVI model..."
    )

    model_scanvi = (
        scvi.model.SCANVI.from_scvi_model(
            model_scvi,
            labels_key=labels_key,
            unlabeled_category=unlabeled_category,
        )
    )

    print("SCANVI initialization completed.")

    return model_scanvi


# =============================================================================
# SCANVI training
# =============================================================================

def train_scanvi(
    model_scanvi,
    max_epochs: int | None,
    accelerator: str,
    devices: int | str | None,
) -> None:
    """
    Train the SCANVI model.
    """
    model_scanvi.view_anndata_setup()

    print("\nTraining SCANVI model...")

    train_kwargs = {}

    if max_epochs is not None:
        train_kwargs["max_epochs"] = max_epochs

    if accelerator != "auto":
        train_kwargs["accelerator"] = accelerator

    if devices is not None:
        train_kwargs["devices"] = devices

    model_scanvi.train(
        **train_kwargs,
    )

    print("SCANVI training completed.")


# =============================================================================
# SCANVI outputs
# =============================================================================

def add_scanvi_outputs(
    adata: ad.AnnData,
    model_scanvi,
) -> None:
    """
    Add the SCANVI latent representation and cell-type predictions.
    """
    print("\nComputing SCANVI latent representation...")

    latent = (
        model_scanvi
        .get_latent_representation()
    )

    adata.obsm["X_scANVI"] = latent

    print(
        f"SCANVI latent representation: "
        f"{latent.shape[0]:,} cells × "
        f"{latent.shape[1]:,} dimensions."
    )

    print(
        "\nPredicting cell types with SCANVI..."
    )

    adata.obs["annotation_scanvi"] = (
        model_scanvi.predict(adata)
    )

    print(
        "\nSCANVI annotation summary:"
    )

    print(
        adata.obs["annotation_scanvi"]
        .value_counts(dropna=False)
    )


# =============================================================================
# Model saving
# =============================================================================

def save_model(
    model,
    output_path: Path,
    overwrite: bool,
    model_name: str,
) -> None:
    """
    Save a scvi-tools model using native serialization.
    """
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        f"\nSaving {model_name} model to: "
        f"{output_path}"
    )

    model.save(
        output_path,
        overwrite=overwrite,
    )

    print(
        f"{model_name} model saved."
    )


# =============================================================================
# Data saving
# =============================================================================

def save_annotated_data(
    adata: ad.AnnData,
    output_h5ad: Path,
) -> None:
    """
    Save the final AnnData object.
    """
    output_h5ad.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        f"\nSaving final AnnData object: "
        f"{output_h5ad}"
    )

    adata.write(
        output_h5ad,
    )

    print("AnnData saved.")


def save_scanvi_annotations(
    adata: ad.AnnData,
    output_annotations: Path,
) -> None:
    """
    Save SCANVI cell-level predictions.
    """
    output_annotations.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    annotation_table = adata.obs[
        ["annotation_scanvi"]
    ].copy()

    annotation_table.to_csv(
        output_annotations,
        sep="\t",
        header=True,
        index=True,
    )

    print(
        f"SCANVI annotations saved to: "
        f"{output_annotations}"
    )


# =============================================================================
# Main workflow
# =============================================================================

def run_analysis(
    input_h5ad: Path,
    output_h5ad: Path,
    scvi_model_path: Path,
    scanvi_model_path: Path,
    output_annotations: Path,
    sample_key: str,
    chip_key: str,
    labels_key: str,
    unlabeled_category: str,
    pct_mt_key: str,
    total_counts_key: str,
    counts_layer: str,
    target_sum: float,
    n_hvgs: int,
    min_mean: float,
    max_mean: float,
    min_disp: float,
    span: float,
    scvi_max_epochs: int | None,
    scanvi_max_epochs: int | None,
    accelerator: str,
    devices: int | str | None,
    library_size: float,
    seed: int,
    overwrite: bool,
) -> None:
    """
    Execute the complete SCVI/SCANVI workflow.
    """
    # -------------------------------------------------------------------------
    # 1. Set reproducibility seed
    # -------------------------------------------------------------------------
    scvi.settings.seed = seed

    print(
        f"Random seed: {seed}"
    )

    # -------------------------------------------------------------------------
    # 2. Read input data
    # -------------------------------------------------------------------------
    print(
        f"\nReading AnnData: {input_h5ad}"
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
    # 3. Validate input
    # -------------------------------------------------------------------------
    validate_input(
        adata=adata,
        sample_key=sample_key,
        chip_key=chip_key,
        labels_key=labels_key,
        pct_mt_key=pct_mt_key,
        total_counts_key=total_counts_key,
        counts_layer=counts_layer,
        unlabeled_category=unlabeled_category,
    )

    # -------------------------------------------------------------------------
    # 4. Normalize and identify HVGs
    # -------------------------------------------------------------------------
    #
    # This follows the original workflow. Importantly, the SoupX count
    # layer remains unchanged and continues to contain the count matrix
    # used by SCVI.
    #
    preprocess_data(
        adata=adata,
        sample_key=sample_key,
        target_sum=target_sum,
        n_hvgs=n_hvgs,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        span=span,
    )

    # -------------------------------------------------------------------------
    # 5. Configure SCVI
    # -------------------------------------------------------------------------
    setup_scvi(
        adata=adata,
        counts_layer=counts_layer,
        sample_key=sample_key,
        chip_key=chip_key,
        pct_mt_key=pct_mt_key,
        total_counts_key=total_counts_key,
    )

    # -------------------------------------------------------------------------
    # 6. Train SCVI
    # -------------------------------------------------------------------------
    model_scvi = train_scvi(
        adata=adata,
        max_epochs=scvi_max_epochs,
        accelerator=accelerator,
        devices=devices,
    )

    # -------------------------------------------------------------------------
    # 7. Generate SCVI latent representation and normalized expression
    # -------------------------------------------------------------------------
    add_scvi_outputs(
        adata=adata,
        model_scvi=model_scvi,
        library_size=library_size,
    )

    # -------------------------------------------------------------------------
    # 8. Save SCVI model
    # -------------------------------------------------------------------------
    save_model(
        model=model_scvi,
        output_path=scvi_model_path,
        overwrite=overwrite,
        model_name="SCVI",
    )

    # -------------------------------------------------------------------------
    # 9. Initialize SCANVI from SCVI
    # -------------------------------------------------------------------------
    model_scanvi = initialize_scanvi(
        model_scvi=model_scvi,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category,
    )

    # -------------------------------------------------------------------------
    # 10. Train SCANVI
    # -------------------------------------------------------------------------
    train_scanvi(
        model_scanvi=model_scanvi,
        max_epochs=scanvi_max_epochs,
        accelerator=accelerator,
        devices=devices,
    )

    # -------------------------------------------------------------------------
    # 11. Generate SCANVI representation and annotations
    # -------------------------------------------------------------------------
    add_scanvi_outputs(
        adata=adata,
        model_scanvi=model_scanvi,
    )

    # -------------------------------------------------------------------------
    # 12. Save SCANVI model
    # -------------------------------------------------------------------------
    save_model(
        model=model_scanvi,
        output_path=scanvi_model_path,
        overwrite=overwrite,
        model_name="SCANVI",
    )

    # -------------------------------------------------------------------------
    # 13. Save final AnnData
    # -------------------------------------------------------------------------
    save_annotated_data(
        adata=adata,
        output_h5ad=output_h5ad,
    )

    # -------------------------------------------------------------------------
    # 14. Save SCANVI predictions
    # -------------------------------------------------------------------------
    save_scanvi_annotations(
        adata=adata,
        output_annotations=output_annotations,
    )

    print(
        "\nSCVI/SCANVI workflow completed successfully."
    )


# =============================================================================
# Command-line interface
# =============================================================================

def parse_devices(
    value: str,
) -> int | str | None:
    """
    Convert the --devices argument to an integer, string, or None.
    """
    if value.lower() == "none":
        return None

    try:
        return int(value)
    except ValueError:
        return value


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Train SCVI and SCANVI models for multi-sample "
            "single-cell RNA-seq integration and annotation."
        )
    )

    # -------------------------------------------------------------------------
    # Input/output
    # -------------------------------------------------------------------------

    parser.add_argument(
        "--input-h5ad",
        type=Path,
        required=True,
        help="Input AnnData file.",
    )

    parser.add_argument(
        "--output-h5ad",
        type=Path,
        required=True,
        help="Final AnnData file containing SCVI/SCANVI results.",
    )

    parser.add_argument(
        "--scvi-model",
        type=Path,
        required=True,
        help="Output directory for the trained SCVI model.",
    )

    parser.add_argument(
        "--scanvi-model",
        type=Path,
        required=True,
        help="Output directory for the trained SCANVI model.",
    )

    parser.add_argument(
        "--output-annotations",
        type=Path,
        required=True,
        help="Output table containing SCANVI cell-type predictions.",
    )

    # -------------------------------------------------------------------------
    # AnnData fields
    # -------------------------------------------------------------------------

    parser.add_argument(
        "--sample-key",
        default=DEFAULT_SAMPLE_KEY,
        help=(
            "Column containing sample identifiers. "
            f"Default: {DEFAULT_SAMPLE_KEY}"
        ),
    )

    parser.add_argument(
        "--chip-key",
        default=DEFAULT_CHIP_KEY,
        help=(
            "Column containing 10x chip identifiers. "
            f"Default: {DEFAULT_CHIP_KEY}"
        ),
    )

    parser.add_argument(
        "--labels-key",
        default=DEFAULT_LABELS_KEY,
        help=(
            "Column containing cell-type labels. "
            f"Default: {DEFAULT_LABELS_KEY}"
        ),
    )

    parser.add_argument(
        "--unlabeled-category",
        default=DEFAULT_UNLABELED_CATEGORY,
        help=(
            "Category indicating unlabeled cells. "
            f"Default: {DEFAULT_UNLABELED_CATEGORY}"
        ),
    )

    parser.add_argument(
        "--pct-mt-key",
        default=DEFAULT_PCT_MT_KEY,
        help=(
            "Column containing mitochondrial percentage. "
            f"Default: {DEFAULT_PCT_MT_KEY}"
        ),
    )

    parser.add_argument(
        "--total-counts-key",
        default=DEFAULT_TOTAL_COUNTS_KEY,
        help=(
            "Column containing total counts per cell. "
            f"Default: {DEFAULT_TOTAL_COUNTS_KEY}"
        ),
    )

    parser.add_argument(
        "--counts-layer",
        default=DEFAULT_COUNTS_LAYER,
        help=(
            "Layer containing raw/SoupX-corrected counts used by SCVI. "
            f"Default: {DEFAULT_COUNTS_LAYER}"
        ),
    )

    # -------------------------------------------------------------------------
    # HVG selection
    # -------------------------------------------------------------------------

    parser.add_argument(
        "--target-sum",
        type=float,
        default=DEFAULT_TARGET_SUM,
        help="Target count for normalization. Default: 10000.",
    )

    parser.add_argument(
        "--n-hvgs",
        type=int,
        default=DEFAULT_N_HVGS,
        help="Number of highly variable genes. Default: 3000.",
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
    # Training
    # -------------------------------------------------------------------------

    parser.add_argument(
        "--scvi-max-epochs",
        type=int,
        default=None,
        help=(
            "Maximum SCVI training epochs. "
            "If omitted, scvi-tools determines the training duration."
        ),
    )

    parser.add_argument(
        "--scanvi-max-epochs",
        type=int,
        default=None,
        help=(
            "Maximum SCANVI training epochs. "
            "If omitted, scvi-tools determines the training duration."
        ),
    )

    parser.add_argument(
        "--accelerator",
        default="auto",
        help=(
            "Training accelerator, e.g. 'auto', 'cpu', or 'gpu'. "
            "Default: auto."
        ),
    )

    parser.add_argument(
        "--devices",
        default="none",
        help=(
            "Number of devices for training. "
            "Use 'none' to let scvi-tools determine the setting. "
            "Default: none."
        ),
    )

    parser.add_argument(
        "--library-size",
        type=float,
        default=1e4,
        help=(
            "Library size used when generating SCVI normalized expression. "
            "Default: 10000."
        ),
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed. Default: 0.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Overwrite existing SCVI/SCANVI model directories."
        ),
    )

    return parser.parse_args()


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    args = parse_arguments()

    run_analysis(
        input_h5ad=args.input_h5ad,
        output_h5ad=args.output_h5ad,
        scvi_model_path=args.scvi_model,
        scanvi_model_path=args.scanvi_model,
        output_annotations=args.output_annotations,
        sample_key=args.sample_key,
        chip_key=args.chip_key,
        labels_key=args.labels_key,
        unlabeled_category=args.unlabeled_category,
        pct_mt_key=args.pct_mt_key,
        total_counts_key=args.total_counts_key,
        counts_layer=args.counts_layer,
        target_sum=args.target_sum,
        n_hvgs=args.n_hvgs,
        min_mean=args.min_mean,
        max_mean=args.max_mean,
        min_disp=args.min_disp,
        span=args.span,
        scvi_max_epochs=args.scvi_max_epochs,
        scanvi_max_epochs=args.scanvi_max_epochs,
        accelerator=args.accelerator,
        devices=parse_devices(args.devices),
        library_size=args.library_size,
        seed=args.seed,
        overwrite=args.overwrite,
    )

"""Comprehensive GWAS visualization suite.

This module provides functions to generate complete sets of GWAS plots.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def _get_p_value_column(df: Any) -> str:
    """Find the p-value column name in a DataFrame.

    Checks common column name variants in priority order.

    Args:
        df: pandas DataFrame with GWAS results

    Returns:
        The column name containing p-values

    Raises:
        KeyError: If no p-value column is found
    """
    candidates = ["p_value", "P", "pvalue", "p_val", "PVALUE", "P_VALUE"]
    for col in candidates:
        if col in df.columns:
            return col
    raise KeyError(f"No p-value column found. Available columns: {list(df.columns)}")


def _get_chrom_column(df: Any) -> str:
    """Find the chromosome column name in a DataFrame."""
    candidates = ["CHROM", "CHR", "chrom", "chr", "chromosome", "CHROMOSOME"]
    for col in candidates:
        if col in df.columns:
            return col
    raise KeyError(f"No chromosome column found. Available columns: {list(df.columns)}")


def _get_pos_column(df: Any) -> str:
    """Find the position column name in a DataFrame."""
    candidates = ["POS", "BP", "pos", "bp", "position", "POSITION"]
    for col in candidates:
        if col in df.columns:
            return col
    raise KeyError(f"No position column found. Available columns: {list(df.columns)}")


def generate_all_plots(
    association_results: Union[Path, List[Any]],
    output_dir: Path,
    pca_file: Optional[Path] = None,
    kinship_file: Optional[Path] = None,
    vcf_file: Optional[Path] = None,
    significance_threshold: float = 5e-8,
) -> Dict[str, Any]:
    """Generate a comprehensive set of GWAS plots.

    Args:
        association_results: Path to association results file, or a list of results
        output_dir: Output directory for plots
        pca_file: Optional PCA file for population structure plots
        kinship_file: Optional kinship matrix file
        vcf_file: Optional VCF file for variant plots
        significance_threshold: P-value threshold for significance

    Returns:
        Dictionary with paths to generated plots, summary statistics,
        'status' ('completed' or 'failed'), and 'num_plots_generated' count.

    Example:
        >>> results = generate_all_plots(results_file, output_dir)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results: Dict[str, Any] = {"plots": {}, "statistics": {}, "files_created": []}

    # Handle empty list input gracefully
    if isinstance(association_results, list):
        if len(association_results) == 0:
            results["status"] = "completed"
            results["num_plots_generated"] = 0
            return results
        # If list of dicts, convert to DataFrame
        try:
            import pandas as pd

            df = pd.DataFrame(association_results)
        except Exception as e:
            logger.error(f"Failed to convert list to DataFrame: {e}")
            results["status"] = "failed"
            results["num_plots_generated"] = 0
            return results
    else:
        association_results = Path(association_results)
        try:
            # Load association results from file
            if association_results.suffix.lower() == ".csv":
                import pandas as pd

                df = pd.read_csv(association_results)
            elif association_results.suffix.lower() == ".tsv":
                import pandas as pd

                df = pd.read_csv(association_results, sep="\t")
            else:
                logger.error(f"Unsupported file format: {association_results.suffix}")
                results["status"] = "failed"
                results["num_plots_generated"] = 0
                return results
        except Exception as e:
            logger.error(f"Failed to load association results: {e}")
            results["status"] = "failed"
            results["num_plots_generated"] = 0
            return results

    try:
        # Resolve flexible column names
        p_col = _get_p_value_column(df)
        chrom_col = _get_chrom_column(df)
        pos_col = _get_pos_column(df)

        # Convert DataFrame rows to list-of-dicts format for general.py plot functions
        records = []
        for _, row in df.iterrows():
            record: Dict[str, Any] = {
                "chrom": str(row[chrom_col]),
                "pos": int(row[pos_col]),
                "p_value": float(row[p_col]),
            }
            # Include optional columns if present
            if "BETA" in df.columns:
                record["beta"] = float(row["BETA"])
            if "SE" in df.columns:
                record["se"] = float(row["SE"])
            records.append(record)

        # Generate Manhattan plot
        manhattan_file = output_dir / "manhattan_plot.png"
        try:
            from .general import manhattan_plot

            plot = manhattan_plot(records, output_path=manhattan_file, significance_threshold=significance_threshold)
            if plot:
                results["plots"]["manhattan"] = str(manhattan_file)
                results["files_created"].append(str(manhattan_file))
        except Exception as e:
            logger.warning(f"Failed to create Manhattan plot: {e}")

        # Generate Q-Q plot
        qq_file = output_dir / "qq_plot.png"
        try:
            from .general import qq_plot

            p_values = df[p_col].dropna().values
            plot = qq_plot(p_values.tolist(), output_path=qq_file)
            if plot:
                results["plots"]["qq"] = str(qq_file)
                results["files_created"].append(str(qq_file))
        except Exception as e:
            logger.warning(f"Failed to create Q-Q plot: {e}")

        # Generate regional plots for top hits
        try:
            from .general import regional_plot

            top_hits = df.nsmallest(5, p_col)

            for i, (_, hit) in enumerate(top_hits.iterrows()):
                chrom = str(hit[chrom_col])
                pos = int(hit[pos_col])
                regional_file = output_dir / f"regional_plot_chr{chrom}_{pos}.png"

                plot = regional_plot(records, chrom, pos - 500000, pos + 500000, output_path=regional_file)
                if plot:
                    results["plots"][f"regional_{i+1}"] = str(regional_file)
                    results["files_created"].append(str(regional_file))
        except Exception as e:
            logger.warning(f"Failed to create regional plots: {e}")

        # Generate effect size plot if BETA column present
        if "BETA" in df.columns:
            effect_file = output_dir / "effect_size_plot.png"
            try:
                from .general import effect_size_plot

                plot = effect_size_plot(records, output_path=effect_file)
                if plot:
                    results["plots"]["effect_size"] = str(effect_file)
                    results["files_created"].append(str(effect_file))
            except Exception as e:
                logger.warning(f"Failed to create effect size plot: {e}")

        # Generate PCA plot if PCA data available
        if pca_file and Path(pca_file).exists():
            try:
                pca_data = np.load(pca_file)
                pca_plot_file = output_dir / "pca_plot.png"

                from .general import pca_plot

                plot = pca_plot(pca_data, output_path=pca_plot_file)
                if plot:
                    results["plots"]["pca"] = str(pca_plot_file)
                    results["files_created"].append(str(pca_plot_file))
            except Exception as e:
                logger.warning(f"Failed to create PCA plot: {e}")

        # Generate kinship heatmap if kinship data available
        if kinship_file and Path(kinship_file).exists():
            try:
                kinship_data = np.load(kinship_file)
                kinship_plot_file = output_dir / "kinship_heatmap.png"

                from .general import kinship_heatmap

                plot = kinship_heatmap(kinship_data, output_path=kinship_plot_file)
                if plot:
                    results["plots"]["kinship"] = str(kinship_plot_file)
                    results["files_created"].append(str(kinship_plot_file))
            except Exception as e:
                logger.warning(f"Failed to create kinship heatmap: {e}")

        # Calculate basic statistics
        try:
            p_values = df[p_col].dropna().astype(float).values
            results["statistics"]["total_variants"] = len(df)
            results["statistics"]["significant_variants"] = int((p_values < significance_threshold).sum())
            if len(p_values) > 0:
                results["statistics"]["lambda_gc"] = float(np.median(p_values) / 0.456)
                results["statistics"]["min_p"] = float(np.min(p_values))
        except Exception as e:
            logger.warning(f"Failed to calculate statistics: {e}")

        num_plots = len(results["files_created"])
        logger.info(f"Generated {num_plots} plots in {output_dir}")
        results["status"] = "completed"
        results["num_plots_generated"] = num_plots

    except Exception as e:
        logger.error(f"Failed to generate plots: {e}")
        results["status"] = "failed"
        results["num_plots_generated"] = len(results["files_created"])

    return results

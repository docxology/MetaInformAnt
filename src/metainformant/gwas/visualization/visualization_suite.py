"""Comprehensive GWAS visualization suite.

This module provides functions to generate complete sets of GWAS plots.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


def generate_all_plots(
    association_results: Path,
    output_dir: Path,
    pca_file: Optional[Path] = None,
    kinship_file: Optional[Path] = None,
    vcf_file: Optional[Path] = None,
    significance_threshold: float = 5e-8,
) -> Dict[str, Any]:
    """Generate a comprehensive set of GWAS plots.

    Args:
        association_results: Path to association results file
        output_dir: Output directory for plots
        pca_file: Optional PCA file for population structure plots
        kinship_file: Optional kinship matrix file
        vcf_file: Optional VCF file for variant plots
        significance_threshold: P-value threshold for significance

    Returns:
        Dictionary with paths to generated plots and summary statistics

    Example:
        >>> results = generate_all_plots(results_file, output_dir)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {"plots": {}, "statistics": {}, "files_created": []}

    try:
        # Load association results
        if association_results.suffix.lower() == ".csv":
            import pandas as pd

            df = pd.read_csv(association_results)
        elif association_results.suffix.lower() == ".tsv":
            import pandas as pd

            df = pd.read_csv(association_results, sep="\t")
        else:
            logger.error(f"Unsupported file format: {association_results.suffix}")
            return results

        # Generate Manhattan plot
        manhattan_file = output_dir / "manhattan_plot.png"
        try:
            from .visualization import manhattan_plot

            plot = manhattan_plot(df, output_file=manhattan_file, significance_threshold=significance_threshold)
            if plot:
                results["plots"]["manhattan"] = str(manhattan_file)
                results["files_created"].append(str(manhattan_file))
        except Exception as e:
            logger.warning(f"Failed to create Manhattan plot: {e}")

        # Generate Q-Q plot
        qq_file = output_dir / "qq_plot.png"
        try:
            from .visualization import qq_plot

            p_values = df["P"].dropna().values
            plot = qq_plot(p_values, output_file=qq_file)
            if plot:
                results["plots"]["qq"] = str(qq_file)
                results["files_created"].append(str(qq_file))
        except Exception as e:
            logger.warning(f"Failed to create Q-Q plot: {e}")

        # Generate regional plots for top hits
        try:
            from .visualization import regional_plot

            top_hits = df.nsmallest(5, "P")  # Top 5 hits

            for i, (_, hit) in enumerate(top_hits.iterrows()):
                chrom = hit["CHR"]
                pos = hit["BP"]
                regional_file = output_dir / f"regional_plot_chr{chrom}_{pos}.png"

                plot = regional_plot(df, chrom, pos - 500000, pos + 500000, output_file=regional_file)
                if plot:
                    results["plots"][f"regional_{i+1}"] = str(regional_file)
                    results["files_created"].append(str(regional_file))
        except Exception as e:
            logger.warning(f"Failed to create regional plots: {e}")

        # Generate PCA plot if PCA data available
        if pca_file and pca_file.exists():
            try:
                import numpy as np

                pca_data = np.load(pca_file)
                pca_plot_file = output_dir / "pca_plot.png"

                from .visualization_population import pca_plot

                plot = pca_plot(pca_data, output_file=pca_plot_file)
                if plot:
                    results["plots"]["pca"] = str(pca_plot_file)
                    results["files_created"].append(str(pca_plot_file))
            except Exception as e:
                logger.warning(f"Failed to create PCA plot: {e}")

        # Calculate basic statistics
        try:
            p_values = df["P"].dropna().values
            results["statistics"]["total_variants"] = len(df)
            results["statistics"]["significant_variants"] = (p_values < significance_threshold).sum()
            results["statistics"]["lambda_gc"] = float(np.median(p_values) / 0.456)
            results["statistics"]["min_p"] = float(np.min(p_values))
        except Exception as e:
            logger.warning(f"Failed to calculate statistics: {e}")

        logger.info(f"Generated {len(results['files_created'])} plots in {output_dir}")

    except Exception as e:
        logger.error(f"Failed to generate plots: {e}")

    return results

"""GWAS visualization utilities.

This module provides functions for creating GWAS visualization plots,
including Manhattan plots, Q-Q plots, and regional association plots.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union, cast

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Import matplotlib with graceful fallback
try:
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    logger.warning("matplotlib not available, visualization functions will return None")

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    logger.warning("numpy not available, some visualizations may not work")

try:
    from scipy import stats as _scipy_stats

    HAS_SCIPY = True
except ImportError:  # pragma: no cover - used only in lean environments
    _scipy_stats = None
    HAS_SCIPY = False

EXPECTED_MEDIAN_CHI2_1DF = 0.454936423119572


def _coerce_float(value: Any, default: Optional[float] = None) -> Optional[float]:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    if not math.isfinite(result):
        return default
    return result


def _chrom_sort_key(chrom: Any) -> tuple:
    """Sort chromosome/contig names naturally across numeric and accession IDs."""
    text = str(chrom).strip()
    lowered = text.lower()
    if lowered.startswith("chr"):
        lowered = lowered[3:]

    if lowered.isdigit():
        return (0, int(lowered), text)

    special = {"x": 23, "y": 24, "m": 25, "mt": 25}
    if lowered in special:
        return (1, special[lowered], text)

    numbers = [int(part) for part in re.findall(r"\d+", lowered)]
    if numbers:
        return (2, numbers, lowered, text)
    return (3, lowered, text)


def _neg_log10_p(p_value: Any, cap: float = 300.0) -> float:
    p = _coerce_float(p_value)
    if p is None:
        return 0.0
    if p <= 0:
        return cap
    if p > 1:
        return 0.0
    return min(-math.log10(max(p, 10 ** (-cap))), cap)


def _normalise_gwas_results(results: Union[List[Dict[str, Any]], Dict[str, Any]]) -> List[Dict[str, Any]]:
    if isinstance(results, dict):
        source_rows: Iterable[Any] = [results]
    else:
        source_rows = list(results or [])

    rows: List[Dict[str, Any]] = []
    for original_index, result in enumerate(source_rows):
        if not isinstance(result, dict):
            continue
        chrom = str(result.get("chrom", result.get("chromosome", "1")))
        pos = _coerce_float(result.get("pos", result.get("position", 0)), 0.0)
        p_val = result.get("p_value", result.get("pval", result.get("pvalue", 1.0)))
        rows.append(
            {
                "chrom": chrom,
                "pos": pos if pos is not None else 0.0,
                "p_value": p_val,
                "neg_log_p": _neg_log10_p(p_val),
                "original_index": original_index,
                "source": result,
            }
        )
    rows.sort(key=lambda row: (_chrom_sort_key(row["chrom"]), row["pos"], row["original_index"]))
    return rows


def _compute_genome_axis(rows: List[Dict[str, Any]]) -> tuple[Dict[str, float], Dict[str, float], float]:
    """Assign cumulative x positions using observed contig extents."""
    if not rows:
        return {}, {}, 0.0

    chrom_order: List[str] = []
    chrom_bounds: Dict[str, tuple[float, float]] = {}
    for row in rows:
        chrom = row["chrom"]
        pos = float(row["pos"])
        if chrom not in chrom_bounds:
            chrom_order.append(chrom)
            chrom_bounds[chrom] = (pos, pos)
        else:
            lo, hi = chrom_bounds[chrom]
            chrom_bounds[chrom] = (min(lo, pos), max(hi, pos))

    largest_extent = max(max(hi - lo, hi, 1.0) for lo, hi in chrom_bounds.values())
    gap = max(1_000.0, min(largest_extent * 0.01, 5_000_000.0))

    offsets: Dict[str, float] = {}
    centers: Dict[str, float] = {}
    current = 0.0
    for chrom in chrom_order:
        lo, hi = chrom_bounds[chrom]
        offsets[chrom] = current
        centers[chrom] = current + (lo + hi) / 2.0
        current += max(hi, lo + 1.0) + gap

    for row in rows:
        row["global_pos"] = offsets[row["chrom"]] + float(row["pos"])

    return offsets, centers, current


def _lambda_gc_from_pvalues(p_values: Sequence[float]) -> Optional[float]:
    """Compute λ_GC via the canonical genomic-control helper (p → χ²(1) → median)."""
    from metainformant.gwas.analysis.correction import lambda_gc_from_p_values

    return lambda_gc_from_p_values([float(p) for p in p_values if math.isfinite(float(p)) and 0 < float(p) <= 1])


def _qq_confidence_band(n: int, alpha: float = 0.05) -> tuple[Any, Any, Any]:
    ranks = np.arange(1, n + 1)
    expected = (ranks - 0.5) / n
    if HAS_SCIPY and _scipy_stats is not None:
        lower_p = _scipy_stats.beta.ppf(alpha / 2.0, ranks, n - ranks + 1)
        upper_p = _scipy_stats.beta.ppf(1.0 - alpha / 2.0, ranks, n - ranks + 1)
    else:
        mean = ranks / (n + 1.0)
        sd = np.sqrt((ranks * (n - ranks + 1.0)) / (((n + 1.0) ** 2) * (n + 2.0)))
        lower_p = mean - 1.96 * sd
        upper_p = mean + 1.96 * sd
    lower_p = np.clip(lower_p, 1e-300, 1.0)
    upper_p = np.clip(upper_p, 1e-300, 1.0)
    return -np.log10(expected), -np.log10(upper_p), -np.log10(lower_p)


def regional_plot(
    results: List[Dict[str, Any]], chrom: str, start: int, end: int, output_path: Optional[Union[str, Path]] = None
) -> Any:
    """Create a regional association plot.

    Args:
        results: GWAS results for the region
        chrom: Chromosome
        start: Start position
        end: End position
        output_path: Path to save the plot (optional)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create regional plot")
        return None

    logger.info(f"Creating regional plot for {chrom}:{start}-{end}")

    # Filter results to the specified region
    region_results = []
    for result in results:
        r_chrom = str(result.get("chrom", result.get("chromosome", "")))
        r_pos = result.get("pos", result.get("position", 0))
        if r_chrom == str(chrom) and start <= r_pos <= end:
            region_results.append(result)

    if not region_results:
        logger.warning(f"No results found in region {chrom}:{start}-{end}")
        return None

    # Extract positions and p-values
    positions = []
    p_values = []
    for result in region_results:
        pos = result.get("pos", result.get("position", 0))
        p_val = result.get("p_value", result.get("pval", 1.0))
        positions.append(pos)
        if p_val > 0:
            p_values.append(-math.log10(p_val))
        else:
            p_values.append(50)  # Cap very small p-values

    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot points
    ax.scatter(positions, p_values, s=20, alpha=0.7, color="blue")

    # Add significance threshold line
    threshold_line = -math.log10(5e-8)
    ax.axhline(y=threshold_line, color="red", linestyle="--", alpha=0.7, label="Genome-wide significance")

    # Labels and title
    ax.set_xlabel(f"Position on chromosome {chrom}")
    ax.set_ylabel("-log₁₀(p-value)")
    ax.set_title(f"Regional Association Plot: {chrom}:{start:,}-{end:,}")
    ax.set_xlim(start, end)
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved regional plot to {output_path}")

    return fig


def pca_plot(
    pca_result: tuple, output_path: Optional[Union[str, Path]] = None, explained_var: Optional[List[float]] = None
) -> Any:
    """Create PCA scatter plot.

    Args:
        pca_result: PCA results tuple (components, variance, loadings)
        output_path: Path to save the plot (optional)
        explained_var: Explained variance ratios

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create PCA plot")
        return None

    logger.info("Creating PCA plot")

    try:
        components, variance, loadings = pca_result

        if len(components) < 2:
            logger.warning("Need at least 2 PCA components for plotting")
            return None

        # Create 2D scatter plot of first two components
        fig, ax = plt.subplots(figsize=(10, 8))

        # Plot points
        ax.scatter(components[0], components[1], s=2, alpha=0.6, color="blue")

        # Labels and title
        pc1_var = explained_var[0] * 100 if explained_var and len(explained_var) > 0 else 0
        pc2_var = explained_var[1] * 100 if explained_var and len(explained_var) > 1 else 0

        ax.set_xlabel(f"PC1 ({pc1_var:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({pc2_var:.1f}% variance)")
        ax.set_title("PCA Scatter Plot")
        ax.grid(True, alpha=0.3)

        # Add explained variance text if available
        if explained_var and len(explained_var) >= 2:
            var_text = f"PC1: {pc1_var:.1f}%, PC2: {pc2_var:.1f}%"
            ax.text(
                0.02,
                0.98,
                var_text,
                transform=ax.transAxes,
                verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
            )

        plt.tight_layout()

        # Save if output path provided
        if output_path:
            output_path = Path(output_path)
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            logger.info(f"Saved PCA plot to {output_path}")

        return fig

    except (ValueError, IndexError, TypeError) as e:
        logger.error(f"Error creating PCA plot: {e}")
        return None


def kinship_heatmap(
    kinship_matrix: Union[np.ndarray, List[List[float]]], output_path: Optional[Union[str, Path]] = None
) -> Any:
    """Create kinship matrix heatmap.

    Args:
        kinship_matrix: Kinship matrix
        output_path: Path to save the plot (optional)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create kinship heatmap")
        return None

    logger.info("Creating kinship heatmap")

    try:
        # Convert to numpy array if needed
        if isinstance(kinship_matrix, list):
            kinship_matrix = np.array(kinship_matrix)

        # Create heatmap
        fig, ax = plt.subplots(figsize=(10, 8))

        # Plot heatmap
        im = ax.imshow(kinship_matrix, cmap="viridis", aspect="equal")

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Kinship coefficient")

        # Labels and title
        ax.set_title("Kinship Matrix Heatmap")
        ax.set_xlabel("Sample")
        ax.set_ylabel("Sample")

        plt.tight_layout()

        # Save if output path provided
        if output_path:
            output_path = Path(output_path)
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            logger.info(f"Saved kinship heatmap to {output_path}")

        return fig

    except Exception as e:
        logger.error(f"Error creating kinship heatmap: {e}")
        return None


def effect_size_plot(results: List[Dict[str, Any]], output_path: Optional[Union[str, Path]] = None) -> Any:
    """Create effect size distribution plot.

    Args:
        results: GWAS results with 'beta' field for effect sizes
        output_path: Path to save the plot (optional)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create effect size plot")
        return None

    logger.info("Creating effect size plot")

    # Extract effect sizes (beta values)
    effect_sizes = []
    for result in results:
        beta = result.get("beta", result.get("effect_size", None))
        if beta is not None:
            effect_sizes.append(float(beta))

    if not effect_sizes:
        logger.warning("No effect sizes found in results")
        return None

    effect_sizes = np.array(effect_sizes)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Histogram of effect sizes
    ax1.hist(effect_sizes, bins=50, edgecolor="black", alpha=0.7, color="steelblue")
    ax1.axvline(x=0, color="red", linestyle="--", alpha=0.7, label="Null effect")
    ax1.axvline(
        x=np.mean(effect_sizes), color="green", linestyle="-", alpha=0.7, label=f"Mean: {np.mean(effect_sizes):.4f}"
    )
    ax1.set_xlabel("Effect Size (Beta)")
    ax1.set_ylabel("Frequency")
    ax1.set_title("Effect Size Distribution")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Box plot
    ax2.boxplot(effect_sizes, vert=True)
    ax2.set_ylabel("Effect Size (Beta)")
    ax2.set_title("Effect Size Box Plot")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved effect size plot to {output_path}")

    return fig


def generate_all_plots(
    association_results: Union[str, Path, List[Dict[str, Any]]],
    output_dir: Union[str, Path],
    pca_file: Optional[Union[str, Path]] = None,
    kinship_file: Optional[Union[str, Path]] = None,
    vcf_file: Optional[Union[str, Path]] = None,
    significance_threshold: float = 5e-8,
) -> Dict[str, Any]:
    """Generate all GWAS visualization plots.

    Args:
        association_results: Path to association results or results data
        output_dir: Output directory for plots
        pca_file: Path to PCA results file
        kinship_file: Path to kinship matrix file
        vcf_file: Path to VCF file
        significance_threshold: Significance threshold

    Returns:
        Dictionary with plot file paths and metadata
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Generating all GWAS plots in {output_dir}")

    plots_generated: Dict[str, Any] = {}

    # Load association results from file or use directly
    results_data: List[Dict[str, Any]] = []
    association_path = Path(association_results) if isinstance(association_results, (str, Path)) else None

    if association_path and association_path.exists():
        import csv
        import json

        suffix = association_path.suffix.lower()
        if suffix == ".json":
            with open(association_path) as fh:
                loaded = json.load(fh)
                results_data = loaded if isinstance(loaded, list) else loaded.get("results", [])
        elif suffix in (".tsv", ".txt"):
            with open(association_path, newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    entry: Dict[str, Any] = {}
                    for k, v in row.items():
                        try:
                            entry[k] = float(v)
                        except (ValueError, TypeError):
                            entry[k] = v
                    results_data.append(entry)
        elif suffix == ".csv":
            with open(association_path, newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    entry = {}
                    for k, v in row.items():
                        try:
                            entry[k] = float(v)
                        except (ValueError, TypeError):
                            entry[k] = v
                    results_data.append(entry)

    # Manhattan plot
    if results_data:
        try:
            manhattan_path = output_dir / "manhattan_plot.png"
            fig = manhattan_plot(  # noqa: F405
                results_data, output_path=manhattan_path, significance_threshold=significance_threshold
            )
            if fig is not None:
                plots_generated["manhattan"] = str(manhattan_path)
                plt.close(fig)
        except Exception as e:
            logger.warning(f"Manhattan plot failed: {e}")

        # Q-Q plot
        try:
            p_vals = cast(
                List[float],
                [
                    r.get("p_value", r.get("pval", r.get("pvalue")))
                    for r in results_data
                    if r.get("p_value", r.get("pval", r.get("pvalue"))) is not None
                ],
            )
            if p_vals:
                qq_path = output_dir / "qq_plot.png"
                fig = qq_plot(p_vals, output_path=qq_path)  # noqa: F405
                if fig is not None:
                    plots_generated["qq"] = str(qq_path)
                    plt.close(fig)
        except Exception as e:
            logger.warning(f"Q-Q plot failed: {e}")

    # PCA plot
    if pca_file:
        try:
            import json as _json

            pca_path_obj = Path(pca_file)
            if pca_path_obj.exists():
                with open(pca_path_obj) as fh:
                    pca_data = _json.load(fh)
                components = pca_data.get("components", [])
                variance = pca_data.get("variance", [])
                loadings = pca_data.get("loadings", [])
                explained_var = pca_data.get("explained_variance", [])
                pca_output = output_dir / "pca_plot.png"
                fig = pca_plot((components, variance, loadings), output_path=pca_output, explained_var=explained_var)
                if fig is not None:
                    plots_generated["pca"] = str(pca_output)
                    plt.close(fig)
        except Exception as e:
            logger.warning(f"PCA plot failed: {e}")

    # Kinship heatmap
    if kinship_file:
        try:
            import json as _json

            kinship_path_obj = Path(kinship_file)
            if kinship_path_obj.exists():
                with open(kinship_path_obj) as fh:
                    kinship_data = _json.load(fh)
                matrix = kinship_data if isinstance(kinship_data, list) else kinship_data.get("matrix", [])
                kinship_output = output_dir / "kinship_plot.png"
                fig = kinship_heatmap(matrix, output_path=kinship_output)
                if fig is not None:
                    plots_generated["kinship"] = str(kinship_output)
                    plt.close(fig)
        except Exception as e:
            logger.warning(f"Kinship heatmap failed: {e}")

    logger.info(f"Generated {len(plots_generated)} plots: {list(plots_generated.keys())}")
    return plots_generated


def missingness_plot(vcf_data: Dict[str, Any], output_path: Optional[Union[str, Path]] = None) -> Any:
    """Create missingness visualization showing per-sample and per-variant missingness.

    Args:
        vcf_data: VCF data dictionary with 'variants' and 'genotypes' keys
        output_path: Path to save the plot (optional)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        logger.warning("matplotlib or numpy not available, cannot create missingness plot")
        return None

    logger.info("Creating missingness plot")

    genotypes = vcf_data.get("genotypes", [])
    if not genotypes:
        logger.warning("No genotype data found for missingness plot")
        return None

    genotypes = np.array(genotypes)
    n_variants, n_samples = genotypes.shape

    # Calculate missingness (assuming -1 or negative values indicate missing)
    sample_missingness = np.mean(genotypes < 0, axis=0) * 100
    variant_missingness = np.mean(genotypes < 0, axis=1) * 100

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Per-sample missingness
    ax1.bar(range(n_samples), sample_missingness, color="steelblue", alpha=0.7)
    ax1.axhline(y=5, color="red", linestyle="--", alpha=0.7, label="5% threshold")
    ax1.set_xlabel("Sample Index")
    ax1.set_ylabel("Missing Rate (%)")
    ax1.set_title(f"Per-Sample Missingness (n={n_samples})")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Per-variant missingness histogram
    ax2.hist(variant_missingness, bins=50, edgecolor="black", alpha=0.7, color="steelblue")
    ax2.axvline(x=5, color="red", linestyle="--", alpha=0.7, label="5% threshold")
    ax2.set_xlabel("Missing Rate (%)")
    ax2.set_ylabel("Number of Variants")
    ax2.set_title(f"Per-Variant Missingness Distribution (n={n_variants})")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved missingness plot to {output_path}")

    return fig


def functional_enrichment_plot(
    results: List[Dict[str, Any]], gff_path: Union[str, Path], output_path: Optional[Union[str, Path]] = None
) -> Any:
    """Create functional enrichment plot showing enrichment of significant variants in functional categories.

    Args:
        results: GWAS results with 'chrom', 'pos', and 'p_value' fields
        gff_path: Path to GFF annotation file
        output_path: Path to save the plot (optional)

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create functional enrichment plot")
        return None

    logger.info("Creating functional enrichment plot")

    # Count significant variants
    significant = [r for r in results if r.get("p_value", r.get("pval", 1.0)) < 5e-8]

    if not significant:
        logger.warning("No significant variants found for enrichment analysis")
        # Create a simple summary plot instead
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No genome-wide significant variants\n(p < 5e-8)", ha="center", va="center", fontsize=14)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.set_title("Functional Enrichment Analysis")

        if output_path:
            output_path = Path(output_path)
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            logger.info(f"Saved functional enrichment plot to {output_path}")

        return fig

    # Create summary plot of p-value distribution by chromosome
    chrom_counts: Dict[str, int] = {}
    for r in significant:
        chrom = str(r.get("chrom", r.get("chromosome", "unknown")))
        chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

    fig, ax = plt.subplots(figsize=(10, 6))

    chroms = sorted(chrom_counts.keys(), key=lambda x: int(x) if x.isdigit() else 999)
    counts = [chrom_counts[c] for c in chroms]

    ax.bar(range(len(chroms)), counts, color="steelblue", alpha=0.7)
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(chroms)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Number of Significant Variants")
    ax.set_title(f"Significant Variants by Chromosome (n={len(significant)})")
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved functional enrichment plot to {output_path}")

    return fig


# This module historically contained a second copy of the plotting functions.
# Keep the import path stable while ensuring every public function is the same
# object as the canonical split implementation.
from metainformant.gwas.visualization._general_impl import *  # noqa: E402,F401,F403

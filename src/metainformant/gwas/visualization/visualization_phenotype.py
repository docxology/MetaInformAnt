"""Phenotype visualization for GWAS.

This module provides plots for visualizing phenotype distributions,
genotype-phenotype relationships, and phenotype correlations with
principal components.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    plt = None

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None


def phenotype_distribution(
    phenotype_values: List[float],
    trait_name: str = "Trait",
    output_file: Optional[Union[str, Path]] = None,
    by_population: Optional[Dict[str, str]] = None,
    sample_ids: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Create a histogram with KDE overlay of phenotype values.

    Displays the distribution of a continuous phenotype trait, optionally
    stratified by population. Includes mean and median reference lines
    and basic descriptive statistics.

    Args:
        phenotype_values: List of numeric phenotype measurements.
        trait_name: Human-readable name of the trait for axis labels.
        output_file: Optional path to save the plot image.
        by_population: Optional mapping of sample_id -> population label
            for stratified distributions.
        sample_ids: Optional list of sample identifiers corresponding to
            phenotype_values (required when by_population is provided).

    Returns:
        Dictionary with status, descriptive statistics, and output path.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for phenotype distribution plot")
        return {"status": "skipped", "output_path": None, "reason": "matplotlib not available"}

    if not HAS_NUMPY:
        logger.warning("numpy not available for phenotype distribution plot")
        return {"status": "skipped", "output_path": None, "reason": "numpy not available"}

    try:
        values = np.array(phenotype_values, dtype=np.float64)
        if len(values) == 0:
            return {"status": "failed", "error": "No phenotype values provided", "output_path": None}

        # Basic statistics
        mean_val = float(np.mean(values))
        median_val = float(np.median(values))
        std_val = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
        n_samples = len(values)

        fig, ax = plt.subplots(figsize=(10, 7))

        if by_population is not None and sample_ids is not None:
            # Stratified distribution by population
            pop_data: Dict[str, List[float]] = {}
            for sid, val in zip(sample_ids, phenotype_values):
                pop = by_population.get(sid, "Unknown")
                if pop not in pop_data:
                    pop_data[pop] = []
                pop_data[pop].append(val)

            colors = plt.cm.tab10(np.linspace(0, 1, max(len(pop_data), 1)))
            for i, (pop_name, pop_vals) in enumerate(sorted(pop_data.items())):
                pop_arr = np.array(pop_vals)
                ax.hist(
                    pop_arr,
                    bins=max(10, int(np.sqrt(len(pop_arr)))),
                    alpha=0.5,
                    color=colors[i],
                    label=f"{pop_name} (n={len(pop_arr)})",
                    density=True,
                )
                # KDE overlay per population
                if len(pop_arr) > 1:
                    kde_x = np.linspace(float(pop_arr.min()), float(pop_arr.max()), 200)
                    bandwidth = 1.06 * float(np.std(pop_arr, ddof=1)) * len(pop_arr) ** (-1.0 / 5.0)
                    if bandwidth > 0:
                        kde_y = np.zeros_like(kde_x)
                        for v in pop_arr:
                            kde_y += np.exp(-0.5 * ((kde_x - v) / bandwidth) ** 2)
                        kde_y /= len(pop_arr) * bandwidth * np.sqrt(2 * np.pi)
                        ax.plot(kde_x, kde_y, color=colors[i], linewidth=2)

            ax.legend(title="Population", loc="upper right")
        else:
            # Single distribution
            n_bins = max(10, int(np.sqrt(len(values))))
            ax.hist(values, bins=n_bins, alpha=0.7, color="steelblue", edgecolor="black", density=True)

            # KDE overlay
            if len(values) > 1:
                kde_x = np.linspace(float(values.min()), float(values.max()), 200)
                bandwidth = 1.06 * std_val * n_samples ** (-1.0 / 5.0)
                if bandwidth > 0:
                    kde_y = np.zeros_like(kde_x)
                    for v in values:
                        kde_y += np.exp(-0.5 * ((kde_x - v) / bandwidth) ** 2)
                    kde_y /= n_samples * bandwidth * np.sqrt(2 * np.pi)
                    ax.plot(kde_x, kde_y, color="darkblue", linewidth=2, label="KDE")

        # Mean and median lines
        ax.axvline(x=mean_val, color="red", linestyle="--", linewidth=1.5, label=f"Mean: {mean_val:.3f}")
        ax.axvline(x=median_val, color="green", linestyle="-.", linewidth=1.5, label=f"Median: {median_val:.3f}")

        ax.set_xlabel(trait_name, fontsize=12)
        ax.set_ylabel("Density", fontsize=12)
        ax.set_title(f"{trait_name} Distribution (n={n_samples})", fontsize=14, pad=20)
        ax.legend(loc="upper right")
        ax.grid(True, alpha=0.3)

        # Stats annotation
        stats_text = f"Mean: {mean_val:.3f}\nStd: {std_val:.3f}\nn: {n_samples}"
        ax.text(
            0.02,
            0.98,
            stats_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
        )

        plt.tight_layout()

        if output_file:
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            logger.info(f"Saved phenotype distribution plot to {output_file}")

        plt.close(fig)
        return {
            "status": "success",
            "n_samples": n_samples,
            "mean": mean_val,
            "median": median_val,
            "std": std_val,
            "output_path": str(output_file) if output_file else None,
        }

    except Exception as e:
        logger.error(f"Error creating phenotype distribution plot: {e}")
        return {"status": "failed", "error": str(e), "output_path": None}


def phenotype_correlation_matrix(
    phenotypes_dict: Dict[str, List[float]],
    output_file: Optional[Union[str, Path]] = None,
    method: str = "pearson",
) -> Dict[str, Any]:
    """Create a correlation heatmap across multiple phenotype traits.

    Computes pairwise correlation between all provided traits and renders
    a heatmap with annotated correlation coefficients.

    Args:
        phenotypes_dict: Mapping of trait_name -> list of numeric values.
            All value lists must have the same length.
        output_file: Optional path to save the plot image.
        method: Correlation method, currently supports "pearson" (default).

    Returns:
        Dictionary with status, correlation matrix, and output path.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for correlation matrix plot")
        return {"status": "skipped", "output_path": None, "reason": "matplotlib not available"}

    if not HAS_NUMPY:
        logger.warning("numpy not available for correlation matrix plot")
        return {"status": "skipped", "output_path": None, "reason": "numpy not available"}

    try:
        if not phenotypes_dict or len(phenotypes_dict) < 2:
            return {"status": "failed", "error": "At least two traits are required", "output_path": None}

        trait_names = list(phenotypes_dict.keys())
        n_traits = len(trait_names)

        # Validate equal lengths
        lengths = [len(phenotypes_dict[t]) for t in trait_names]
        if len(set(lengths)) != 1:
            return {"status": "failed", "error": "All trait value lists must have the same length", "output_path": None}

        # Build data matrix: rows = samples, cols = traits
        data_matrix = np.column_stack([np.array(phenotypes_dict[t], dtype=np.float64) for t in trait_names])

        # Compute correlation matrix
        corr_matrix = np.corrcoef(data_matrix, rowvar=False)

        fig, ax = plt.subplots(figsize=(max(8, n_traits * 1.5), max(6, n_traits * 1.2)))

        # Heatmap
        im = ax.imshow(corr_matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("Correlation coefficient (r)", fontsize=11)

        # Annotate cells with r values
        for i in range(n_traits):
            for j in range(n_traits):
                r_val = corr_matrix[i, j]
                text_color = "white" if abs(r_val) > 0.6 else "black"
                ax.text(
                    j, i, f"{r_val:.2f}", ha="center", va="center", fontsize=10, color=text_color, fontweight="bold"
                )

        # Axis labels
        ax.set_xticks(range(n_traits))
        ax.set_xticklabels(trait_names, rotation=45, ha="right", fontsize=10)
        ax.set_yticks(range(n_traits))
        ax.set_yticklabels(trait_names, fontsize=10)
        ax.set_title("Phenotype Correlation Matrix", fontsize=14, pad=20)

        plt.tight_layout()

        if output_file:
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            logger.info(f"Saved correlation matrix plot to {output_file}")

        plt.close(fig)

        # Convert correlation matrix to serializable format
        corr_dict: Dict[str, Dict[str, float]] = {}
        for i, t1 in enumerate(trait_names):
            corr_dict[t1] = {}
            for j, t2 in enumerate(trait_names):
                corr_dict[t1][t2] = float(corr_matrix[i, j])

        return {
            "status": "success",
            "n_traits": n_traits,
            "trait_names": trait_names,
            "correlation_matrix": corr_dict,
            "output_path": str(output_file) if output_file else None,
        }

    except Exception as e:
        logger.error(f"Error creating correlation matrix plot: {e}")
        return {"status": "failed", "error": str(e), "output_path": None}


def genotype_phenotype_boxplot(
    genotypes: List[int],
    phenotypes: List[float],
    variant_id: str = "Variant",
    output_file: Optional[Union[str, Path]] = None,
    metadata: Optional[Dict] = None,
) -> Dict[str, Any]:
    """Create a box/violin plot of phenotype values by genotype class.

    Genotype classes are 0 (Ref/Ref), 1 (Ref/Alt), and 2 (Alt/Alt).
    Includes sample count per group and a manual F-test p-value computed
    from sum of squares between and within groups (no scipy dependency).

    Args:
        genotypes: List of integer genotype values (0, 1, or 2).
        phenotypes: List of numeric phenotype values.
        variant_id: Identifier for the variant being plotted.
        output_file: Optional path to save the plot image.
        metadata: Optional dictionary of additional variant metadata.

    Returns:
        Dictionary with status, group statistics, F-test results, and output path.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for genotype-phenotype boxplot")
        return {"status": "skipped", "output_path": None, "reason": "matplotlib not available"}

    if not HAS_NUMPY:
        logger.warning("numpy not available for genotype-phenotype boxplot")
        return {"status": "skipped", "output_path": None, "reason": "numpy not available"}

    try:
        if len(genotypes) != len(phenotypes):
            return {
                "status": "failed",
                "error": "genotypes and phenotypes must have the same length",
                "output_path": None,
            }

        if len(genotypes) == 0:
            return {"status": "failed", "error": "No data provided", "output_path": None}

        genotype_arr = np.array(genotypes, dtype=int)
        phenotype_arr = np.array(phenotypes, dtype=np.float64)

        genotype_labels = {0: "Ref/Ref (0)", 1: "Ref/Alt (1)", 2: "Alt/Alt (2)"}

        # Group phenotypes by genotype
        groups: Dict[int, np.ndarray] = {}
        for geno in [0, 1, 2]:
            mask = genotype_arr == geno
            if np.any(mask):
                groups[geno] = phenotype_arr[mask]

        if not groups:
            return {"status": "failed", "error": "No valid genotype groups found", "output_path": None}

        # Manual one-way ANOVA F-test
        grand_mean = float(np.mean(phenotype_arr))
        n_total = len(phenotype_arr)
        k = len(groups)  # number of groups

        ss_between = 0.0
        ss_within = 0.0
        for geno, group_vals in groups.items():
            group_mean = float(np.mean(group_vals))
            n_g = len(group_vals)
            ss_between += n_g * (group_mean - grand_mean) ** 2
            ss_within += float(np.sum((group_vals - group_mean) ** 2))

        df_between = k - 1
        df_within = n_total - k

        f_statistic = None
        p_value = None
        if df_between > 0 and df_within > 0 and ss_within > 0:
            ms_between = ss_between / df_between
            ms_within = ss_within / df_within
            f_statistic = ms_between / ms_within

            # Approximate p-value using the F-distribution via the regularized
            # incomplete beta function. We use a numerical approximation.
            # For F(d1, d2) with observed f_stat:
            # p = 1 - I_x(d1/2, d2/2), where x = d1*f / (d1*f + d2)
            x = (df_between * f_statistic) / (df_between * f_statistic + df_within)
            p_value = 1.0 - _regularized_incomplete_beta(x, df_between / 2.0, df_within / 2.0)

        fig, ax = plt.subplots(figsize=(8, 6))

        # Prepare box plot data
        box_data = []
        box_labels = []
        box_positions = []
        group_stats: Dict[str, Dict[str, Any]] = {}

        for pos, geno in enumerate([0, 1, 2]):
            if geno in groups:
                vals = groups[geno]
                box_data.append(vals)
                label = genotype_labels[geno]
                box_labels.append(f"{label}\nn={len(vals)}")
                box_positions.append(pos)
                group_stats[genotype_labels[geno]] = {
                    "n": len(vals),
                    "mean": float(np.mean(vals)),
                    "median": float(np.median(vals)),
                    "std": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                }

        # Box plot with individual points
        bp = ax.boxplot(box_data, positions=box_positions, widths=0.6, patch_artist=True, showmeans=True)

        colors_box = ["#4ECDC4", "#FFD93D", "#FF6B6B"]
        for i, patch in enumerate(bp["boxes"]):
            geno_idx = [0, 1, 2][i] if i < 3 else i
            patch.set_facecolor(colors_box[geno_idx % 3])
            patch.set_alpha(0.7)

        # Overlay individual points with jitter
        for i, (pos, vals) in enumerate(zip(box_positions, box_data)):
            jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(vals))
            ax.scatter(
                pos + jitter,
                vals,
                alpha=0.4,
                s=15,
                color="black",
                zorder=3,
            )

        ax.set_xticks(box_positions)
        ax.set_xticklabels(box_labels, fontsize=10)
        ax.set_ylabel("Phenotype Value", fontsize=12)
        ax.set_xlabel("Genotype", fontsize=12)

        title_str = f"Genotype-Phenotype: {variant_id}"
        if f_statistic is not None and p_value is not None:
            title_str += f"\nF={f_statistic:.2f}, p={p_value:.4g}"
        ax.set_title(title_str, fontsize=13, pad=15)
        ax.grid(True, alpha=0.3, axis="y")

        plt.tight_layout()

        if output_file:
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            logger.info(f"Saved genotype-phenotype boxplot to {output_file}")

        plt.close(fig)
        return {
            "status": "success",
            "variant_id": variant_id,
            "n_samples": n_total,
            "n_groups": k,
            "group_stats": group_stats,
            "f_statistic": f_statistic,
            "p_value": p_value,
            "output_path": str(output_file) if output_file else None,
        }

    except Exception as e:
        logger.error(f"Error creating genotype-phenotype boxplot: {e}")
        return {"status": "failed", "error": str(e), "output_path": None}


def _regularized_incomplete_beta(x: float, a: float, b: float, n_iter: int = 200) -> float:
    """Compute the regularized incomplete beta function I_x(a, b).

    Uses a continued fraction expansion (Lentz's method) for numerical
    stability. This avoids a dependency on scipy for the F-test p-value.

    Args:
        x: Upper limit of integration (0 <= x <= 1).
        a: First shape parameter (> 0).
        b: Second shape parameter (> 0).
        n_iter: Maximum number of continued fraction iterations.

    Returns:
        Value of the regularized incomplete beta function.
    """
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0

    # Use the symmetry relation for better convergence
    if x > (a + 1.0) / (a + b + 2.0):
        return 1.0 - _regularized_incomplete_beta(1.0 - x, b, a, n_iter)

    # Log-beta function via lgamma
    import math

    ln_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)

    # Front factor
    ln_front = a * math.log(x) + b * math.log(1.0 - x) - ln_beta - math.log(a)

    # Continued fraction (Lentz's algorithm)
    tiny = 1e-30
    f = tiny
    c = tiny
    d = 0.0

    for m in range(1, n_iter + 1):
        # Even step: m2 = 2*m
        numerator = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
        d = 1.0 + numerator * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + numerator / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        f *= c * d

        # Odd step
        numerator = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + numerator * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + numerator / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = c * d
        f *= delta

        if abs(delta - 1.0) < 1e-10:
            break

    return math.exp(ln_front) * f


def top_hits_genotype_phenotype(
    assoc_results: List[Dict[str, Any]],
    genotype_matrix: List[List[int]],
    phenotypes: List[float],
    output_dir: Union[str, Path],
    n_top: int = 5,
) -> Dict[str, Any]:
    """Auto-generate genotype-phenotype boxplots for top associated variants.

    Selects the top N variants by p-value from association results and
    generates individual boxplots for each.

    Args:
        assoc_results: List of association result dictionaries. Each must
            contain at minimum a "p_value" key, and optionally "variant_id".
        genotype_matrix: List of genotype vectors, one per variant
            (same order as assoc_results). Each vector has one entry per sample.
        phenotypes: List of phenotype values, one per sample.
        output_dir: Directory to write output plot files.
        n_top: Number of top hits to plot (default 5).

    Returns:
        Dictionary with status, list of generated file paths, and summary.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for top hits plots")
        return {"status": "skipped", "output_path": None, "reason": "matplotlib not available"}

    if not HAS_NUMPY:
        logger.warning("numpy not available for top hits plots")
        return {"status": "skipped", "output_path": None, "reason": "numpy not available"}

    try:
        if not assoc_results:
            return {"status": "failed", "error": "No association results provided", "output_path": None}

        if len(assoc_results) != len(genotype_matrix):
            return {
                "status": "failed",
                "error": "assoc_results and genotype_matrix must have the same length",
                "output_path": None,
            }

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Sort by p-value, take top N
        indexed_results = [(i, r) for i, r in enumerate(assoc_results)]
        indexed_results.sort(key=lambda x: x[1].get("p_value", 1.0))
        top_results = indexed_results[: min(n_top, len(indexed_results))]

        generated_files: List[str] = []
        variant_summaries: List[Dict[str, Any]] = []

        for rank, (orig_idx, result) in enumerate(top_results, start=1):
            variant_id = result.get("variant_id", f"variant_{orig_idx}")
            p_val = result.get("p_value", 1.0)
            geno_vector = genotype_matrix[orig_idx]

            out_file = output_dir / f"top{rank}_{variant_id}.png"

            box_result = genotype_phenotype_boxplot(
                genotypes=geno_vector,
                phenotypes=phenotypes,
                variant_id=f"{variant_id} (p={p_val:.2e})",
                output_file=out_file,
                metadata=result,
            )

            if box_result["status"] == "success":
                generated_files.append(str(out_file))
                variant_summaries.append(
                    {
                        "rank": rank,
                        "variant_id": variant_id,
                        "p_value": p_val,
                        "output_path": str(out_file),
                    }
                )

        return {
            "status": "success",
            "n_plotted": len(generated_files),
            "generated_files": generated_files,
            "variant_summaries": variant_summaries,
            "output_path": str(output_dir),
        }

    except Exception as e:
        logger.error(f"Error creating top hits plots: {e}")
        return {"status": "failed", "error": str(e), "output_path": None}


def phenotype_pca_correlation(
    pcs: List[List[float]],
    phenotype_values: List[float],
    output_file: Optional[Union[str, Path]] = None,
    trait_name: str = "Trait",
) -> Dict[str, Any]:
    """Create a 1x2 panel showing phenotype vs PC1 and phenotype vs PC2.

    Each panel is a scatter plot colored by phenotype intensity with
    a Pearson correlation coefficient displayed in the panel title.

    Args:
        pcs: List of PC vectors per sample, where each inner list
            contains [PC1, PC2, ...] values.
        phenotype_values: List of phenotype values, one per sample.
        output_file: Optional path to save the plot image.
        trait_name: Human-readable name of the trait for labels.

    Returns:
        Dictionary with status, correlation values, and output path.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for phenotype PCA correlation plot")
        return {"status": "skipped", "output_path": None, "reason": "matplotlib not available"}

    if not HAS_NUMPY:
        logger.warning("numpy not available for phenotype PCA correlation plot")
        return {"status": "skipped", "output_path": None, "reason": "numpy not available"}

    try:
        pcs_arr = np.array(pcs, dtype=np.float64)
        pheno_arr = np.array(phenotype_values, dtype=np.float64)

        if len(pcs_arr) != len(pheno_arr):
            return {
                "status": "failed",
                "error": "pcs and phenotype_values must have the same length",
                "output_path": None,
            }

        if pcs_arr.ndim != 2 or pcs_arr.shape[1] < 2:
            return {"status": "failed", "error": "pcs must have at least 2 components per sample", "output_path": None}

        if len(pheno_arr) < 2:
            return {"status": "failed", "error": "At least 2 samples required", "output_path": None}

        pc1 = pcs_arr[:, 0]
        pc2 = pcs_arr[:, 1]

        # Pearson correlations
        r_pc1 = float(np.corrcoef(pc1, pheno_arr)[0, 1])
        r_pc2 = float(np.corrcoef(pc2, pheno_arr)[0, 1])

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Panel 1: Phenotype vs PC1
        sc1 = ax1.scatter(
            pc1, pheno_arr, c=pheno_arr, cmap="viridis", alpha=0.7, s=40, edgecolors="gray", linewidths=0.3
        )
        ax1.set_xlabel("PC1", fontsize=12)
        ax1.set_ylabel(trait_name, fontsize=12)
        ax1.set_title(f"{trait_name} vs PC1 (r={r_pc1:.3f})", fontsize=13)
        ax1.grid(True, alpha=0.3)
        plt.colorbar(sc1, ax=ax1, label=trait_name, shrink=0.8)

        # Panel 2: Phenotype vs PC2
        sc2 = ax2.scatter(
            pc2, pheno_arr, c=pheno_arr, cmap="viridis", alpha=0.7, s=40, edgecolors="gray", linewidths=0.3
        )
        ax2.set_xlabel("PC2", fontsize=12)
        ax2.set_ylabel(trait_name, fontsize=12)
        ax2.set_title(f"{trait_name} vs PC2 (r={r_pc2:.3f})", fontsize=13)
        ax2.grid(True, alpha=0.3)
        plt.colorbar(sc2, ax=ax2, label=trait_name, shrink=0.8)

        fig.suptitle(f"Phenotype-PCA Correlation: {trait_name}", fontsize=15, y=1.02)
        plt.tight_layout()

        if output_file:
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            logger.info(f"Saved phenotype PCA correlation plot to {output_file}")

        plt.close(fig)
        return {
            "status": "success",
            "n_samples": len(pheno_arr),
            "r_pc1": r_pc1,
            "r_pc2": r_pc2,
            "output_path": str(output_file) if output_file else None,
        }

    except Exception as e:
        logger.error(f"Error creating phenotype PCA correlation plot: {e}")
        return {"status": "failed", "error": str(e), "output_path": None}

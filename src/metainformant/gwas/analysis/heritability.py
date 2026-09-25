"""SNP heritability estimation and partitioning.

Estimates narrow-sense SNP heritability (h2_SNP) using REML variance component
estimation on genomic relationship matrices. Supports per-chromosome partitioning
and visualization of heritability contributions.

The core model:
    y = mu + g + e
    where Var(g) = sigma_g^2 * K, Var(e) = sigma_e^2 * I
    h2 = sigma_g^2 / (sigma_g^2 + sigma_e^2)

Reference: Yang et al. (2011) Nature Genetics 43:519-525 (GCTA-GREML).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging
from metainformant.gwas.heritability.estimation import greml_simple

logger = logging.get_logger(__name__)

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None


def estimate_heritability(
    kinship_matrix: Any,
    phenotypes: List[float],
    method: str = "reml",
) -> Dict[str, Any]:
    """Estimate SNP heritability (h2_SNP) using REML variance component estimation.

    Delegates to the shared REML core in
    :func:`metainformant.gwas.heritability.estimation.greml_simple`, which
    eigendecomposes the kinship matrix, rotates the phenotypes into the
    eigenspace, and maximizes the restricted log-likelihood by grid search
    followed by golden-section refinement. Keeping a single REML core avoids
    divergent duplicate implementations.

    The core model:
        y = mu + g + e
        where Var(g) = sigma_g^2 * K, Var(e) = sigma_e^2 * I
        h2 = sigma_g^2 / (sigma_g^2 + sigma_e^2)

    Reference: Yang et al. (2011) Nature Genetics 43:519-525 (GCTA-GREML).

    Args:
        kinship_matrix: Kinship/GRM matrix (n x n), numpy array or list of lists.
        phenotypes: Phenotype values for each sample.
        method: Estimation method label. Only REML (GREML) is implemented;
            the label is reported back in the result dictionary.

    Returns:
        Dictionary with status, h2 estimate, standard error, variance components,
        log-likelihood, sample size, and method used.
    """
    result = greml_simple(kinship_matrix, phenotypes)
    if result.get("status") == "success":
        result["method"] = method
    return result


def partition_heritability_by_chromosome(
    kinship_matrices: Dict[int, Any],
    phenotypes: List[float],
) -> Dict[str, Any]:
    """Partition SNP heritability by chromosome.

    Estimates per-chromosome h2 by fitting each chromosome's kinship matrix
    independently (one-at-a-time approach). The total h2 is the sum of
    per-chromosome estimates.

    Args:
        kinship_matrices: Maps chromosome number to kinship matrix for
            variants on that chromosome.
        phenotypes: Phenotype values for each sample.

    Returns:
        Dictionary with per-chromosome h2 estimates, total h2, and
        number of chromosomes analyzed.
    """
    if not HAS_NUMPY:
        return {
            "status": "error",
            "message": "numpy is required for heritability partitioning",
        }

    n = len(phenotypes)
    if n < 3:
        return {
            "status": "error",
            "message": f"Need at least 3 samples, got {n}",
        }

    if not kinship_matrices:
        return {"status": "error", "message": "No kinship matrices provided"}

    logger.info(f"Partitioning heritability across {len(kinship_matrices)} chromosomes")

    per_chromosome: Dict[str, Dict[str, float]] = {}
    total_h2 = 0.0

    for chrom, K_chr in sorted(kinship_matrices.items()):
        result = estimate_heritability(K_chr, phenotypes)

        if result["status"] == "success":
            chr_h2 = result["h2"]
            chr_se = result["h2_se"]
            per_chromosome[str(chrom)] = {"h2": chr_h2, "h2_se": chr_se}
            total_h2 += chr_h2
        else:
            logger.warning(
                f"Chromosome {chrom} estimation failed: {result.get('message', 'unknown')}"
            )
            per_chromosome[str(chrom)] = {"h2": 0.0, "h2_se": 0.0}

    # Cap total h2 at 1.0 (sum of independent estimates can exceed 1)
    total_h2 = min(total_h2, 1.0)

    return {
        "status": "success",
        "per_chromosome": per_chromosome,
        "total_h2": float(total_h2),
        "n_chromosomes": len(kinship_matrices),
    }


def heritability_bar_chart(
    h2_data: Dict[str, Any],
    output_file: Optional[Union[str, Path]] = None,
    title: str = "SNP Heritability by Chromosome",
) -> Dict[str, Any]:
    """Create a bar chart of per-chromosome heritability estimates.

    Plots h2 per chromosome with error bars (h2_se), draws a horizontal
    line at the total h2, and colors bars by relative contribution.

    Args:
        h2_data: Output from partition_heritability_by_chromosome containing
            per_chromosome and total_h2 keys.
        output_file: Path to save the figure. If None, the figure is not saved.
        title: Chart title.

    Returns:
        Dictionary with status and output path (if saved).
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return {
            "status": "skipped",
            "output_path": None,
            "message": "matplotlib not available",
        }

    per_chromosome = h2_data.get("per_chromosome", {})
    total_h2 = h2_data.get("total_h2", 0.0)

    if not per_chromosome:
        return {
            "status": "failed",
            "output_path": None,
            "message": "No per-chromosome data",
        }

    # Sort chromosomes numerically where possible
    def _chr_sort_key(c: str) -> Tuple[int, str]:
        try:
            return (0, str(int(c)).zfill(5))
        except ValueError:
            return (1, c)

    sorted_chroms = sorted(per_chromosome.keys(), key=_chr_sort_key)
    h2_values = [per_chromosome[c]["h2"] for c in sorted_chroms]
    h2_se_values = [per_chromosome[c]["h2_se"] for c in sorted_chroms]

    try:
        fig, ax = plt.subplots(figsize=(max(8, len(sorted_chroms) * 0.6), 5))

        # Color bars by relative contribution
        max_h2 = max(h2_values) if max(h2_values) > 0 else 1.0
        colors = plt.cm.YlOrRd([v / max_h2 for v in h2_values])

        x_positions = range(len(sorted_chroms))
        ax.bar(
            x_positions,
            h2_values,
            yerr=h2_se_values,
            capsize=3,
            color=colors,
            edgecolor="gray",
            linewidth=0.5,
        )

        # Horizontal line at total h2
        ax.axhline(
            y=total_h2,
            color="steelblue",
            linestyle="--",
            linewidth=1.5,
            label=f"Total h2 = {total_h2:.3f}",
        )

        ax.set_xlabel("Chromosome")
        ax.set_ylabel("Heritability (h2)")
        ax.set_title(title)
        ax.set_xticks(list(x_positions))
        ax.set_xticklabels(sorted_chroms, rotation=45 if len(sorted_chroms) > 10 else 0)
        ax.legend(loc="upper right")
        ax.set_ylim(bottom=0)

        plt.tight_layout()

        output_path_str: Optional[str] = None
        if output_file is not None:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=150, bbox_inches="tight")
            output_path_str = str(output_path)
            logger.info(f"Heritability bar chart saved to {output_path}")

        plt.close(fig)

        return {"status": "success", "output_path": output_path_str}

    except Exception as e:
        logger.warning(f"Heritability bar chart failed: {e}")
        return {"status": "failed", "output_path": None, "message": str(e)}

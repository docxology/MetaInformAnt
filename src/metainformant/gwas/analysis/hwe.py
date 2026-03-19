"""Hardy-Weinberg Equilibrium (HWE) analysis for GWAS QC.

Provides per-variant HWE chi-square tests, genome-wide distribution
reporting, and MAF-stratified flagging. Detects genotyping errors,
selection signals, and batch effects.

Key references:
  - Wigginton et al. (2005) Am J Hum Genet 76:887-893 (exact HWE test)
  - Laurie et al. (2010) Genet Epidemiol — stratified HWE QC
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple, Union
from pathlib import Path

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def hwe_chi2_test(
    n_hom_ref: int,
    n_het: int,
    n_hom_alt: int,
) -> Tuple[float, float]:
    """Chi-square test for Hardy-Weinberg equilibrium.

    Tests H₀: genotype frequencies follow HWE (p²AA + 2p(1-p)Aa + (1-p)²aa).

    Args:
        n_hom_ref: Count of homozygous reference genotypes (0/0)
        n_het: Count of heterozygous genotypes (0/1)
        n_hom_alt: Count of homozygous alt genotypes (1/1)

    Returns:
        Tuple of (chi2_statistic, p_value). Returns (0.0, 1.0) for monomorphic.
    """
    n = n_hom_ref + n_het + n_hom_alt
    if n == 0:
        return 0.0, 1.0

    # Allele frequencies
    p_alt = (2 * n_hom_alt + n_het) / (2 * n)
    p_ref = 1.0 - p_alt

    # Expected counts under HWE
    e_hom_ref = n * p_ref ** 2
    e_het = n * 2 * p_ref * p_alt
    e_hom_alt = n * p_alt ** 2

    # Check for monomorphic (all expected in one class)
    if e_hom_ref < 0.01 or e_het < 0.01 or e_hom_alt < 0.01:
        return 0.0, 1.0

    # Pearson chi-square with 1 df (after imposing N and allele freq constraints)
    chi2 = (
        (n_hom_ref - e_hom_ref) ** 2 / e_hom_ref
        + (n_het - e_het) ** 2 / e_het
        + (n_hom_alt - e_hom_alt) ** 2 / e_hom_alt
    )

    # p-value from chi-square CDF with 1 df
    p_value = _chi2_sf(chi2, df=1)
    return chi2, p_value


def hwe_flag_variants(
    vcf_data: Dict[str, Any],
    threshold: float = 1e-6,
) -> List[Dict[str, Any]]:
    """Compute HWE statistics for all variants in a parsed VCF.

    Args:
        vcf_data: Parsed VCF dict (from parse_vcf_full); needs 'variants'
            and 'genotypes' keys.
        threshold: P-value threshold below which variants are flagged
            (default: 1e-6, commonly 1e-4 to 1e-6 in GWAS QC).

    Returns:
        List of dicts:
            chrom, pos, snp_id, maf, n_hom_ref, n_het, n_hom_alt,
            chi2, p_hwe, flagged (bool), obs_het, exp_het
    """
    variants = vcf_data.get("variants", [])
    genotypes_raw = vcf_data.get("genotypes", {})
    samples = vcf_data.get("samples", [])
    n_samples = len(samples)

    # genotypes may be:
    #   list[list] — sample-major (n_samples x n_variants), from parse_vcf_full
    #   dict[sample_id, list] — keyed by sample ID
    _DOSAGE_TO_GT = {0: "0/0", 1: "0/1", 2: "1/1", "0": "0/0", "1": "0/1", "2": "1/1"}

    if isinstance(genotypes_raw, list):
        def _gt(sample_idx: int, var_idx: int) -> str:
            row = genotypes_raw[sample_idx] if sample_idx < len(genotypes_raw) else []
            raw = row[var_idx] if var_idx < len(row) else None
            return _DOSAGE_TO_GT.get(raw, str(raw)) if raw is not None else "./."
        sample_indices = list(range(n_samples))
    else:
        def _gt(sample_idx: int, var_idx: int) -> str:
            sample_id = samples[sample_idx] if sample_idx < len(samples) else ""
            gt_list = genotypes_raw.get(sample_id, [])
            raw = gt_list[var_idx] if var_idx < len(gt_list) else None
            return _DOSAGE_TO_GT.get(raw, str(raw)) if raw is not None else "./."
        sample_indices = list(range(n_samples))

    results = []
    for i, var in enumerate(variants):
        snp_id = var.get("id", f"var_{i}")
        chrom = var.get("chrom", "")
        pos = var.get("pos", 0)

        # Count genotypes using format-agnostic accessor
        n_hom_ref = n_het = n_hom_alt = n_miss = 0
        for s_idx in sample_indices:
            gt = _gt(s_idx, i)
            if gt in ("0/0", "0|0"):
                n_hom_ref += 1
            elif gt in ("0/1", "1/0", "0|1", "1|0"):
                n_het += 1
            elif gt in ("1/1", "1|1"):
                n_hom_alt += 1
            else:
                n_miss += 1

        n_called = n_hom_ref + n_het + n_hom_alt
        if n_called == 0:
            continue

        # Allele freq
        alt_alleles = 2 * n_hom_alt + n_het
        total_alleles = 2 * n_called
        maf = min(alt_alleles / total_alleles, 1 - alt_alleles / total_alleles) if total_alleles > 0 else 0.0

        chi2, p_val = hwe_chi2_test(n_hom_ref, n_het, n_hom_alt)

        p_alt = (2 * n_hom_alt + n_het) / (2 * n_called)
        p_ref = 1.0 - p_alt
        exp_het = 2 * p_ref * p_alt
        obs_het = n_het / n_called if n_called > 0 else 0.0

        results.append({
            "chrom": chrom,
            "pos": pos,
            "snp_id": snp_id,
            "maf": round(maf, 4),
            "n_hom_ref": n_hom_ref,
            "n_het": n_het,
            "n_hom_alt": n_hom_alt,
            "n_missing": n_miss,
            "chi2_hwe": round(chi2, 4),
            "p_hwe": p_val,
            "flagged": p_val < threshold,
            "obs_het": round(obs_het, 4),
            "exp_het": round(exp_het, 4),
        })

    logger.info(f"HWE test: {len(results)} variants tested, "
                f"{sum(1 for r in results if r['flagged'])} flagged at p<{threshold:.0e}")
    return results


def hwe_distribution_plot(
    hwe_results: List[Dict[str, Any]],
    output_path: Optional[Union[str, Path]] = None,
    threshold: float = 1e-6,
) -> Any:
    """Plot genome-wide HWE p-value distribution and obs vs expected heterozygosity.

    Two-panel figure:
      Left:  −log10(p_HWE) distribution histogram with threshold line
      Right: Observed vs expected heterozygosity scatter, colored by MAF

    Args:
        hwe_results: Output of hwe_flag_variants()
        output_path: Where to save the figure
        threshold: HWE flagging threshold (draws vertical dashed line)

    Returns:
        matplotlib Figure, or None if matplotlib unavailable.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for HWE plot")
        return None

    if not hwe_results:
        logger.warning("No HWE results to plot")
        return None

    p_vals = [r["p_hwe"] for r in hwe_results if r["p_hwe"] > 0]
    obs_hets = [r["obs_het"] for r in hwe_results]
    exp_hets = [r["exp_het"] for r in hwe_results]
    mafs = [r["maf"] for r in hwe_results]
    n_flagged = sum(1 for r in hwe_results if r["flagged"])

    log_p = [-math.log10(max(p, 1e-300)) for p in p_vals]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f"Hardy-Weinberg Equilibrium Diagnostic\n"
        f"N={len(hwe_results):,} variants, {n_flagged} flagged (p<{threshold:.0e})",
        fontsize=12, fontweight="bold",
    )

    # Left panel: HWE -log10(p) histogram
    ax1.hist(log_p, bins=min(60, max(10, len(log_p) // 10)),
             color="#4C72B0", edgecolor="white", linewidth=0.5, alpha=0.85)
    ax1.axvline(-math.log10(threshold), color="#C44E52", linewidth=2,
                linestyle="--", label=f"Threshold (p={threshold:.0e})")
    ax1.set_xlabel("−log₁₀(HWE p-value)", fontsize=11)
    ax1.set_ylabel("Number of variants", fontsize=11)
    ax1.set_title("HWE p-value distribution", fontsize=10)
    ax1.legend(fontsize=9)

    # Add annotation
    ax1.text(0.97, 0.95, f"{n_flagged} flagged\n({100*n_flagged/len(hwe_results):.1f}%)",
             transform=ax1.transAxes, ha="right", va="top",
             fontsize=9, color="#C44E52",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # Right panel: obs vs expected het scatter, colored by MAF
    if HAS_NUMPY:
        mafs_arr = [r["maf"] for r in hwe_results]
        scatter = ax2.scatter(exp_hets, obs_hets, c=mafs_arr, cmap="viridis",
                              s=4, alpha=0.4, linewidths=0)
        plt.colorbar(scatter, ax=ax2, label="MAF", shrink=0.8)
    else:
        ax2.scatter(exp_hets, obs_hets, s=4, alpha=0.3, color="#4C72B0")

    lim = max(max(obs_hets, default=1), max(exp_hets, default=1)) * 1.05
    ax2.plot([0, lim], [0, lim], "r--", linewidth=1.5, label="obs = exp (HWE)")
    ax2.set_xlabel("Expected heterozygosity", fontsize=11)
    ax2.set_ylabel("Observed heterozygosity", fontsize=11)
    ax2.set_title("Obs vs Expected Heterozygosity", fontsize=10)
    ax2.legend(fontsize=9)
    ax2.set_xlim(0, lim)
    ax2.set_ylim(0, lim)

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"HWE distribution plot saved to {output_path}")
    return fig


# ── Pure-Python chi-square survival function ──────────────────────────────────

def _chi2_sf(x: float, df: int = 1) -> float:
    """Survival function of chi-square distribution P(X > x).

    Uses regularized incomplete gamma function for df=1 (simplest case).
    For df=1: P(X > x) = erfc(sqrt(x/2)).
    """
    if x <= 0:
        return 1.0
    if df == 1:
        return math.erfc(math.sqrt(x / 2.0))
    # For df=2: P(X > x) = exp(-x/2)
    if df == 2:
        return math.exp(-x / 2.0)
    # General case via regularized incomplete gamma (series expansion)
    return _regularized_gamma_q(df / 2.0, x / 2.0)


def _regularized_gamma_q(a: float, x: float, max_iter: int = 200, tol: float = 1e-12) -> float:
    """Q(a, x) = 1 - P(a, x): upper regularized incomplete gamma function.

    Uses continued fraction representation (Lentz method).
    """
    if x < 0:
        return 1.0
    if x == 0:
        return 1.0

    # For x < a + 1: use series expansion for P(a, x), then Q = 1 - P
    if x < a + 1:
        # Series expansion for P(a, x)
        term = 1.0 / a
        sum_p = term
        for n in range(1, max_iter):
            term *= x / (a + n)
            sum_p += term
            if abs(term) < tol * abs(sum_p):
                break
        # P(a, x) = exp(-x) * x^a * sum_p / Gamma(a)
        log_p = -x + a * math.log(x) - _log_gamma(a) + math.log(sum_p)
        p = math.exp(log_p) if log_p < 0 else 1.0
        return max(0.0, min(1.0, 1.0 - p))

    # Continued fraction for Q(a, x)
    f = x + 1 - a
    c = 1.0 / 1e-300
    d = 1.0 / f
    h = d
    for i in range(1, max_iter):
        an = -i * (i - a)
        f += 2
        d = an * d + f
        if abs(d) < 1e-300:
            d = 1e-300
        c = f + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < tol:
            break
    log_q = -x + a * math.log(x) - _log_gamma(a) + math.log(h)
    return max(0.0, min(1.0, math.exp(log_q)))


def _log_gamma(x: float) -> float:
    """Natural log of Gamma function via Lanczos approximation."""
    g = 7
    c = [
        0.99999999999980993, 676.5203681218851, -1259.1392167224028,
        771.32342877765313, -176.61502916214059, 12.507343278686905,
        -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7,
    ]
    if x < 0.5:
        return math.log(math.pi / math.sin(math.pi * x)) - _log_gamma(1 - x)
    x -= 1
    a = c[0]
    t = x + g + 0.5
    for i in range(1, g + 2):
        a += c[i] / (x + i)
    return 0.5 * math.log(2 * math.pi) + (x + 0.5) * math.log(t) - t + math.log(a)

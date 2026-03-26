"""GWAS summary statistics output utilities.

This module provides functions for writing GWAS results in standard formats,
including full summary statistics TSV, significant hits tables, and JSON summaries.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def write_summary_statistics(
    results: List[Dict[str, Any]],
    variant_info: List[Dict[str, Any]],
    output_path: Union[str, Path],
) -> Path:
    """Write full GWAS summary statistics to TSV file.

    Output format follows standard GWAS summary statistics convention:
    CHR  POS  SNP  REF  ALT  BETA  SE  P  N  MAF

    Args:
        results: List of association test result dictionaries
        variant_info: List of variant info dictionaries with chrom, pos, id, ref, alt
        output_path: Path to output TSV file

    Returns:
        Path to written file

    Raises:
        ValueError: If results and variant_info have different lengths
    """
    if len(results) != len(variant_info):
        raise ValueError(
            f"Results ({len(results)}) and variant_info ({len(variant_info)}) "
            "must have the same length"
        )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(
        f"Writing summary statistics for {len(results)} variants to {output_path}"
    )

    header = "CHR\tPOS\tSNP\tREF\tALT\tBETA\tSE\tP\tN\tMAF\n"

    with open(output_path, "w") as f:
        f.write(header)
        for result, vinfo in zip(results, variant_info):
            chrom = vinfo.get("chrom", ".")
            pos = vinfo.get("pos", 0)
            snp_id = vinfo.get("id", ".")
            ref = vinfo.get("ref", ".")
            alt = vinfo.get("alt", ".")
            if isinstance(alt, list):
                alt = ",".join(alt)

            beta = result.get("beta", 0.0)
            se = result.get("se", 0.0)
            p_value = result.get("p_value", 1.0)
            n_samples = result.get("n_samples", 0)
            maf = result.get("maf", 0.0)

            f.write(
                f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\t"
                f"{beta:.6g}\t{se:.6g}\t{p_value:.6e}\t{n_samples}\t{maf:.4f}\n"
            )

    logger.info(f"Summary statistics written to {output_path}")
    return output_path


def write_significant_hits(
    results: List[Dict[str, Any]],
    variant_info: List[Dict[str, Any]],
    output_path: Union[str, Path],
    threshold: float = 5e-8,
) -> Path:
    """Write genome-wide significant hits to TSV file.

    Args:
        results: List of association test result dictionaries
        variant_info: List of variant info dictionaries
        output_path: Path to output TSV file
        threshold: Significance threshold (default: 5e-8 genome-wide)

    Returns:
        Path to written file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    significant = [
        (result, vinfo)
        for result, vinfo in zip(results, variant_info)
        if result.get("p_value", 1.0) < threshold
    ]

    # Sort by p-value
    significant.sort(key=lambda x: x[0].get("p_value", 1.0))

    logger.info(
        f"Writing {len(significant)} significant hits (threshold={threshold:.2e}) to {output_path}"
    )

    with open(output_path, "w") as f:
        f.write("CHR\tPOS\tSNP\tREF\tALT\tBETA\tSE\tP\tN\tMAF\n")
        for result, vinfo in significant:
            chrom = vinfo.get("chrom", ".")
            pos = vinfo.get("pos", 0)
            snp_id = vinfo.get("id", ".")
            ref = vinfo.get("ref", ".")
            alt = vinfo.get("alt", ".")
            if isinstance(alt, list):
                alt = ",".join(alt)

            beta = result.get("beta", 0.0)
            se = result.get("se", 0.0)
            p_value = result.get("p_value", 1.0)
            n_samples = result.get("n_samples", 0)
            maf = result.get("maf", 0.0)

            f.write(
                f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\t"
                f"{beta:.6g}\t{se:.6g}\t{p_value:.6e}\t{n_samples}\t{maf:.4f}\n"
            )

    return output_path


def create_results_summary(
    results: List[Dict[str, Any]],
    output_path: Union[str, Path],
    threshold: float = 5e-8,
) -> Dict[str, Any]:
    """Create JSON summary of GWAS results.

    Args:
        results: List of association test result dictionaries
        output_path: Path to output JSON file
        threshold: Genome-wide significance threshold

    Returns:
        Summary dictionary with lambda_gc, n_significant, top_hits
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    p_values = [r.get("p_value", 1.0) for r in results]
    n_tests = len(p_values)

    # Compute genomic inflation factor (lambda_gc)
    lambda_gc = _compute_lambda_gc(p_values)

    # Count significant hits
    n_significant = sum(1 for p in p_values if p < threshold)
    n_suggestive = sum(1 for p in p_values if p < 1e-5)

    # Top hits (sorted by p-value)
    indexed_results = [(i, r) for i, r in enumerate(results)]
    indexed_results.sort(key=lambda x: x[1].get("p_value", 1.0))
    top_hits = []
    for idx, result in indexed_results[:20]:
        hit = {
            "variant_index": idx,
            "p_value": result.get("p_value", 1.0),
            "beta": result.get("beta", 0.0),
            "se": result.get("se", 0.0),
        }
        if "variant_id" in result:
            hit["variant_id"] = result["variant_id"]
        if "chrom" in result:
            hit["chrom"] = result["chrom"]
        if "pos" in result:
            hit["pos"] = result["pos"]
        top_hits.append(hit)

    summary = {
        "n_variants_tested": n_tests,
        "lambda_gc": lambda_gc,
        "n_significant_genome_wide": n_significant,
        "n_suggestive": n_suggestive,
        "significance_threshold": threshold,
        "top_hits": top_hits,
    }

    with open(output_path, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(
        f"Results summary: {n_tests} variants, lambda_gc={lambda_gc:.3f}, "
        f"{n_significant} genome-wide significant"
    )
    return summary


# ── Public analysis functions (extracted from orchestrator inline logic) ──


def normal_cdf(x: float) -> float:
    """Standard normal CDF via Abramowitz & Stegun (7.1.26) approximation.

    Pure-Python — no scipy dependency.
    """
    import math as _math

    sign = 1 if x >= 0 else -1
    x = abs(x) / _math.sqrt(2)
    t = 1.0 / (1.0 + 0.3275911 * x)
    y = 1.0 - (
        (((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t
        + 0.254829592
    ) * t * _math.exp(-x * x)
    return 0.5 * (1.0 + sign * y)


def calculate_hwe_pvalue(obs_hets: int, obs_hom1: int, obs_hom2: int) -> float:
    """Calculate Hardy-Weinberg equilibrium p-value using chi-squared approximation.

    Args:
        obs_hets: Number of heterozygotes (0/1)
        obs_hom1: Number of reference homozygotes (0/0)
        obs_hom2: Number of alternate homozygotes (1/1)

    Returns:
        p-value for HWE deviation test
    """
    n_samples = obs_hets + obs_hom1 + obs_hom2
    if n_samples == 0:
        return 1.0

    n_alleles = 2 * n_samples
    n_ref = 2 * obs_hom1 + obs_hets
    n_alt = 2 * obs_hom2 + obs_hets

    p = n_ref / n_alleles
    q = n_alt / n_alleles

    exp_hom1 = (p * p) * n_samples
    exp_hets = (2 * p * q) * n_samples
    exp_hom2 = (q * q) * n_samples

    # Chi-square statistic
    chi2 = 0.0
    for obs, exp in [(obs_hom1, exp_hom1), (obs_hets, exp_hets), (obs_hom2, exp_hom2)]:
        if exp > 0:
            chi2 += ((obs - exp) ** 2) / exp

    # p-value from chi2 with 1 dof
    try:
        from scipy import stats

        return float(stats.chi2.sf(chi2, 1))
    except ImportError:
        # Fallback to standard normal squared approximation
        return float(2 * (1 - normal_cdf(math.sqrt(chi2))))


def calculate_titv_ratio(variant_info: List[Dict[str, Any]]) -> float:
    """Calculate the Transition/Transversion (Ti/Tv) ratio for a set of variants.

    Args:
        variant_info: List of dictionary entries with 'ref' and 'alt' alleles

    Returns:
        Ti/Tv ratio float
    """
    transitions = 0
    transversions = 0

    ti_pairs = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    for vinfo in variant_info:
        ref = vinfo.get("ref", "").upper()
        alt = vinfo.get("alt", "")
        if isinstance(alt, list):
            alt = alt[0].upper() if alt else ""
        elif isinstance(alt, str):
            alt = alt.upper()
        else:
            continue

        if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
            if (ref, alt) in ti_pairs:
                transitions += 1
            else:
                transversions += 1

    if transversions == 0:
        return float(transitions) if transitions > 0 else 0.0
    return transitions / transversions


def calculate_ld_decay(
    genotype_matrix: List[List[int]],
    variants_info: List[Dict[str, Any]],
    max_distance: int = 1000000,
    bins: int = 50,
) -> Dict[str, Any]:
    """Calculate Linkage Disequilibrium (r^2) decay across genomic distances.

    Args:
        genotype_matrix: List of genotypes (variants x samples)
        variants_info: List of variant metadata
        max_distance: Maximum base-pair distance to compute LD
        bins: Number of distance bins for decay curve

    Returns:
        Dictionary with 'distances' (bin centers) and 'r2_means' (mean r^2 in each bin)
    """
    if not HAS_NUMPY:  # noqa: F821
        raise ImportError("numpy is required for LD decay calculation")

    logger.info(f"Computing LD decay up to {max_distance}bp across {bins} bins...")

    # Pre-filter for valid variants (no multiallelics etc if needed)
    G = np.array(genotype_matrix, dtype=float)  # noqa: F821

    # Calculate minor allele frequencies
    # If G has -1 for missing, we should mask them. For simplicity assume complete here or imputed
    np.mean(G, axis=1)  # noqa: F821

    # Basic LD logic (this scales at O(N^2), so we sample random pairs if there are many)
    n_vars = G.shape[0]
    sample_pairs = min(10000, n_vars * (n_vars - 1) // 2)

    r2_vals = []
    distances = []

    import random

    pairs = set()
    attempts = 0
    while len(pairs) < sample_pairs and attempts < sample_pairs * 5:
        attempts += 1
        i = random.randint(0, n_vars - 1)
        j = random.randint(0, n_vars - 1)
        if i == j:
            continue

        if (j, i) in pairs or (i, j) in pairs:
            continue

        chrom_i = str(variants_info[i].get("chrom", "1"))
        chrom_j = str(variants_info[j].get("chrom", "2"))

        if chrom_i != chrom_j:
            continue

        pos_i = int(variants_info[i].get("pos", 0))
        pos_j = int(variants_info[j].get("pos", 0))
        dist = abs(pos_i - pos_j)

        if dist > max_distance:
            continue

        pairs.add((min(i, j), max(i, j)))

        # Calculate r^2
        g_i = G[i]
        g_j = G[j]

        cov = np.cov(g_i, g_j)[0, 1]  # noqa: F821
        var_i = np.var(g_i)  # noqa: F821
        var_j = np.var(g_j)  # noqa: F821

        if var_i > 0 and var_j > 0:
            r2 = (cov**2) / (var_i * var_j)
            r2_vals.append(r2)
            distances.append(dist)

    # Binning
    bin_edges = np.linspace(0, max_distance, bins + 1)  # noqa: F821
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    r2_means = np.zeros(bins)  # noqa: F821

    if distances:
        dist_array = np.array(distances)  # noqa: F821
        r2_array = np.array(r2_vals)  # noqa: F821
        for b in range(bins):
            mask = (dist_array >= bin_edges[b]) & (dist_array < bin_edges[b + 1])
            if np.any(mask):  # noqa: F821
                r2_means[b] = np.mean(r2_array[mask])  # noqa: F821

    return {
        "distances": bin_centers.tolist(),
        "r2_means": r2_means.tolist(),
        "n_pairs": len(distances),
    }


def compute_comprehensive_summary(
    association_results: List[Dict[str, Any]],
    *,
    significance_threshold: float = 5e-8,
) -> Dict[str, Any]:
    """Compute a comprehensive GWAS summary statistics dictionary.

    Produces MAF, effect-size, SE, z-score distributions; calibration
    quantiles; significance counts; and top-10 hits.  Extracted from
    script ``08_post_gwas.py`` step 8f.

    Args:
        association_results: List of result dicts, each containing at least
            ``p_value``, ``beta``, ``se``, ``MAF``, ``chrom``, ``pos``, ``snp``.
        significance_threshold: Genome-wide threshold.

    Returns:
        Nested summary dictionary.
    """
    import statistics as _stats

    p_values = [r["p_value"] for r in association_results]
    betas = [r["beta"] for r in association_results]
    ses = [r["se"] for r in association_results]
    mafs = [r["MAF"] for r in association_results]
    abs_betas = [abs(b) for b in betas]
    z_scores = [b / s if s > 0 else 0.0 for b, s in zip(betas, ses)]

    # λ_GC
    lambda_gc = _compute_lambda_gc(p_values)

    # Calibration quantiles
    q_levels = [0.01, 0.05, 0.25, 0.50, 0.75]
    sorted_p = sorted(p_values)
    n_p = len(sorted_p)
    calibration_quantiles: Dict[str, Dict[str, Any]] = {}
    for q in q_levels:
        idx = min(int(q * n_p), n_p - 1)
        obs = sorted_p[idx]
        exp = q
        calibration_quantiles[f"q{int(q * 100):02d}"] = {
            "observed": round(obs, 6),
            "expected": round(exp, 6),
            "ratio_obs_exp": round(obs / exp, 4) if exp > 0 else float("nan"),
        }

    # Significance counts
    sig_gw = sum(1 for p in p_values if p < significance_threshold)
    sig_suggestive = sum(1 for p in p_values if p < 1e-5)
    sig_nominal = sum(1 for p in p_values if p < 0.05)

    # Per-chromosome variant counts
    chrom_counts: Dict[str, int] = {}
    for r in association_results:
        c = r["chrom"]
        chrom_counts[c] = chrom_counts.get(c, 0) + 1

    # Top 10 hits
    sorted_by_p = sorted(association_results, key=lambda r: r["p_value"])
    top_hits = [
        {
            "snp": r["snp"],
            "chrom": r["chrom"],
            "pos": r["pos"],
            "beta": round(r["beta"], 4),
            "se": round(r["se"], 4),
            "p_value": r["p_value"],
            "MAF": r["MAF"],
            "z_score": round(r["beta"] / r["se"], 4) if r["se"] > 0 else 0.0,
        }
        for r in sorted_by_p[:10]
    ]

    return {
        "n_variants": len(association_results),
        "n_chromosomes": len(chrom_counts),
        "maf_distribution": {
            "mean": round(_stats.mean(mafs), 6),
            "median": round(_stats.median(mafs), 6),
            "min": round(min(mafs), 6),
            "max": round(max(mafs), 6),
            "sd": round(_stats.stdev(mafs), 6) if len(mafs) > 1 else 0.0,
        },
        "effect_size_distribution": {
            "mean_abs_beta": round(_stats.mean(abs_betas), 6),
            "median_abs_beta": round(_stats.median(abs_betas), 6),
            "sd_beta": round(_stats.stdev(betas), 6) if len(betas) > 1 else 0.0,
            "range_min": round(min(betas), 6),
            "range_max": round(max(betas), 6),
        },
        "se_distribution": {
            "mean_se": round(_stats.mean(ses), 6),
            "median_se": round(_stats.median(ses), 6),
            "min_se": round(min(ses), 6),
        },
        "z_score_distribution": {
            "mean_abs_z": round(_stats.mean([abs(z) for z in z_scores]), 6),
            "max_abs_z": round(max([abs(z) for z in z_scores]), 6),
        },
        "p_value_calibration": {
            "lambda_gc": round(lambda_gc, 6),
            "median_p": round(_stats.median(p_values), 6),
            "min_p": min(p_values),
            "calibration_quantiles": calibration_quantiles,
        },
        "significance_counts": {
            "genome_wide_5e-8": sig_gw,
            "suggestive_1e-5": sig_suggestive,
            "nominal_0.05": sig_nominal,
        },
        "variants_per_chromosome": chrom_counts,
        "top_10_hits": top_hits,
    }


def sign_test(
    association_results: List[Dict[str, Any]],
    *,
    significance_threshold: float = 5e-8,
) -> Dict[str, Any]:
    """Binomial sign test on effect direction (positive vs negative betas).

    Uses the normal approximation for n > 30.  Extracted from script
    ``08_post_gwas.py`` step 8h.

    Args:
        association_results: List of dicts with ``beta`` and ``p_value`` keys.
        significance_threshold: GW threshold for sub-analysis.

    Returns:
        Dict with n_positive, n_negative, z_statistic, p_value_two_sided,
        interpretation, and GW-hit breakdowns.
    """
    betas = [r["beta"] for r in association_results]
    n_pos = sum(1 for b in betas if b > 0)
    n_neg = sum(1 for b in betas if b < 0)
    n_total = n_pos + n_neg

    if n_total == 0:
        return {"n_positive_beta": 0, "n_negative_beta": 0, "n_total": 0}

    p_hat = n_pos / n_total
    import math

    se_hat = math.sqrt(0.25 / n_total)
    z_sign = (p_hat - 0.5) / se_hat
    p_sign = 2.0 * min(normal_cdf(z_sign), 1.0 - normal_cdf(z_sign))

    gw_hits = [r for r in association_results if r["p_value"] < significance_threshold]
    n_gw_pos = sum(1 for r in gw_hits if r["beta"] > 0)

    interpretation = (
        "Significant directional enrichment (positive effects dominant)"
        if p_sign < 0.05 and n_pos > n_neg
        else "Significant directional enrichment (negative effects dominant)"
        if p_sign < 0.05 and n_neg > n_pos
        else "No significant directional enrichment (balanced +/- effects)"
    )

    return {
        "n_positive_beta": n_pos,
        "n_negative_beta": n_neg,
        "n_total": n_total,
        "p_hat_positive": round(p_hat, 4),
        "z_statistic": round(z_sign, 4),
        "p_value_two_sided": round(p_sign, 6),
        "interpretation": interpretation,
        "gw_significant_hits": len(gw_hits),
        "gw_positive_beta": n_gw_pos,
    }


def per_chromosome_summary(
    association_results: List[Dict[str, Any]],
    *,
    significance_threshold: float = 5e-8,
) -> List[Dict[str, Any]]:
    """Group-by-chromosome variant summary table.

    Extracted from script ``08_post_gwas.py`` step 8j.

    Args:
        association_results: Result dicts with ``chrom``, ``p_value``,
            ``beta``, ``MAF``.
        significance_threshold: GW threshold.

    Returns:
        List of per-chromosome summary dicts sorted by chromosome.
    """
    chrom_data: Dict[str, List[Dict[str, Any]]] = {}
    for r in association_results:
        chrom_data.setdefault(r["chrom"], []).append(r)

    rows = []
    for chrom in sorted(chrom_data.keys()):
        variants = chrom_data[chrom]
        ps = [v["p_value"] for v in variants]
        bs = [v["beta"] for v in variants]
        ms = [v["MAF"] for v in variants]
        n_gw = sum(1 for p in ps if p < significance_threshold)
        n_sug = sum(1 for p in ps if p < 1e-5)
        rows.append(
            {
                "CHR": chrom,
                "N_VARIANTS": len(variants),
                "MIN_P": min(ps),
                "MEAN_ABS_BETA": round(sum(abs(b) for b in bs) / len(bs), 4),
                "MEAN_MAF": round(sum(ms) / len(ms), 4),
                "N_GW_SIG": n_gw,
                "N_SUGGESTIVE": n_sug,
            }
        )

    return rows


def _compute_lambda_gc(p_values: List[float]) -> float:
    """Compute genomic inflation factor (lambda_gc).

    lambda_gc = median(chi2_observed) / median(chi2_expected)
    where chi2_expected under null with 1 df has median = 0.4549.

    Args:
        p_values: List of p-values

    Returns:
        Genomic inflation factor
    """
    if not p_values:
        return 1.0

    # Convert p-values to chi-squared statistics (1 df)
    chi2_values = []
    for p in p_values:
        if p <= 0 or p >= 1:
            continue
        # chi2 = qchisq(1-p, df=1); for p close to 0, use log transform
        try:
            # Use inverse chi-squared: for df=1, chi2 = (z_score)^2
            # z = Phi^-1(1-p/2); chi2 = z^2
            # Approximate using erfinv
            z = _inverse_normal_cdf(1 - p / 2)
            chi2_values.append(z * z)
        except (ValueError, OverflowError):
            continue

    if not chi2_values:
        return 1.0

    chi2_values.sort()
    n = len(chi2_values)
    median_observed = (
        chi2_values[n // 2]
        if n % 2 == 1
        else (chi2_values[n // 2 - 1] + chi2_values[n // 2]) / 2
    )

    # Expected median of chi-squared with df=1
    expected_median = 0.4549364

    return median_observed / expected_median if expected_median > 0 else 1.0


def _inverse_normal_cdf(p: float) -> float:
    """Approximate inverse normal CDF (probit function).

    Uses rational approximation from Abramowitz and Stegun.

    Args:
        p: Probability in (0, 1)

    Returns:
        z-score
    """
    if p <= 0 or p >= 1:
        raise ValueError(f"p must be in (0, 1), got {p}")

    if p < 0.5:
        return -_inverse_normal_cdf(1 - p)

    # Rational approximation for p >= 0.5
    t = math.sqrt(-2 * math.log(1 - p))

    c0 = 2.515517
    c1 = 0.802853
    c2 = 0.010328
    d1 = 1.432788
    d2 = 0.189269
    d3 = 0.001308

    z = t - (c0 + c1 * t + c2 * t * t) / (1 + d1 * t + d2 * t * t + d3 * t * t * t)
    return z

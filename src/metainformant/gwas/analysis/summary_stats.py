"""GWAS summary statistics output utilities.

This module provides functions for writing GWAS results in standard formats,
including full summary statistics TSV, significant hits tables, and JSON summaries.

Missing values are written as ``NA`` in the numeric columns (the marker used by
GWAS-SSF, PLINK and GCTA, and read as missing by pandas and R), never as a
fabricated number: an absent q-value or allele frequency is not 1 or 0.
"""

from __future__ import annotations

import json
import math
import random
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Union

from metainformant.core.utils import logging
from metainformant.gwas.analysis.correction import fdr_correction, lambda_gc_from_p_values

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    if not TYPE_CHECKING:  # runtime-only fallback; type checkers keep the real module type
        np = None

logger = logging.get_logger(__name__)

#: Marker written for a numeric value that is absent or not finite.
MISSING_VALUE = "NA"
#: Marker written for an absent text field (chromosome, id, alleles), as in VCF.
MISSING_TEXT = "."

_SUMMARY_HEADER = "CHR\tPOS\tSNP\tREF\tALT\tBETA\tSE\tP\tQ_FDR\tN\tMAF\n"
_Q_VALUE_KEYS = ("q_value", "q_fdr", "Q_FDR")
_N_KEYS = ("n_samples", "N")
_MAF_KEYS = ("maf", "MAF")


def _first_present(record: Dict[str, Any], keys: Sequence[str]) -> Any:
    """Return the first value among ``keys`` that is present and not None."""
    for key in keys:
        value = record.get(key)
        if value is not None:
            return value
    return None


def _finite_float(value: Any) -> Optional[float]:
    """Return ``value`` as a finite float, or None if it is missing or not a finite number."""
    if value is None or isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _first_finite(record: Dict[str, Any], keys: Sequence[str]) -> Optional[float]:
    """Return the first value among ``keys`` that is a finite number, or None."""
    for key in keys:
        number = _finite_float(record.get(key))
        if number is not None:
            return number
    return None


def _tested_results(results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """The results whose test produced a p-value.

    A test that failed carries NaN (or None) for its p-value: it is missing, not p = 1, so the summaries
    below leave it out of their counts, calibration and top hits instead of counting it as a null test.
    """
    return [result for result in results if _finite_float(result.get("p_value")) is not None]


def _format_number(value: Any, spec: str) -> str:
    """Format ``value`` with ``spec``; a missing or non-finite value becomes ``NA``."""
    number = _finite_float(value)
    return MISSING_VALUE if number is None else format(number, spec)


def _format_count(value: Any) -> str:
    """Format an integer-valued column (position, sample size); missing becomes ``NA``."""
    number = _finite_float(value)
    if number is None:
        return MISSING_VALUE
    return str(int(number)) if number.is_integer() else format(number, ".6g")


def _format_text(value: Any) -> str:
    """Format a text column; a missing or empty value becomes ``.``."""
    if isinstance(value, (list, tuple)):
        value = ",".join(str(item) for item in value)
    return MISSING_TEXT if value is None or value == "" else str(value)


def _benjamini_hochberg_q_values(results: List[Dict[str, Any]]) -> List[Optional[float]]:
    """Benjamini-Hochberg q-values over the p-values of ``results``.

    The family of tests is every result with a valid p-value (finite, in [0, 1]);
    results without one get None.
    """
    p_values: List[Optional[float]] = []
    for result in results:
        p_value = _finite_float(result.get("p_value"))
        p_values.append(p_value if p_value is not None and 0.0 <= p_value <= 1.0 else None)

    valid = [p for p in p_values if p is not None]
    if not valid:
        return [None] * len(results)

    adjusted = iter(fdr_correction(valid, method="bh")["adjusted_p_values"])
    return [next(adjusted) if p is not None else None for p in p_values]


def minor_allele_frequency(genotype_row: Sequence[Any]) -> Optional[float]:
    """Minor allele frequency from one variant's 0/1/2 alternate-allele dosages.

    Any other code (such as -1 or NaN) is treated as missing, as in
    ``quality.filter_by_maf``. Returns None when no genotype is called.
    """
    alt_alleles = 0
    n_called = 0
    for genotype in genotype_row:
        if genotype in (0, 1, 2):
            alt_alleles += int(genotype)
            n_called += 1
    if n_called == 0:
        return None
    alt_frequency = alt_alleles / (2 * n_called)
    return min(alt_frequency, 1.0 - alt_frequency)


def _format_row(
    result: Dict[str, Any],
    variant: Dict[str, Any],
    computed_q_value: Optional[float],
    computed_maf: Optional[float],
) -> str:
    """Format one summary-statistics TSV row.

    The result's own finite q-value and MAF win; ``computed_*`` fill in only when the
    result carries none. Anything still unknown is written as ``NA``.
    """
    q_value = _first_finite(result, _Q_VALUE_KEYS)
    if q_value is None:
        q_value = computed_q_value
    maf = _first_finite(result, _MAF_KEYS)
    if maf is None:
        maf = computed_maf

    return (
        f"{_format_text(variant.get('chrom'))}\t{_format_count(variant.get('pos'))}\t"
        f"{_format_text(variant.get('id'))}\t{_format_text(variant.get('ref'))}\t{_format_text(variant.get('alt'))}\t"
        f"{_format_number(result.get('beta'), '.6g')}\t{_format_number(result.get('se'), '.6g')}\t"
        f"{_format_number(result.get('p_value'), '.6e')}\t{_format_number(q_value, '.6e')}\t"
        f"{_format_count(_first_present(result, _N_KEYS))}\t{_format_number(maf, '.4f')}\n"
    )


def _derived_columns(
    results: List[Dict[str, Any]],
    genotypes: Optional[Sequence[Sequence[Any]]],
) -> tuple[List[Optional[float]], List[Optional[float]]]:
    """Compute the per-result q-values and MAFs that results do not already carry.

    Q_FDR is Benjamini-Hochberg over the p-values of ``results``; MAF comes from the
    genotype row of the same index. Both are None where the inputs do not exist.
    """
    if genotypes is not None and len(genotypes) != len(results):
        raise ValueError(f"Genotypes ({len(genotypes)}) and results ({len(results)}) must have the same length")

    needs_q = any(_first_finite(result, _Q_VALUE_KEYS) is None for result in results)
    q_values: List[Optional[float]] = [None] * len(results)
    if needs_q:
        q_values = _benjamini_hochberg_q_values(results)

    mafs: List[Optional[float]] = [None] * len(results)
    if genotypes is not None:
        for i, result in enumerate(results):
            if _first_finite(result, _MAF_KEYS) is None:
                mafs[i] = minor_allele_frequency(genotypes[i])
    return q_values, mafs


def write_summary_statistics(
    results: List[Dict[str, Any]],
    variant_info: List[Dict[str, Any]],
    output_path: Union[str, Path],
    *,
    genotypes: Optional[Sequence[Sequence[Any]]] = None,
) -> Path:
    """Write full GWAS summary statistics to TSV file.

    Output format follows standard GWAS summary statistics convention, with
    FDR-adjusted q-values:
    CHR  POS  SNP  REF  ALT  BETA  SE  P  Q_FDR  N  MAF

    Nothing is invented for values the results do not have. Q_FDR is the
    result's own ``q_value``/``q_fdr``/``Q_FDR`` when present, otherwise the
    Benjamini-Hochberg q-value over the p-values of ``results`` (so pass the full
    set of tests, not a subset). MAF is the result's own ``maf``/``MAF`` when
    present, otherwise computed from ``genotypes`` when given. Any other missing
    or non-finite numeric value is written as ``NA``, and a missing text field as
    ``.``.

    Args:
        results: List of association test result dictionaries
        variant_info: List of variant info dictionaries with chrom, pos, id, ref, alt
        output_path: Path to output TSV file
        genotypes: Optional genotype matrix (variants x samples) of 0/1/2
            alternate-allele dosages, with any other code treated as missing,
            used to compute MAF for results that lack one.

    Returns:
        Path to written file

    Raises:
        ValueError: If results and variant_info (or genotypes) have different lengths
    """
    if len(results) != len(variant_info):
        raise ValueError(
            f"Results ({len(results)}) and variant_info ({len(variant_info)}) " "must have the same length"
        )

    q_values, mafs = _derived_columns(results, genotypes)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Writing summary statistics for {len(results)} variants to {output_path}")

    with open(output_path, "w") as f:
        f.write(_SUMMARY_HEADER)
        for result, vinfo, q_value, maf in zip(results, variant_info, q_values, mafs):
            f.write(_format_row(result, vinfo, q_value, maf))

    logger.info(f"Summary statistics written to {output_path}")
    return output_path


def write_significant_hits(
    results: List[Dict[str, Any]],
    variant_info: List[Dict[str, Any]],
    output_path: Union[str, Path],
    threshold: float = 5e-8,
    *,
    genotypes: Optional[Sequence[Sequence[Any]]] = None,
) -> Path:
    """Write genome-wide significant hits to TSV file.

    Columns and missing-value handling are those of :func:`write_summary_statistics`.
    Q_FDR and MAF are derived over the full ``results`` before the significant
    hits are selected, so a q-value reflects every test performed.

    Args:
        results: List of association test result dictionaries
        variant_info: List of variant info dictionaries
        output_path: Path to output TSV file
        threshold: Significance threshold (default: 5e-8 genome-wide)
        genotypes: Optional genotype matrix (variants x samples), as for
            :func:`write_summary_statistics`.

    Returns:
        Path to written file

    Raises:
        ValueError: If genotypes and results have different lengths
    """
    q_values, mafs = _derived_columns(results, genotypes)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    significant = []
    for result, vinfo, q_value, maf in zip(results, variant_info, q_values, mafs):
        p_value = _finite_float(result.get("p_value"))
        if p_value is not None and p_value < threshold:
            significant.append((p_value, result, vinfo, q_value, maf))

    # Sort by p-value
    significant.sort(key=lambda hit: hit[0])

    logger.info(f"Writing {len(significant)} significant hits (threshold={threshold:.2e}) to {output_path}")

    with open(output_path, "w") as f:
        f.write(_SUMMARY_HEADER)
        for _, result, vinfo, q_value, maf in significant:
            f.write(_format_row(result, vinfo, q_value, maf))

    return output_path


def create_results_summary(
    results: List[Dict[str, Any]],
    output_path: Union[str, Path],
    threshold: float = 5e-8,
) -> Dict[str, Any]:
    """Create JSON summary of GWAS results.

    A test that failed has NaN (or None) for its p-value. It is left out of ``n_variants_tested``, the
    inflation factor, the significance counts and the top hits, and is reported in ``n_tests_failed``.

    Args:
        results: List of association test result dictionaries
        output_path: Path to output JSON file
        threshold: Genome-wide significance threshold

    Returns:
        Summary dictionary with n_variants_tested, n_tests_failed, lambda_gc, n_significant, top_hits
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    tested = _tested_results(results)
    p_values = [float(r["p_value"]) for r in tested]
    n_tests = len(p_values)
    n_failed = len(results) - n_tests

    # Compute genomic inflation factor via the canonical chi-square conversion
    lambda_gc = lambda_gc_from_p_values(p_values)
    if lambda_gc is None:
        lambda_gc = 1.0

    # Count significant hits
    n_significant = sum(1 for p in p_values if p < threshold)
    n_suggestive = sum(1 for p in p_values if p < 1e-5)

    # Top hits (sorted by p-value; variant_index is the position in ``results``)
    indexed_results = [(i, r) for i, r in enumerate(results) if _finite_float(r.get("p_value")) is not None]
    indexed_results.sort(key=lambda x: float(x[1]["p_value"]))
    top_hits = []
    for idx, result in indexed_results[:20]:
        hit = {
            "variant_index": idx,
            "p_value": float(result["p_value"]),
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
        "n_tests_failed": n_failed,
        "lambda_gc": lambda_gc,
        "n_significant_genome_wide": n_significant,
        "n_suggestive": n_suggestive,
        "significance_threshold": threshold,
        "top_hits": top_hits,
    }

    with open(output_path, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(
        f"Results summary: {n_tests} variants, lambda_gc={lambda_gc:.3f}, " f"{n_significant} genome-wide significant"
    )
    if n_failed:
        logger.warning(f"Results summary: {n_failed} failed test(s) left out of the counts and top hits")
    return summary


# ── Public analysis functions (extracted from orchestrator inline logic) ──


def normal_cdf(x: float) -> float:
    """Standard normal CDF via Abramowitz & Stegun (7.1.26) approximation.

    Pure-Python — no scipy dependency.
    """
    sign = 1 if x >= 0 else -1
    x = abs(x) / math.sqrt(2)
    t = 1.0 / (1.0 + 0.3275911 * x)
    y = 1.0 - (
        (((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592
    ) * t * math.exp(-x * x)
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
    seed: int = 42,
) -> Dict[str, Any]:
    """Calculate Linkage Disequilibrium (r^2) decay across genomic distances.

    Variant pairs are drawn at random, so the result depends on ``seed``. The draw uses its own generator:
    the same inputs and seed always give the same curve, and the global ``random`` state is neither read nor
    advanced.

    Args:
        genotype_matrix: List of genotypes (variants x samples)
        variants_info: List of variant metadata
        max_distance: Maximum base-pair distance to compute LD
        bins: Number of distance bins for decay curve
        seed: Seed of the generator that draws the variant pairs

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

    # Basic LD logic (this scales at O(N^2), so we sample random pairs if there are many)
    n_vars = G.shape[0]
    sample_pairs = min(10000, n_vars * (n_vars - 1) // 2)

    r2_vals = []
    distances = []

    rng = random.Random(seed)
    pairs: set[tuple[int, int]] = set()
    attempts = 0
    while len(pairs) < sample_pairs and attempts < sample_pairs * 5:
        attempts += 1
        i = rng.randint(0, n_vars - 1)
        j = rng.randint(0, n_vars - 1)
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

    A test that failed has NaN (or None) for its p-value and statistics. It is left out of every
    distribution and count here and reported in ``n_tests_failed``.

    Args:
        association_results: List of result dicts, each containing at least
            ``p_value``, ``beta``, ``se``, ``MAF``, ``chrom``, ``pos``, ``snp``.
        significance_threshold: Genome-wide threshold.

    Returns:
        Nested summary dictionary.
    """
    import statistics as _stats

    n_failed = len(association_results) - len(_tested_results(association_results))
    association_results = _tested_results(association_results)

    p_values = [r["p_value"] for r in association_results]
    betas = [r["beta"] for r in association_results]
    ses = [r["se"] for r in association_results]
    mafs = [r["MAF"] for r in association_results]
    abs_betas = [abs(b) for b in betas]
    z_scores = [b / s if s > 0 else 0.0 for b, s in zip(betas, ses)]

    # λ_GC via the canonical chi-square conversion
    lambda_gc = lambda_gc_from_p_values(p_values)
    if lambda_gc is None:
        lambda_gc = 1.0

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
        "n_tests_failed": n_failed,
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

    A test that failed has NaN (or None) for its p-value and effect; it takes no part in the test.

    Args:
        association_results: List of dicts with ``beta`` and ``p_value`` keys.
        significance_threshold: GW threshold for sub-analysis.

    Returns:
        Dict with n_positive, n_negative, z_statistic, p_value_two_sided,
        interpretation, and GW-hit breakdowns.
    """
    association_results = _tested_results(association_results)
    betas = [r["beta"] for r in association_results]
    n_pos = sum(1 for b in betas if b > 0)
    n_neg = sum(1 for b in betas if b < 0)
    n_total = n_pos + n_neg

    if n_total == 0:
        return {"n_positive_beta": 0, "n_negative_beta": 0, "n_total": 0}

    p_hat = n_pos / n_total
    se_hat = math.sqrt(0.25 / n_total)
    z_sign = (p_hat - 0.5) / se_hat
    p_sign = 2.0 * min(normal_cdf(z_sign), 1.0 - normal_cdf(z_sign))

    gw_hits = [r for r in association_results if r["p_value"] < significance_threshold]
    n_gw_pos = sum(1 for r in gw_hits if r["beta"] > 0)

    interpretation = (
        "Significant directional enrichment (positive effects dominant)"
        if p_sign < 0.05 and n_pos > n_neg
        else (
            "Significant directional enrichment (negative effects dominant)"
            if p_sign < 0.05 and n_neg > n_pos
            else "No significant directional enrichment (balanced +/- effects)"
        )
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

    A test that failed has NaN (or None) for its p-value and statistics; it is left out of the table.

    Returns:
        List of per-chromosome summary dicts sorted by chromosome.
    """
    chrom_data: Dict[str, List[Dict[str, Any]]] = {}
    for r in _tested_results(association_results):
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

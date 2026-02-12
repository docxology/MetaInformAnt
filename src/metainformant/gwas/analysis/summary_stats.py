"""GWAS summary statistics output utilities.

This module provides functions for writing GWAS results in standard formats,
including full summary statistics TSV, significant hits tables, and JSON summaries.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

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
            f"Results ({len(results)}) and variant_info ({len(variant_info)}) " "must have the same length"
        )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Writing summary statistics for {len(results)} variants to {output_path}")

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
        (result, vinfo) for result, vinfo in zip(results, variant_info) if result.get("p_value", 1.0) < threshold
    ]

    # Sort by p-value
    significant.sort(key=lambda x: x[0].get("p_value", 1.0))

    logger.info(f"Writing {len(significant)} significant hits (threshold={threshold:.2e}) to {output_path}")

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
        f"Results summary: {n_tests} variants, lambda_gc={lambda_gc:.3f}, " f"{n_significant} genome-wide significant"
    )
    return summary


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
    median_observed = chi2_values[n // 2] if n % 2 == 1 else (chi2_values[n // 2 - 1] + chi2_values[n // 2]) / 2

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

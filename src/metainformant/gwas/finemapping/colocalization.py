"""Multi-trait colocalization analysis for GWAS fine-mapping.

Extends pairwise colocalization to multiple traits, provides specialized
GWAS-eQTL colocalization with gene-level summaries, implements the CLPP
(colocalization posterior probability) method from individual PIPs, and
supports region-based colocalization with automatic window selection.

Evidence is built from Wakefield approximate Bayes factors, which need the
scale of each trait's effect estimates: pass standard errors, or sample sizes
(with an optional MAF and case fractions), as described in
:func:`metainformant.gwas.finemapping.credible_sets.compute_log_bayes_factors`.
Traits given neither are scored on the z-score scale. All accumulation is in log
space, so strong signals (|z| in the tens or hundreds) cannot overflow.

References:
    - Giambartolomei et al. (2014) PLoS Genetics 10:e1004383 (coloc)
    - Hormozdiari et al. (2016) AJHG 98:667-680 (eCAVIAR / CLPP)
    - Wallace (2020) PLoS Genetics 16:e1008720 (coloc extensions)
    - Wakefield (2009) Genetic Epidemiology 33:79-86 (approximate Bayes factors)
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from numbers import Real
from typing import Any, cast

from metainformant.core.utils.logging import get_logger
from metainformant.gwas.finemapping.credible_sets import (
    colocalization as _pairwise_coloc,
    compute_log_bayes_factors,
    log_sum_exp,
)

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def multi_trait_coloc(
    z_scores: dict,
    ld_matrix: list[list[float]] | None = None,
    *,
    prior_variances: Mapping[str, float] | None = None,
    standard_errors: Mapping[str, float | Sequence[float]] | None = None,
    sample_sizes: Mapping[str, float | Sequence[float]] | None = None,
    maf: float | Sequence[float] | None = None,
    case_fractions: Mapping[str, float] | None = None,
) -> dict:
    """Multi-trait colocalization extending pairwise analysis to N traits.

    For each pair of traits, runs pairwise colocalization and then
    synthesizes results into an overall summary. Also computes a
    multi-trait posterior probability of a single shared causal variant
    using the product-of-Bayes-factors approach. Bayes factors are summed
    in log space, so the result stays finite however strong the signals are.

    The per-trait mappings below are keyed by the trait names of ``z_scores``;
    a trait absent from every mapping is scored on the z-score scale. See
    :func:`metainformant.gwas.finemapping.credible_sets.compute_log_bayes_factors`
    for the Bayes-factor model and its defaults.

    Args:
        z_scores: Dictionary mapping trait names to lists of Z-scores.
            All Z-score lists must have the same length (same variant set).
            Example: {"height": [1.2, -0.5, ...], "bmi": [0.8, 2.1, ...]}
        ld_matrix: Optional LD matrix shared across all traits.
        prior_variances: Optional prior variance W of the true effect per trait.
        standard_errors: Optional standard errors per trait (scalar or one per
            variant).
        sample_sizes: Optional sample size per trait (scalar or one per variant).
        maf: Optional allele frequency (scalar or one per variant) shared by all
            traits, used with ``sample_sizes``.
        case_fractions: Optional case fraction per case-control trait; traits
            left out are treated as quantitative.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - pairwise: Dict of pairwise coloc results keyed by trait pair
            - multi_trait_pp_shared: Posterior probability all traits share
              a single causal variant
            - n_traits: Number of traits analyzed
            - n_variants: Number of variants in the region
            - summary: Human-readable summary
    """
    if not z_scores:
        return {"status": "error", "message": "No Z-scores provided"}

    trait_names = list(z_scores.keys())
    n_traits = len(trait_names)

    if n_traits < 2:
        return {
            "status": "error",
            "message": f"Need at least 2 traits for colocalization, got {n_traits}",
        }

    # Validate all Z-score lists have the same length
    lengths = {name: len(zs) for name, zs in z_scores.items()}
    unique_lengths = set(lengths.values())
    if len(unique_lengths) > 1:
        return {
            "status": "error",
            "message": f"Z-score lengths differ across traits: {lengths}",
        }

    n_variants = list(lengths.values())[0]
    if n_variants == 0:
        return {"status": "error", "message": "Z-score lists are empty"}

    per_trait_inputs = {
        "prior_variances": prior_variances,
        "standard_errors": standard_errors,
        "sample_sizes": sample_sizes,
        "case_fractions": case_fractions,
    }
    for label, mapping in per_trait_inputs.items():
        unknown = sorted(str(name) for name in (mapping or {}) if name not in z_scores)
        if unknown:
            return {"status": "error", "message": f"{label} names unknown traits: {unknown}"}

    logger.info(f"Running multi-trait colocalization: {n_traits} traits, {n_variants} variants")

    def _scale_kwargs(name: str) -> dict[str, Any]:
        return {
            "prior_variance": (prior_variances or {}).get(name),
            "standard_errors": (standard_errors or {}).get(name),
            "sample_size": (sample_sizes or {}).get(name),
            "maf": maf,
            "case_fraction": (case_fractions or {}).get(name),
        }

    # Per-trait log Bayes factors (also validates the scale inputs of every trait)
    all_log_bfs: dict[str, list[float]] = {}
    try:
        for name in trait_names:
            all_log_bfs[name] = compute_log_bayes_factors(z_scores[name], **_scale_kwargs(name))
    except ValueError as exc:
        return {"status": "error", "message": str(exc)}

    # Step 1: Run all pairwise colocalization tests
    pairwise_results: dict[str, dict] = {}
    for i in range(n_traits):
        for j in range(i + 1, n_traits):
            trait_a = trait_names[i]
            trait_b = trait_names[j]
            pair_key = f"{trait_a}_vs_{trait_b}"
            kw_a = _scale_kwargs(trait_a)
            kw_b = _scale_kwargs(trait_b)

            result = _pairwise_coloc(
                z_scores_1=z_scores[trait_a],
                z_scores_2=z_scores[trait_b],
                ld_matrix=ld_matrix,
                prior_variance_1=kw_a["prior_variance"],
                prior_variance_2=kw_b["prior_variance"],
                standard_errors_1=kw_a["standard_errors"],
                standard_errors_2=kw_b["standard_errors"],
                sample_size_1=kw_a["sample_size"],
                sample_size_2=kw_b["sample_size"],
                maf=maf,
                case_fraction_1=kw_a["case_fraction"],
                case_fraction_2=kw_b["case_fraction"],
            )
            pairwise_results[pair_key] = result

    # Step 2: Compute multi-trait shared causal variant probability.
    # For each variant the product of BFs across traits is the sum of log BFs.
    log_product_bfs = [sum(all_log_bfs[name][v_idx] for name in trait_names) for v_idx in range(n_variants)]
    log_sum_product = log_sum_exp(log_product_bfs)

    # Null BF for multi-trait (all traits have no association)
    log_null = math.log(n_variants)

    # Individual trait BFs (at least one trait has association)
    # Use sum of individual Bayes factors as a rough alternative
    log_sum_individual = log_sum_exp(lbf for log_bfs in all_log_bfs.values() for lbf in log_bfs)

    # Normalize
    log_total_evidence = log_sum_exp([log_null, log_sum_individual, log_sum_product])
    pp_shared = math.exp(log_sum_product - log_total_evidence)

    # Step 3: Summarize pairwise PP(H4) values
    pp_h4_values = []
    for pair_key, result in pairwise_results.items():
        if result.get("status") == "success":
            pp_h4_values.append(result.get("PP_H4", 0.0))

    mean_pp_h4 = sum(pp_h4_values) / len(pp_h4_values) if pp_h4_values else 0.0
    min_pp_h4 = min(pp_h4_values) if pp_h4_values else 0.0

    # Determine summary
    if pp_shared > 0.8:
        verdict = "Strong evidence for a shared causal variant across all traits"
    elif pp_shared > 0.5:
        verdict = "Moderate evidence for a shared causal variant"
    elif mean_pp_h4 > 0.5:
        verdict = "Pairwise colocalization supports some shared signals"
    else:
        verdict = "Limited evidence for shared causal variants"

    summary = (
        f"{verdict}. Multi-trait PP(shared)={pp_shared:.4f}, "
        f"mean pairwise PP(H4)={mean_pp_h4:.4f}, "
        f"min pairwise PP(H4)={min_pp_h4:.4f}"
    )

    return {
        "status": "success",
        "pairwise": pairwise_results,
        "multi_trait_pp_shared": float(pp_shared),
        "n_traits": n_traits,
        "n_variants": n_variants,
        "summary": summary,
    }


def eqtl_coloc(
    gwas_z: list[float],
    eqtl_z: list[float],
    gene_id: str,
    ld_matrix: list[list[float]] | None = None,
    *,
    gwas_prior_variance: float | None = None,
    eqtl_prior_variance: float | None = None,
    gwas_standard_errors: float | Sequence[float] | None = None,
    eqtl_standard_errors: float | Sequence[float] | None = None,
    gwas_sample_size: float | Sequence[float] | None = None,
    eqtl_sample_size: float | Sequence[float] | None = None,
    maf: float | Sequence[float] | None = None,
    gwas_case_fraction: float | None = None,
) -> dict:
    """GWAS-eQTL colocalization with gene-level summary.

    Specialized wrapper around pairwise colocalization for the common use
    case of testing whether a GWAS signal colocalizes with an eQTL signal
    for a specific gene. Provides gene-level interpretation and variant
    annotation.

    The Bayes factors need the scale of each study's effect estimates: pass
    standard errors, or sample sizes (with ``maf`` and, for a case-control GWAS,
    ``gwas_case_fraction``). The eQTL trait is quantitative. A study given neither
    is scored on the z-score scale. See
    :func:`metainformant.gwas.finemapping.credible_sets.compute_log_bayes_factors`
    for the model and the default priors.

    Args:
        gwas_z: Z-scores from GWAS summary statistics.
        eqtl_z: Z-scores from eQTL analysis for the target gene.
        gene_id: Gene identifier (e.g., ENSG ID or gene symbol) for
            annotation and reporting.
        ld_matrix: Optional LD matrix for the region.
        gwas_prior_variance: Prior variance W of the true GWAS effect.
        eqtl_prior_variance: Prior variance W of the true eQTL effect.
        gwas_standard_errors: GWAS standard errors (scalar or one per variant).
        eqtl_standard_errors: eQTL standard errors (scalar or one per variant).
        gwas_sample_size: GWAS sample size (scalar or one per variant).
        eqtl_sample_size: eQTL sample size (scalar or one per variant).
        maf: Allele frequency shared by both studies, for the sample-size model.
        gwas_case_fraction: Fraction of cases when the GWAS is case-control.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - gene_id: The gene identifier
            - coloc_result: Full colocalization result
            - pp_h4: Posterior probability of shared causal variant
            - is_colocalized: Boolean (PP_H4 > 0.8)
            - lead_variant_gwas: Index of lead GWAS variant
            - lead_variant_eqtl: Index of lead eQTL variant
            - shared_variant: Index of most likely shared variant (if H4)
            - interpretation: Human-readable interpretation
    """
    if not gwas_z or not eqtl_z:
        return {"status": "error", "message": "GWAS and eQTL Z-scores are required"}

    if len(gwas_z) != len(eqtl_z):
        return {
            "status": "error",
            "message": (f"Z-score lengths differ: GWAS={len(gwas_z)}, eQTL={len(eqtl_z)}"),
        }

    n_variants = len(gwas_z)
    logger.info(f"Running GWAS-eQTL colocalization for gene {gene_id} ({n_variants} variants)")

    # Run pairwise colocalization
    coloc_result = _pairwise_coloc(
        z_scores_1=gwas_z,
        z_scores_2=eqtl_z,
        ld_matrix=ld_matrix,
        prior_variance_1=gwas_prior_variance,
        prior_variance_2=eqtl_prior_variance,
        standard_errors_1=gwas_standard_errors,
        standard_errors_2=eqtl_standard_errors,
        sample_size_1=gwas_sample_size,
        sample_size_2=eqtl_sample_size,
        maf=maf,
        case_fraction_1=gwas_case_fraction,
    )

    if coloc_result.get("status") != "success":
        return {
            "status": "error",
            "message": f"Colocalization failed: {coloc_result.get('message', 'unknown')}",
            "gene_id": gene_id,
        }

    pp_h4 = coloc_result.get("PP_H4", 0.0)
    is_colocalized = pp_h4 > 0.8

    # Find lead variants for each trait
    lead_gwas = max(range(n_variants), key=lambda i: abs(gwas_z[i]))
    lead_eqtl = max(range(n_variants), key=lambda i: abs(eqtl_z[i]))

    # Most likely shared variant: highest product of Bayes factors, i.e. the
    # largest posterior of being the shared causal variant given H4
    snp_pp_h4 = coloc_result["snp_pp_h4"]
    shared_variant = max(range(n_variants), key=lambda i: snp_pp_h4[i])

    # Build interpretation
    if is_colocalized:
        interp = (
            f"Strong colocalization (PP.H4={pp_h4:.3f}) between GWAS signal "
            f"and eQTL for {gene_id}. The most likely shared causal variant "
            f"is at index {shared_variant}, suggesting this gene mediates "
            f"the GWAS association."
        )
    elif pp_h4 > 0.5:
        interp = (
            f"Moderate evidence for colocalization (PP.H4={pp_h4:.3f}) between "
            f"GWAS signal and eQTL for {gene_id}. The shared causal variant "
            f"hypothesis is plausible but not definitive."
        )
    elif coloc_result.get("PP_H3", 0.0) > 0.5:
        interp = (
            f"Evidence suggests distinct causal variants (PP.H3="
            f"{coloc_result['PP_H3']:.3f}) for GWAS and eQTL signals at "
            f"{gene_id}. Lead GWAS variant (index {lead_gwas}) differs "
            f"from lead eQTL variant (index {lead_eqtl})."
        )
    else:
        pp_h1 = coloc_result.get("PP_H1", 0.0)
        pp_h2 = coloc_result.get("PP_H2", 0.0)
        interp = (
            f"Limited colocalization evidence for {gene_id}. "
            f"PP.H4={pp_h4:.3f}, PP.H3={coloc_result.get('PP_H3', 0):.3f}. "
            f"GWAS-only PP.H1={pp_h1:.3f}, eQTL-only PP.H2={pp_h2:.3f}."
        )

    return {
        "status": "success",
        "gene_id": gene_id,
        "coloc_result": coloc_result,
        "pp_h4": float(pp_h4),
        "is_colocalized": is_colocalized,
        "lead_variant_gwas": lead_gwas,
        "lead_variant_eqtl": lead_eqtl,
        "shared_variant": shared_variant,
        "interpretation": interp,
    }


def compute_clpp(
    pip_1: list[float],
    pip_2: list[float],
) -> dict:
    """Compute colocalization posterior probability from individual PIPs (CLPP).

    The CLPP method (Hormozdiari et al. 2016) estimates the probability
    of a shared causal variant using per-variant posterior inclusion
    probabilities from two independent fine-mapping analyses. The CLPP
    for a shared causal variant at position i is PIP1_i * PIP2_i, and
    the overall CLPP is the sum across all variants.

    Args:
        pip_1: Posterior inclusion probabilities for trait 1.
        pip_2: Posterior inclusion probabilities for trait 2.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - clpp: Overall colocalization posterior probability
            - per_variant_clpp: Per-variant CLPP values
            - max_clpp_index: Index of variant with highest CLPP
            - max_clpp_value: Highest per-variant CLPP
            - is_colocalized: Boolean (CLPP > 0.01 threshold)
            - n_variants: Number of variants
    """
    if not pip_1 or not pip_2:
        return {"status": "error", "message": "Both PIP lists are required"}

    if len(pip_1) != len(pip_2):
        return {
            "status": "error",
            "message": (f"PIP lengths differ: trait 1 has {len(pip_1)}, " f"trait 2 has {len(pip_2)}"),
        }

    n_variants = len(pip_1)
    logger.debug(f"Computing CLPP for {n_variants} variants")

    # Per-variant CLPP: PIP1_i * PIP2_i
    per_variant_clpp = [p1 * p2 for p1, p2 in zip(pip_1, pip_2)]

    # Overall CLPP: sum of per-variant CLPPs
    clpp = sum(per_variant_clpp)

    # Find variant with highest CLPP
    max_idx = 0
    max_val = per_variant_clpp[0]
    for i in range(1, n_variants):
        if per_variant_clpp[i] > max_val:
            max_val = per_variant_clpp[i]
            max_idx = i

    # Colocalization threshold: CLPP > 0.01 (Hormozdiari et al.)
    is_colocalized = clpp > 0.01

    return {
        "status": "success",
        "clpp": float(clpp),
        "per_variant_clpp": per_variant_clpp,
        "max_clpp_index": max_idx,
        "max_clpp_value": float(max_val),
        "is_colocalized": is_colocalized,
        "n_variants": n_variants,
    }


def regional_coloc(
    z_scores_dict: dict,
    region: dict,
    ld_matrix: list[list[float]] | None = None,
    *,
    prior_variances: Mapping[str, float] | None = None,
    standard_errors: Mapping[str, float | Sequence[float]] | None = None,
    sample_sizes: Mapping[str, float | Sequence[float]] | None = None,
    maf: float | Sequence[float] | None = None,
    case_fractions: Mapping[str, float] | None = None,
) -> dict:
    """Run colocalization across a genomic region with automatic window selection.

    Takes Z-scores for multiple traits across a genomic region and
    automatically selects sub-windows centered on the strongest signals.
    Runs colocalization within each window and returns results for the
    best window as well as per-window summaries.

    Args:
        z_scores_dict: Dictionary mapping trait names to lists of Z-scores.
            All lists must have the same length.
        region: Dictionary describing the genomic region with keys:
            - chrom: Chromosome name/number
            - start: Start position (bp)
            - end: End position (bp)
            - positions: Optional list of variant positions (bp) matching
              the Z-score indices. If not provided, variants are assumed
              to be uniformly spaced.
        ld_matrix: Optional LD matrix for the region.
        prior_variances: Optional prior variance W of the true effect per trait.
        standard_errors: Optional standard errors per trait: a scalar, or one
            value per variant of the whole region (sliced with each window).
        sample_sizes: Optional sample size per trait, scalar or per variant.
        maf: Optional allele frequency shared by all traits, scalar or per variant.
        case_fractions: Optional case fraction per case-control trait.
            These scale inputs are forwarded to :func:`multi_trait_coloc`.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - best_window: Coloc result for the best sub-window
            - all_windows: List of per-window coloc results
            - n_windows: Number of windows tested
            - region: The input region specification
    """
    if not z_scores_dict:
        return {"status": "error", "message": "No Z-scores provided"}

    trait_names = list(z_scores_dict.keys())
    n_traits = len(trait_names)

    if n_traits < 2:
        return {
            "status": "error",
            "message": f"Need at least 2 traits, got {n_traits}",
        }

    # Validate lengths
    lengths = {name: len(zs) for name, zs in z_scores_dict.items()}
    unique_lengths = set(lengths.values())
    if len(unique_lengths) > 1:
        return {
            "status": "error",
            "message": f"Z-score lengths differ: {lengths}",
        }

    n_variants = list(lengths.values())[0]
    if n_variants == 0:
        return {"status": "error", "message": "Z-score lists are empty"}

    chrom = region.get("chrom", "unknown")
    start = region.get("start", 0)
    end = region.get("end", n_variants)
    positions = region.get("positions")

    if positions is None:
        # Generate uniform positions
        if n_variants > 1:
            step = (end - start) / (n_variants - 1) if n_variants > 1 else 1
            positions = [start + int(i * step) for i in range(n_variants)]
        else:
            positions = [start]

    if len(positions) != n_variants:
        return {
            "status": "error",
            "message": (f"Positions length ({len(positions)}) does not match " f"Z-score length ({n_variants})"),
        }

    # Per-variant scale inputs must cover the whole region so windows can slice them
    per_variant_inputs: list[tuple[str, float | Sequence[float] | None]] = [("maf", maf)]
    per_variant_inputs += [(f"standard_errors[{name!r}]", v) for name, v in (standard_errors or {}).items()]
    per_variant_inputs += [(f"sample_sizes[{name!r}]", v) for name, v in (sample_sizes or {}).items()]
    for label, values in per_variant_inputs:
        count = _variant_count(values)
        if count is not None and count != n_variants:
            return {
                "status": "error",
                "message": f"{label} has {count} values but the region has {n_variants} variants",
            }

    logger.info(
        f"Running regional colocalization: {n_traits} traits, " f"{n_variants} variants, chr{chrom}:{start}-{end}"
    )

    # Determine window size and windows
    region_size = end - start
    window_size = _select_window_size(region_size, n_variants)
    windows = _generate_windows(positions, window_size, n_variants)

    all_window_results: list[dict] = []

    for w_idx, (w_start_idx, w_end_idx) in enumerate(windows):
        # Slice Z-scores for this window
        window_z: dict[str, list[float]] = {}
        for name in trait_names:
            window_z[name] = z_scores_dict[name][w_start_idx:w_end_idx]

        # Slice LD matrix if provided
        window_ld: list[list[float]] | None = None
        if ld_matrix is not None:
            window_ld = [row[w_start_idx:w_end_idx] for row in ld_matrix[w_start_idx:w_end_idx]]

        # Run multi-trait coloc for this window, slicing per-variant scale inputs
        window_result = multi_trait_coloc(
            z_scores=window_z,
            ld_matrix=window_ld,
            prior_variances=prior_variances,
            standard_errors=_slice_trait_values(standard_errors, w_start_idx, w_end_idx),
            sample_sizes=_slice_trait_values(sample_sizes, w_start_idx, w_end_idx),
            maf=_slice_variant_values(maf, w_start_idx, w_end_idx),
            case_fractions=case_fractions,
        )

        window_result["window_index"] = w_idx
        window_result["window_start_idx"] = w_start_idx
        window_result["window_end_idx"] = w_end_idx
        window_result["window_start_bp"] = positions[w_start_idx]
        window_result["window_end_bp"] = positions[min(w_end_idx - 1, n_variants - 1)]

        all_window_results.append(window_result)

    # Select best window by highest multi-trait PP(shared)
    best_window = None
    best_pp = -1.0
    for wr in all_window_results:
        pp = wr.get("multi_trait_pp_shared", 0.0)
        if pp > best_pp:
            best_pp = pp
            best_window = wr

    return {
        "status": "success",
        "best_window": best_window,
        "all_windows": all_window_results,
        "n_windows": len(all_window_results),
        "region": {
            "chrom": chrom,
            "start": start,
            "end": end,
            "n_variants": n_variants,
        },
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _variant_count(values: float | Sequence[float] | None) -> int | None:
    """Return how many per-variant values ``values`` holds, or None for a scalar/absent value."""
    if values is None or isinstance(values, Real):
        return None
    return len(cast("Sequence[float]", values))


def _slice_variant_values(values: float | Sequence[float] | None, start: int, end: int) -> Any:
    """Slice a per-variant sequence to ``[start, end)``; scalars and None pass through."""
    if _variant_count(values) is None:
        return values
    return list(cast("Sequence[float]", values))[start:end]


def _slice_trait_values(
    values: Mapping[str, float | Sequence[float]] | None, start: int, end: int
) -> dict[str, Any] | None:
    """Slice every per-variant sequence of a per-trait mapping to ``[start, end)``."""
    if values is None:
        return None
    return {name: _slice_variant_values(trait_values, start, end) for name, trait_values in values.items()}


def _select_window_size(region_size: int, n_variants: int) -> int:
    """Select an appropriate window size for regional colocalization.

    Uses heuristics based on region size and variant density.

    Args:
        region_size: Size of the region in base pairs.
        n_variants: Number of variants in the region.

    Returns:
        Window size in number of variants.
    """
    # Aim for windows of ~250kb or ~100 variants, whichever is smaller
    if n_variants <= 100:
        return n_variants  # Single window

    # Target ~100 variants per window
    target_size = min(100, n_variants)

    # But ensure at least 20 variants per window
    return max(target_size, 20)


def _generate_windows(
    positions: list[int],
    window_size: int,
    n_variants: int,
) -> list[tuple[int, int]]:
    """Generate overlapping windows across the region.

    Windows overlap by 50% to avoid missing signals at window boundaries.

    Args:
        positions: Variant positions in base pairs.
        window_size: Number of variants per window.
        n_variants: Total number of variants.

    Returns:
        List of (start_index, end_index) tuples.
    """
    if window_size >= n_variants:
        return [(0, n_variants)]

    step = max(window_size // 2, 1)  # 50% overlap
    windows: list[tuple[int, int]] = []

    start = 0
    while start < n_variants:
        end = min(start + window_size, n_variants)
        windows.append((start, end))
        if end >= n_variants:
            break
        start += step

    return windows

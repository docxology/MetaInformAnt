"""Multi-trait colocalization analysis for GWAS fine-mapping.

Extends pairwise colocalization to multiple traits, provides specialized
GWAS-eQTL colocalization with gene-level summaries, implements the CLPP
(colocalization posterior probability) method from individual PIPs, and
supports region-based colocalization over sliding windows.

Evidence is built from Wakefield approximate Bayes factors, which need the
scale of each trait's effect estimates: pass standard errors, or sample sizes
(with an optional MAF and case fractions), as described in
:func:`metainformant.gwas.finemapping.credible_sets.compute_log_bayes_factors`.
Traits given neither are scored on the z-score scale. All accumulation is in log
space, so strong signals (|z| in the tens or hundreds) cannot overflow.

The probability that ``N`` traits share one causal variant
(:func:`multi_trait_coloc`) is the exact posterior of coloc's model generalised to ``N``
traits (a partition model in the manner of MOLOC and HyPrColoc): each trait has no causal
variant or one, traits with the same causal variant form a group, and a group of ``k``
traits at one variant has prior ``p1 * (p12 / p1) ** (k - 1)``. For two traits it is coloc's
``PP.H4``. It was a heuristic, ``sum_v prod_t BF / (n + sum_t sum_v BF + sum_v prod_t BF)``,
whose flat priors and missing "distinct variants" hypothesis made it call a null trait
"shared" whenever the other traits shared a signal (see :func:`multi_trait_coloc`).

References:
    - Giambartolomei et al. (2014) PLoS Genetics 10:e1004383 (coloc)
    - Giambartolomei et al. (2018) Bioinformatics 34:2538-2545 (MOLOC)
    - Foley et al. (2021) Nature Communications 12:764 (HyPrColoc)
    - Hormozdiari et al. (2016) AJHG 98:667-680 (eCAVIAR / CLPP)
    - Wallace (2020) PLoS Genetics 16:e1008720 (coloc extensions)
    - Wakefield (2009) Genetic Epidemiology 33:79-86 (approximate Bayes factors)
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from numbers import Real
from typing import Any, cast

import numpy as np

from metainformant.core.utils.logging import get_logger
from metainformant.gwas.finemapping.credible_sets import colocalization as _pairwise_coloc, compute_log_bayes_factors

logger = get_logger(__name__)

#: PP above which a shared causal variant is called strong evidence (coloc's conventional PP.H4 cut-off).
COLOC_STRONG_THRESHOLD = 0.8
#: PP above which a shared causal variant is called moderate evidence.
COLOC_MODERATE_THRESHOLD = 0.5
#: Most traits :func:`multi_trait_coloc` evaluates exactly: the cost grows as ``2 ** N`` per variant.
MAX_EXACT_TRAITS = 10

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
    prior_p1: float = 1e-4,
    prior_p12: float = 1e-5,
) -> dict:
    """Multi-trait colocalization extending pairwise analysis to N traits.

    For each pair of traits, runs pairwise colocalization and then
    synthesizes results into an overall summary. ``multi_trait_pp_shared`` is the
    posterior probability that **all** traits have the same causal variant.

    The model is coloc's, generalised to ``N`` traits (see the module docstring). Each trait has
    either no causal variant in the region or exactly one; traits with the same causal variant
    form a group. A set of ``k`` traits sharing one specific variant has prior
    ``p1 * (p12 / p1) ** (k - 1)`` (``p1`` for one trait, ``p12`` for two, ``p12 ** 2 / p1`` for
    three, ...), and the Bayes factors of the traits multiply. The probability of "all share" is
    the weight of the configurations in which every trait has the same variant, over the weight of
    every configuration, computed exactly (a dynamic programme over variants and trait subsets, in
    log space, so the result stays finite however strong the signals are). For two traits it is
    coloc's ``PP.H4``, and the pairwise results use the same priors. The traits are taken as
    independent studies (no sample overlap), and each has at most one causal variant: LD is not
    modelled, so two distinct causal variants in strong LD can still look shared, as they do in
    coloc. Cost is ``O(n_variants * 2 ** n_traits)`` and at most :data:`MAX_EXACT_TRAITS` traits are
    accepted; analyse subsets of more traits.

    The previous statistic was a heuristic that gave the null hypothesis, "some trait
    associated" and "all share" equal prior weight and had no hypothesis of distinct causal
    variants. In simulations it called a null third trait "shared" in 99% of regions where the other
    two shared a signal, and "shared" in 87% where the third had a distinct signal, at the 0.8
    cut-off; the posterior here flags about 1% of those.

    The per-trait mappings below are keyed by the trait names of ``z_scores``;
    a trait absent from every mapping is scored on the z-score scale. See
    :func:`metainformant.gwas.finemapping.credible_sets.compute_log_bayes_factors`
    for the Bayes-factor model and its defaults.

    Args:
        z_scores: Dictionary mapping trait names to lists of Z-scores.
            All Z-score lists must have the same length (same variant set).
            Example: {"height": [1.2, -0.5, ...], "bmi": [0.8, 2.1, ...]}
        ld_matrix: Accepted for API compatibility and **not used**: the Bayes-factor model does not
            model LD, so the result is the same with or without it (a warning is logged if given).
        prior_variances: Optional prior variance W of the true effect per trait.
        standard_errors: Optional standard errors per trait (scalar or one per
            variant).
        sample_sizes: Optional sample size per trait (scalar or one per variant).
        maf: Optional allele frequency (scalar or one per variant) shared by all
            traits, used with ``sample_sizes``.
        case_fractions: Optional case fraction per case-control trait; traits
            left out are treated as quantitative.
        prior_p1: Prior probability that a given variant is the causal variant of a given trait.
        prior_p12: Prior probability that a given variant is the causal variant of two given traits
            (coloc's ``p12``); larger groups follow ``p1 * (p12 / p1) ** (k - 1)``.

    Returns:
        Dictionary with keys:
            - status: "success" or "error"
            - pairwise: Dict of pairwise coloc results keyed by trait pair
            - multi_trait_pp_shared: Posterior probability all traits share
              a single causal variant (``PP.H4`` for two traits)
            - multi_trait_pp_all_associated: Posterior probability that every trait has a
              causal variant in the region, shared or not (an upper bound on the previous key; the
              gap between the two says the traits are associated but not shared)
            - snp_pp_shared: Per variant, the probability of being the shared causal variant
              given that all traits share one (sums to 1)
            - priors: The priors used: ``p1``, ``p12`` and ``p0`` (the weight of "no trait associated")
            - n_traits: Number of traits analyzed
            - n_variants: Number of variants in the region
            - summary: Human-readable summary
    """
    if ld_matrix is not None:
        logger.warning("multi_trait_coloc: ld_matrix is not used by the Wakefield Bayes-factor model and is ignored")

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

    if n_traits > MAX_EXACT_TRAITS:
        return {
            "status": "error",
            "message": (
                f"{n_traits} traits exceed the {MAX_EXACT_TRAITS} for which the posterior is computed exactly "
                "(its cost doubles with every trait); analyse subsets of the traits"
            ),
        }

    n_variants = list(lengths.values())[0]
    if n_variants == 0:
        return {"status": "error", "message": "Z-score lists are empty"}

    try:
        log_prior_p0 = _log_null_weight(n_traits, prior_p1, prior_p12)
    except ValueError as exc:
        return {"status": "error", "message": str(exc)}

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
                ld_matrix=None,
                prior_p1=prior_p1,
                prior_p2=prior_p1,
                prior_p12=prior_p12,
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
            if result.get("status") != "success":
                reason = result.get("message", "unknown error")
                return {"status": "error", "message": f"pairwise colocalization of {pair_key} failed: {reason}"}
            pairwise_results[pair_key] = result

    # Step 2: the exact posterior that every trait shares one causal variant
    posterior = _shared_causal_posterior([all_log_bfs[name] for name in trait_names], prior_p1, prior_p12, log_prior_p0)
    pp_shared = posterior["pp_shared"]

    # Step 3: Summarize pairwise PP(H4) values
    pp_h4_values = [result["PP_H4"] for result in pairwise_results.values()]
    mean_pp_h4 = sum(pp_h4_values) / len(pp_h4_values)
    min_pp_h4 = min(pp_h4_values)

    # Determine summary
    if pp_shared > COLOC_STRONG_THRESHOLD:
        verdict = "Strong evidence for a shared causal variant across all traits"
    elif pp_shared > COLOC_MODERATE_THRESHOLD:
        verdict = "Moderate evidence for a shared causal variant"
    elif mean_pp_h4 > COLOC_MODERATE_THRESHOLD:
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
        "multi_trait_pp_all_associated": float(posterior["pp_all_associated"]),
        "snp_pp_shared": posterior["snp_pp_shared"],
        "priors": {"p1": prior_p1, "p12": prior_p12, "p0": math.exp(log_prior_p0)},
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
        ld_matrix: Accepted for API compatibility and **not used**: the Bayes-factor model does not
            model LD, so the result is the same with or without it (a warning is logged if given).
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
    if ld_matrix is not None:
        logger.warning("eqtl_coloc: ld_matrix is not used by the Wakefield Bayes-factor model and is ignored")

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
        ld_matrix=None,
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

    pp_h4 = coloc_result["PP_H4"]
    is_colocalized = pp_h4 > COLOC_STRONG_THRESHOLD

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
    elif pp_h4 > COLOC_MODERATE_THRESHOLD:
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
    """Run colocalization across a genomic region over overlapping windows.

    Takes Z-scores for multiple traits across a genomic region and cuts the
    variants into windows of up to 100 variants that overlap by half (one window for
    a region of 100 variants or fewer); the windows follow the order of the variants, not
    the signals. Runs :func:`multi_trait_coloc` within each window and returns the window
    with the highest ``multi_trait_pp_shared`` as well as per-window summaries.

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
        ld_matrix: Accepted for API compatibility and **not used** (see :func:`multi_trait_coloc`);
            a warning is logged if given.
        prior_variances: Optional prior variance W of the true effect per trait.
        standard_errors: Optional standard errors per trait: a scalar, or one
            value per variant of the whole region (sliced with each window).
        sample_sizes: Optional sample size per trait, scalar or per variant.
        maf: Optional allele frequency shared by all traits, scalar or per variant.
        case_fractions: Optional case fraction per case-control trait.
            These scale inputs are forwarded to :func:`multi_trait_coloc`.

    Returns:
        Dictionary with keys:
            - status: "success", or "error" when the input is invalid or no
              window could be analysed (every window returned an error)
            - best_window: Coloc result for the best sub-window, chosen among
              the windows that succeeded (absent on error)
            - all_windows: List of per-window coloc results, including any
              window that failed, which carries its own ``status`` and ``message``
            - n_windows: Number of windows tested
            - n_failed_windows: Number of windows whose analysis returned an error
            - region: The input region specification
            - message: Why the analysis failed (only when ``status`` is "error")
    """
    if ld_matrix is not None:
        logger.warning("regional_coloc: ld_matrix is not used by the Wakefield Bayes-factor model and is ignored")

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
    window_size = _select_window_size(n_variants)
    windows = _generate_windows(positions, window_size, n_variants)

    all_window_results: list[dict] = []

    for w_idx, (w_start_idx, w_end_idx) in enumerate(windows):
        # Slice Z-scores for this window
        window_z: dict[str, list[float]] = {}
        for name in trait_names:
            window_z[name] = z_scores_dict[name][w_start_idx:w_end_idx]

        # Run multi-trait coloc for this window, slicing per-variant scale inputs
        window_result = multi_trait_coloc(
            z_scores=window_z,
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

    # A window that errored has no PP(shared): it can neither be the best window nor count as a result.
    # Its own ``status``/``message`` stay in ``all_windows``.
    region_summary = {
        "chrom": chrom,
        "start": start,
        "end": end,
        "n_variants": n_variants,
    }
    successful_windows = [wr for wr in all_window_results if wr.get("status") == "success"]
    n_failed_windows = len(all_window_results) - len(successful_windows)

    if not successful_windows:
        reasons = list(
            dict.fromkeys(str(wr.get("message", "unknown error")) for wr in all_window_results)  # distinct, in order
        )
        return {
            "status": "error",
            "message": f"Colocalization failed in all {len(all_window_results)} windows: " + "; ".join(reasons),
            "all_windows": all_window_results,
            "n_windows": len(all_window_results),
            "n_failed_windows": n_failed_windows,
            "region": region_summary,
        }

    # Select best window by highest multi-trait PP(shared); the earliest window wins a tie
    best_window = successful_windows[0]
    for wr in successful_windows[1:]:
        if wr["multi_trait_pp_shared"] > best_window["multi_trait_pp_shared"]:
            best_window = wr

    return {
        "status": "success",
        "best_window": best_window,
        "all_windows": all_window_results,
        "n_windows": len(all_window_results),
        "n_failed_windows": n_failed_windows,
        "region": region_summary,
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _log_group_prior(size: int, prior_p1: float, prior_p12: float) -> float:
    """Log prior that one specific variant is the causal variant of one specific set of ``size`` traits.

    ``p1 * (p12 / p1) ** (size - 1)``: ``p1`` for one trait and ``p12`` for two, as in coloc.
    """
    return math.log(prior_p1) + (size - 1) * (math.log(prior_p12) - math.log(prior_p1))


def _log_null_weight(n_traits: int, prior_p1: float, prior_p12: float) -> float:
    """Log weight of "no trait has a causal variant": one minus the prior mass of every non-empty trait set.

    For two traits this is ``1 - p1 - p2 - p12``, the weight :func:`...credible_sets.colocalization`
    gives its H0.

    Raises:
        ValueError: If a prior is not a positive finite number, or the priors sum to more than 1.
    """
    candidates: tuple[tuple[str, Any], ...] = (("prior_p1", prior_p1), ("prior_p12", prior_p12))
    for name, value in candidates:
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be a positive finite number, got {value!r}")
    mass = sum(
        math.comb(n_traits, size) * math.exp(_log_group_prior(size, prior_p1, prior_p12))
        for size in range(1, n_traits + 1)
    )
    if mass > 1.0:
        raise ValueError(
            f"Prior probabilities sum to more than 1 for {n_traits} traits (p1={prior_p1}, p12={prior_p12})"
        )
    return math.log1p(-mass) if mass < 1.0 else -math.inf


def _shared_causal_posterior(
    log_bfs: Sequence[Sequence[float]], prior_p1: float, prior_p12: float, log_prior_p0: float
) -> dict[str, Any]:
    """Exact posterior that every trait shares one causal variant, and that every trait has one.

    ``log_bfs[t][v]`` is the natural-log Bayes factor of trait ``t`` at variant ``v``. A configuration
    gives every trait no causal variant or one; the traits that have the same one form a group, and its
    weight is the product of the group priors (:func:`_log_group_prior`) and of the Bayes factors of
    its traits. The sum over all configurations is a dynamic programme over variants: ``state[U]`` is the log
    weight of the configurations over the variants seen so far in which exactly the traits in the bitmask
    ``U`` have a causal variant, and each variant hosts at most one group, of traits not yet placed. Every
    term is positive, so nothing cancels, and everything is in log space.

    Returns:
        ``pp_shared``, ``pp_all_associated`` and ``snp_pp_shared`` (see :func:`multi_trait_coloc`).
    """
    log_bf = np.asarray(log_bfs, dtype=float)
    n_traits, n_variants = log_bf.shape
    n_sets = 1 << n_traits
    full = n_sets - 1

    # log_group[S, v]: log prior of the trait set S sharing a variant, plus its summed log Bayes factors at v
    summed = np.zeros((n_sets, n_variants))
    log_group = np.full((n_sets, n_variants), -np.inf)
    for mask in range(1, n_sets):
        lowest = mask & -mask
        summed[mask] = summed[mask ^ lowest] + log_bf[lowest.bit_length() - 1]
        log_group[mask] = _log_group_prior(mask.bit_count(), prior_p1, prior_p12) + summed[mask]

    everything = np.arange(n_sets)
    unplaced = [None] + [everything[(everything & mask) == 0] for mask in range(1, n_sets)]

    state = np.full(n_sets, -np.inf)
    state[0] = 0.0  # no variant seen: no trait has a causal variant yet (weight 1, rescaled by the null weight)
    for variant in range(n_variants):
        updated = state.copy()  # this variant is nobody's causal variant
        for mask in range(1, n_sets):
            rest = unplaced[mask]
            updated[rest | mask] = np.logaddexp(updated[rest | mask], state[rest] + log_group[mask, variant])
        state = updated

    log_total = np.logaddexp(log_prior_p0, np.logaddexp.reduce(state[1:]))
    log_shared = np.logaddexp.reduce(log_group[full])
    snp_pp_shared = np.exp(summed[full] - np.logaddexp.reduce(summed[full]))
    return {
        "pp_shared": float(min(1.0, math.exp(log_shared - log_total))),
        "pp_all_associated": float(min(1.0, math.exp(state[full] - log_total))),
        "snp_pp_shared": [float(value) for value in snp_pp_shared],
    }


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


def _select_window_size(n_variants: int) -> int:
    """Window size, in variants, for regional colocalization: 100, or all variants if there are fewer.

    The size depends only on the number of variants; the base-pair length of the region plays no part.

    Args:
        n_variants: Number of variants in the region.

    Returns:
        Window size in number of variants.
    """
    return min(100, n_variants)


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

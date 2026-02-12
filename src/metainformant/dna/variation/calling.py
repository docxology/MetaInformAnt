"""Variant calling utilities for DNA sequence analysis.

This module provides tools for calling variants from pileup data using
Bayesian genotyping models, filtering variant calls, merging calls from
multiple callers, computing variant statistics, and annotating variant
sequence context for mutation signature analysis.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Valid DNA bases
BASES = ("A", "C", "G", "T")

# Transition pairs (purine-purine or pyrimidine-pyrimidine substitutions)
TRANSITIONS = {
    ("A", "G"),
    ("G", "A"),
    ("C", "T"),
    ("T", "C"),
}

# All possible single-nucleotide substitution types
SUBSTITUTION_TYPES = {(ref, alt) for ref in BASES for alt in BASES if ref != alt}


def call_variants_pileup(
    pileup_data: list[dict],
    min_depth: int = 10,
    min_qual: float = 20.0,
    min_alt_freq: float = 0.2,
) -> list[dict]:
    """Call SNP variants from pileup data using a simple Bayesian genotyping model.

    For each position in the pileup, computes the posterior probability of each
    possible genotype (homozygous reference, heterozygous, homozygous alternate)
    given the observed base counts and base qualities. Positions meeting the
    minimum depth, quality, and alternate allele frequency thresholds are
    reported as variant calls.

    Args:
        pileup_data: List of pileup position dicts, each containing:
            - chrom: Chromosome/contig name
            - pos: 1-based position
            - ref: Reference base
            - depth: Read depth at this position
            - bases: Dict mapping bases to counts (e.g., {"A": 50, "G": 12})
            - quals: Optional dict mapping bases to average quality scores
        min_depth: Minimum read depth to consider (default 10).
        min_qual: Minimum variant quality score in Phred scale (default 20.0).
        min_alt_freq: Minimum alternate allele frequency (default 0.2).

    Returns:
        List of variant call dicts, each containing:
            - chrom: Chromosome name
            - pos: 1-based position
            - ref: Reference allele
            - alt: Alternate allele
            - qual: Variant quality (Phred-scaled)
            - depth: Total depth
            - alt_depth: Alternate allele depth
            - alt_freq: Alternate allele frequency
            - genotype: Predicted genotype ("0/0", "0/1", or "1/1")
            - genotype_quality: Genotype quality (Phred-scaled)
            - filter: "PASS" or list of filter reasons

    Raises:
        ValueError: If pileup_data is empty or contains invalid entries.
    """
    if not pileup_data:
        raise ValueError("Pileup data must not be empty")

    variants: list[dict] = []

    for pileup_pos in pileup_data:
        chrom = pileup_pos.get("chrom", "unknown")
        pos = pileup_pos.get("pos", 0)
        ref = pileup_pos.get("ref", "N").upper()
        depth = pileup_pos.get("depth", 0)
        bases = pileup_pos.get("bases", {})
        quals = pileup_pos.get("quals", {})

        # Skip low-depth positions
        if depth < min_depth:
            continue

        # Normalize base counts
        base_counts: Dict[str, int] = {}
        for base in BASES:
            base_counts[base] = bases.get(base, 0)

        # Find alternate alleles
        ref_count = base_counts.get(ref, 0)
        total_count = sum(base_counts.values())

        if total_count == 0:
            continue

        for alt_base in BASES:
            if alt_base == ref:
                continue

            alt_count = base_counts.get(alt_base, 0)
            if alt_count == 0:
                continue

            alt_freq = alt_count / total_count

            if alt_freq < min_alt_freq:
                continue

            # Bayesian genotyping
            avg_qual = quals.get(alt_base, 30.0)
            error_rate = 10 ** (-avg_qual / 10.0)

            genotype_result = genotype_variants(
                allele_counts={"ref": ref_count, "alt": alt_count},
                ploidy=2,
                error_rate=error_rate,
            )

            gt = genotype_result.get("genotype", "0/1")
            gq = genotype_result.get("genotype_quality", 0.0)

            # Compute variant quality (Phred-scaled probability that site is NOT reference)
            prob_ref_hom = genotype_result.get("prob_hom_ref", 1.0)
            if prob_ref_hom < 1.0:
                variant_qual = -10.0 * math.log10(max(prob_ref_hom, 1e-300))
            else:
                variant_qual = 0.0

            if variant_qual < min_qual:
                continue

            # Apply filters
            filters: list[str] = []
            if depth < min_depth:
                filters.append("LowDepth")
            if alt_freq < min_alt_freq:
                filters.append("LowAF")
            if variant_qual < min_qual:
                filters.append("LowQual")

            filter_str = "PASS" if not filters else ";".join(filters)

            variants.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt_base,
                    "qual": round(variant_qual, 2),
                    "depth": depth,
                    "alt_depth": alt_count,
                    "alt_freq": round(alt_freq, 4),
                    "genotype": gt,
                    "genotype_quality": round(gq, 2),
                    "filter": filter_str,
                }
            )

    # Sort by chromosome and position
    variants.sort(key=lambda v: (v["chrom"], v["pos"]))
    logger.debug(f"Called {len(variants)} variants from {len(pileup_data)} pileup positions")
    return variants


def genotype_variants(
    allele_counts: dict,
    ploidy: int = 2,
    error_rate: float = 0.01,
) -> dict:
    """Compute maximum likelihood genotype from allele counts.

    Uses a simple binomial/multinomial likelihood model to compute the
    posterior probability of each possible genotype given the observed
    allele counts and sequencing error rate.

    Args:
        allele_counts: Dict with "ref" and "alt" integer counts.
        ploidy: Organism ploidy level (default 2 for diploid).
        error_rate: Per-base sequencing error rate (default 0.01).

    Returns:
        Dict containing:
            - genotype: Predicted genotype string ("0/0", "0/1", or "1/1")
            - genotype_quality: Phred-scaled quality of the genotype call
            - prob_hom_ref: Posterior probability of homozygous reference
            - prob_het: Posterior probability of heterozygous
            - prob_hom_alt: Posterior probability of homozygous alternate
            - status: "success"

    Raises:
        ValueError: If allele_counts is empty or ploidy < 1.
    """
    if not allele_counts:
        raise ValueError("Allele counts must not be empty")
    if ploidy < 1:
        raise ValueError("Ploidy must be >= 1")

    ref_count = allele_counts.get("ref", 0)
    alt_count = allele_counts.get("alt", 0)
    total = ref_count + alt_count

    if total == 0:
        return {
            "genotype": "./.",
            "genotype_quality": 0.0,
            "prob_hom_ref": 0.0,
            "prob_het": 0.0,
            "prob_hom_alt": 0.0,
            "status": "success",
        }

    # Clamp error rate
    error_rate = max(1e-10, min(0.5, error_rate))

    # For diploid genotyping, compute log-likelihoods of each genotype
    # P(data | genotype) using binomial model
    # Genotype 0/0: expected alt freq = error_rate
    # Genotype 0/1: expected alt freq = 0.5
    # Genotype 1/1: expected alt freq = 1.0 - error_rate

    if ploidy == 2:
        genotypes = [
            ("0/0", error_rate),
            ("0/1", 0.5),
            ("1/1", 1.0 - error_rate),
        ]
    else:
        # For other ploidies, create genotypes for each dosage
        genotypes = []
        for dosage in range(ploidy + 1):
            expected_freq = dosage / ploidy
            # Adjust for error rate
            adj_freq = expected_freq * (1.0 - error_rate) + (1.0 - expected_freq) * error_rate
            gt_str = "/".join(["0"] * (ploidy - dosage) + ["1"] * dosage)
            genotypes.append((gt_str, adj_freq))

    # Compute log-likelihoods
    log_likelihoods: list[Tuple[str, float]] = []
    for gt_label, expected_alt_freq in genotypes:
        # Binomial log-likelihood
        p = max(1e-300, min(1.0 - 1e-300, expected_alt_freq))
        ll = alt_count * math.log(p) + ref_count * math.log(1.0 - p)
        log_likelihoods.append((gt_label, ll))

    # Convert to posterior probabilities (uniform prior)
    max_ll = max(ll for _, ll in log_likelihoods)
    scaled = [(gt, math.exp(ll - max_ll)) for gt, ll in log_likelihoods]
    total_prob = sum(prob for _, prob in scaled)

    posteriors = {gt: prob / total_prob for gt, prob in scaled}

    # Find maximum likelihood genotype
    best_gt = max(posteriors, key=posteriors.get)  # type: ignore[arg-type]
    best_prob = posteriors[best_gt]

    # Genotype quality (Phred-scaled probability of wrong genotype)
    wrong_prob = 1.0 - best_prob
    gq = -10.0 * math.log10(max(wrong_prob, 1e-300)) if wrong_prob > 0 else 99.0

    result = {
        "genotype": best_gt,
        "genotype_quality": round(min(99.0, gq), 2),
        "prob_hom_ref": round(posteriors.get("0/0", 0.0), 6),
        "prob_het": round(posteriors.get("0/1", 0.0), 6),
        "prob_hom_alt": round(posteriors.get("1/1", 0.0), 6),
        "status": "success",
    }

    logger.debug(f"Genotyped: ref={ref_count}, alt={alt_count} -> {best_gt} (GQ={gq:.1f})")
    return result


def filter_variants(
    variants: list[dict],
    min_depth: int = 10,
    min_qual: float = 30.0,
    max_strand_bias: float = 0.01,
) -> list[dict]:
    """Apply multi-criteria filters to variant calls.

    Filters variants based on depth, quality, and optionally strand bias.
    Variants that fail any filter have the corresponding reason appended
    to their filter field; passing variants are marked "PASS".

    Args:
        variants: List of variant call dicts (from call_variants_pileup or similar).
        min_depth: Minimum read depth (default 10).
        min_qual: Minimum variant quality in Phred scale (default 30.0).
        max_strand_bias: Maximum strand bias p-value threshold (default 0.01).
            If a variant has a "strand_bias" field with p-value below this
            threshold, it is flagged. Ignored if field is absent.

    Returns:
        List of filtered variant dicts. Each dict contains all original fields
        plus an updated "filter" field. Only variants that PASS all filters
        are included in the output.

    Raises:
        ValueError: If variants list is None.
    """
    if variants is None:
        raise ValueError("Variants list must not be None")

    filtered: list[dict] = []

    for variant in variants:
        result = dict(variant)
        fail_reasons: list[str] = []

        # Depth filter
        depth = variant.get("depth", 0)
        if depth < min_depth:
            fail_reasons.append(f"LowDepth({depth}<{min_depth})")

        # Quality filter
        qual = variant.get("qual", 0.0)
        if qual < min_qual:
            fail_reasons.append(f"LowQual({qual:.1f}<{min_qual:.1f})")

        # Strand bias filter
        strand_bias = variant.get("strand_bias")
        if strand_bias is not None and strand_bias < max_strand_bias:
            fail_reasons.append(f"StrandBias(p={strand_bias:.4f})")

        # Alternate allele frequency sanity check
        alt_freq = variant.get("alt_freq", 0.0)
        if alt_freq <= 0.0:
            fail_reasons.append("NoAltReads")

        if not fail_reasons:
            result["filter"] = "PASS"
            filtered.append(result)
        else:
            result["filter"] = ";".join(fail_reasons)
            # Do not include failed variants in output
            logger.debug(
                f"Filtered out variant at {variant.get('chrom', '?')}:{variant.get('pos', '?')} "
                f"- {result['filter']}"
            )

    logger.debug(
        f"Filtered {len(variants)} variants -> {len(filtered)} passed " f"(min_depth={min_depth}, min_qual={min_qual})"
    )
    return filtered


def merge_variant_calls(
    call_sets: list[list[dict]],
    min_callers: int = 2,
) -> list[dict]:
    """Merge variant calls from multiple callers requiring minimum agreement.

    Groups variants by position (chrom + pos + ref + alt) across all call
    sets, then retains only those called by at least `min_callers` callers.
    For merged variants, quality and depth are taken from the maximum across
    callers, and allele frequencies are averaged.

    Args:
        call_sets: List of call sets, where each call set is a list of
            variant dicts (from call_variants_pileup, filter_variants, etc.).
        min_callers: Minimum number of callers that must agree for a variant
            to be included in the merged output (default 2).

    Returns:
        List of merged variant dicts, each containing:
            - All standard variant fields (chrom, pos, ref, alt, etc.)
            - num_callers: Number of callers that called this variant
            - caller_indices: List of 0-based caller indices that agreed
            - max_qual: Maximum quality across callers
            - avg_alt_freq: Average alternate allele frequency across callers
            - status: "success" for the overall merge

    Raises:
        ValueError: If call_sets is empty or min_callers < 1.
    """
    if not call_sets:
        raise ValueError("Call sets must not be empty")
    if min_callers < 1:
        raise ValueError("min_callers must be >= 1")
    if min_callers > len(call_sets):
        logger.warning(
            f"min_callers ({min_callers}) > number of call sets ({len(call_sets)}); " "no variants will pass"
        )

    # Index variants by position key
    variant_registry: Dict[str, list[Tuple[int, dict]]] = defaultdict(list)

    for caller_idx, call_set in enumerate(call_sets):
        for variant in call_set:
            chrom = variant.get("chrom", "unknown")
            pos = variant.get("pos", 0)
            ref = variant.get("ref", "")
            alt = variant.get("alt", "")
            key = f"{chrom}:{pos}:{ref}:{alt}"
            variant_registry[key].append((caller_idx, variant))

    # Merge variants meeting minimum caller threshold
    merged: list[dict] = []

    for key, caller_variants in variant_registry.items():
        num_callers = len(caller_variants)
        if num_callers < min_callers:
            continue

        # Aggregate across callers
        caller_indices = [idx for idx, _ in caller_variants]
        quals = [v.get("qual", 0.0) for _, v in caller_variants]
        depths = [v.get("depth", 0) for _, v in caller_variants]
        alt_freqs = [v.get("alt_freq", 0.0) for _, v in caller_variants]
        alt_depths = [v.get("alt_depth", 0) for _, v in caller_variants]

        # Use the first variant as a template
        _, template = caller_variants[0]

        merged_variant = {
            "chrom": template.get("chrom", "unknown"),
            "pos": template.get("pos", 0),
            "ref": template.get("ref", ""),
            "alt": template.get("alt", ""),
            "qual": round(max(quals), 2),
            "depth": max(depths),
            "alt_depth": max(alt_depths),
            "alt_freq": round(sum(alt_freqs) / len(alt_freqs), 4),
            "genotype": template.get("genotype", "./."),
            "genotype_quality": round(max(v.get("genotype_quality", 0.0) for _, v in caller_variants), 2),
            "filter": "PASS",
            "num_callers": num_callers,
            "caller_indices": caller_indices,
            "max_qual": round(max(quals), 2),
            "avg_alt_freq": round(sum(alt_freqs) / len(alt_freqs), 4),
        }

        merged.append(merged_variant)

    # Sort by position
    merged.sort(key=lambda v: (v["chrom"], v["pos"]))

    logger.debug(
        f"Merged {sum(len(cs) for cs in call_sets)} variants from "
        f"{len(call_sets)} callers -> {len(merged)} consensus variants "
        f"(min_callers={min_callers})"
    )
    return merged


def compute_variant_stats(
    variants: list[dict],
) -> dict:
    """Compute summary statistics for a set of variant calls.

    Calculates transition/transversion ratio, variant density, allele
    frequency spectrum, heterozygosity/homozygosity ratio, and per-chromosome
    variant counts.

    Args:
        variants: List of variant call dicts.

    Returns:
        Dict containing:
            - total_variants: Total number of variants
            - ti_tv_ratio: Transition/transversion ratio (None if no transversions)
            - transitions: Number of transitions
            - transversions: Number of transversions
            - variant_density: Variants per position (if positions span a range)
            - allele_frequency_spectrum: Dict mapping AF bins to counts
            - het_hom_ratio: Heterozygous/homozygous ratio (None if no hom)
            - het_count: Number of heterozygous calls
            - hom_count: Number of homozygous alternate calls
            - snp_count: Number of single-nucleotide variants
            - indel_count: Number of insertions/deletions
            - per_chrom: Dict mapping chromosome to variant count
            - substitution_spectrum: Dict mapping "X>Y" to count
            - status: "success"

    Raises:
        ValueError: If variants is None.
    """
    if variants is None:
        raise ValueError("Variants must not be None")

    if not variants:
        return {
            "total_variants": 0,
            "ti_tv_ratio": None,
            "transitions": 0,
            "transversions": 0,
            "variant_density": 0.0,
            "allele_frequency_spectrum": {},
            "het_hom_ratio": None,
            "het_count": 0,
            "hom_count": 0,
            "snp_count": 0,
            "indel_count": 0,
            "per_chrom": {},
            "substitution_spectrum": {},
            "status": "success",
        }

    transitions = 0
    transversions = 0
    het_count = 0
    hom_alt_count = 0
    snp_count = 0
    indel_count = 0
    per_chrom: Dict[str, int] = defaultdict(int)
    substitution_spectrum: Dict[str, int] = defaultdict(int)
    allele_freqs: list[float] = []
    positions: list[int] = []

    for variant in variants:
        ref = variant.get("ref", "").upper()
        alt = variant.get("alt", "").upper()
        chrom = variant.get("chrom", "unknown")
        pos = variant.get("pos", 0)
        gt = variant.get("genotype", "")
        af = variant.get("alt_freq", 0.0)

        per_chrom[chrom] += 1
        positions.append(pos)
        allele_freqs.append(af)

        # SNP vs indel classification
        if len(ref) == 1 and len(alt) == 1:
            snp_count += 1

            # Ti/Tv classification
            pair = (ref, alt)
            if pair in TRANSITIONS:
                transitions += 1
            elif pair in SUBSTITUTION_TYPES:
                transversions += 1

            # Substitution spectrum
            sub_key = f"{ref}>{alt}"
            substitution_spectrum[sub_key] += 1
        else:
            indel_count += 1

        # Het/hom classification
        if gt in ("0/1", "1/0", "0|1", "1|0"):
            het_count += 1
        elif gt in ("1/1", "1|1"):
            hom_alt_count += 1

    # Ti/Tv ratio
    ti_tv_ratio = transitions / transversions if transversions > 0 else None

    # Het/hom ratio
    het_hom_ratio = het_count / hom_alt_count if hom_alt_count > 0 else None

    # Variant density (variants per base pair, based on position span)
    if len(positions) >= 2:
        span = max(positions) - min(positions)
        variant_density = len(variants) / span if span > 0 else 0.0
    else:
        variant_density = 0.0

    # Allele frequency spectrum (binned)
    af_spectrum: Dict[str, int] = defaultdict(int)
    bins = [
        (0.0, 0.01, "singleton"),
        (0.01, 0.05, "rare"),
        (0.05, 0.10, "low_frequency"),
        (0.10, 0.25, "intermediate"),
        (0.25, 0.50, "common"),
        (0.50, 1.01, "high_frequency"),
    ]
    for af in allele_freqs:
        for low, high, label in bins:
            if low <= af < high:
                af_spectrum[label] += 1
                break

    logger.debug(
        f"Variant stats: {len(variants)} total, Ti/Tv={ti_tv_ratio}, "
        f"Het/Hom={het_hom_ratio}, SNPs={snp_count}, indels={indel_count}"
    )

    return {
        "total_variants": len(variants),
        "ti_tv_ratio": round(ti_tv_ratio, 4) if ti_tv_ratio is not None else None,
        "transitions": transitions,
        "transversions": transversions,
        "variant_density": round(variant_density, 6),
        "allele_frequency_spectrum": dict(af_spectrum),
        "het_hom_ratio": round(het_hom_ratio, 4) if het_hom_ratio is not None else None,
        "het_count": het_count,
        "hom_count": hom_alt_count,
        "snp_count": snp_count,
        "indel_count": indel_count,
        "per_chrom": dict(per_chrom),
        "substitution_spectrum": dict(substitution_spectrum),
        "status": "success",
    }


def annotate_variant_context(
    variants: list[dict],
    reference: str,
    window: int = 5,
) -> list[dict]:
    """Add sequence context and trinucleotide context to variant calls.

    For each variant, extracts the surrounding reference sequence context
    and the trinucleotide context (for mutation signature analysis using
    the standard 96-channel SBS framework).

    Args:
        variants: List of variant dicts with "pos" (0-based) and "ref"/"alt" fields.
        reference: Reference DNA sequence.
        window: Number of flanking bases on each side for context (default 5).

    Returns:
        List of annotated variant dicts, each containing all original fields plus:
            - upstream_context: Sequence upstream of the variant (window bp)
            - downstream_context: Sequence downstream of the variant (window bp)
            - full_context: Complete context string (upstream + [ref/alt] + downstream)
            - trinucleotide_context: 3-bp context (1 bp upstream + ref + 1 bp downstream)
            - trinucleotide_mutation: Trinucleotide mutation label (e.g., "ACA>ATA")
            - sbs_channel: SBS96 mutation channel in pyrimidine context

    Raises:
        ValueError: If reference is empty or variants is None.
    """
    if not reference:
        raise ValueError("Reference sequence must not be empty")
    if variants is None:
        raise ValueError("Variants list must not be None")

    ref_upper = reference.upper()
    annotated: list[dict] = []

    # Pyrimidine bases for SBS96 standardization
    complement = {"A": "T", "T": "A", "C": "C", "G": "G", "N": "N"}
    # Actually we need full complement for reverse complement
    full_complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

    for variant in variants:
        result = dict(variant)
        pos = variant.get("pos", 0)
        ref_allele = variant.get("ref", "").upper()
        alt_allele = variant.get("alt", "").upper()

        # Extract context windows
        upstream_start = max(0, pos - window)
        downstream_end = min(len(ref_upper), pos + len(ref_allele) + window)

        upstream_context = ref_upper[upstream_start:pos]
        downstream_context = ref_upper[pos + len(ref_allele) : downstream_end]
        full_context = f"{upstream_context}[{ref_allele}/{alt_allele}]{downstream_context}"

        result["upstream_context"] = upstream_context
        result["downstream_context"] = downstream_context
        result["full_context"] = full_context

        # Trinucleotide context (for SNPs only)
        if len(ref_allele) == 1 and len(alt_allele) == 1:
            tri_start = max(0, pos - 1)
            tri_end = min(len(ref_upper), pos + 2)

            if pos > 0 and pos < len(ref_upper) - 1:
                trinuc = ref_upper[pos - 1 : pos + 2]
                trinuc_mut = f"{trinuc}>{ref_upper[pos - 1]}{alt_allele}{ref_upper[pos + 1]}"
            elif pos == 0:
                trinuc = f"N{ref_allele}{ref_upper[pos + 1] if pos + 1 < len(ref_upper) else 'N'}"
                trinuc_mut = f"{trinuc}>N{alt_allele}{ref_upper[pos + 1] if pos + 1 < len(ref_upper) else 'N'}"
            else:
                trinuc = f"{ref_upper[pos - 1]}{ref_allele}N"
                trinuc_mut = f"{trinuc}>{ref_upper[pos - 1]}{alt_allele}N"

            result["trinucleotide_context"] = trinuc
            result["trinucleotide_mutation"] = trinuc_mut

            # SBS96 channel (standardize to pyrimidine context: C or T as reference)
            sbs_channel = _to_sbs96_channel(ref_allele, alt_allele, trinuc, full_complement)
            result["sbs_channel"] = sbs_channel
        else:
            result["trinucleotide_context"] = None
            result["trinucleotide_mutation"] = None
            result["sbs_channel"] = None

        annotated.append(result)

    logger.debug(f"Annotated context for {len(annotated)} variants (window={window})")
    return annotated


def _to_sbs96_channel(
    ref: str,
    alt: str,
    trinucleotide: str,
    complement: Dict[str, str],
) -> str:
    """Convert a mutation to the standard SBS96 pyrimidine context.

    The SBS96 framework represents all single-base substitutions in the
    context of a pyrimidine (C or T) reference base. If the reference is
    a purine (A or G), the reverse complement is used.

    Args:
        ref: Reference base (single character).
        alt: Alternate base (single character).
        trinucleotide: 3-bp trinucleotide context string.
        complement: Dict mapping bases to their complements.

    Returns:
        SBS96 channel string in format "X[R>A]Y" where R is C or T.
    """
    if len(trinucleotide) != 3:
        return f"?[{ref}>{alt}]?"

    if ref in ("C", "T"):
        # Already in pyrimidine context
        return f"{trinucleotide[0]}[{ref}>{alt}]{trinucleotide[2]}"
    else:
        # Reverse complement to get pyrimidine context
        rc_ref = complement.get(ref, "N")
        rc_alt = complement.get(alt, "N")
        rc_5prime = complement.get(trinucleotide[2], "N")
        rc_3prime = complement.get(trinucleotide[0], "N")
        return f"{rc_5prime}[{rc_ref}>{rc_alt}]{rc_3prime}"

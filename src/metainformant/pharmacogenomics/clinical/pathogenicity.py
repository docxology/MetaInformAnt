"""ACMG variant classification and ClinVar pathogenicity assessment.

Implements the ACMG/AMP (American College of Medical Genetics and Genomics /
Association for Molecular Pathology) 5-tier variant classification system
based on the 2015 Richards et al. guidelines. Evaluates individual evidence
criteria (PVS, PS, PM, PP, BA, BS, BP) and combines them into a final
classification: Pathogenic, Likely Pathogenic, VUS, Likely Benign, or Benign.

Also provides ClinVar data querying (from local parsed data) and gnomAD
population frequency checks for the BA1/BS1 benign criteria.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


class ACMGClassification(str, Enum):
    """ACMG 5-tier variant classification.

    Per Richards et al., 2015 (Genet Med 17:405-424).
    """

    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely Pathogenic"
    VUS = "Uncertain Significance"
    LIKELY_BENIGN = "Likely Benign"
    BENIGN = "Benign"


class ACMGCriteria(str, Enum):
    """Individual ACMG evidence criteria.

    Pathogenic criteria:
        PVS1: Null variant in gene where LOF is a known disease mechanism
        PS1-4: Strong pathogenic evidence
        PM1-6: Moderate pathogenic evidence
        PP1-5: Supporting pathogenic evidence

    Benign criteria:
        BA1: Stand-alone benign (allele frequency > 5%)
        BS1-4: Strong benign evidence
        BP1-7: Supporting benign evidence
    """

    # Very strong pathogenic
    PVS1 = "PVS1"

    # Strong pathogenic
    PS1 = "PS1"
    PS2 = "PS2"
    PS3 = "PS3"
    PS4 = "PS4"

    # Moderate pathogenic
    PM1 = "PM1"
    PM2 = "PM2"
    PM3 = "PM3"
    PM4 = "PM4"
    PM5 = "PM5"
    PM6 = "PM6"

    # Supporting pathogenic
    PP1 = "PP1"
    PP2 = "PP2"
    PP3 = "PP3"
    PP4 = "PP4"
    PP5 = "PP5"

    # Stand-alone benign
    BA1 = "BA1"

    # Strong benign
    BS1 = "BS1"
    BS2 = "BS2"
    BS3 = "BS3"
    BS4 = "BS4"

    # Supporting benign
    BP1 = "BP1"
    BP2 = "BP2"
    BP3 = "BP3"
    BP4 = "BP4"
    BP5 = "BP5"
    BP6 = "BP6"
    BP7 = "BP7"


# ── ACMG criteria descriptions ────────────────────────────────────────────────

_CRITERIA_DESCRIPTIONS: dict[str, str] = {
    "PVS1": "Null variant (nonsense, frameshift, canonical +/-1 or 2 splice sites, initiation codon, "
    "single or multi-exon deletion) in gene where LOF is a known mechanism of disease",
    "PS1": "Same amino acid change as a previously established pathogenic variant regardless of nucleotide change",
    "PS2": "De novo (both maternity and paternity confirmed) in patient with disease and no family history",
    "PS3": "Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene product",
    "PS4": "Prevalence of variant in affected individuals significantly increased compared to controls",
    "PM1": "Located in a mutational hot spot and/or critical/well-established functional domain",
    "PM2": "Absent from controls (or at extremely low frequency in gnomAD)",
    "PM3": "For recessive disorders, detected in trans with a pathogenic variant",
    "PM4": "Protein length changes due to in-frame deletions/insertions or stop-loss variants",
    "PM5": "Novel missense change at an amino acid residue where a different pathogenic missense has been seen before",
    "PM6": "Assumed de novo, but without confirmation of paternity and maternity",
    "PP1": "Cosegregation with disease in multiple affected family members",
    "PP2": "Missense in gene with low rate of benign missense and where missense are common mechanism of disease",
    "PP3": "Multiple lines of computational evidence support deleterious effect (CADD, REVEL, SIFT, PolyPhen)",
    "PP4": "Patient phenotype or family history highly specific for disease with single genetic etiology",
    "PP5": "Reputable source recently reports variant as pathogenic (may need independent evaluation)",
    "BA1": "Allele frequency > 5% in any general continental population dataset (gnomAD, ExAC, 1000G)",
    "BS1": "Allele frequency greater than expected for disorder (disease-specific threshold)",
    "BS2": "Observed in a healthy adult individual for a fully penetrant early-onset condition",
    "BS3": "Well-established in vitro or in vivo functional studies show no damaging effect",
    "BS4": "Lack of segregation in affected members of a family",
    "BP1": "Missense variant in gene for which primarily truncating variants are known to cause disease",
    "BP2": "Observed in trans with a pathogenic variant for a fully penetrant dominant disorder, "
    "or observed in cis with a pathogenic variant",
    "BP3": "In-frame deletions/insertions in a repetitive region without known function",
    "BP4": "Multiple lines of computational evidence suggest no impact on gene product",
    "BP5": "Variant found in a case with an alternate molecular basis for disease",
    "BP6": "Reputable source recently reports variant as benign (may need independent evaluation)",
    "BP7": "Synonymous variant with no predicted splice impact and no conservation",
}


# ── Built-in ClinVar-like data ─────────────────────────────────────────────────

_CLINVAR_DATA: dict[str, dict[str, Any]] = {
    "rs3892097": {
        "variant_id": "rs3892097",
        "gene": "CYP2D6",
        "hgvs": "NM_000106.6:c.506-1G>A",
        "classification": "Pathogenic",
        "review_status": "criteria provided, multiple submitters, no conflicts",
        "stars": 3,
        "condition": "CYP2D6-related pharmacogenomic condition",
        "last_evaluated": "2023-06",
        "submitter_count": 5,
    },
    "rs4244285": {
        "variant_id": "rs4244285",
        "gene": "CYP2C19",
        "hgvs": "NM_000769.4:c.681G>A",
        "classification": "Pathogenic",
        "review_status": "criteria provided, multiple submitters, no conflicts",
        "stars": 3,
        "condition": "CYP2C19-related pharmacogenomic condition",
        "last_evaluated": "2023-06",
        "submitter_count": 4,
    },
    "rs1799853": {
        "variant_id": "rs1799853",
        "gene": "CYP2C9",
        "hgvs": "NM_000771.4:c.430C>T",
        "classification": "Pathogenic/Likely pathogenic",
        "review_status": "criteria provided, multiple submitters",
        "stars": 2,
        "condition": "Warfarin sensitivity",
        "last_evaluated": "2023-04",
        "submitter_count": 3,
    },
    "rs1057910": {
        "variant_id": "rs1057910",
        "gene": "CYP2C9",
        "hgvs": "NM_000771.4:c.1075A>C",
        "classification": "Pathogenic",
        "review_status": "criteria provided, multiple submitters, no conflicts",
        "stars": 3,
        "condition": "Warfarin sensitivity",
        "last_evaluated": "2023-06",
        "submitter_count": 5,
    },
    "rs3918290": {
        "variant_id": "rs3918290",
        "gene": "DPYD",
        "hgvs": "NM_000110.4:c.1905+1G>A",
        "classification": "Pathogenic",
        "review_status": "criteria provided, multiple submitters, no conflicts",
        "stars": 3,
        "condition": "Dihydropyrimidine dehydrogenase deficiency",
        "last_evaluated": "2023-08",
        "submitter_count": 6,
    },
    "rs4149056": {
        "variant_id": "rs4149056",
        "gene": "SLCO1B1",
        "hgvs": "NM_006446.5:c.521T>C",
        "classification": "Pathogenic/Likely pathogenic",
        "review_status": "criteria provided, conflicting interpretations",
        "stars": 1,
        "condition": "Simvastatin-related myopathy",
        "last_evaluated": "2023-05",
        "submitter_count": 3,
    },
}


def classify_variant_acmg(
    variant: dict[str, Any],
    evidence: dict[str, bool] | None = None,
) -> dict[str, Any]:
    """Classify a variant using the ACMG 5-tier system.

    Takes variant information and either pre-evaluated evidence criteria or
    auto-evaluates criteria from the variant data, then aggregates into a
    final classification.

    Args:
        variant: Variant data dictionary with keys like:
            - "rsid": rsID (e.g., "rs3892097")
            - "gene": Gene symbol
            - "consequence": Variant consequence (e.g., "missense", "frameshift")
            - "allele_frequency": Population allele frequency dict or float
            - "computational_predictions": Dict of tool -> prediction
            - "functional_data": Functional study evidence
            - "segregation_data": Family segregation data
            - "de_novo": Whether variant is de novo
        evidence: Optional pre-evaluated criteria dict mapping ACMGCriteria name
            to True/False. If None, criteria are auto-evaluated from variant data.

    Returns:
        Dictionary with:
            - "classification": ACMGClassification value
            - "criteria_met": List of met ACMG criteria
            - "criteria_details": Dict of criterion -> description for met criteria
            - "score_summary": Count of criteria at each evidence level
            - "confidence": Assessment confidence ("high", "moderate", "low")
    """
    if evidence is None:
        evidence = apply_acmg_criteria(variant)

    # Collect met criteria
    criteria_met: list[str] = [k for k, v in evidence.items() if v]
    criteria_details: dict[str, str] = {c: _CRITERIA_DESCRIPTIONS.get(c, "No description") for c in criteria_met}

    # Aggregate into classification
    classification = aggregate_evidence(evidence)

    # Determine confidence
    met_count = len(criteria_met)
    if met_count >= 3:
        confidence = "high"
    elif met_count >= 2:
        confidence = "moderate"
    else:
        confidence = "low"

    # Count criteria by strength level
    score_summary = {
        "very_strong_pathogenic": sum(1 for c in criteria_met if c.startswith("PVS")),
        "strong_pathogenic": sum(1 for c in criteria_met if c.startswith("PS")),
        "moderate_pathogenic": sum(1 for c in criteria_met if c.startswith("PM")),
        "supporting_pathogenic": sum(1 for c in criteria_met if c.startswith("PP")),
        "standalone_benign": sum(1 for c in criteria_met if c.startswith("BA")),
        "strong_benign": sum(1 for c in criteria_met if c.startswith("BS")),
        "supporting_benign": sum(1 for c in criteria_met if c.startswith("BP")),
    }

    result = {
        "classification": classification,
        "criteria_met": sorted(criteria_met),
        "criteria_details": criteria_details,
        "score_summary": score_summary,
        "confidence": confidence,
    }

    logger.info(
        "ACMG classification for %s: %s (%d criteria met, confidence: %s)",
        variant.get("rsid", variant.get("hgvs", "unknown")),
        classification.value,
        met_count,
        confidence,
    )

    return result


def apply_acmg_criteria(
    variant_data: dict[str, Any],
) -> dict[str, bool]:
    """Evaluate individual ACMG criteria from variant data.

    Systematically evaluates each ACMG criterion based on available variant
    information. Each criterion is assessed independently.

    Args:
        variant_data: Variant information dictionary (see classify_variant_acmg)

    Returns:
        Dictionary mapping criterion name (str) to met/not-met (bool)
    """
    criteria: dict[str, bool] = {}

    consequence = variant_data.get("consequence", "").lower()
    gene = variant_data.get("gene", "").upper()

    # Get allele frequency (handle both dict and float)
    af_data = variant_data.get("allele_frequency", 0.0)
    if isinstance(af_data, dict):
        max_af = max(af_data.values()) if af_data else 0.0
    else:
        max_af = float(af_data) if af_data else 0.0

    comp_preds = variant_data.get("computational_predictions", {})
    functional = variant_data.get("functional_data", {})
    segregation = variant_data.get("segregation_data", {})

    # ── PVS1: Null variant in LOF-mechanism gene ──────────────────────────
    lof_consequences = {
        "frameshift",
        "nonsense",
        "stop_gained",
        "splice_donor",
        "splice_acceptor",
        "start_lost",
        "initiator_codon",
        "exon_deletion",
    }
    criteria["PVS1"] = consequence in lof_consequences

    # ── PS1: Same amino acid change as known pathogenic ───────────────────
    criteria["PS1"] = bool(variant_data.get("same_aa_pathogenic", False))

    # ── PS2: Confirmed de novo ────────────────────────────────────────────
    de_novo = variant_data.get("de_novo", False)
    paternity_confirmed = variant_data.get("paternity_confirmed", False)
    criteria["PS2"] = bool(de_novo and paternity_confirmed)

    # ── PS3: Well-established functional studies ──────────────────────────
    functional_result = functional.get("result", "")
    functional_quality = functional.get("quality", "")
    criteria["PS3"] = functional_result.lower() in (
        "damaging",
        "loss_of_function",
        "deleterious",
    ) and functional_quality.lower() in ("strong", "well-established")

    # ── PS4: Significantly increased prevalence in affected vs controls ───
    case_count = variant_data.get("case_count", 0)
    control_count = variant_data.get("control_count", 0)
    odds_ratio = variant_data.get("odds_ratio", 0.0)
    criteria["PS4"] = odds_ratio > 5.0 and case_count >= 5

    # ── PM1: Mutational hotspot or critical functional domain ─────────────
    criteria["PM1"] = bool(variant_data.get("in_hotspot", False) or variant_data.get("in_functional_domain", False))

    # ── PM2: Absent from controls ─────────────────────────────────────────
    pm2_threshold = variant_data.get("pm2_threshold", 0.0001)
    criteria["PM2"] = max_af < pm2_threshold

    # ── PM3: In trans with pathogenic (recessive) ─────────────────────────
    criteria["PM3"] = bool(variant_data.get("in_trans_pathogenic", False))

    # ── PM4: Protein length change (in-frame) ────────────────────────────
    inframe_consequences = {"inframe_deletion", "inframe_insertion", "stop_lost"}
    criteria["PM4"] = consequence in inframe_consequences

    # ── PM5: Novel missense at known pathogenic residue ──────────────────
    criteria["PM5"] = bool(variant_data.get("different_aa_pathogenic_residue", False))

    # ── PM6: Assumed de novo (unconfirmed) ────────────────────────────────
    criteria["PM6"] = bool(de_novo and not paternity_confirmed)

    # ── PP1: Cosegregation ────────────────────────────────────────────────
    seg_lod = segregation.get("lod_score", 0.0)
    criteria["PP1"] = seg_lod >= 1.5

    # ── PP2: Missense in constrained gene ─────────────────────────────────
    criteria["PP2"] = "missense" in consequence and bool(variant_data.get("gene_missense_constrained", False))

    # ── PP3: Computational evidence (deleterious) ─────────────────────────
    deleterious_count = 0
    total_preds = 0
    for tool, prediction in comp_preds.items():
        total_preds += 1
        pred_lower = str(prediction).lower()
        if pred_lower in ("damaging", "deleterious", "probably_damaging", "disease_causing", "pathogenic", "high"):
            deleterious_count += 1
    criteria["PP3"] = total_preds >= 2 and deleterious_count > total_preds / 2

    # ── PP4: Phenotype specificity ────────────────────────────────────────
    criteria["PP4"] = bool(variant_data.get("phenotype_specific", False))

    # ── PP5: Reputable source reports pathogenic ──────────────────────────
    criteria["PP5"] = bool(variant_data.get("reputable_source_pathogenic", False))

    # ── BA1: Stand-alone benign (AF > 5%) ─────────────────────────────────
    criteria["BA1"] = max_af > 0.05

    # ── BS1: Allele frequency greater than expected for disease ───────────
    disease_af_threshold = variant_data.get("disease_af_threshold", 0.01)
    criteria["BS1"] = max_af > disease_af_threshold and not criteria["BA1"]

    # ── BS2: Observed in healthy individual (penetrant condition) ─────────
    criteria["BS2"] = bool(variant_data.get("observed_healthy_adult", False))

    # ── BS3: Functional studies show no damaging effect ───────────────────
    criteria["BS3"] = functional_result.lower() in (
        "benign",
        "normal",
        "no_effect",
        "neutral",
    ) and functional_quality.lower() in ("strong", "well-established")

    # ── BS4: Lack of segregation ──────────────────────────────────────────
    criteria["BS4"] = bool(segregation.get("does_not_segregate", False))

    # ── BP1: Missense in truncating-variant disease gene ──────────────────
    criteria["BP1"] = "missense" in consequence and bool(variant_data.get("gene_truncating_mechanism", False))

    # ── BP2: Observed in trans/cis with pathogenic (dominant) ─────────────
    criteria["BP2"] = bool(variant_data.get("in_trans_dominant_pathogenic", False))

    # ── BP3: In-frame in repetitive region ────────────────────────────────
    criteria["BP3"] = consequence in inframe_consequences and bool(variant_data.get("in_repetitive_region", False))

    # ── BP4: Computational evidence (benign) ──────────────────────────────
    benign_count = 0
    for tool, prediction in comp_preds.items():
        pred_lower = str(prediction).lower()
        if pred_lower in ("benign", "tolerated", "neutral", "polymorphism", "low"):
            benign_count += 1
    criteria["BP4"] = total_preds >= 2 and benign_count > total_preds / 2

    # ── BP5: Alternate molecular basis found ──────────────────────────────
    criteria["BP5"] = bool(variant_data.get("alternate_molecular_basis", False))

    # ── BP6: Reputable source reports benign ──────────────────────────────
    criteria["BP6"] = bool(variant_data.get("reputable_source_benign", False))

    # ── BP7: Synonymous with no splice impact ─────────────────────────────
    criteria["BP7"] = "synonymous" in consequence and not bool(variant_data.get("splice_impact", False))

    met = [k for k, v in criteria.items() if v]
    logger.debug("ACMG criteria evaluation: %d met (%s)", len(met), ", ".join(sorted(met)))

    return criteria


def aggregate_evidence(
    criteria: dict[str, bool],
) -> ACMGClassification:
    """Combine ACMG criteria into a final 5-tier classification.

    Implements the combining rules from Richards et al. 2015 (Table 5):

    Pathogenic:
        (i)   1 Very Strong + >=1 Strong
        (ii)  1 Very Strong + >=2 Moderate
        (iii) 1 Very Strong + 1 Moderate + 1 Supporting
        (iv)  1 Very Strong + >=2 Supporting
        (v)   >=2 Strong
        (vi)  1 Strong + >=3 Moderate
        (vii) 1 Strong + 2 Moderate + >=2 Supporting
        (viii) 1 Strong + 1 Moderate + >=4 Supporting

    Likely Pathogenic:
        (i)   1 Very Strong + 1 Moderate
        (ii)  1 Strong + 1-2 Moderate
        (iii) 1 Strong + >=2 Supporting
        (iv)  >=3 Moderate
        (v)   2 Moderate + >=2 Supporting
        (vi)  1 Moderate + >=4 Supporting

    Benign:
        (i) 1 Stand-alone (BA1)
        (ii) >=2 Strong

    Likely Benign:
        (i) 1 Strong + 1 Supporting
        (ii) >=2 Supporting

    VUS: Everything else

    Args:
        criteria: Dictionary of criterion name -> True/False

    Returns:
        ACMGClassification enum value
    """
    met = {k for k, v in criteria.items() if v}

    # Count criteria by category
    pvs = sum(1 for c in met if c.startswith("PVS"))
    ps = sum(1 for c in met if c.startswith("PS"))
    pm = sum(1 for c in met if c.startswith("PM"))
    pp = sum(1 for c in met if c.startswith("PP"))
    ba = sum(1 for c in met if c.startswith("BA"))
    bs = sum(1 for c in met if c.startswith("BS"))
    bp = sum(1 for c in met if c.startswith("BP"))

    # ── Benign rules (check first - BA1 is standalone) ────────────────────
    if ba >= 1:
        return ACMGClassification.BENIGN
    if bs >= 2:
        return ACMGClassification.BENIGN

    # ── Likely Benign rules ───────────────────────────────────────────────
    if bs >= 1 and bp >= 1:
        return ACMGClassification.LIKELY_BENIGN
    if bp >= 2:
        return ACMGClassification.LIKELY_BENIGN

    # ── Pathogenic rules ──────────────────────────────────────────────────
    if pvs >= 1:
        if ps >= 1:
            return ACMGClassification.PATHOGENIC
        if pm >= 2:
            return ACMGClassification.PATHOGENIC
        if pm >= 1 and pp >= 1:
            return ACMGClassification.PATHOGENIC
        if pp >= 2:
            return ACMGClassification.PATHOGENIC

    if ps >= 2:
        return ACMGClassification.PATHOGENIC
    if ps >= 1:
        if pm >= 3:
            return ACMGClassification.PATHOGENIC
        if pm >= 2 and pp >= 2:
            return ACMGClassification.PATHOGENIC
        if pm >= 1 and pp >= 4:
            return ACMGClassification.PATHOGENIC

    # ── Likely Pathogenic rules ───────────────────────────────────────────
    if pvs >= 1 and pm >= 1:
        return ACMGClassification.LIKELY_PATHOGENIC

    if ps >= 1:
        if 1 <= pm <= 2:
            return ACMGClassification.LIKELY_PATHOGENIC
        if pp >= 2:
            return ACMGClassification.LIKELY_PATHOGENIC

    if pm >= 3:
        return ACMGClassification.LIKELY_PATHOGENIC
    if pm >= 2 and pp >= 2:
        return ACMGClassification.LIKELY_PATHOGENIC
    if pm >= 1 and pp >= 4:
        return ACMGClassification.LIKELY_PATHOGENIC

    # ── Default: VUS ──────────────────────────────────────────────────────
    return ACMGClassification.VUS


def query_clinvar(
    variant_id: str,
) -> dict[str, Any] | None:
    """Look up ClinVar classification for a variant.

    Queries the local ClinVar data cache for a variant's classification,
    review status, and associated conditions.

    Args:
        variant_id: Variant identifier (rsID, ClinVar variation ID, or HGVS)

    Returns:
        ClinVar record dictionary, or None if not found.
        Keys: variant_id, gene, hgvs, classification, review_status,
        stars, condition, last_evaluated, submitter_count
    """
    variant_normalized = variant_id.strip().lower()

    # Search by rsID
    for key, record in _CLINVAR_DATA.items():
        if key.lower() == variant_normalized:
            logger.info("ClinVar lookup for %s: %s", variant_id, record["classification"])
            return dict(record)

    # Search by HGVS
    for key, record in _CLINVAR_DATA.items():
        if record.get("hgvs", "").lower() == variant_normalized:
            logger.info("ClinVar lookup for %s (HGVS): %s", variant_id, record["classification"])
            return dict(record)

    logger.debug("No ClinVar record found for %s", variant_id)
    return None


def check_gnomad_frequency(
    variant: dict[str, Any],
    population: str | None = None,
    threshold: float = 0.05,
) -> dict[str, Any]:
    """Check variant population frequency against gnomAD thresholds.

    Evaluates whether a variant exceeds frequency thresholds relevant to
    ACMG benign criteria (BA1 at 5%, BS1 at disease-specific threshold).

    Args:
        variant: Variant dictionary with "allele_frequency" key containing
            either a float or a dict of population -> frequency
        population: Optional specific population to check (e.g., "European",
            "East_Asian"). If None, checks maximum across all populations.
        threshold: Frequency threshold for classification (default 0.05 for BA1)

    Returns:
        Dictionary with:
            - "max_frequency": Maximum allele frequency observed
            - "population": Population with highest frequency
            - "exceeds_threshold": Whether threshold is exceeded
            - "ba1_triggered": Whether BA1 benign criterion is met (AF > 5%)
            - "bs1_triggered": Whether BS1 benign criterion is met (AF > disease threshold)
            - "frequencies": All population frequencies if available
    """
    af_data = variant.get("allele_frequency", 0.0)

    if isinstance(af_data, dict):
        frequencies = {k: float(v) for k, v in af_data.items()}
    elif isinstance(af_data, (int, float)):
        frequencies = {"Global": float(af_data)}
    else:
        frequencies = {"Global": 0.0}

    if population:
        check_freq = frequencies.get(population, 0.0)
        max_pop = population
    else:
        if frequencies:
            max_pop = max(frequencies, key=lambda k: frequencies[k])
            check_freq = frequencies[max_pop]
        else:
            max_pop = "Unknown"
            check_freq = 0.0

    ba1 = check_freq > 0.05
    bs1 = check_freq > threshold and not ba1

    result = {
        "max_frequency": check_freq,
        "population": max_pop,
        "exceeds_threshold": check_freq > threshold,
        "ba1_triggered": ba1,
        "bs1_triggered": bs1,
        "frequencies": frequencies,
    }

    if ba1:
        logger.info(
            "BA1 triggered for variant (AF=%.4f in %s, threshold=%.4f)",
            check_freq,
            max_pop,
            threshold,
        )
    elif bs1:
        logger.info(
            "BS1 triggered for variant (AF=%.4f in %s, disease threshold=%.4f)",
            check_freq,
            max_pop,
            threshold,
        )

    return result

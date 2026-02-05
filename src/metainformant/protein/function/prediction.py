"""Protein function prediction utilities.

This module provides tools for predicting protein function from domain
composition, subcellular localization, solubility, physicochemical properties,
intrinsically disordered regions, active sites, and post-translational
modification sites. All predictions use pure Python implementations with
biologically grounded algorithms and scoring matrices.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Amino acid property tables
# ---------------------------------------------------------------------------

# Kyte-Doolittle hydrophobicity scale
_KD_HYDROPHOBICITY: Dict[str, float] = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
}

# Amino acid molecular weights (residue masses, Da)
_AA_MW: Dict[str, float] = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "Q": 128.1307,
    "E": 129.1155,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
    "U": 150.0388,
    "O": 237.3018,
}

# pKa values for isoelectric point calculation
_PKA_NTERM = 9.69
_PKA_CTERM = 2.34
_PKA_SIDE: Dict[str, float] = {
    "D": 3.65,
    "E": 4.25,
    "C": 8.18,
    "Y": 10.07,
    "H": 6.00,
    "K": 10.53,
    "R": 12.48,
}

# IUPred-like disorder propensity scores (simplified)
# Higher values = more disordered
_DISORDER_PROPENSITY: Dict[str, float] = {
    "A": 0.06,
    "R": 0.18,
    "N": 0.18,
    "D": 0.19,
    "C": -0.20,
    "Q": 0.20,
    "E": 0.30,
    "G": 0.17,
    "H": 0.05,
    "I": -0.49,
    "L": -0.34,
    "K": 0.26,
    "M": -0.19,
    "F": -0.42,
    "P": 0.41,
    "S": 0.14,
    "T": 0.04,
    "W": -0.49,
    "Y": -0.27,
    "V": -0.38,
}


# ---------------------------------------------------------------------------
# Domain-to-function mapping
# ---------------------------------------------------------------------------

_DOMAIN_FUNCTION_MAP: Dict[str, Dict[str, Any]] = {
    "Protein kinase domain": {
        "molecular_function": "protein kinase activity",
        "biological_process": "protein phosphorylation",
        "go_terms": ["GO:0004672", "GO:0006468"],
        "ec_number": "2.7.x.x",
    },
    "Trypsin-like serine protease": {
        "molecular_function": "serine-type endopeptidase activity",
        "biological_process": "proteolysis",
        "go_terms": ["GO:0004252", "GO:0006508"],
        "ec_number": "3.4.21.x",
    },
    "Ras family GTPase": {
        "molecular_function": "GTPase activity",
        "biological_process": "signal transduction",
        "go_terms": ["GO:0003924", "GO:0007264"],
        "ec_number": "3.6.5.x",
    },
    "RNA recognition motif (RRM)": {
        "molecular_function": "RNA binding",
        "biological_process": "mRNA processing",
        "go_terms": ["GO:0003723", "GO:0006397"],
        "ec_number": None,
    },
    "C2H2 zinc finger": {
        "molecular_function": "DNA binding",
        "biological_process": "transcription regulation",
        "go_terms": ["GO:0003677", "GO:0006355"],
        "ec_number": None,
    },
    "Homeodomain": {
        "molecular_function": "sequence-specific DNA binding",
        "biological_process": "embryonic development",
        "go_terms": ["GO:0043565", "GO:0009790"],
        "ec_number": None,
    },
    "WD40 repeat": {
        "molecular_function": "protein binding",
        "biological_process": "signal transduction",
        "go_terms": ["GO:0005515", "GO:0007165"],
        "ec_number": None,
    },
    "Helicase C-terminal domain": {
        "molecular_function": "ATP-dependent helicase activity",
        "biological_process": "RNA unwinding",
        "go_terms": ["GO:0004386", "GO:0008026"],
        "ec_number": "3.6.4.x",
    },
}


def predict_function_from_domains(domains: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Predict protein function from domain composition using rule-based mapping.

    Maps identified protein domains to Gene Ontology molecular functions,
    biological processes, and EC numbers using a curated domain-to-function
    knowledge base.

    Args:
        domains: List of domain hit dicts, each containing at minimum
            a ``name`` key with the domain name string.

    Returns:
        Dict with keys:
            - molecular_functions: list of predicted molecular functions
            - biological_processes: list of predicted biological processes
            - go_terms: list of associated GO term identifiers
            - ec_numbers: list of EC numbers (if enzymatic)
            - confidence: overall prediction confidence (0.0 to 1.0)
            - domain_contributions: dict mapping each domain to its
              functional annotation
            - is_enzyme: bool indicating predicted enzymatic activity
    """
    if not domains:
        return {
            "molecular_functions": [],
            "biological_processes": [],
            "go_terms": [],
            "ec_numbers": [],
            "confidence": 0.0,
            "domain_contributions": {},
            "is_enzyme": False,
        }

    mol_functions: List[str] = []
    bio_processes: List[str] = []
    go_terms: List[str] = []
    ec_numbers: List[str] = []
    contributions: Dict[str, Dict[str, Any]] = {}
    matched_count = 0

    for domain in domains:
        name = domain.get("name", "")
        if name in _DOMAIN_FUNCTION_MAP:
            func_info = _DOMAIN_FUNCTION_MAP[name]
            matched_count += 1

            mf = func_info["molecular_function"]
            bp = func_info["biological_process"]
            if mf not in mol_functions:
                mol_functions.append(mf)
            if bp not in bio_processes:
                bio_processes.append(bp)

            for gt in func_info.get("go_terms", []):
                if gt not in go_terms:
                    go_terms.append(gt)

            ec = func_info.get("ec_number")
            if ec and ec not in ec_numbers:
                ec_numbers.append(ec)

            contributions[name] = {
                "molecular_function": mf,
                "biological_process": bp,
                "go_terms": func_info.get("go_terms", []),
            }

    confidence = matched_count / len(domains) if domains else 0.0
    is_enzyme = len(ec_numbers) > 0

    logger.debug(
        "Function prediction: %d functions from %d/%d domains",
        len(mol_functions),
        matched_count,
        len(domains),
    )

    return {
        "molecular_functions": mol_functions,
        "biological_processes": bio_processes,
        "go_terms": go_terms,
        "ec_numbers": ec_numbers,
        "confidence": round(min(confidence, 1.0), 4),
        "domain_contributions": contributions,
        "is_enzyme": is_enzyme,
    }


# ---------------------------------------------------------------------------
# Subcellular localization prediction
# ---------------------------------------------------------------------------

# Localization signal patterns
_NLS_PATTERNS: List[Dict[str, Any]] = [
    {"name": "monopartite_NLS", "pattern": r"[KR]{4,6}", "weight": 1.0},
    {"name": "bipartite_NLS", "pattern": r"[KR]{2}.{9,12}[KR]{3}", "weight": 1.5},
    {"name": "SV40_NLS", "pattern": r"P[KR]{3,5}", "weight": 1.2},
]

_MTS_PATTERN = r"^M[LFIVMA]{1,3}[RK].{0,4}[LFIVMA]{2,4}[^DE]{3,8}[RK]"

_ER_RETENTION_PATTERNS = [
    r"KDEL$",  # ER retention signal (C-terminal)
    r"KKXX$",  # ER retrieval signal (C-terminal dilysine)
    r"[KR]DEL$",
]

_PEROXISOMAL_SIGNALS = [
    r"[SA]KL$",  # PTS1 (C-terminal)
    r"^M.{0,20}R[LI].{5}HL",  # PTS2 (N-terminal)
]


def predict_localization(sequence: str) -> Dict[str, Any]:
    """Predict subcellular localization from sequence features.

    Analyses the protein sequence for targeting signals including nuclear
    localization signals (NLS), mitochondrial targeting sequences (MTS),
    ER retention signals, peroxisomal targeting, and general sequence
    composition features. Amino acid composition is used as a secondary
    classifier.

    Args:
        sequence: Protein amino acid sequence string.

    Returns:
        Dict with keys:
            - predicted_location: most likely compartment string
            - scores: dict mapping compartment names to confidence scores
            - signals_found: list of detected targeting signals
            - composition_features: amino acid composition features used
            - confidence: overall prediction confidence (0.0 to 1.0)
    """
    if not sequence:
        return {
            "predicted_location": "unknown",
            "scores": {},
            "signals_found": [],
            "composition_features": {},
            "confidence": 0.0,
        }

    seq = sequence.upper()
    signals: List[Dict[str, Any]] = []
    scores: Dict[str, float] = {
        "nucleus": 0.0,
        "mitochondria": 0.0,
        "endoplasmic_reticulum": 0.0,
        "cytoplasm": 0.2,  # default baseline
        "secreted": 0.0,
        "membrane": 0.0,
        "peroxisome": 0.0,
    }

    # Check for nuclear localization signals
    for nls_def in _NLS_PATTERNS:
        for match in re.finditer(nls_def["pattern"], seq):
            signals.append(
                {
                    "type": "NLS",
                    "name": nls_def["name"],
                    "position": match.start(),
                    "sequence": match.group(),
                }
            )
            scores["nucleus"] += nls_def["weight"]

    # Check for mitochondrial targeting sequence (N-terminal)
    mts_match = re.search(_MTS_PATTERN, seq[:60])
    if mts_match:
        signals.append(
            {
                "type": "MTS",
                "name": "mitochondrial_targeting",
                "position": 0,
                "sequence": mts_match.group(),
            }
        )
        scores["mitochondria"] += 1.5

    # Additional MTS check: amphipathic helix at N-terminus
    if len(seq) >= 20:
        n_term_20 = seq[:20]
        basic_count = sum(1 for aa in n_term_20 if aa in ("R", "K"))
        acidic_count = sum(1 for aa in n_term_20 if aa in ("D", "E"))
        hydrophobic = sum(1 for aa in n_term_20 if aa in ("L", "I", "V", "F", "A"))
        if basic_count >= 3 and acidic_count == 0 and hydrophobic >= 5:
            scores["mitochondria"] += 0.8

    # Check for ER retention signals
    for pattern in _ER_RETENTION_PATTERNS:
        if re.search(pattern, seq):
            signals.append(
                {
                    "type": "ER_retention",
                    "name": "ER_retention_signal",
                    "position": len(seq) - 4,
                    "sequence": seq[-4:],
                }
            )
            scores["endoplasmic_reticulum"] += 1.5
            break

    # Check for peroxisomal targeting signals
    for pattern in _PEROXISOMAL_SIGNALS:
        match = re.search(pattern, seq)
        if match:
            signals.append(
                {
                    "type": "PTS",
                    "name": "peroxisomal_targeting",
                    "position": match.start(),
                    "sequence": match.group(),
                }
            )
            scores["peroxisome"] += 1.5

    # Signal peptide heuristic (N-terminal hydrophobic)
    if len(seq) >= 25:
        n_term = seq[:25]
        hydro_score = sum(_KD_HYDROPHOBICITY.get(aa, 0.0) for aa in n_term[5:20]) / 15
        if hydro_score > 1.0:
            scores["secreted"] += 0.8
            signals.append(
                {
                    "type": "signal_peptide",
                    "name": "secretory_signal",
                    "position": 0,
                    "sequence": n_term,
                }
            )

    # Transmembrane heuristic
    for i in range(0, max(1, len(seq) - 19)):
        window = seq[i : i + 19]
        if len(window) < 19:
            break
        avg_hydro = sum(_KD_HYDROPHOBICITY.get(aa, 0.0) for aa in window) / 19
        if avg_hydro > 1.6:
            scores["membrane"] += 0.5
            break

    # Composition-based features
    comp = _compute_composition(seq)
    composition_features = {
        "basic_fraction": comp.get("K", 0.0) + comp.get("R", 0.0),
        "acidic_fraction": comp.get("D", 0.0) + comp.get("E", 0.0),
        "hydrophobic_fraction": sum(comp.get(aa, 0.0) for aa in "IVLFMA"),
        "charged_fraction": sum(comp.get(aa, 0.0) for aa in "KRDEH"),
    }

    # Composition adjustments
    if composition_features["basic_fraction"] > 0.15:
        scores["nucleus"] += 0.3
    if composition_features["hydrophobic_fraction"] > 0.45:
        scores["membrane"] += 0.3

    # Normalize scores to probabilities
    total = sum(scores.values())
    if total > 0:
        scores = {k: round(v / total, 4) for k, v in scores.items()}

    predicted = max(scores, key=scores.get)  # type: ignore[arg-type]
    confidence = scores[predicted]

    logger.debug(
        "Localization prediction: %s (confidence=%.2f, %d signals)",
        predicted,
        confidence,
        len(signals),
    )

    return {
        "predicted_location": predicted,
        "scores": scores,
        "signals_found": signals,
        "composition_features": composition_features,
        "confidence": round(confidence, 4),
    }


def _compute_composition(seq: str) -> Dict[str, float]:
    """Compute fractional amino acid composition.

    Args:
        seq: Amino acid sequence.

    Returns:
        Dict mapping single-letter AA codes to their fraction (0.0-1.0).
    """
    if not seq:
        return {}
    counts = Counter(seq)
    total = len(seq)
    return {aa: count / total for aa, count in counts.items()}


# ---------------------------------------------------------------------------
# Solubility prediction
# ---------------------------------------------------------------------------


def predict_solubility(sequence: str) -> Dict[str, Any]:
    """Predict protein solubility from sequence features.

    Uses a combination of charge, hydrophobicity, molecular weight,
    and amino acid composition to estimate whether a recombinantly
    expressed protein is likely to be soluble. The model is based on
    the Wilkinson-Harrison solubility model and the charge-hydrophobicity
    relationship.

    Args:
        sequence: Protein amino acid sequence string.

    Returns:
        Dict with keys:
            - is_soluble: bool prediction
            - solubility_score: numerical score (higher = more soluble)
            - contributing_factors: dict with individual feature scores
            - charge_at_ph7: net charge at pH 7.0
            - mean_hydrophobicity: average Kyte-Doolittle score
            - molecular_weight: MW in Da
            - recommendations: list of suggestions to improve solubility
    """
    if not sequence:
        return {
            "is_soluble": False,
            "solubility_score": 0.0,
            "contributing_factors": {},
            "charge_at_ph7": 0.0,
            "mean_hydrophobicity": 0.0,
            "molecular_weight": 0.0,
            "recommendations": [],
        }

    seq = sequence.upper()

    # Compute features
    charge = _charge_at_ph(seq, 7.0)
    abs_charge = abs(charge)
    mean_hydro = sum(_KD_HYDROPHOBICITY.get(aa, 0.0) for aa in seq) / len(seq)
    mw = sum(_AA_MW.get(aa, 110.0) for aa in seq) + 18.015  # +water

    # Amino acid composition features
    comp = _compute_composition(seq)
    fraction_charged = sum(comp.get(aa, 0.0) for aa in "KRDEH")
    fraction_aromatic = sum(comp.get(aa, 0.0) for aa in "FWY")
    fraction_polar = sum(comp.get(aa, 0.0) for aa in "STNQ")
    fraction_cys = comp.get("C", 0.0)
    fraction_pro = comp.get("P", 0.0)

    # Wilkinson-Harrison-like scoring
    factors: Dict[str, float] = {}

    # Charge contributes positively to solubility
    charge_score = min(abs_charge / len(seq) * 50, 1.0)
    factors["charge"] = round(charge_score, 4)

    # High hydrophobicity hurts solubility
    hydro_score = max(0.0, 1.0 - (mean_hydro + 0.5) / 3.0)
    factors["hydrophobicity"] = round(hydro_score, 4)

    # Charged residues help
    charged_score = min(fraction_charged * 4.0, 1.0)
    factors["charged_residues"] = round(charged_score, 4)

    # Aromatic residues can reduce solubility at high fractions
    aromatic_score = max(0.0, 1.0 - fraction_aromatic * 5.0)
    factors["aromatic_content"] = round(aromatic_score, 4)

    # Cysteine content (disulfide bonds can either help or hurt)
    cys_score = 0.5 if fraction_cys < 0.03 else max(0.0, 1.0 - fraction_cys * 10.0)
    factors["cysteine_content"] = round(cys_score, 4)

    # Molecular weight penalty (very large proteins are harder to express)
    if mw < 30000:
        mw_score = 0.8
    elif mw < 60000:
        mw_score = 0.6
    elif mw < 100000:
        mw_score = 0.4
    else:
        mw_score = 0.2
    factors["molecular_weight"] = round(mw_score, 4)

    # Proline content (helps prevent aggregation)
    pro_score = min(fraction_pro * 10.0, 0.8)
    factors["proline_content"] = round(pro_score, 4)

    # Polar residue content
    polar_score = min(fraction_polar * 3.0, 1.0)
    factors["polar_residues"] = round(polar_score, 4)

    # Overall solubility score (weighted average)
    weights = {
        "charge": 0.15,
        "hydrophobicity": 0.25,
        "charged_residues": 0.15,
        "aromatic_content": 0.10,
        "cysteine_content": 0.05,
        "molecular_weight": 0.10,
        "proline_content": 0.05,
        "polar_residues": 0.15,
    }

    solubility_score = sum(factors[k] * weights[k] for k in weights)
    is_soluble = solubility_score >= 0.45

    # Generate recommendations
    recommendations: List[str] = []
    if mean_hydro > 0.0:
        recommendations.append("High average hydrophobicity; consider fusion with solubility tag (MBP, SUMO)")
    if fraction_charged < 0.15:
        recommendations.append("Low charged residue content; surface charge engineering may improve solubility")
    if mw > 80000:
        recommendations.append("Large protein; consider expressing as domains or truncations")
    if fraction_cys > 0.05:
        recommendations.append("High cysteine content; consider reducing conditions or cysteine mutations")
    if fraction_aromatic > 0.15:
        recommendations.append("High aromatic content; consider lower expression temperature")
    if not recommendations:
        recommendations.append("Sequence features are favorable for soluble expression")

    return {
        "is_soluble": is_soluble,
        "solubility_score": round(solubility_score, 4),
        "contributing_factors": factors,
        "charge_at_ph7": round(charge, 4),
        "mean_hydrophobicity": round(mean_hydro, 4),
        "molecular_weight": round(mw, 2),
        "recommendations": recommendations,
    }


# ---------------------------------------------------------------------------
# Comprehensive physicochemical properties
# ---------------------------------------------------------------------------


def compute_physicochemical(sequence: str) -> Dict[str, Any]:
    """Compute comprehensive physicochemical properties of a protein.

    Calculates molecular weight, isoelectric point, molar extinction
    coefficient, instability index, aliphatic index, GRAVY score, and
    net charge at pH 7. All calculations use established bioinformatics
    formulas.

    Args:
        sequence: Protein amino acid sequence string.

    Returns:
        Dict with keys:
            - molecular_weight: MW in Daltons
            - isoelectric_point: pI value
            - extinction_coefficient: dict with reduced and oxidized values
            - instability_index: Guruprasad et al. instability index
            - aliphatic_index: Ikai aliphatic index
            - gravy: Grand Average of Hydropathy
            - charge_at_ph7: net charge at pH 7.0
            - length: sequence length
            - amino_acid_composition: dict of AA fractions
            - formula_weight: formula weight including water

    Raises:
        ValueError: If sequence is empty.
    """
    if not sequence:
        raise ValueError("Sequence must be non-empty")

    seq = sequence.upper()
    n = len(seq)
    comp = Counter(seq)

    # Molecular weight
    mw = sum(_AA_MW.get(aa, 110.0) * count for aa, count in comp.items())
    mw += 18.015  # water molecule

    # Isoelectric point (bisection method)
    pi = _compute_isoelectric_point(seq)

    # Extinction coefficient at 280 nm (Pace method)
    n_trp = comp.get("W", 0)
    n_tyr = comp.get("Y", 0)
    n_cys = comp.get("C", 0)
    ext_reduced = n_trp * 5500 + n_tyr * 1490
    ext_oxidized = ext_reduced + (n_cys // 2) * 125

    # Instability index (Guruprasad et al. 1990)
    instability = _compute_instability_index(seq)

    # Aliphatic index (Ikai 1980)
    frac_a = comp.get("A", 0) / n
    frac_v = comp.get("V", 0) / n
    frac_i = comp.get("I", 0) / n
    frac_l = comp.get("L", 0) / n
    aliphatic = 100 * (frac_a + 2.9 * frac_v + 3.9 * (frac_i + frac_l))

    # GRAVY
    gravy_score = sum(_KD_HYDROPHOBICITY.get(aa, 0.0) for aa in seq) / n

    # Charge at pH 7
    charge = _charge_at_ph(seq, 7.0)

    # Amino acid composition (fractions)
    aa_comp = {aa: count / n for aa, count in comp.items()}

    return {
        "molecular_weight": round(mw, 2),
        "isoelectric_point": round(pi, 2),
        "extinction_coefficient": {
            "reduced": float(ext_reduced),
            "oxidized": float(ext_oxidized),
            "n_trp": n_trp,
            "n_tyr": n_tyr,
            "n_cys": n_cys,
        },
        "instability_index": round(instability, 2),
        "is_stable": instability < 40.0,
        "aliphatic_index": round(aliphatic, 2),
        "gravy": round(gravy_score, 4),
        "charge_at_ph7": round(charge, 4),
        "length": n,
        "amino_acid_composition": {k: round(v, 4) for k, v in sorted(aa_comp.items())},
        "formula_weight": round(mw, 2),
    }


def _compute_isoelectric_point(seq: str) -> float:
    """Compute isoelectric point using bisection method.

    Args:
        seq: Uppercase amino acid sequence.

    Returns:
        Isoelectric point (pI).
    """
    low, high = 0.0, 14.0
    for _ in range(100):
        mid = (low + high) / 2.0
        charge = _charge_at_ph(seq, mid)
        if charge > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def _charge_at_ph(seq: str, ph: float) -> float:
    """Calculate net charge at a given pH.

    Args:
        seq: Uppercase amino acid sequence.
        ph: pH value.

    Returns:
        Net charge.
    """
    charge = 0.0
    # N-terminus (positive)
    charge += 1.0 / (1.0 + 10 ** (ph - _PKA_NTERM))
    # C-terminus (negative)
    charge -= 1.0 / (1.0 + 10 ** (_PKA_CTERM - ph))
    # Side chains
    for aa in seq:
        if aa in ("D", "E", "C", "Y"):
            charge -= 1.0 / (1.0 + 10 ** (_PKA_SIDE[aa] - ph))
        elif aa in ("H", "K", "R"):
            charge += 1.0 / (1.0 + 10 ** (ph - _PKA_SIDE[aa]))
    return charge


def _compute_instability_index(seq: str) -> float:
    """Compute instability index (Guruprasad et al. 1990).

    Args:
        seq: Uppercase amino acid sequence.

    Returns:
        Instability index value. Values > 40 indicate likely unstable in vitro.
    """
    if len(seq) < 2:
        return 0.0

    # DIWV weight values for common dipeptides
    diwv: Dict[str, float] = {
        "WW": 1.0,
        "WC": 1.0,
        "WM": 24.68,
        "WH": 24.68,
        "CW": -14.03,
        "CH": 33.60,
        "CC": 1.0,
        "CF": -14.03,
        "FY": 33.60,
        "FF": 1.0,
        "FK": -14.03,
        "FL": 1.0,
        "GG": 13.34,
        "GE": -14.03,
        "GA": -7.49,
        "GR": 1.0,
        "EE": 33.60,
        "ED": -14.03,
        "EK": 1.0,
        "EA": 11.0,
        "KK": 1.0,
        "KD": 1.0,
        "KE": 33.60,
        "KR": 33.60,
        "DD": 1.0,
        "DG": 1.0,
        "DE": 1.0,
        "DR": -6.54,
        "AA": 1.0,
        "AE": 1.0,
        "AG": 1.0,
        "AV": 1.0,
        "RR": 58.28,
        "RK": 44.94,
        "RL": 1.0,
        "RA": 1.0,
        "VV": 1.0,
        "VA": 1.0,
        "VL": 1.0,
        "VK": -7.49,
        "LL": 1.0,
        "LA": 1.0,
        "LK": -7.49,
        "LR": 1.0,
        "SS": 1.0,
        "SA": 1.0,
        "SG": 1.0,
        "SR": 44.94,
        "PP": 20.26,
        "PG": 1.0,
        "PA": 20.26,
        "PV": 20.26,
        "II": 1.0,
        "IA": 1.0,
        "IL": 20.26,
        "IV": -7.49,
        "TT": 1.0,
        "TA": 1.0,
        "TG": -7.49,
        "TR": 1.0,
        "YY": 33.60,
        "YF": 33.60,
        "YW": -9.37,
        "YA": 24.68,
        "HH": 1.0,
        "HW": -1.88,
        "HR": 1.0,
        "HA": 1.0,
        "QQ": 20.26,
        "QR": 1.0,
        "QE": 33.60,
        "QK": 1.0,
        "NN": 1.0,
        "NR": 1.0,
        "NG": -14.03,
        "NK": 24.68,
        "MM": -1.88,
        "MK": 1.0,
        "ML": 1.0,
        "MA": 1.0,
    }

    score = 0.0
    for i in range(len(seq) - 1):
        dipeptide = seq[i : i + 2]
        score += diwv.get(dipeptide, 1.0)

    return (10.0 / len(seq)) * score


# ---------------------------------------------------------------------------
# Intrinsically disordered region prediction
# ---------------------------------------------------------------------------


def predict_disordered_regions(
    sequence: str,
    method: str = "iupred_like",
    window: int = 21,
    threshold: float = 0.5,
) -> List[Dict[str, Any]]:
    """Predict intrinsically disordered regions using amino acid propensity.

    Implements a simplified IUPred-like algorithm that scores each residue
    for disorder propensity based on its amino acid type and local sequence
    context. Contiguous regions above the threshold are reported as
    intrinsically disordered regions (IDRs).

    Args:
        sequence: Protein amino acid sequence string.
        method: Prediction method. Currently supports ``"iupred_like"``
            (amino acid propensity with sliding window).
        window: Sliding window size for smoothing (must be odd).
        threshold: Minimum score to classify a region as disordered.

    Returns:
        List of dicts, each with:
            - start: 0-based start position
            - end: 0-based end position (exclusive)
            - length: region length
            - mean_score: average disorder score in the region
            - max_score: maximum disorder score in the region
            - sequence: subsequence of the disordered region

    Raises:
        ValueError: If method is not supported.
    """
    if method not in ("iupred_like",):
        raise ValueError(f"Unsupported disorder prediction method: {method}")

    if not sequence or len(sequence) < window:
        return []

    seq = sequence.upper()
    n = len(seq)

    # Compute raw per-residue disorder propensity
    raw_scores = [_DISORDER_PROPENSITY.get(aa, 0.0) for aa in seq]

    # Smooth with sliding window
    half_w = window // 2
    smoothed: List[float] = []
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        window_vals = raw_scores[start:end]
        smoothed.append(sum(window_vals) / len(window_vals))

    # Normalize to 0-1 range using sigmoid-like transformation
    normalized: List[float] = []
    for s in smoothed:
        # Map propensity scores to probability-like values
        prob = 1.0 / (1.0 + math.exp(-5.0 * s))
        normalized.append(prob)

    # Find contiguous regions above threshold
    regions: List[Dict[str, Any]] = []
    i = 0
    while i < n:
        if normalized[i] >= threshold:
            region_start = i
            region_scores: List[float] = []
            while i < n and normalized[i] >= threshold:
                region_scores.append(normalized[i])
                i += 1
            region_end = i

            # Require minimum length of 5 residues
            if region_end - region_start >= 5:
                regions.append(
                    {
                        "start": region_start,
                        "end": region_end,
                        "length": region_end - region_start,
                        "mean_score": round(sum(region_scores) / len(region_scores), 4),
                        "max_score": round(max(region_scores), 4),
                        "sequence": seq[region_start:region_end],
                    }
                )
        else:
            i += 1

    logger.debug("Predicted %d disordered regions", len(regions))
    return regions


# ---------------------------------------------------------------------------
# Active site prediction
# ---------------------------------------------------------------------------

# Known catalytic motifs within specific domain types
_ACTIVE_SITE_PATTERNS: List[Dict[str, Any]] = [
    {
        "name": "Serine protease catalytic triad",
        "domains": {"Trypsin-like serine protease"},
        "patterns": [
            {"residue": "H", "pattern": r"[GSAC]H[ILVFYWM]", "role": "catalytic histidine"},
            {"residue": "D", "pattern": r"D[ILVFA][ASTG]", "role": "catalytic aspartate"},
            {"residue": "S", "pattern": r"G[DNSG]SG[GAS]", "role": "catalytic serine"},
        ],
    },
    {
        "name": "Protein kinase active site",
        "domains": {"Protein kinase domain"},
        "patterns": [
            {"residue": "K", "pattern": r"V[AI]K", "role": "ATP binding lysine"},
            {"residue": "D", "pattern": r"[HY]RD[ILVMA]", "role": "catalytic aspartate"},
            {"residue": "D", "pattern": r"DFG", "role": "DFG motif aspartate"},
        ],
    },
    {
        "name": "GTPase catalytic site",
        "domains": {"Ras family GTPase"},
        "patterns": [
            {"residue": "T", "pattern": r"DTAGQ", "role": "catalytic threonine"},
            {"residue": "G", "pattern": r"G.{4}GK[ST]", "role": "P-loop glycine"},
        ],
    },
    {
        "name": "Metalloprotease zinc binding",
        "domains": set(),  # generic, no domain restriction
        "patterns": [
            {"residue": "H", "pattern": r"HE..H", "role": "zinc binding histidine"},
            {"residue": "E", "pattern": r"HE..H", "role": "catalytic glutamate"},
        ],
    },
    {
        "name": "Cysteine protease active site",
        "domains": set(),
        "patterns": [
            {"residue": "C", "pattern": r"Q.C[GSAW]", "role": "catalytic cysteine"},
            {"residue": "H", "pattern": r"[ILVFA]H[AG]", "role": "catalytic histidine"},
        ],
    },
]


def find_active_sites(
    sequence: str,
    domains: List[Dict[str, Any]] | None = None,
) -> List[Dict[str, Any]]:
    """Identify potential catalytic residues from conserved motif patterns.

    Scans the protein sequence for known catalytic site motifs. When domain
    annotations are provided, the search is focused to domain-specific
    patterns for higher specificity. Without domain information, all generic
    catalytic motif patterns are tested.

    Args:
        sequence: Protein amino acid sequence string.
        domains: Optional list of domain hit dicts with ``name`` keys.
            Used to restrict pattern search to domain-appropriate motifs.

    Returns:
        List of dicts, each with:
            - site_name: name of the catalytic site type
            - residue: the catalytic amino acid
            - position: 0-based position in the sequence
            - role: functional role of the residue
            - motif_matched: the regex pattern that matched
            - context: surrounding sequence context (5 residues each side)
            - confidence: confidence score (0.0 to 1.0)
    """
    if not sequence:
        return []

    seq = sequence.upper()
    domain_names: set[str] = set()
    if domains:
        domain_names = {d.get("name", "") for d in domains}

    results: List[Dict[str, Any]] = []

    for site_def in _ACTIVE_SITE_PATTERNS:
        # Check if this pattern applies (domain match or generic)
        required_domains = site_def.get("domains", set())
        if required_domains and not (required_domains & domain_names):
            # No matching domain, but still check if no domains provided
            if domains is not None:
                continue

        # Confidence is higher when domain context matches
        base_confidence = 0.8 if (required_domains & domain_names) else 0.4
        if not required_domains:
            base_confidence = 0.5

        for pattern_def in site_def["patterns"]:
            target_residue = pattern_def["residue"]
            pattern = pattern_def["pattern"]
            role = pattern_def["role"]

            try:
                for match in re.finditer(pattern, seq):
                    # Find the target residue within the match
                    matched_str = match.group()
                    for offset, aa in enumerate(matched_str):
                        if aa == target_residue:
                            abs_pos = match.start() + offset
                            # Get context
                            ctx_start = max(0, abs_pos - 5)
                            ctx_end = min(len(seq), abs_pos + 6)
                            context = seq[ctx_start:ctx_end]

                            results.append(
                                {
                                    "site_name": site_def["name"],
                                    "residue": target_residue,
                                    "position": abs_pos,
                                    "role": role,
                                    "motif_matched": pattern,
                                    "context": context,
                                    "confidence": round(base_confidence, 4),
                                }
                            )
                            break  # only first target residue per match

            except re.error as exc:
                logger.warning("Regex error in active site pattern: %s", exc)

    # Deduplicate by position
    seen_positions: set[int] = set()
    deduped: List[Dict[str, Any]] = []
    for result in sorted(results, key=lambda r: (-r["confidence"], r["position"])):
        if result["position"] not in seen_positions:
            deduped.append(result)
            seen_positions.add(result["position"])

    deduped.sort(key=lambda r: r["position"])

    logger.debug("Found %d potential active site residues", len(deduped))
    return deduped


# ---------------------------------------------------------------------------
# Post-translational modification site prediction
# ---------------------------------------------------------------------------

# PTM motif definitions
_PTM_PATTERNS: List[Dict[str, Any]] = [
    # Phosphorylation
    {
        "type": "phosphorylation",
        "subtype": "Ser/Thr kinase (basophilic)",
        "pattern": r"[RK]{2,3}.{0,2}[ST]",
        "target_residue": "ST",
        "description": "Basophilic kinase substrate motif (PKA/PKC-like)",
    },
    {
        "type": "phosphorylation",
        "subtype": "Ser/Thr kinase (proline-directed)",
        "pattern": r"[ST]P",
        "target_residue": "ST",
        "description": "Proline-directed kinase substrate (CDK/MAPK-like)",
    },
    {
        "type": "phosphorylation",
        "subtype": "Tyr kinase",
        "pattern": r"[EDQN].{0,2}Y.{0,2}[ILVFM]",
        "target_residue": "Y",
        "description": "Tyrosine kinase substrate motif",
    },
    {
        "type": "phosphorylation",
        "subtype": "CK2 substrate",
        "pattern": r"[ST].{2}[DE]",
        "target_residue": "ST",
        "description": "Casein kinase 2 substrate motif",
    },
    # N-linked glycosylation
    {
        "type": "glycosylation",
        "subtype": "N-linked",
        "pattern": r"N[^P][ST][^P]",
        "target_residue": "N",
        "description": "N-X-S/T sequon for N-linked glycosylation",
    },
    # O-linked glycosylation (simplified)
    {
        "type": "glycosylation",
        "subtype": "O-linked",
        "pattern": r"[ST]P.P",
        "target_residue": "ST",
        "description": "O-linked glycosylation motif (mucin-type)",
    },
    # Ubiquitination
    {
        "type": "ubiquitination",
        "subtype": "ubiquitin conjugation",
        "pattern": r"[IVLMA].K.[DERED]",
        "target_residue": "K",
        "description": "Ubiquitination site (lysine in degradation context)",
    },
    # SUMOylation
    {
        "type": "sumoylation",
        "subtype": "SUMO conjugation",
        "pattern": r"[IVLMF]K.E",
        "target_residue": "K",
        "description": "SUMOylation consensus motif (psi-K-X-E)",
    },
    # Acetylation
    {
        "type": "acetylation",
        "subtype": "lysine acetylation",
        "pattern": r"[GAS]K[RKPS]",
        "target_residue": "K",
        "description": "Lysine acetylation motif",
    },
    # Myristoylation (N-terminal)
    {
        "type": "myristoylation",
        "subtype": "N-myristoylation",
        "pattern": r"^MG[^EDRKHPFYW].[STAGCN][^P]",
        "target_residue": "G",
        "description": "N-myristoylation motif (position 2 glycine)",
    },
    # GPI anchor signal
    {
        "type": "gpi_anchor",
        "subtype": "GPI attachment",
        "pattern": r"[GANDSTC][GASV][GANDSTC].{4,12}[LIVMFWG]{5,}$",
        "target_residue": "GANDSTC",
        "description": "GPI anchor signal (C-terminal)",
    },
]


def predict_post_translational_mods(sequence: str) -> List[Dict[str, Any]]:
    """Predict post-translational modification sites from sequence motifs.

    Scans for phosphorylation sites (Ser/Thr/Tyr with kinase-specific
    motifs), N-linked glycosylation (N-X-S/T sequon), O-linked
    glycosylation, ubiquitination, SUMOylation, acetylation,
    myristoylation, and GPI anchor signals.

    Args:
        sequence: Protein amino acid sequence string.

    Returns:
        List of dicts sorted by position, each with:
            - type: PTM type (phosphorylation, glycosylation, etc.)
            - subtype: specific PTM subtype
            - position: 0-based position of the modified residue
            - residue: the modified amino acid
            - motif: the matched motif sequence
            - description: human-readable description
            - confidence: prediction confidence (0.0 to 1.0)
    """
    if not sequence:
        return []

    seq = sequence.upper()
    results: List[Dict[str, Any]] = []

    for ptm_def in _PTM_PATTERNS:
        try:
            for match in re.finditer(ptm_def["pattern"], seq):
                matched_str = match.group()
                target_chars = set(ptm_def["target_residue"])

                # Find the target residue within the match
                for offset, aa in enumerate(matched_str):
                    if aa in target_chars:
                        abs_pos = match.start() + offset

                        # Confidence based on motif specificity
                        confidence = 0.6
                        motif_len = len(matched_str)
                        if motif_len >= 5:
                            confidence = 0.75
                        if motif_len >= 7:
                            confidence = 0.85

                        # Surface accessibility heuristic: charged/polar
                        # neighbors increase confidence
                        window_start = max(0, abs_pos - 3)
                        window_end = min(len(seq), abs_pos + 4)
                        local = seq[window_start:window_end]
                        polar_count = sum(1 for c in local if c in "STNQKRDEH")
                        if polar_count >= 3:
                            confidence = min(confidence + 0.1, 1.0)

                        results.append(
                            {
                                "type": ptm_def["type"],
                                "subtype": ptm_def["subtype"],
                                "position": abs_pos,
                                "residue": aa,
                                "motif": matched_str,
                                "description": ptm_def["description"],
                                "confidence": round(confidence, 4),
                            }
                        )
                        break  # only first target per match

        except re.error as exc:
            logger.warning("Regex error in PTM pattern: %s", exc)

    # Sort by position and remove duplicate positions for same PTM type
    results.sort(key=lambda r: (r["position"], r["type"]))

    seen: set[Tuple[int, str]] = set()
    deduped: List[Dict[str, Any]] = []
    for r in results:
        key = (r["position"], r["type"])
        if key not in seen:
            deduped.append(r)
            seen.add(key)

    logger.debug("Predicted %d PTM sites", len(deduped))
    return deduped

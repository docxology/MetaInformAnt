"""Protein sequence processing and analysis utilities.

This module provides core functionality for reading, writing, and analyzing
protein sequences in FASTA format.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

from metainformant.core import io
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Canonical Kyte-Doolittle hydropathy scale (Kyte & Doolittle 1982).
# Single source of truth for the protein domain; do not re-declare local copies.
KD_HYDROPHOBICITY: Dict[str, float] = {
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

# Dipeptide instability weight values (DIWV) of Guruprasad, Reddy & Pandit, Protein Eng 4:155-161 (1990): the weight
# of every ordered pair of the twenty amino acids ("AC" is an Ala followed by a Cys), all 400 of them; 247 are the
# neutral weight 1.0. The values were checked against the table of Biopython (Bio.SeqUtils.ProtParamData.DIWV).
# fmt: off
DIWV: Dict[str, float] = {
    "AA": 1.0, "AC": 44.94, "AD": -7.49, "AE": 1.0, "AF": 1.0, "AG": 1.0, "AH": -7.49,
    "AI": 1.0, "AK": 1.0, "AL": 1.0, "AM": 1.0, "AN": 1.0, "AP": 20.26, "AQ": 1.0,
    "AR": 1.0, "AS": 1.0, "AT": 1.0, "AV": 1.0, "AW": 1.0, "AY": 1.0, "CA": 1.0,
    "CC": 1.0, "CD": 20.26, "CE": 1.0, "CF": 1.0, "CG": 1.0, "CH": 33.6, "CI": 1.0,
    "CK": 1.0, "CL": 20.26, "CM": 33.6, "CN": 1.0, "CP": 20.26, "CQ": -6.54, "CR": 1.0,
    "CS": 1.0, "CT": 33.6, "CV": -6.54, "CW": 24.68, "CY": 1.0, "DA": 1.0, "DC": 1.0,
    "DD": 1.0, "DE": 1.0, "DF": -6.54, "DG": 1.0, "DH": 1.0, "DI": 1.0, "DK": -7.49,
    "DL": 1.0, "DM": 1.0, "DN": 1.0, "DP": 1.0, "DQ": 1.0, "DR": -6.54, "DS": 20.26,
    "DT": -14.03, "DV": 1.0, "DW": 1.0, "DY": 1.0, "EA": 1.0, "EC": 44.94, "ED": 20.26,
    "EE": 33.6, "EF": 1.0, "EG": 1.0, "EH": -6.54, "EI": 20.26, "EK": 1.0, "EL": 1.0,
    "EM": 1.0, "EN": 1.0, "EP": 20.26, "EQ": 20.26, "ER": 1.0, "ES": 20.26, "ET": 1.0,
    "EV": 1.0, "EW": -14.03, "EY": 1.0, "FA": 1.0, "FC": 1.0, "FD": 13.34, "FE": 1.0,
    "FF": 1.0, "FG": 1.0, "FH": 1.0, "FI": 1.0, "FK": -14.03, "FL": 1.0, "FM": 1.0,
    "FN": 1.0, "FP": 20.26, "FQ": 1.0, "FR": 1.0, "FS": 1.0, "FT": 1.0, "FV": 1.0,
    "FW": 1.0, "FY": 33.601, "GA": -7.49, "GC": 1.0, "GD": 1.0, "GE": -6.54, "GF": 1.0,
    "GG": 13.34, "GH": 1.0, "GI": -7.49, "GK": -7.49, "GL": 1.0, "GM": 1.0, "GN": -7.49,
    "GP": 1.0, "GQ": 1.0, "GR": 1.0, "GS": 1.0, "GT": -7.49, "GV": 1.0, "GW": 13.34,
    "GY": -7.49, "HA": 1.0, "HC": 1.0, "HD": 1.0, "HE": 1.0, "HF": -9.37, "HG": -9.37,
    "HH": 1.0, "HI": 44.94, "HK": 24.68, "HL": 1.0, "HM": 1.0, "HN": 24.68, "HP": -1.88,
    "HQ": 1.0, "HR": 1.0, "HS": 1.0, "HT": -6.54, "HV": 1.0, "HW": -1.88, "HY": 44.94,
    "IA": 1.0, "IC": 1.0, "ID": 1.0, "IE": 44.94, "IF": 1.0, "IG": 1.0, "IH": 13.34,
    "II": 1.0, "IK": -7.49, "IL": 20.26, "IM": 1.0, "IN": 1.0, "IP": -1.88, "IQ": 1.0,
    "IR": 1.0, "IS": 1.0, "IT": 1.0, "IV": -7.49, "IW": 1.0, "IY": 1.0, "KA": 1.0,
    "KC": 1.0, "KD": 1.0, "KE": 1.0, "KF": 1.0, "KG": -7.49, "KH": 1.0, "KI": -7.49,
    "KK": 1.0, "KL": -7.49, "KM": 33.6, "KN": 1.0, "KP": -6.54, "KQ": 24.64, "KR": 33.6,
    "KS": 1.0, "KT": 1.0, "KV": -7.49, "KW": 1.0, "KY": 1.0, "LA": 1.0, "LC": 1.0,
    "LD": 1.0, "LE": 1.0, "LF": 1.0, "LG": 1.0, "LH": 1.0, "LI": 1.0, "LK": -7.49,
    "LL": 1.0, "LM": 1.0, "LN": 1.0, "LP": 20.26, "LQ": 33.6, "LR": 20.26, "LS": 1.0,
    "LT": 1.0, "LV": 1.0, "LW": 24.68, "LY": 1.0, "MA": 13.34, "MC": 1.0, "MD": 1.0,
    "ME": 1.0, "MF": 1.0, "MG": 1.0, "MH": 58.28, "MI": 1.0, "MK": 1.0, "ML": 1.0,
    "MM": -1.88, "MN": 1.0, "MP": 44.94, "MQ": -6.54, "MR": -6.54, "MS": 44.94, "MT": -1.88,
    "MV": 1.0, "MW": 1.0, "MY": 24.68, "NA": 1.0, "NC": -1.88, "ND": 1.0, "NE": 1.0,
    "NF": -14.03, "NG": -14.03, "NH": 1.0, "NI": 44.94, "NK": 24.68, "NL": 1.0, "NM": 1.0,
    "NN": 1.0, "NP": -1.88, "NQ": -6.54, "NR": 1.0, "NS": 1.0, "NT": -7.49, "NV": 1.0,
    "NW": -9.37, "NY": 1.0, "PA": 20.26, "PC": -6.54, "PD": -6.54, "PE": 18.38, "PF": 20.26,
    "PG": 1.0, "PH": 1.0, "PI": 1.0, "PK": 1.0, "PL": 1.0, "PM": -6.54, "PN": 1.0,
    "PP": 20.26, "PQ": 20.26, "PR": -6.54, "PS": 20.26, "PT": 1.0, "PV": 20.26, "PW": -1.88,
    "PY": 1.0, "QA": 1.0, "QC": -6.54, "QD": 20.26, "QE": 20.26, "QF": -6.54, "QG": 1.0,
    "QH": 1.0, "QI": 1.0, "QK": 1.0, "QL": 1.0, "QM": 1.0, "QN": 1.0, "QP": 20.26,
    "QQ": 20.26, "QR": 1.0, "QS": 44.94, "QT": 1.0, "QV": -6.54, "QW": 1.0, "QY": -6.54,
    "RA": 1.0, "RC": 1.0, "RD": 1.0, "RE": 1.0, "RF": 1.0, "RG": -7.49, "RH": 20.26,
    "RI": 1.0, "RK": 1.0, "RL": 1.0, "RM": 1.0, "RN": 13.34, "RP": 20.26, "RQ": 20.26,
    "RR": 58.28, "RS": 44.94, "RT": 1.0, "RV": 1.0, "RW": 58.28, "RY": -6.54, "SA": 1.0,
    "SC": 33.6, "SD": 1.0, "SE": 20.26, "SF": 1.0, "SG": 1.0, "SH": 1.0, "SI": 1.0,
    "SK": 1.0, "SL": 1.0, "SM": 1.0, "SN": 1.0, "SP": 44.94, "SQ": 20.26, "SR": 20.26,
    "SS": 20.26, "ST": 1.0, "SV": 1.0, "SW": 1.0, "SY": 1.0, "TA": 1.0, "TC": 1.0,
    "TD": 1.0, "TE": 20.26, "TF": 13.34, "TG": -7.49, "TH": 1.0, "TI": 1.0, "TK": 1.0,
    "TL": 1.0, "TM": 1.0, "TN": -14.03, "TP": 1.0, "TQ": -6.54, "TR": 1.0, "TS": 1.0,
    "TT": 1.0, "TV": 1.0, "TW": -14.03, "TY": 1.0, "VA": 1.0, "VC": 1.0, "VD": -14.03,
    "VE": 1.0, "VF": 1.0, "VG": -7.49, "VH": 1.0, "VI": 1.0, "VK": -1.88, "VL": 1.0,
    "VM": 1.0, "VN": 1.0, "VP": 20.26, "VQ": 1.0, "VR": 1.0, "VS": 1.0, "VT": -7.49,
    "VV": 1.0, "VW": 1.0, "VY": -6.54, "WA": -14.03, "WC": 1.0, "WD": 1.0, "WE": 1.0,
    "WF": 1.0, "WG": -9.37, "WH": 24.68, "WI": 1.0, "WK": 1.0, "WL": 13.34, "WM": 24.68,
    "WN": 13.34, "WP": 1.0, "WQ": 1.0, "WR": 1.0, "WS": 1.0, "WT": -14.03, "WV": -7.49,
    "WW": 1.0, "WY": 1.0, "YA": 24.68, "YC": 1.0, "YD": 24.68, "YE": -6.54, "YF": 1.0,
    "YG": -7.49, "YH": 13.34, "YI": 1.0, "YK": 1.0, "YL": 1.0, "YM": 44.94, "YN": 1.0,
    "YP": 13.34, "YQ": 1.0, "YR": -15.91, "YS": 1.0, "YT": -7.49, "YV": 1.0, "YW": -9.37,
    "YY": 13.34,
}
# fmt: on

#: Mass in Da of the water molecule that a polypeptide gains on top of the sum of its residue masses
WATER_MASS = 18.01528


def _store_fasta_record(sequences: Dict[str, str], record_id: str, sequence: str) -> None:
    """Add a record to ``sequences``; a repeated ID is an error because the first record would be lost."""
    if record_id in sequences:
        raise ValueError(f"Duplicate sequence ID {record_id!r}")
    sequences[record_id] = sequence


def read_fasta(path: Union[str, Path]) -> Dict[str, str]:
    """Read protein sequences from a FASTA file.

    Args:
        path: Path to the FASTA file

    Returns:
        Dictionary mapping sequence IDs (the first word of each header) to upper-case sequence strings

    Raises:
        FileNotFoundError: If the FASTA file doesn't exist
        ValueError: If the FASTA file format is invalid, it has no sequences, or two records share an ID (the
            dictionary could keep only one of them)
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    sequences: Dict[str, str] = {}
    current_id = None
    current_seq: List[str] = []

    try:
        with io.open_text_auto(path) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    # Save previous sequence
                    if current_id is not None:
                        _store_fasta_record(sequences, current_id, "".join(current_seq))

                    # Start new sequence
                    current_id = line[1:].split()[0]  # Take first word as ID
                    current_seq = []
                else:
                    # Sequence line
                    if current_id is None:
                        raise ValueError(f"FASTA file must start with '>' header line at {path}:{line_num}")
                    current_seq.append(line.upper())

        # Save last sequence
        if current_id is not None:
            _store_fasta_record(sequences, current_id, "".join(current_seq))

    except Exception as e:
        raise ValueError(f"Error reading FASTA file {path}: {e}")

    if not sequences:
        raise ValueError(f"No sequences found in FASTA file: {path}")

    logger.debug(f"Read {len(sequences)} protein sequences from {path}")
    return sequences


def write_fasta(sequences: Dict[str, str], path: Union[str, Path], line_width: int = 60) -> None:
    """Write protein sequences to a FASTA file.

    Args:
        sequences: Dictionary mapping sequence IDs to sequence strings
        path: Output file path
        line_width: Maximum characters per line (0 for no wrapping)
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with io.open_text_auto(path, "w") as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n")

            if line_width > 0:
                # Wrap sequence lines
                for i in range(0, len(sequence), line_width):
                    f.write(sequence[i : i + line_width] + "\n")
            else:
                # Single line per sequence
                f.write(sequence + "\n")

    logger.debug(f"Wrote {len(sequences)} protein sequences to {path}")


def validate_protein_sequence(seq: str) -> bool:
    """Validate that a sequence contains only valid amino acid characters.

    Args:
        seq: Sequence to validate

    Returns:
        True if sequence is valid protein, False otherwise
    """
    if not seq:
        return False

    # Standard amino acid codes (including ambiguous)
    valid_chars = set("ACDEFGHIKLMNPQRSTVWYBXZJUO-")
    return all(c in valid_chars for c in seq.upper())


def sequence_length(seq: str) -> int:
    """Get the length of a protein sequence.

    Args:
        seq: Protein sequence string

    Returns:
        Sequence length
    """
    return len(seq)


def molecular_weight(seq: str) -> float:
    """Calculate molecular weight of a protein sequence.

    The weight is the sum of the average residue masses plus one water molecule (``WATER_MASS``, 18.015 Da), which
    the free termini add to the residues (so a single glycine is 75.07 Da, the mass of the free amino acid). Letters
    that are not an amino acid with a known residue mass are skipped with a warning, and their mass is not counted.

    Args:
        seq: Protein sequence string

    Returns:
        Molecular weight in Da (0.0 when the sequence holds no residue of known mass)
    """
    if not seq:
        return 0.0

    # Average residue weights (including water loss)
    residue_weights = {
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
        "U": 150.0388,  # Selenocysteine
        "O": 237.3018,  # Pyrrolysine
    }

    weight = 0.0
    n_known = 0
    for aa in seq.upper():
        if aa in residue_weights:
            weight += residue_weights[aa]
            n_known += 1
        else:
            logger.warning(f"Unknown amino acid: {aa}")

    return weight + WATER_MASS if n_known else 0.0


def isoelectric_point(seq: str) -> float:
    """Calculate isoelectric point using Henderson-Hasselbalch bisection method.

    Uses pKa values for amino acid side chains, N-terminus, and C-terminus
    to find the pH at which the net charge is zero.

    Args:
        seq: Protein sequence string

    Returns:
        Isoelectric point (pI)
    """
    if not seq:
        return 0.0

    seq = seq.upper()

    # pKa values
    pka_nterm = 9.69
    pka_cterm = 2.34
    pka_side = {
        "D": 3.65,
        "E": 4.25,  # acidic
        "C": 8.18,
        "Y": 10.07,  # weakly acidic
        "H": 6.00,
        "K": 10.53,
        "R": 12.48,  # basic
    }

    def _charge_at_ph(ph: float) -> float:
        """Calculate net charge at given pH."""
        charge = 0.0
        # N-terminus (positive)
        charge += 1.0 / (1.0 + 10 ** (ph - pka_nterm))
        # C-terminus (negative)
        charge -= 1.0 / (1.0 + 10 ** (pka_cterm - ph))
        # Side chains
        for aa in seq:
            if aa in ("D", "E", "C", "Y"):
                charge -= 1.0 / (1.0 + 10 ** (pka_side[aa] - ph))
            elif aa in ("H", "K", "R"):
                charge += 1.0 / (1.0 + 10 ** (ph - pka_side[aa]))
        return charge

    # Bisection method
    low, high = 0.0, 14.0
    for _ in range(100):
        mid = (low + high) / 2.0
        charge = _charge_at_ph(mid)
        if charge > 0:
            low = mid
        else:
            high = mid
    return round((low + high) / 2.0, 2)


def find_motifs(seq: str, motif_patterns: List[str]) -> Dict[str, List[int]]:
    """Find occurrences of motifs in a protein sequence.

    Args:
        seq: Protein sequence string
        motif_patterns: List of motif patterns (can include regex)

    Returns:
        Dictionary mapping motif patterns to lists of start positions
    """
    results = {}

    for pattern in motif_patterns:
        positions = []
        try:
            # Try as regex first
            for match in re.finditer(pattern.upper(), seq.upper()):
                positions.append(match.start())
        except re.error:
            # Fall back to simple string matching
            pattern_upper = pattern.upper()
            start = 0
            while True:
                pos = seq.upper().find(pattern_upper, start)
                if pos == -1:
                    break
                positions.append(pos)
                start = pos + 1

        results[pattern] = positions

    return results


def hydropathy_score(seq: str, window_size: int = 19) -> List[float]:
    """Calculate hydropathy scores using sliding window.

    Args:
        seq: Protein sequence string
        window_size: Size of sliding window

    Returns:
        List of hydropathy scores
    """
    # Kyte-Doolittle hydropathy scale (canonical table above)
    scores = []
    seq_upper = seq.upper()

    for i in range(len(seq_upper) - window_size + 1):
        window = seq_upper[i : i + window_size]
        window_score = sum(KD_HYDROPHOBICITY.get(aa, 0.0) for aa in window) / window_size
        scores.append(window_score)

    return scores


def transmembrane_regions(seq: str, threshold: float = 1.6) -> List[Tuple[int, int]]:
    """Predict transmembrane regions using hydropathy analysis.

    Args:
        seq: Protein sequence string
        threshold: Hydropathy threshold for TM regions

    Returns:
        List of (start, end) tuples for predicted TM regions
    """
    scores = hydropathy_score(seq, window_size=19)
    regions = []

    i = 0
    while i < len(scores):
        if scores[i] >= threshold:
            # Start of potential TM region
            start = i
            while i < len(scores) and scores[i] >= threshold:
                i += 1
            end = i - 1

            # Check if region is long enough (typically 20-25 residues)
            if end - start >= 15:
                regions.append((start, end))
        else:
            i += 1

    return regions


def amino_acid_composition(seq: str) -> Dict[str, float]:
    """Calculate amino acid composition of a protein sequence.

    Args:
        seq: Protein sequence string

    Returns:
        Dictionary mapping amino acids to percentages
    """
    if not seq:
        return {}

    from collections import Counter

    seq_upper = seq.upper()
    counts = Counter(seq_upper)
    total = sum(counts.values())

    composition = {}
    for aa, count in counts.items():
        composition[aa] = (count / total) * 100

    return composition


def extinction_coefficient(seq: str, disulfide_bonds: int = 0) -> Dict[str, float]:
    """Calculate molar extinction coefficient at 280 nm.

    Uses the Pace method based on Trp, Tyr, and Cys content.

    Args:
        seq: Protein sequence string
        disulfide_bonds: Number of disulfide bonds (cystines)

    Returns:
        Dictionary with extinction coefficients (reduced and oxidized)
    """
    if not seq:
        return {"reduced": 0.0, "oxidized": 0.0}

    seq = seq.upper()
    n_trp = seq.count("W")
    n_tyr = seq.count("Y")
    n_cys = seq.count("C")

    # Pace et al. (1995) coefficients
    ext_reduced = n_trp * 5500 + n_tyr * 1490
    ext_oxidized = ext_reduced + disulfide_bonds * 125

    return {
        "reduced": float(ext_reduced),
        "oxidized": float(ext_oxidized),
        "n_trp": n_trp,
        "n_tyr": n_tyr,
        "n_cys": n_cys,
    }


def gravy(seq: str) -> float:
    """Calculate Grand Average of Hydropathy (GRAVY) score.

    GRAVY is the average Kyte-Doolittle hydropathy value across
    the entire sequence. Positive = hydrophobic, negative = hydrophilic.

    Args:
        seq: Protein sequence string

    Returns:
        GRAVY score
    """
    if not seq:
        return 0.0

    total = sum(KD_HYDROPHOBICITY.get(aa, 0.0) for aa in seq.upper())
    return total / len(seq)


def instability_index(seq: str) -> float:
    """Calculate protein instability index (Guruprasad et al., 1990).

    ``II = (10 / L) * sum(DIWV[x_i x_i+1])`` over the ``L - 1`` neighbouring pairs of residues, with the weights of
    :data:`DIWV`. Values above 40 predict a protein that is unstable in a test tube. A pair that contains a letter
    other than the twenty standard amino acids has no weight in the table and counts with the neutral weight 1.0.

    Args:
        seq: Protein sequence string

    Returns:
        Instability index value (0.0 for a sequence of fewer than two residues)

    Example:
        >>> instability_index("WM")  # (10 / 2) * DIWV["WM"], and DIWV["WM"] is 24.68
        123.4
    """
    if len(seq) < 2:
        return 0.0

    seq = seq.upper()
    score = sum(DIWV.get(seq[i : i + 2], 1.0) for i in range(len(seq) - 1))
    return (10.0 / len(seq)) * score


def aromaticity(seq: str) -> float:
    """Calculate aromaticity index (frequency of aromatic amino acids).

    Aromatic amino acids: Phe (F), Trp (W), Tyr (Y).

    Args:
        seq: Protein sequence string

    Returns:
        Aromaticity value (0.0 to 1.0)
    """
    if not seq:
        return 0.0

    seq = seq.upper()
    aromatic_count = seq.count("F") + seq.count("W") + seq.count("Y")
    return aromatic_count / len(seq)


def charge_at_ph(seq: str, ph: float = 7.0) -> float:
    """Calculate net charge of protein at a given pH.

    Args:
        seq: Protein sequence string
        ph: pH value (default 7.0)

    Returns:
        Net charge at the given pH
    """
    if not seq:
        return 0.0

    seq = seq.upper()

    pka_nterm = 9.69
    pka_cterm = 2.34
    pka_side = {
        "D": 3.65,
        "E": 4.25,
        "C": 8.18,
        "Y": 10.07,
        "H": 6.00,
        "K": 10.53,
        "R": 12.48,
    }

    charge = 0.0
    charge += 1.0 / (1.0 + 10 ** (ph - pka_nterm))
    charge -= 1.0 / (1.0 + 10 ** (pka_cterm - ph))

    for aa in seq:
        if aa in ("D", "E", "C", "Y"):
            charge -= 1.0 / (1.0 + 10 ** (pka_side[aa] - ph))
        elif aa in ("H", "K", "R"):
            charge += 1.0 / (1.0 + 10 ** (ph - pka_side[aa]))

    return round(charge, 4)


def sequence_summary(seq: str) -> Dict[str, Any]:
    """Generate a comprehensive summary of protein sequence properties.

    Args:
        seq: Protein sequence string

    Returns:
        Dictionary with all computed sequence properties
    """
    if not seq:
        return {}

    return {
        "length": sequence_length(seq),
        "molecular_weight": molecular_weight(seq),
        "isoelectric_point": isoelectric_point(seq),
        "gravy": gravy(seq),
        "aromaticity": aromaticity(seq),
        "instability_index": instability_index(seq),
        "extinction_coefficient": extinction_coefficient(seq),
        "charge_at_ph7": charge_at_ph(seq, 7.0),
        "amino_acid_composition": amino_acid_composition(seq),
        "is_valid": validate_protein_sequence(seq),
    }


# Aliases for backward compatibility
parse_fasta = read_fasta
is_valid_protein_sequence = validate_protein_sequence


def calculate_aa_composition(seq: str) -> Dict[str, float]:
    """Calculate amino acid composition as fractions (0-1).

    Args:
        seq: Protein sequence

    Returns:
        Dictionary mapping amino acids to their fractional abundance
    """
    comp = amino_acid_composition(seq)
    return {aa: pct / 100.0 for aa, pct in comp.items()}


def kmer_frequencies(seq: str, k: int = 2) -> Dict[str, int]:
    """Calculate k-mer frequencies in a protein sequence.

    Args:
        seq: Protein sequence
        k: Length of k-mers

    Returns:
        Dictionary mapping k-mers to their counts
    """
    if k < 1 or k > len(seq):
        return {}

    counts: Dict[str, int] = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        counts[kmer] = counts.get(kmer, 0) + 1

    return counts

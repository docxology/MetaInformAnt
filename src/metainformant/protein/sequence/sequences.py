"""Protein sequence processing and analysis utilities.

This module provides core functionality for reading, writing, and analyzing
protein sequences in FASTA format.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, Iterator, List, Tuple, Union

from metainformant.core import io, logging

logger = logging.get_logger(__name__)


def read_fasta(path: Union[str, Path]) -> Dict[str, str]:
    """Read protein sequences from a FASTA file.

    Args:
        path: Path to the FASTA file

    Returns:
        Dictionary mapping sequence IDs to sequence strings

    Raises:
        FileNotFoundError: If the FASTA file doesn't exist
        ValueError: If the FASTA file format is invalid
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    sequences = {}
    current_id = None
    current_seq = []

    try:
        with io.open_text_auto(path) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    # Save previous sequence
                    if current_id is not None:
                        sequences[current_id] = "".join(current_seq)

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
            sequences[current_id] = "".join(current_seq)

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

    Args:
        seq: Protein sequence string

    Returns:
        Molecular weight in Da
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
    for aa in seq.upper():
        if aa in residue_weights:
            weight += residue_weights[aa]
        else:
            logger.warning(f"Unknown amino acid: {aa}")

    return weight


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
    # Kyte-Doolittle hydropathy scale
    kd_scale = {
        "I": 4.5,
        "V": 4.2,
        "L": 3.8,
        "P": 1.6,
        "A": 1.8,
        "W": -0.9,
        "M": 1.9,
        "H": -3.2,
        "T": -0.7,
        "F": 2.8,
        "C": 2.5,
        "N": -3.5,
        "Q": -3.5,
        "Y": -1.3,
        "D": -3.5,
        "E": -3.5,
        "K": -3.9,
        "R": -4.5,
        "S": -0.8,
        "G": -0.4,
    }

    scores = []
    seq_upper = seq.upper()

    for i in range(len(seq_upper) - window_size + 1):
        window = seq_upper[i : i + window_size]
        window_score = sum(kd_scale.get(aa, 0.0) for aa in window) / window_size
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

    kd_scale = {
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

    total = sum(kd_scale.get(aa, 0.0) for aa in seq.upper())
    return total / len(seq)


def instability_index(seq: str) -> float:
    """Calculate protein instability index (Guruprasad et al., 1990).

    Values above 40 indicate the protein is likely unstable in vitro.

    Args:
        seq: Protein sequence string

    Returns:
        Instability index value
    """
    if len(seq) < 2:
        return 0.0

    # DIWV weight values for dipeptides (subset of most impactful pairs)
    diwv = {
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

    seq = seq.upper()
    score = 0.0
    for i in range(len(seq) - 1):
        dipeptide = seq[i : i + 2]
        score += diwv.get(dipeptide, 1.0)

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

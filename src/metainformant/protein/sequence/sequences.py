"""Protein sequence processing and analysis utilities.

This module provides core functionality for reading, writing, and analyzing
protein sequences in FASTA format.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterator, List, Tuple, Union

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
    """Calculate isoelectric point of a protein sequence.

    Args:
        seq: Protein sequence string

    Returns:
        Isoelectric point (pI)
    """
    # Simplified pI calculation
    # In practice, would use more sophisticated algorithms

    acidic_count = seq.upper().count("D") + seq.upper().count("E")
    basic_count = seq.upper().count("R") + seq.upper().count("H") + seq.upper().count("K")

    if acidic_count == 0 and basic_count == 0:
        return 7.0  # Neutral

    # Rough approximation
    if acidic_count > basic_count:
        return 4.0  # Acidic
    elif basic_count > acidic_count:
        return 10.0  # Basic
    else:
        return 7.0  # Neutral


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

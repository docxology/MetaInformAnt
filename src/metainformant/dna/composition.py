"""DNA sequence composition analysis utilities.

This module provides functions for analyzing nucleotide composition,
GC content, skew calculations, and other sequence composition metrics.
"""

from __future__ import annotations

from typing import List

from metainformant.core import logging

logger = logging.get_logger(__name__)


def gc_content(seq: str) -> float:
    """Calculate GC content of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        GC content as a fraction (0.0 to 1.0)
    """
    if not seq:
        return 0.0

    seq_upper = seq.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    total_count = len(seq_upper)

    return gc_count / total_count if total_count > 0 else 0.0


def gc_skew(seq: str) -> float:
    """Calculate GC skew of a DNA sequence.

    GC skew = (G - C) / (G + C)

    Args:
        seq: DNA sequence string

    Returns:
        GC skew value
    """
    if not seq:
        return 0.0

    seq_upper = seq.upper()
    g_count = seq_upper.count('G')
    c_count = seq_upper.count('C')

    if g_count + c_count == 0:
        return 0.0

    return (g_count - c_count) / (g_count + c_count)


def cumulative_gc_skew(seq: str) -> List[float]:
    """Calculate cumulative GC skew along a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        List of cumulative GC skew values
    """
    if not seq:
        return []

    seq_upper = seq.upper()
    skews = []

    g_count = 0
    c_count = 0

    for base in seq_upper:
        if base == 'G':
            g_count += 1
        elif base == 'C':
            c_count += 1

        if g_count + c_count > 0:
            skew = (g_count - c_count) / (g_count + c_count)
        else:
            skew = 0.0

        skews.append(skew)

    return skews


def at_skew(seq: str) -> float:
    """Calculate AT skew of a DNA sequence.

    AT skew = (A - T) / (A + T)

    Args:
        seq: DNA sequence string

    Returns:
        AT skew value
    """
    if not seq:
        return 0.0

    seq_upper = seq.upper()
    a_count = seq_upper.count('A')
    t_count = seq_upper.count('T')

    if a_count + t_count == 0:
        return 0.0

    return (a_count - t_count) / (a_count + t_count)


def melting_temperature(seq: str) -> float:
    """Calculate DNA melting temperature using Wallace rule.

    Tm = 4(G+C) + 2(A+T)

    Args:
        seq: DNA sequence string

    Returns:
        Melting temperature in Celsius
    """
    if not seq:
        return 0.0

    seq_upper = seq.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    at_count = seq_upper.count('A') + seq_upper.count('T')

    return 4 * gc_count + 2 * at_count


def nucleotide_frequencies(seq: str) -> dict[str, float]:
    """Calculate nucleotide frequencies in a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        Dictionary mapping nucleotides to frequencies
    """
    if not seq:
        return {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0}

    seq_upper = seq.upper()
    length = len(seq_upper)

    frequencies = {
        'A': seq_upper.count('A') / length,
        'T': seq_upper.count('T') / length,
        'G': seq_upper.count('G') / length,
        'C': seq_upper.count('C') / length,
    }

    return frequencies


def dinucleotide_frequencies(seq: str) -> dict[str, float]:
    """Calculate dinucleotide frequencies in a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        Dictionary mapping dinucleotides to frequencies
    """
    if len(seq) < 2:
        return {}

    from collections import Counter

    dinucs = []
    seq_upper = seq.upper()

    for i in range(len(seq_upper) - 1):
        dinuc = seq_upper[i:i + 2]
        dinucs.append(dinuc)

    counts = Counter(dinucs)
    total = sum(counts.values())

    frequencies = {}
    for dinuc, count in counts.items():
        frequencies[dinuc] = count / total if total > 0 else 0.0

    return frequencies


def codon_frequencies(seq: str) -> dict[str, float]:
    """Calculate codon frequencies in a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        Dictionary mapping codons to frequencies
    """
    if len(seq) < 3:
        return {}

    from collections import Counter

    codons = []
    seq_upper = seq.upper()

    for i in range(0, len(seq_upper) - 2, 3):
        codon = seq_upper[i:i + 3]
        codons.append(codon)

    counts = Counter(codons)
    total = sum(counts.values())

    frequencies = {}
    for codon, count in counts.items():
        frequencies[codon] = count / total if total > 0 else 0.0

    return frequencies


def amino_acid_frequencies(seq: str) -> dict[str, float]:
    """Calculate amino acid frequencies from DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        Dictionary mapping amino acids to frequencies
    """
    if not seq:
        return {}

    # Translate DNA to protein
    from metainformant.dna.transcription import translate_dna

    try:
        protein = translate_dna(seq)
    except Exception:
        return {}

    from collections import Counter

    # Count amino acids (exclude stop codons)
    aa_counts = Counter(protein.replace('*', ''))

    total = sum(aa_counts.values())
    frequencies = {}

    for aa, count in aa_counts.items():
        frequencies[aa] = count / total if total > 0 else 0.0

    return frequencies


def purine_pyrimidine_ratio(seq: str) -> float:
    """Calculate purine/pyrimidine ratio.

    Purines: A, G
    Pyrimidines: C, T

    Args:
        seq: DNA sequence string

    Returns:
        Purine/pyrimidine ratio
    """
    if not seq:
        return 1.0

    seq_upper = seq.upper()
    purines = seq_upper.count('A') + seq_upper.count('G')
    pyrimidines = seq_upper.count('C') + seq_upper.count('T')

    if pyrimidines == 0:
        return float('inf') if purines > 0 else 1.0

    return purines / pyrimidines


def base_pair_stacking_energy(seq: str) -> float:
    """Calculate approximate base pair stacking energy.

    This is a simplified calculation based on nearest-neighbor
    thermodynamic parameters.

    Args:
        seq: DNA sequence string

    Returns:
        Approximate stacking energy contribution
    """
    if len(seq) < 2:
        return 0.0

    # Simplified stacking energies (kcal/mol at 37°C)
    # Based on SantaLucia (1998) nearest-neighbor parameters
    stacking_energies = {
        'AA': -1.0, 'AT': -0.9, 'AG': -1.3, 'AC': -1.5,
        'TA': -0.6, 'TT': -1.0, 'TG': -1.4, 'TC': -1.3,
        'GA': -1.1, 'GT': -1.5, 'GG': -1.6, 'GC': -2.1,
        'CA': -1.7, 'CT': -1.8, 'CG': -2.7, 'CC': -1.5,
    }

    seq_upper = seq.upper()
    energy = 0.0

    for i in range(len(seq_upper) - 1):
        dinuc = seq_upper[i:i + 2]
        energy += stacking_energies.get(dinuc, 0.0)

    return energy


def sequence_complexity(seq: str, window_size: int = 10) -> List[float]:
    """Calculate sequence complexity in sliding windows.

    Args:
        seq: DNA sequence string
        window_size: Size of sliding window

    Returns:
        List of complexity scores for each window
    """
    if len(seq) < window_size:
        return []

    complexities = []

    for i in range(len(seq) - window_size + 1):
        window = seq[i:i + window_size]

        # Calculate Shannon entropy
        from collections import Counter
        import math

        counts = Counter(window.upper())
        total = sum(counts.values())

        entropy = 0.0
        for count in counts.values():
            if count > 0:
                p = count / total
                entropy -= p * math.log2(p)

        # Normalize by maximum entropy (2 bits for DNA)
        complexity = entropy / 2.0
        complexities.append(complexity)

    return complexities


def cpg_islands(seq: str, min_length: int = 200, min_gc: float = 0.5,
                min_cpg_ratio: float = 0.6) -> List[tuple[int, int]]:
    """Identify CpG islands in DNA sequence.

    Args:
        seq: DNA sequence string
        min_length: Minimum island length
        min_gc: Minimum GC content threshold
        min_cpg_ratio: Minimum CpG ratio threshold

    Returns:
        List of tuples (start_pos, end_pos) for CpG islands
    """
    islands = []
    seq_upper = seq.upper()

    i = 0
    while i < len(seq_upper):
        # Look for regions with high GC and CpG content
        start = i
        gc_count = 0
        cpg_count = 0
        total_bases = 0

        while i < len(seq_upper) and total_bases < min_length:
            if seq_upper[i] in 'GC':
                gc_count += 1

            # Check for CpG dinucleotide
            if i < len(seq_upper) - 1 and seq_upper[i:i + 2] == 'CG':
                cpg_count += 1

            total_bases += 1
            i += 1

        # Check if this region meets criteria
        if total_bases >= min_length:
            gc_ratio = gc_count / total_bases
            expected_cpg = (seq_upper[start:i].count('C') * seq_upper[start:i].count('G')) / total_bases
            observed_cpg = cpg_count

            if expected_cpg > 0:
                cpg_ratio = observed_cpg / expected_cpg
            else:
                cpg_ratio = 0.0

            if gc_ratio >= min_gc and cpg_ratio >= min_cpg_ratio:
                islands.append((start, i - 1))

        i = start + 1  # Slide window

    return islands



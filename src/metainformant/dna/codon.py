"""Codon usage analysis and codon optimization utilities.

This module provides tools for analyzing codon usage bias, calculating codon
adaptation index (CAI), and performing codon optimization for heterologous expression.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)


# Standard genetic code (NCBI translation table 1)
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'G', 'GCC': 'G', 'GCA': 'G', 'GCG': 'G',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def codon_usage(seq: str) -> Dict[str, float]:
    """Calculate codon usage frequencies from a DNA sequence.

    Args:
        seq: DNA sequence (must be divisible by 3)

    Returns:
        Dictionary mapping codons to their frequencies

    Raises:
        ValueError: If sequence length is not divisible by 3

    Example:
        >>> seq = "ATGGCCATTGTAATGGGCC"
        >>> usage = codon_usage(seq)
        >>> "ATG" in usage
        True
    """
    if len(seq) % 3 != 0:
        raise ValueError("Sequence length must be divisible by 3")

    if not seq:
        return {}

    # Count codons
    codon_counts = {}
    total_codons = 0

    for i in range(0, len(seq), 3):
        codon = seq[i:i+3].upper()
        if len(codon) == 3 and codon in GENETIC_CODE:
            codon_counts[codon] = codon_counts.get(codon, 0) + 1
            total_codons += 1

    # Convert to frequencies
    codon_frequencies = {}
    for codon, count in codon_counts.items():
        codon_frequencies[codon] = count / total_codons if total_codons > 0 else 0.0

    return codon_frequencies


def cai(sequence: str, reference_usage: Optional[Dict[str, float]] = None) -> float:
    """Calculate Codon Adaptation Index (CAI) for a sequence.

    CAI measures how well codon usage matches a reference set (typically
    highly expressed genes).

    Args:
        sequence: DNA sequence
        reference_usage: Reference codon usage frequencies (optional)

    Returns:
        CAI value (0.0 to 1.0)

    Example:
        >>> seq = "ATGGCCATTGTAATGGGCC"
        >>> cai_value = cai(seq)
        >>> 0.0 <= cai_value <= 1.0
        True
    """
    if not sequence or len(sequence) % 3 != 0:
        return 0.0

    # Use default reference usage if not provided
    if reference_usage is None:
        reference_usage = _get_default_reference_usage()

    # Calculate CAI
    cai_values = []

    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3].upper()
        if len(codon) == 3 and codon in GENETIC_CODE:
            # Skip stop codons
            if GENETIC_CODE[codon] == '*':
                continue

            # Get reference frequency for this codon
            ref_freq = reference_usage.get(codon, 0.0)

            # Get maximum frequency for synonymous codons
            max_freq = _get_max_synonymous_frequency(codon, reference_usage)

            if max_freq > 0:
                cai_values.append(ref_freq / max_freq)

    return (sum(cai_values) / len(cai_values)) ** (1.0 / len(cai_values)) if cai_values else 0.0


def _get_default_reference_usage() -> Dict[str, float]:
    """Get default reference codon usage (E. coli highly expressed genes)."""
    # Simplified E. coli reference usage (subset)
    return {
        'ATG': 0.47, 'TTA': 0.07, 'TTG': 0.13, 'CTT': 0.12, 'CTC': 0.11, 'CTA': 0.04, 'CTG': 0.47,
        'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11, 'GTT': 0.28, 'GTC': 0.15, 'GTA': 0.11, 'GTG': 0.46,
        'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.12, 'TCG': 0.14, 'AGT': 0.14, 'AGC': 0.25,
        'CCT': 0.18, 'CCC': 0.12, 'CCA': 0.20, 'CCG': 0.50, 'ACT': 0.19, 'ACC': 0.40, 'ACA': 0.16, 'ACG': 0.25,
        'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33, 'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15,
        'TAT': 0.59, 'TAC': 0.41, 'CAT': 0.57, 'CAC': 0.43, 'CAA': 0.34, 'CAG': 0.66,
        'AAT': 0.49, 'AAC': 0.51, 'AAA': 0.74, 'AAG': 0.26, 'GAT': 0.63, 'GAC': 0.37, 'GAA': 0.68, 'GAG': 0.32,
        'TGT': 0.46, 'TGC': 0.54, 'TGG': 1.00, 'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11,
        'AGA': 0.04, 'AGG': 0.04
    }


def _get_max_synonymous_frequency(codon: str, reference_usage: Dict[str, float]) -> float:
    """Get maximum frequency for codons encoding the same amino acid."""
    if codon not in GENETIC_CODE:
        return 0.0

    amino_acid = GENETIC_CODE[codon]

    # Find all codons for this amino acid
    synonymous_codons = [c for c, aa in GENETIC_CODE.items() if aa == amino_acid]

    # Get their frequencies
    frequencies = [reference_usage.get(c, 0.0) for c in synonymous_codons]

    return max(frequencies) if frequencies else 0.0


def gc_content_codon_positions(seq: str) -> Dict[str, float]:
    """Calculate GC content at each codon position.

    Args:
        seq: DNA sequence

    Returns:
        Dictionary with GC content for each position (1, 2, 3)

    Example:
        >>> seq = "ATGGCCATTGTA"
        >>> gc_pos = gc_content_codon_positions(seq)
        >>> len(gc_pos) == 3
        True
    """
    if len(seq) % 3 != 0:
        seq = seq[:len(seq) - (len(seq) % 3)]  # Trim incomplete codons

    if not seq:
        return {'position_1': 0.0, 'position_2': 0.0, 'position_3': 0.0}

    gc_counts = {'position_1': 0, 'position_2': 0, 'position_3': 0}
    total_counts = {'position_1': 0, 'position_2': 0, 'position_3': 0}

    for i in range(0, len(seq), 3):
        codon = seq[i:i+3].upper()
        if len(codon) == 3:
            for pos in range(3):
                nucleotide = codon[pos]
                pos_key = f'position_{pos+1}'
                total_counts[pos_key] += 1
                if nucleotide in 'GC':
                    gc_counts[pos_key] += 1

    # Calculate percentages
    gc_content = {}
    for pos in range(1, 4):
        pos_key = f'position_{pos}'
        total = total_counts[pos_key]
        gc_content[pos_key] = (gc_counts[pos_key] / total * 100) if total > 0 else 0.0

    return gc_content


def synonymous_codons(amino_acid: str) -> List[str]:
    """Get all codons that encode a given amino acid.

    Args:
        amino_acid: Single-letter amino acid code

    Returns:
        List of codons encoding the amino acid

    Example:
        >>> codons = synonymous_codons("L")
        >>> len(codons) == 6  # Leucine has 6 codons
        True
    """
    return [codon for codon, aa in GENETIC_CODE.items() if aa == amino_acid]


def optimize_codons(sequence: str, target_usage: Dict[str, float]) -> str:
    """Optimize codon usage for a target organism.

    Args:
        sequence: Input DNA sequence
        target_usage: Target codon usage frequencies

    Returns:
        Codon-optimized sequence

    Example:
        >>> seq = "ATGGCCATTGTA"
        >>> optimized = optimize_codons(seq, {"GCC": 0.8, "GCT": 0.2})
        >>> len(optimized) == len(seq)
        True
    """
    if len(sequence) % 3 != 0:
        raise ValueError("Sequence length must be divisible by 3")

    if not sequence:
        return sequence

    optimized_sequence = []

    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3].upper()

        if codon not in GENETIC_CODE:
            # Keep unknown codons as-is
            optimized_sequence.append(codon)
            continue

        amino_acid = GENETIC_CODE[codon]

        # Skip stop codons
        if amino_acid == '*':
            optimized_sequence.append(codon)
            continue

        # Find best codon for this amino acid
        synonymous = synonymous_codons(amino_acid)
        if len(synonymous) <= 1:
            # No choice for this amino acid
            optimized_sequence.append(codon)
            continue

        # Choose codon with highest frequency in target usage
        best_codon = max(synonymous, key=lambda c: target_usage.get(c, 0.0))
        optimized_sequence.append(best_codon)

    return ''.join(optimized_sequence)


def calculate_enc(sequence: str, reference_usage: Optional[Dict[str, float]] = None) -> float:
    """Calculate Effective Number of Codons (ENC).

    ENC measures codon usage bias, with lower values indicating stronger bias.

    Args:
        sequence: DNA sequence
        reference_usage: Reference codon usage (optional)

    Returns:
        ENC value (20-61, where 20 = extreme bias, 61 = no bias)

    Example:
        >>> seq = "ATGGCCATTGTAATGGGCC"
        >>> enc = calculate_enc(seq)
        >>> 20 <= enc <= 61
        True
    """
    codon_freq = codon_usage(sequence)
    if not codon_freq:
        return 61.0  # Maximum ENC (no bias)

    # Use reference usage if provided, otherwise use observed
    usage = reference_usage or codon_freq

    # Group codons by amino acid
    aa_groups = {}
    for codon, freq in usage.items():
        if codon in GENETIC_CODE:
            aa = GENETIC_CODE[codon]
            if aa not in aa_groups:
                aa_groups[aa] = {}
            aa_groups[aa][codon] = freq

    # Calculate ENC
    total_enc = 0.0
    total_aa = 0

    for aa, codons in aa_groups.items():
        if aa == '*' or len(codons) <= 1:
            continue

        # Calculate F values for this amino acid
        codon_list = list(codons.keys())
        n = len(codon_list)

        if n > 1:
            # Calculate homozygosity
            f_values = []
            total_freq = sum(codons.values())

            if total_freq > 0:
                for codon in codon_list:
                    freq = codons[codon] / total_freq
                    f_values.append(freq ** 2)

                # ENC formula for this amino acid
                f_avg = sum(f_values)
                if f_avg > 0:
                    aa_enc = 1.0 / f_avg
                    total_enc += min(aa_enc, n)  # Cap at maximum for this AA
                    total_aa += 1

    return total_enc / total_aa if total_aa > 0 else 61.0


def detect_codon_bias(sequence: str) -> Dict[str, float]:
    """Detect codon usage bias patterns in a sequence.

    Args:
        sequence: DNA sequence

    Returns:
        Dictionary with bias metrics

    Example:
        >>> seq = "ATGGCCATTGTAATGGGCC"
        >>> bias = detect_codon_bias(seq)
        >>> "cai" in bias
        True
    """
    codon_freq = codon_usage(sequence)

    metrics = {
        'cai': cai(sequence),
        'enc': calculate_enc(sequence),
        'gc_positions': gc_content_codon_positions(sequence),
        'total_codons': sum(codon_freq.values()),
        'unique_codons': len(codon_freq)
    }

    return metrics


def back_translate(protein_seq: str, codon_preferences: Optional[Dict[str, str]] = None) -> str:
    """Back-translate a protein sequence to DNA using preferred codons.

    Args:
        protein_seq: Protein sequence (single-letter codes)
        codon_preferences: Dictionary mapping amino acids to preferred codons

    Returns:
        DNA sequence

    Example:
        >>> protein = "MA"
        >>> dna = back_translate(protein)
        >>> len(dna) == 6  # 2 amino acids * 3 nucleotides
        True
    """
    if not codon_preferences:
        # Use default preferences (most common codons)
        codon_preferences = {
            'A': 'GCU', 'R': 'CGU', 'N': 'AAU', 'D': 'GAU', 'C': 'UGU',
            'Q': 'CAA', 'E': 'GAA', 'G': 'GGU', 'H': 'CAU', 'I': 'AUU',
            'L': 'UUA', 'K': 'AAA', 'M': 'AUG', 'F': 'UUU', 'P': 'CCU',
            'S': 'UCU', 'T': 'ACU', 'W': 'UGG', 'Y': 'UAU', 'V': 'GUU',
            '*': 'UAA'  # Stop codon
        }

    dna_sequence = []

    for aa in protein_seq.upper():
        if aa in codon_preferences:
            dna_sequence.append(codon_preferences[aa])
        else:
            # Unknown amino acid
            dna_sequence.append('NNN')

    return ''.join(dna_sequence)








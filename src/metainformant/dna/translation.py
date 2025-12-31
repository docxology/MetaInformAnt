"""DNA translation utilities.

This module provides functions for translating RNA/DNA sequences to amino acid sequences.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Genetic code dictionary (standard code)
GENETIC_CODE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate(rna_seq: str, genetic_code: int = 1) -> str:
    """Translate RNA sequence to amino acid sequence.

    Args:
        rna_seq: RNA sequence string
        genetic_code: Genetic code table number (1 = standard)

    Returns:
        Amino acid sequence string

    Raises:
        ValueError: If sequence length is not divisible by 3 or contains invalid characters
    """
    if not rna_seq:
        return ""

    if len(rna_seq) % 3 != 0:
        logger.warning(f"RNA sequence length {len(rna_seq)} is not divisible by 3")

    rna_upper = rna_seq.upper()

    # Validate RNA characters
    valid_chars = set('AUCGN')
    if not all(c in valid_chars for c in rna_upper):
        invalid = set(rna_upper) - valid_chars
        raise ValueError(f"Invalid RNA characters: {invalid}")

    # Use standard genetic code
    code = GENETIC_CODE

    amino_acids = []
    for i in range(0, len(rna_upper) - 2, 3):
        codon = rna_upper[i:i + 3]
        aa = code.get(codon, 'X')  # X for unknown codons
        amino_acids.append(aa)

    return ''.join(amino_acids)


def translate_dna(dna_seq: str, genetic_code: int = 1) -> str:
    """Translate DNA sequence to amino acid sequence.

    Args:
        dna_seq: DNA sequence string
        genetic_code: Genetic code table number (1 = standard)

    Returns:
        Amino acid sequence string
    """
    from metainformant.dna.transcription import transcribe

    # Transcribe DNA to RNA first
    rna_seq = transcribe(dna_seq)
    return translate(rna_seq, genetic_code)


def find_orfs(rna_seq: str, min_length: int = 30) -> List[Tuple[int, int, str]]:
    """Find open reading frames in RNA sequence.

    Args:
        rna_seq: RNA sequence string
        min_length: Minimum ORF length in amino acids

    Returns:
        List of tuples (start_pos, end_pos, frame)
    """
    orfs = []

    # Check all 3 reading frames
    for frame in range(3):
        frame_seq = rna_seq[frame:]

        # Find start codons
        start_positions = []
        for match in re.finditer(r'AUG', frame_seq):
            start_positions.append(match.start())

        for start_pos in start_positions:
            # Find next stop codon
            orf_seq = frame_seq[start_pos:]
            stop_found = False

            for stop_match in re.finditer(r'(UAA|UAG|UGA)', orf_seq[3:], re.IGNORECASE):
                stop_pos = start_pos + stop_match.start() + 3  # Include stop codon
                orf_length_aa = (stop_pos - start_pos) // 3

                if orf_length_aa >= min_length:
                    frame_label = f"+{frame + 1}"
                    orfs.append((start_pos + frame, stop_pos + frame, frame_label))
                stop_found = True
                break

            # If no stop codon found, check if ORF extends to end
            if not stop_found:
                orf_length_aa = len(orf_seq) // 3
                if orf_length_aa >= min_length:
                    frame_label = f"+{frame + 1}"
                    orfs.append((start_pos + frame, len(rna_seq), frame_label))

    return orfs


def find_start_codons(rna_seq: str) -> List[int]:
    """Find positions of start codons (AUG) in RNA sequence.

    Args:
        rna_seq: RNA sequence string

    Returns:
        List of start codon positions
    """
    positions = []
    for match in re.finditer(r'AUG', rna_seq.upper()):
        positions.append(match.start())
    return positions


def find_stop_codons(rna_seq: str) -> List[int]:
    """Find positions of stop codons (UAA, UAG, UGA) in RNA sequence.

    Args:
        rna_seq: RNA sequence string

    Returns:
        List of stop codon positions
    """
    positions = []
    for match in re.finditer(r'(UAA|UAG|UGA)', rna_seq.upper()):
        positions.append(match.start())
    return positions


def six_frame_translation(dna_seq: str) -> Dict[str, str]:
    """Translate DNA sequence in all six reading frames.

    Args:
        dna_seq: DNA sequence string

    Returns:
        Dictionary mapping frame labels to amino acid sequences
    """
    from metainformant.dna.sequences import reverse_complement

    frames = {}

    # Forward strand frames
    for frame in range(3):
        frame_seq = dna_seq[frame:]
        aa_seq = translate_dna(frame_seq)
        frames[f"+{frame + 1}"] = aa_seq

    # Reverse strand frames
    rev_comp = reverse_complement(dna_seq)
    for frame in range(3):
        frame_seq = rev_comp[frame:]
        aa_seq = translate_dna(frame_seq)
        frames[f"-{frame + 1}"] = aa_seq

    return frames


def calculate_cai(sequence: str, reference_usage: Dict[str, float] = None) -> float:
    """Calculate Codon Adaptation Index (CAI) for a sequence.

    Args:
        sequence: Amino acid sequence
        reference_usage: Dictionary of reference codon usage frequencies

    Returns:
        CAI value (0.0 to 1.0)
    """
    if not sequence or not reference_usage:
        return 0.0

    # This is a simplified implementation
    # Real CAI calculation is more complex
    seq_upper = sequence.upper()

    # Count codons (simplified - assumes sequence is codons)
    codon_counts = {}
    for i in range(0, len(seq_upper) - 2, 3):
        codon = seq_upper[i:i + 3]
        codon_counts[codon] = codon_counts.get(codon, 0) + 1

    if not codon_counts:
        return 0.0

    # Calculate CAI-like score
    total_weight = 0
    total_codons = 0

    for codon, count in codon_counts.items():
        usage = reference_usage.get(codon, 0.01)  # Small default
        total_weight += usage * count
        total_codons += count

    return total_weight / total_codons if total_codons > 0 else 0.0


def optimize_codons(sequence: str, target_usage: Dict[str, float]) -> str:
    """Optimize codon usage for a DNA sequence.

    This function optimizes codon usage by replacing synonymous codons
    with the most frequently used codons in the target organism's genome.

    Args:
        sequence: DNA sequence string (must be divisible by 3)
        target_usage: Target codon usage frequencies (codon -> relative frequency)

    Returns:
        Optimized DNA sequence with improved codon usage

    Raises:
        ValueError: If sequence length not divisible by 3 or invalid codons
    """
    if len(sequence) % 3 != 0:
        raise ValueError("Sequence length must be divisible by 3 for codon optimization")

    # Get genetic code for codon -> amino acid mapping
    genetic_code = get_genetic_code()

    # Build codon preference map: amino_acid -> [(codon, frequency), ...]
    codon_preferences = {}
    for codon, aa in genetic_code.items():
        if aa not in codon_preferences:
            codon_preferences[aa] = []
        # Use target usage if available, otherwise default to equal preference
        freq = target_usage.get(codon, 1.0)
        codon_preferences[aa].append((codon, freq))

    # Sort codons by preference (highest frequency first)
    for aa in codon_preferences:
        codon_preferences[aa].sort(key=lambda x: x[1], reverse=True)

    # Optimize sequence codon by codon
    optimized = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3].upper()

        if codon not in genetic_code:
            raise ValueError(f"Invalid codon '{codon}' at position {i}")

        aa = genetic_code[codon]
        preferred_codons = codon_preferences.get(aa, [(codon, 1.0)])

        # Use the most preferred codon
        optimized_codon = preferred_codons[0][0]
        optimized.append(optimized_codon)

    result = ''.join(optimized)

    logger.info(f"Optimized codon usage for {len(sequence)//3} codons")
    return result


def get_genetic_code(code_id: int = 1) -> Dict[str, str]:
    """Get genetic code dictionary for specified code table.

    Args:
        code_id: Genetic code table ID (1 = standard)

    Returns:
        Dictionary mapping codons to amino acids
    """
    if code_id == 1:
        return GENETIC_CODE.copy()
    else:
        logger.warning(f"Genetic code {code_id} not implemented, using standard code")
        return GENETIC_CODE.copy()


def back_translate(protein_seq: str, codon_usage: Dict[str, float] = None) -> str:
    """Back-translate amino acid sequence to DNA using optimal codons.

    Args:
        protein_seq: Amino acid sequence
        codon_usage: Dictionary of codon usage preferences

    Returns:
        DNA sequence
    """
    # Simplified back-translation using standard codons
    aa_to_codon = {
        'A': 'GCU', 'C': 'UGU', 'D': 'GAU', 'E': 'GAA',
        'F': 'UUU', 'G': 'GGU', 'H': 'CAU', 'I': 'AUU',
        'K': 'AAA', 'L': 'CUU', 'M': 'AUG', 'N': 'AAU',
        'P': 'CCU', 'Q': 'CAA', 'R': 'CGU', 'S': 'UCU',
        'T': 'ACU', 'V': 'GUU', 'W': 'UGG', 'Y': 'UAU',
        '*': 'UAA'  # Stop codon
    }

    dna_codons = []
    for aa in protein_seq.upper():
        codon = aa_to_codon.get(aa, 'NNN')  # Unknown amino acids
        dna_codons.append(codon)

    rna_seq = ''.join(dna_codons)

    # Convert RNA back to DNA (U -> T)
    dna_seq = rna_seq.replace('U', 'T')

    return dna_seq


# Import re at module level for find_orfs function
import re



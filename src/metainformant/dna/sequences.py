"""DNA sequence processing and analysis utilities.

This module provides core functionality for reading, writing, and analyzing
DNA sequences in FASTA format, including sequence validation, basic statistics,
and sequence manipulation operations.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterator, List, Tuple, Union

from metainformant.core import io, logging

logger = logging.get_logger(__name__)


def read_fasta(path: Union[str, Path]) -> Dict[str, str]:
    """Read sequences from a FASTA file.

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

                if line.startswith('>'):
                    # Save previous sequence
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)

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
            sequences[current_id] = ''.join(current_seq)

    except Exception as e:
        raise ValueError(f"Error reading FASTA file {path}: {e}")

    if not sequences:
        raise ValueError(f"No sequences found in FASTA file: {path}")

    logger.debug(f"Read {len(sequences)} sequences from {path}")
    return sequences


def write_fasta(sequences: Dict[str, str], path: Union[str, Path], line_width: int = 60) -> None:
    """Write sequences to a FASTA file.

    Args:
        sequences: Dictionary mapping sequence IDs to sequence strings
        path: Output file path
        line_width: Maximum characters per line (0 for no wrapping)
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with io.open_text_auto(path, 'w') as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n")

            if line_width > 0:
                # Wrap sequence lines
                for i in range(0, len(sequence), line_width):
                    f.write(sequence[i:i + line_width] + '\n')
            else:
                # Single line per sequence
                f.write(sequence + '\n')

    logger.debug(f"Wrote {len(sequences)} sequences to {path}")


def reverse_complement(seq: str) -> str:
    """Generate the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        Reverse complement sequence

    Raises:
        ValueError: If sequence contains invalid characters
    """
    if not seq:
        return ""

    # Validate sequence
    if not validate_dna_sequence(seq):
        raise ValueError(f"Invalid DNA sequence: {seq}")

    # Complement mapping
    complement = str.maketrans('ATCGatcg', 'TAGCtagc')

    # Reverse and complement
    return seq.translate(complement)[::-1]


def gc_content(seq: str) -> float:
    """Calculate GC content of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        GC content as a fraction (0.0 to 1.0)

    Raises:
        ValueError: If sequence is empty or invalid
    """
    if not seq:
        raise ValueError("Cannot calculate GC content of empty sequence")

    if not validate_dna_sequence(seq):
        raise ValueError(f"Invalid DNA sequence: {seq}")

    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    total_count = len(seq)

    return gc_count / total_count if total_count > 0 else 0.0


def sequence_length(seq: str) -> int:
    """Get the length of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        Sequence length
    """
    return len(seq)


def validate_dna_sequence(seq: str) -> bool:
    """Validate that a sequence contains only valid DNA characters.

    Args:
        seq: Sequence to validate

    Returns:
        True if sequence is valid DNA, False otherwise
    """
    if not seq:
        return False

    # Allow IUPAC ambiguity codes
    valid_chars = set('ATCGNUWSMKRYBDHVatcgnuwsmkrybdhv-')
    return all(c in valid_chars for c in seq)


def find_motifs(seq: str, motif_patterns: List[str]) -> Dict[str, List[int]]:
    """Find occurrences of motifs in a DNA sequence.

    Args:
        seq: DNA sequence string
        motif_patterns: List of motif patterns (can include IUPAC codes)

    Returns:
        Dictionary mapping motif patterns to lists of start positions
    """
    results = {}

    for pattern in motif_patterns:
        positions = []
        # Convert IUPAC to regex
        regex_pattern = _iupac_to_regex(pattern.upper())
        for match in re.finditer(regex_pattern, seq.upper()):
            positions.append(match.start())
        results[pattern] = positions

    return results


def find_repeats(seq: str, min_length: int = 3) -> Dict[str, List[int]]:
    """Find repeated sequences in DNA.

    Args:
        seq: DNA sequence string
        min_length: Minimum repeat length

    Returns:
        Dictionary mapping repeat sequences to lists of start positions
    """
    repeats = {}

    # Simple repeat finding - look for exact repeats
    seq_upper = seq.upper()
    for i in range(len(seq_upper) - min_length + 1):
        repeat_seq = seq_upper[i:i + min_length]
        if repeat_seq in seq_upper[i + min_length:]:
            if repeat_seq not in repeats:
                repeats[repeat_seq] = []
            repeats[repeat_seq].append(i)

    return repeats


def calculate_sequence_complexity(seq: str) -> float:
    """Calculate sequence complexity using Shannon entropy.

    Args:
        seq: DNA sequence string

    Returns:
        Sequence complexity score (0.0 to 1.0)
    """
    if not seq:
        return 0.0

    from collections import Counter
    import math

    # Count nucleotide frequencies
    counts = Counter(seq.upper())
    total = sum(counts.values())

    # Calculate Shannon entropy
    entropy = 0.0
    for count in counts.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log2(p)

    # Normalize by maximum possible entropy (2 bits for DNA)
    max_entropy = 2.0
    complexity = entropy / max_entropy if max_entropy > 0 else 0.0

    return min(complexity, 1.0)  # Cap at 1.0


def find_orfs(seq: str, min_length: int = 30) -> List[Tuple[int, int, str]]:
    """Find open reading frames in a DNA sequence.

    Args:
        seq: DNA sequence string
        min_length: Minimum ORF length in amino acids

    Returns:
        List of tuples (start_pos, end_pos, frame)
    """
    orfs = []
    seq_upper = seq.upper()

    # Check all 6 reading frames
    for frame in range(6):
        if frame < 3:
            # Forward strand
            search_seq = seq_upper[frame:]
        else:
            # Reverse strand
            search_seq = reverse_complement(seq_upper)[frame - 3:]

        # Find start codons
        start_positions = []
        for match in re.finditer(r'ATG', search_seq):
            start_positions.append(match.start())

        for start_pos in start_positions:
            # Find next stop codon
            orf_seq = search_seq[start_pos:]
            stop_found = False

            for stop_match in re.finditer(r'(TAA|TAG|TGA)', orf_seq[3:], re.IGNORECASE):
                stop_pos = start_pos + stop_match.start() + 3  # Include stop codon
                orf_length_aa = (stop_pos - start_pos) // 3

                if orf_length_aa >= min_length:
                    frame_label = f"{'+' if frame < 3 else '-'}{frame % 3 + 1}"
                    orfs.append((start_pos, stop_pos, frame_label))
                stop_found = True
                break

            # If no stop codon found, check if ORF extends to end
            if not stop_found:
                orf_length_aa = len(orf_seq) // 3
                if orf_length_aa >= min_length:
                    frame_label = f"{'+' if frame < 3 else '-'}{frame % 3 + 1}"
                    orfs.append((start_pos, len(search_seq), frame_label))

    return orfs


def find_start_codons(seq: str) -> List[int]:
    """Find positions of start codons (ATG) in DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        List of start positions
    """
    positions = []
    for match in re.finditer(r'ATG', seq.upper()):
        positions.append(match.start())
    return positions


def find_stop_codons(seq: str) -> List[int]:
    """Find positions of stop codons (TAA, TAG, TGA) in DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        List of start positions
    """
    positions = []
    for match in re.finditer(r'(TAA|TAG|TGA)', seq.upper()):
        positions.append(match.start())
    return positions


def calculate_sequence_entropy(seq: str, k: int = 1) -> float:
    """Calculate sequence entropy using k-mer frequencies.

    Args:
        seq: DNA sequence string
        k: k-mer size

    Returns:
        Sequence entropy score
    """
    if len(seq) < k:
        return 0.0

    from collections import Counter
    import math

    # Count k-mers
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    counts = Counter(kmers)
    total = sum(counts.values())

    # Calculate entropy
    entropy = 0.0
    for count in counts.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log2(p)

    return entropy


def detect_sequence_bias(seq: str) -> Dict[str, float]:
    """Detect nucleotide composition biases in DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        Dictionary with bias metrics for each nucleotide
    """
    if not seq:
        return {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0}

    seq_upper = seq.upper()
    length = len(seq_upper)

    counts = {
        'A': seq_upper.count('A') / length,
        'T': seq_upper.count('T') / length,
        'G': seq_upper.count('G') / length,
        'C': seq_upper.count('C') / length,
    }

    # Expected frequency for no bias
    expected = 0.25

    # Calculate bias as deviation from expected
    biases = {}
    for nuc, freq in counts.items():
        biases[nuc] = abs(freq - expected)

    return biases


def calculate_gc_skew(seq: str) -> float:
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


def calculate_at_skew(seq: str) -> float:
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


def find_palindromes(seq: str, min_length: int = 4) -> List[Tuple[str, int, int]]:
    """Find palindromic sequences in DNA.

    Args:
        seq: DNA sequence string
        min_length: Minimum palindrome length

    Returns:
        List of tuples (palindrome_seq, start_pos, end_pos)
    """
    palindromes = []
    seq_upper = seq.upper()

    for i in range(len(seq_upper) - min_length + 1):
        for j in range(min_length, len(seq_upper) - i + 1):
            substring = seq_upper[i:i + j]
            if _is_palindromic(substring):
                palindromes.append((substring, i, i + j - 1))

    return palindromes


def calculate_melting_temperature(seq: str, method: str = "wallace") -> float:
    """Calculate DNA melting temperature.

    Args:
        seq: DNA sequence string
        method: Calculation method ('wallace' or 'gc')

    Returns:
        Melting temperature in Celsius

    Raises:
        ValueError: If method is not supported
    """
    if not seq:
        raise ValueError("Cannot calculate melting temperature of empty sequence")

    if method.lower() == "wallace":
        # Wallace rule: Tm = 4(G+C) + 2(A+T)
        gc_count = seq.upper().count('G') + seq.upper().count('C')
        at_count = seq.upper().count('A') + seq.upper().count('T')
        return 4 * gc_count + 2 * at_count
    elif method.lower() == "gc":
        # GC content based: Tm = 64.9 + 41*(G+C-16.4)/N
        gc_count = seq.upper().count('G') + seq.upper().count('C')
        length = len(seq)
        if length == 0:
            return 0.0
        return 64.9 + 41 * (gc_count - 16.4) / length
    else:
        raise ValueError(f"Unsupported method: {method}. Use 'wallace' or 'gc'")


def calculate_codon_usage(seq: str) -> Dict[str, float]:
    """Calculate codon usage frequencies.

    Args:
        seq: DNA sequence string (must be divisible by 3)

    Returns:
        Dictionary mapping codons to frequency
    """
    if len(seq) % 3 != 0:
        logger.warning(f"Sequence length {len(seq)} is not divisible by 3")

    from collections import Counter

    codons = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3].upper()
        codons.append(codon)

    counts = Counter(codons)
    total = sum(counts.values())

    frequencies = {}
    for codon, count in counts.items():
        frequencies[codon] = count / total if total > 0 else 0.0

    return frequencies


def dna_complementarity_score(seq1: str, seq2: str) -> float:
    """Calculate complementarity score between two DNA sequences.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence

    Returns:
        Complementarity score (0.0 to 1.0)
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be the same length")

    complement_map = str.maketrans('ATCG', 'TAGC')
    seq2_complement = seq2.upper().translate(complement_map)

    matches = sum(1 for a, b in zip(seq1.upper(), seq2_complement) if a == b)
    return matches / len(seq1) if seq1 else 0.0


def _iupac_to_regex(pattern: str) -> str:
    """Convert IUPAC ambiguity codes to regex pattern.

    Args:
        pattern: IUPAC pattern string

    Returns:
        Regular expression pattern
    """
    iupac_codes = {
        'N': '[ATCG]',
        'U': 'T',  # RNA U becomes T in DNA
        'W': '[AT]',
        'S': '[CG]',
        'M': '[AC]',
        'K': '[GT]',
        'R': '[AG]',
        'Y': '[CT]',
        'B': '[CGT]',
        'D': '[AGT]',
        'H': '[ACT]',
        'V': '[ACG]',
    }

    regex = ''
    for char in pattern:
        regex += iupac_codes.get(char.upper(), char.upper())

    return regex


def _is_palindromic(seq: str) -> bool:
    """Check if a DNA sequence is palindromic.

    Args:
        seq: DNA sequence string

    Returns:
        True if sequence is palindromic
    """
    complement = str.maketrans('ATCG', 'TAGC')
    return seq == seq.translate(complement)[::-1]



from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Dict

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta(path: Path) -> Dict[str, str]:
    """Parse FASTA file and return dictionary of sequence ID to sequence.
    
    Args:
        path: Path to FASTA file
        
    Returns:
        Dictionary mapping sequence IDs to sequences
    """
    records: Dict[str, str] = {}
    current_id: str | None = None
    current_seq_parts: list[str] = []
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                records[current_id] = "".join(current_seq_parts)
            current_id = line[1:].split()[0]
            current_seq_parts = []
        else:
            current_seq_parts.append(line)
    if current_id is not None:
        records[current_id] = "".join(current_seq_parts)
    return records


def is_valid_protein_sequence(seq: str) -> bool:
    """Check if sequence contains only valid amino acid characters.
    
    Args:
        seq: Protein sequence string
        
    Returns:
        True if sequence contains only valid amino acids, False otherwise
    """
    if not seq:
        return False
    return all((c.isalpha() and c.upper() in VALID_AA) for c in seq)


def calculate_aa_composition(seq: str) -> Dict[str, float]:
    """Calculate amino acid composition (frequencies) of a protein sequence.
    
    Args:
        seq: Protein sequence string
        
    Returns:
        Dictionary mapping amino acid codes to their frequencies (0.0 to 1.0)
    """
    seq_upper = [c for c in seq.upper() if c in VALID_AA]
    n = len(seq_upper)
    if n == 0:
        return {aa: 0.0 for aa in sorted(VALID_AA)}
    counts = Counter(seq_upper)
    return {aa: counts.get(aa, 0) / n for aa in sorted(VALID_AA)}


def kmer_frequencies(seq: str, *, k: int) -> Dict[str, int]:
    """Calculate k-mer frequencies in a protein sequence.
    
    Args:
        seq: Protein sequence string
        k: K-mer length (must be positive)
        
    Returns:
        Dictionary mapping k-mers to their counts
        
    Raises:
        ValueError: If k is not positive
    """
    if k <= 0:
        raise ValueError("k must be positive")
    seq = seq.upper()
    freq: Dict[str, int] = {}
    for i in range(0, max(0, len(seq) - k + 1)):
        kmer = seq[i : i + k]
        if any(c not in VALID_AA for c in kmer):
            continue
        freq[kmer] = freq.get(kmer, 0) + 1
    return freq

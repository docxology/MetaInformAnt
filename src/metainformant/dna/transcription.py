from __future__ import annotations

_DNA_TO_RNA = str.maketrans(
    {
        "A": "A",
        "T": "U",
        "G": "G",
        "C": "C",
        "a": "a",
        "t": "u",
        "g": "g",
        "c": "c",
    }
)

_RNA_TO_DNA = str.maketrans(
    {
        "A": "A",
        "U": "T",
        "G": "G",
        "C": "C",
        "a": "a",
        "u": "t",
        "g": "g",
        "c": "c",
    }
)


def transcribe_dna_to_rna(seq: str) -> str:
    """Transcribe DNA to RNA (T->U). Preserves case; non-ATGC characters pass through."""
    return (seq or "").translate(_DNA_TO_RNA)


def reverse_transcribe_rna_to_dna(seq: str) -> str:
    """Reverse transcribe RNA to DNA (U->T). Preserves case; non-AUGC characters pass through."""
    return (seq or "").translate(_RNA_TO_DNA)

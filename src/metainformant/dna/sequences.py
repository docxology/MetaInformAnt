from __future__ import annotations

from typing import Dict

from Bio import SeqIO


def read_fasta(path: str) -> Dict[str, str]:
    """Read a FASTA file into a dictionary of id -> sequence string.

    The sequence IDs are taken from the FASTA record identifiers up to the first whitespace.
    """
    records = SeqIO.parse(path, "fasta")
    seqs: Dict[str, str] = {}
    for rec in records:
        identifier = rec.id
        seqs[identifier] = str(rec.seq)
    return seqs



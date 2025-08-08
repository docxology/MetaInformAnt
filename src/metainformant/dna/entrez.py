from __future__ import annotations

from typing import Any

from Bio import Entrez, SeqIO


def get_genome_from_ncbi(genome_id: str, *, email: str) -> Any:
    Entrez.email = email
    handle = Entrez.efetch(db="nuccore", id=genome_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "fasta")
    return record



from __future__ import annotations

from typing import Any

from Bio import Entrez, SeqIO


def get_genome_from_ncbi(genome_id: str, *, email: str) -> Any:
    """Fetch a sequence record from NCBI nuccore as FASTA.

    Parameters
    - genome_id: e.g., "NC_001422.1" (PhiX174)
    - email: contact email required by NCBI Entrez
    """
    Entrez.email = email
    with Entrez.efetch(db="nuccore", id=genome_id, rettype="fasta", retmode="text") as handle:
        record = SeqIO.read(handle, "fasta")
    return record

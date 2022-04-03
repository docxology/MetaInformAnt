from Bio import Entrez, SeqIO
import logging
from utils import get_db_client, get_entrez_options
from typing import Optional


def get_genome_from_ncbi(genome_id: str):
    entrez_ops = get_entrez_options()
    # Entrez.email =
    handle = Entrez.efetch(
        email=entrez_ops["email"],
        db=entrez_ops["db"],
        id=genome_id,
        rettype=entrez_ops["rettype"],
        retmode=entrez_ops["retmode"],
    )
    logging.debug(handle.read())
    record = SeqIO.read(handle, "fasta")
    return record

"""Accepts a genome ID, fetches it from NCBI and pushes it to the db + file storage"""
import sys
from dataclasses import dataclass
from Bio import Entrez, SeqIO
import logging
from utils import get_db_client
from typing import Optional

logging.getLogger().setLevel(logging.DEBUG)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


@dataclass
class OrganismGenome:
    """Class for tracking the genome object"""

    genome_id: str
    file_path: str
    file_contents: str = None


def get_genome_from_ncbi(genome_id: str, handle):
    record = SeqIO.read(handle, "fasta")
    return record


def check_db_if_exists(genome_id: str, db_cursor) -> (None | OrganismGenome):
    """Checks db for genome id and if it exists calls query to create genome obj."""
    logging.debug("Checking if genome exists...")
    genome_id = genome_id.strip().lower()
    logging.debug(
        f"SELECT count(*) from metainformant.genomes where genome_id={genome_id}"
    )
    db_cursor.execute(
        f"SELECT count(*) from metainformant.genomes where genome_id={genome_id}"
    )
    row_count = db_cursor.fetchone()
    logging.debug(row_count[0])

    if row_count == "0":
        return None

    db_cursor.execute(
        f"SELECT file_path, ncbi_id from metainformant.genomes where genome_id={genome_id} order by created_on desc limit 1"
    )
    # logging.debug(db_cursor.fetchone()[0])
    file_path = db_cursor.fetchone()[0]
    genome_obj = OrganismGenome(genome_id=genome_id, file_path=file_path)
    return genome_obj


def upload_file_to_local_dir(OrganismGenome) -> bool:
    pass


def upload_file_to_gcs(OrganismGenome) -> bool:
    pass


def add_record_to_db(OrganismGenome, postgresdb) -> bool:
    pass


def get_genome(genome_id) -> OrganismGenome:
    Entrez.email = "admin@metainformant.com"  # parameterize?
    # check db first, if exists fetch object
    handle = Entrez.efetch(db="nuccore", id=genome_id, rettype="gb", retmode="text")
    db_cursor = get_db_client()

    genome_obj = check_db_if_exists(genome_id, db_cursor)

    if not genome_obj:
        pass
        # otherwise generate new object
    return genome_obj


get_genome("6469")

"""Accepts a genome ID, fetches it from NCBI and pushes it to the db + file storage"""
import sys
from dataclasses import dataclass

@dataclass
class OrganismGenome:
    """Class for tracking the genome object"""
    genome_id: str
    file_path: str
    file_contents: str

def get_genome_from_ncbi(genome_id: str) -> OrganismGenome:
    pass


def upload_file_to_local_dir(OrganismGenome) -> bool:
    pass

def upload_file_to_gcs(OrganismGenome) -> bool:
    pass

def add_record_to_db(OrganismGenome, postgresdb) -> bool:
    pass
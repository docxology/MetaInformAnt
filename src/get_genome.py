"""Accepts a genome ID, fetches it from NCBI and pushes it to the db + file storage"""
import sys
from dataclasses import dataclass
from Bio import Entrez, SeqIO
 

@dataclass
class OrganismGenome:
    """Class for tracking the genome object"""
    genome_id: str
    file_path: str
    file_contents: str

def get_genome_from_ncbi(genome_id: str):
    Entrez.email =  "admin@metainformant.com" #parameterize?
    handle = Entrez.efetch(db="nuccore", id=genome_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "fasta")
    return record

def check_db_if_exists(genome_id: str) -> (None | OrganismGenome):
    return None

def upload_file_to_local_dir(OrganismGenome) -> bool:
    pass

def upload_file_to_gcs(OrganismGenome) -> bool:
    pass

def add_record_to_db(OrganismGenome, postgresdb) -> bool:
    pass

def get_genome(genome_id) -> OrganismGenome:
    # check db first, if exists fetch object
    genome_obj = check_db_if_exists(genome_id)

    if not genome_obj:
        pass
        # otherwise generate new object
    return genome_obj
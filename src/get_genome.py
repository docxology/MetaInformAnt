"""Accepts a genome ID, fetches it from NCBI and pushes it to the db + file storage"""
import email
import sys
from Bio import Entrez, SeqIO
import logging
from utils import get_db_client, get_entrez_options
from typing import Optional
from ncbi_datasets_helper import (
    download_genome_data_package,
    get_metadata_by_single_accesion,
)
from base_classes import NCBIGenome, OrganismGenome
from utils import save_file_to_local_dir

logging.getLogger().setLevel(logging.DEBUG)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def fetch_genome(search_term, search_term_type):
    # genome_organism = download_genome_data_package([search_term], "test.zip")
    # md = get_metadata_by_single_accesion([search_term])
    new_genome = NCBIGenome(search_term=search_term, search_term_type=search_term_type)
    new_genome.search_router()
    logging.debug(new_genome)
    return new_genome

    # return md


def check_db_if_exists(
    search_term: str, db_cursor, search_field="assembly_accession_id"
) -> (None | OrganismGenome):
    """Checks db for genome id and if it exists calls query to create genome obj."""
    logging.debug("Checking if genome exists...")
    search_term = search_term.strip().lower()
    db_cursor.execute(
        "SELECT count(*) from metainformant.genome where %s=%s",
        (search_field, search_term),
    )

    row_count = db_cursor.fetchone()[0]
    logging.debug("Count: %d", row_count)

    if row_count == 0:
        return None

    db_cursor.execute(
        "SELECT assembly_accession_id, file_path from metainformant.genome where %s=%s order by created_on desc limit 1",
        (search_field, search_term),
    )

    # db_cursor.execute(
    #     f'SELECT genome_id, file_path, ncbi_id from metainformant.genomes where {search_field}="{search_term}" order by created_on desc limit 1'
    # )
    # logging.debug(db_cursor.fetchone()[0])
    genome_metadata = db_cursor.fetchone()
    genome_obj = OrganismGenome(
        genome_id=genome_metadata[0], file_path=genome_metadata[1]
    )
    return genome_obj


def check_for_taxonomic_name(NCBIGenome):
    pass


def add_record_to_db(NCBIGenome) -> bool:
    conn, db_cursor = get_db_client()
    db_cursor.execute(
        "INSERT INTO metainformant.taxonomic_names (tx_id, genus, epithet) VALUES(%s, %s, %s)",
        ("test", "test", "Stest"),
    )
    conn.commit()  # <- We MUST commit to reflect the inserted data
    db_cursor.close()
    conn.close()

    pass


def get_genome(search_term, search_term_type="assembly_accession_id") -> OrganismGenome:
    _, db_cursor = get_db_client()
    # check db first, if exists fetch object
    genome_obj = check_db_if_exists(search_term, db_cursor, search_term_type)
    logging.info(genome_obj)
    if not genome_obj:
        new_genome = fetch_genome(search_term, search_term_type)
        add_record_to_db(new_genome)
        pass
        # otherwise generate new object
    return genome_obj


# get_genome("GCA_000001215.4")

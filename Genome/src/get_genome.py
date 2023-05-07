"""Accepts a genome ID, fetches it from NCBI and pushes it to the db + file storage"""
import email
import json
import sys
from Bio import Entrez, SeqIO
import json
import logging
from utils import get_db_client, get_entrez_options
from typing import Optional
from ncbi_datasets_helper import (
    download_genome_data_package,
    get_metadata_by_single_accesion,
)
from base_classes import NCBIGenome, OrganismGenome
from utils import save_file_to_local_dir


def fetch_genome(search_term, search_term_type):
    new_genome = NCBIGenome(search_term=search_term, search_term_type=search_term_type)
    new_genome.search_router()
    return new_genome

    # return md


def check_db_if_exists(
    search_term: str, db_cursor, search_field="assembly_accession_id"
) -> (None | OrganismGenome):
    """Checks db for genome id and if it exists calls query to create genome obj."""
    logging.debug("Checking if genome exists...")
    db_cursor.execute(
        f"SELECT count(*) from metainformant.genome where {search_field}='{search_term}'",
    )
    logging.debug(db_cursor.query)
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
    # TO DO: recreate NCBI object here
    return 1


def get_taxonomic_name(genome: NCBIGenome, conn, db_cursor):
    # Check for existence
    db_cursor.execute(
        f"SELECT count(*) from metainformant.taxonomic_names where tx_id='{genome.tx_id}'",
    )
    row_count = db_cursor.fetchone()[0]
    if row_count == 0:
        # Insert the taxonomic name
        db_cursor.execute(
            "INSERT INTO metainformant.taxonomic_names (tx_id, genus, epithet) VALUES(%s, %s, %s) ON CONFLICT DO NOTHING",
            (genome.tx_id, genome.taxon_name, genome.epithet),
        )
        conn.commit()
    db_cursor.execute(
        f"SELECT id FROM metainformant.taxonomic_names where tx_id = '{genome.tx_id}'",
    )
    return db_cursor.fetchone()[0]


def add_record_to_db(genome: NCBIGenome) -> bool:
    conn, db_cursor = get_db_client()

    # fetch taxonomic_name id
    taxonomix_row_id = get_taxonomic_name(genome, conn, db_cursor)
    logging.debug(taxonomix_row_id)

    # Insert genome into table
    db_cursor.execute(
        "INSERT INTO metainformant.genome (assembly_accession_id, taxon_name, ncbi_metadata, file_path) VALUES(%s, %s, %s, %s)",
        (
            genome.accesion_id,
            taxonomix_row_id,
            json.dumps(genome.assembly_metadata.to_dict()),
            genome.filename,
        ),
    )
    conn.commit()
    db_cursor.close()
    conn.close()

    pass


def get_genome(search_term, search_term_type="assembly_accession_id") -> OrganismGenome:
    _, db_cursor = get_db_client()
    # check db first, if exists fetch object
    genome_obj = check_db_if_exists(search_term, db_cursor, search_term_type)
    if not genome_obj:
        new_genome = fetch_genome(search_term, search_term_type)
        add_record_to_db(new_genome)
        pass
        # otherwise generate new object
    return genome_obj


# get_genome("GCA_000001215.4")

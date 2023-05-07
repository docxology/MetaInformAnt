""" The command line arguments:
- pass a txt file (eg genome.txt a simple list of genome ids)
- get_genome.py

- acccesion id {id} not found
- 89 files found, 80 downloaded, 17 does not exist for text files

test command: "python src/main.py --get_single_genome_by_id=GCA_000001215.4"
python src/main.py --file test_data/test_multiple_ids.txt
"""
import argparse
from typing import Set
from get_genome import get_genome
import logging as log
from utils import GenomeFetchFailure, load_file_to_set
import sys
from ncbi_datasets_helper import get_accession_by_tax_id

log.getLogger().setLevel(log.INFO)
log.getLogger().addHandler(log.StreamHandler(sys.stdout))

# TODO: Batch requests to ncbi api
def process_ids(genome_ids: Set[str]):
    failed_genomes = []
    succeded_genomes = []
    log.info("Starting to process %s genome ids", len(genome_ids))
    # TODO: if there is more than one i should really batch them
    for idx, genome_id in enumerate(genome_ids):
        try:
            log.info("Processing %s of %s", idx + 1, len(genome_ids))
            get_genome(genome_id.strip())
            succeded_genomes.append(genome_id)
        except GenomeFetchFailure:
            failed_genomes.append(genome_id)
            pass
    log.info(
        "PROCESSING OF ITEMS FINISHED\nOf the %s genomes submitted, %s succeded and %s failed.",
        len(genome_ids),
        len(succeded_genomes),
        len(failed_genomes),
    )
    log.info("The following ids failed to be retrieved from NCBI: %s", failed_genomes)


def fetch_and_handle_ids_for_species(tax_id: str):
    list_genomes = get_accession_by_tax_id(tax_id)
    str_prompt = (
        f"Fetched {len(list_genomes)} accession ids for {tax_id}. Download all? y/n \n"
    )
    if input(str_prompt) != "y":
        exit()
    process_ids(list_genomes)


def route_query(args):
    if args.get_single_genome_by_id:
        process_ids([args.get_single_genome_by_id])
    if args.list_genome_ids_by_tax_id:
        fetch_and_handle_ids_for_species(args.list_genome_ids_by_tax_id)
    if args.file:
        genome_ids = load_file_to_set(args.file)
        process_ids(genome_ids)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a genome")
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--get_single_genome_by_id", help="gets a single genome by its id")
    g.add_argument("--file", help="gets multiple genomes by its id")
    g.add_argument(
        "--list_genome_ids_by_tax_id",
        help="returns a list of genomes for the species id provided",
    )

    args = parser.parse_args()

    route_query(args)

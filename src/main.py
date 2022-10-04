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

log.getLogger().setLevel(log.INFO)
log.getLogger().addHandler(log.StreamHandler(sys.stdout))

# TODO: Batch requests to ncbi api
def process_ids(genome_ids: Set[str]):
    failed_genomes = []
    succeded_genomes = []
    log.info("Starting to process %s genome ids", len(genome_ids))
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


def route_query(args):
    if args.get_single_genome_by_id:
        process_ids([args.get_single_genome_by_id])
    if args.file:
        genome_ids = load_file_to_set(args.file)
        process_ids(genome_ids)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a genome")
    parser.add_argument(
        "--get_single_genome_by_id", help="gets a single genome by its id"
    )
    parser.add_argument("--file", help="gets multiple genomes by its id")

    args = parser.parse_args()

    route_query(args)

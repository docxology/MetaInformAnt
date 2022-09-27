""" The command line arguments:
- pass a txt file (eg genome.txt a simple list of genome ids)
- get_genome.py

- acccesion id {id} not found
- 89 files found, 80 downloaded, 17 does not exist for text files

test command: "python src/main.py --get_single_genome_by_id=GCA_000001215.4"
"""
import argparse
from get_genome import get_genome


def route_query(args):
    print(args.get_single_genome_by_id)
    get_genome((args.get_single_genome_by_id).strip())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a genome")
    parser.add_argument(
        "--get_single_genome_by_id", help="gets a single genome by its id"
    )

    args = parser.parse_args()

    route_query(args)

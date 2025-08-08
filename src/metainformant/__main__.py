import argparse
import sys

from .tests import run_all_tests
from .dna.genomes import is_valid_assembly_accession

def main() -> None:
    parser = argparse.ArgumentParser(prog="metainformant", description="METAINFORMANT CLI")
    subparsers = parser.add_subparsers(dest="command")

    # dna subcommand
    dna_parser = subparsers.add_parser("dna", help="DNA-related operations")
    dna_sub = dna_parser.add_subparsers(dest="dna_cmd")
    fetch = dna_sub.add_parser("fetch", help="Fetch genome by assembly accession")
    fetch.add_argument("--assembly", required=True, help="NCBI assembly accession, e.g., GCF_*")

    # tests subcommand
    tests_parser = subparsers.add_parser("tests", help="Run repository test suite")
    tests_parser.add_argument("pytest_args", nargs=argparse.REMAINDER, help="Arguments passed to pytest")

    args = parser.parse_args()

    if args.command == "dna" and args.dna_cmd == "fetch":
        if not is_valid_assembly_accession(args.assembly):
            print(f"Invalid assembly accession: {args.assembly}")
            sys.exit(2)
        print(f"Validated assembly accession: {args.assembly}")
        # Future: call into actual fetch workflow
        return

    if args.command == "tests":
        exit_code = run_all_tests(args.pytest_args)
        sys.exit(exit_code)

    parser.print_help()


if __name__ == "__main__":
    main()



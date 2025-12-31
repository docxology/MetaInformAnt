"""CLI entry point for METAINFORMANT.

This module provides the command-line interface for the METAINFORMANT
bioinformatics toolkit. Run with --help for usage information.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__


def main() -> int:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="metainformant",
        description="METAINFORMANT: Comprehensive Bioinformatics Toolkit for Multi-Omic Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  metainformant --version
  metainformant --help

For detailed usage of specific modules, import them directly in Python:
  python -c "from metainformant.dna import sequences; help(sequences.read_fasta)"
        """,
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"METAINFORMANT {__version__}",
        help="Show version information and exit",
    )

    parser.add_argument(
        "--modules",
        action="store_true",
        help="List available modules",
    )

    args = parser.parse_args()

    if args.modules:
        _list_modules()
        return 0

    # If no arguments provided, show help
    if len(sys.argv) == 1:
        parser.print_help()
        return 0

    return 0


def _list_modules() -> None:
    """List all available modules."""
    modules = [
        ("core", "Shared utilities and infrastructure"),
        ("dna", "DNA sequence analysis and genomics"),
        ("rna", "RNA-seq workflows and amalgkit integration"),
        ("protein", "Protein sequence and structure analysis"),
        ("gwas", "Genome-wide association studies"),
        ("math", "Mathematical biology and theoretical modeling"),
        ("information", "Information-theoretic analysis"),
        ("life_events", "Life course and temporal analysis"),
        ("visualization", "Plotting and visualization tools"),
        ("networks", "Biological network analysis"),
        ("multiomics", "Cross-omics data integration"),
        ("singlecell", "Single-cell RNA-seq analysis"),
        ("simulation", "Synthetic data generation"),
        ("quality", "Data quality control"),
        ("ml", "Machine learning for biological data"),
        ("ontology", "Gene ontology and functional annotation"),
        ("phenotype", "Phenotypic trait analysis"),
        ("ecology", "Ecological and community analysis"),
        ("epigenome", "Epigenomic data analysis"),
    ]

    print("Available METAINFORMANT modules:")
    print("=" * 50)

    for name, description in modules:
        print(f"    {description}")

    print("\nImport modules in Python:")
    print("  from metainformant import dna, rna, protein  # etc.")


if __name__ == "__main__":
    sys.exit(main())

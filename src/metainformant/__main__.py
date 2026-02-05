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

    subparsers = parser.add_subparsers(dest="command")

    # Protein subcommands
    protein_parser = subparsers.add_parser("protein", help="Protein analysis commands")
    protein_sub = protein_parser.add_subparsers(dest="protein_command")

    # protein taxon-ids
    taxon_parser = protein_sub.add_parser("taxon-ids", help="Read and validate taxon IDs")
    taxon_parser.add_argument("--file", required=True, help="Path to taxon ID file")

    # protein comp
    comp_parser = protein_sub.add_parser("comp", help="Amino acid composition from FASTA")
    comp_parser.add_argument("--fasta", required=True, help="Path to FASTA file")

    # protein rmsd-ca
    rmsd_parser = protein_sub.add_parser("rmsd-ca", help="RMSD between CA atoms of two PDB files")
    rmsd_parser.add_argument("--pdb-a", required=True, help="Path to first PDB file")
    rmsd_parser.add_argument("--pdb-b", required=True, help="Path to second PDB file")

    args = parser.parse_args()

    if args.modules:
        _list_modules()
        return 0

    if args.command == "protein":
        return _handle_protein(args)

    # If no arguments provided, show help
    if len(sys.argv) == 1:
        parser.print_help()
        return 0

    return 0


def _handle_protein(args: argparse.Namespace) -> int:
    """Handle protein subcommands."""
    import numpy as np

    cmd = args.protein_command

    if cmd == "taxon-ids":
        from .protein.sequence.proteomes import read_taxon_ids

        ids = read_taxon_ids(Path(args.file))
        print(" ".join(ids))
        return 0

    elif cmd == "comp":
        from .protein.sequence.sequences import amino_acid_composition, read_fasta

        sequences = read_fasta(Path(args.fasta))
        for name, seq in sequences.items():
            comp = amino_acid_composition(seq)
            parts = [f"{aa}:{frac:.4f}" for aa, frac in sorted(comp.items()) if frac > 0]
            print(f"{name}\t{','.join(parts)}")
        return 0

    elif cmd == "rmsd-ca":
        from .protein.structure.general import compute_rmsd_kabsch
        from .protein.structure.io import read_pdb_ca_coordinates

        ca_a = read_pdb_ca_coordinates(Path(args.pdb_a))
        ca_b = read_pdb_ca_coordinates(Path(args.pdb_b))
        rmsd = compute_rmsd_kabsch(np.array(ca_a), np.array(ca_b))
        print(f"{rmsd:.6f}")
        return 0

    return 1


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
        ("menu", "Interactive menu and discovery system"),
    ]

    print("Available METAINFORMANT modules:")
    print("=" * 50)

    for name, description in modules:
        print(f"  {name:15} - {description}")

    print("\nImport modules in Python:")
    print("  from metainformant import dna, rna, protein  # etc.")


if __name__ == "__main__":
    sys.exit(main())

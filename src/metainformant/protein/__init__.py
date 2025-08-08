"""Protein domain functionality."""

from .proteomes import read_taxon_ids
from .sequences import parse_fasta, is_valid_protein_sequence, calculate_aa_composition, kmer_frequencies
from .pdb import fetch_pdb_structure
from .uniprot import map_ids_uniprot
from .secondary import simple_helix_coil_propensity
from .structure import compute_rmsd_kabsch

__all__ = [
    "read_taxon_ids",
    "parse_fasta",
    "is_valid_protein_sequence",
    "calculate_aa_composition",
    "kmer_frequencies",
    "fetch_pdb_structure",
    "map_ids_uniprot",
    "simple_helix_coil_propensity",
    "compute_rmsd_kabsch",
]

"""Protein domain functionality."""

from .pdb import fetch_pdb_structure
from .proteomes import read_taxon_ids
from .secondary import simple_helix_coil_propensity
from .sequences import calculate_aa_composition, is_valid_protein_sequence, kmer_frequencies, parse_fasta
from .structure import compute_rmsd_kabsch
from .uniprot import map_ids_uniprot

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

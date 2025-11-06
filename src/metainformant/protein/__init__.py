"""Protein domain functionality.

This module provides comprehensive tools for protein sequence and structure analysis,
including database integration, structure prediction, and functional annotation.
"""

from .alphafold import build_alphafold_url, fetch_alphafold_model
from .alignment import needleman_wunsch, pairwise_identity
from .contacts import compute_ca_contact_pairs
from .interpro import fetch_interpro_domains
from .pdb import fetch_pdb_structure
from .proteomes import read_taxon_ids
from .secondary import simple_helix_coil_propensity
from .sequences import calculate_aa_composition, is_valid_protein_sequence, kmer_frequencies, parse_fasta
from .structure import compute_rmsd_kabsch
from .structure_io import read_pdb_ca_coordinates

try:
    from .structure_analysis import (
        analyze_post_translational_modifications,
        analyze_protein_stability,
        identify_domains,
        predict_protein_family,
        predict_secondary_structure,
    )
    _structure_analysis_available = True
except ImportError:
    _structure_analysis_available = False

from .uniprot import fetch_uniprot_fasta, map_ids_uniprot

__all__ = [
    # Proteome utilities
    "read_taxon_ids",
    # Sequence analysis
    "parse_fasta",
    "is_valid_protein_sequence",
    "calculate_aa_composition",
    "kmer_frequencies",
    # Sequence alignment
    "pairwise_identity",
    "needleman_wunsch",
    # Structure retrieval
    "fetch_pdb_structure",
    "read_pdb_ca_coordinates",
    "compute_rmsd_kabsch",
    "compute_ca_contact_pairs",
    # Structure prediction
    "build_alphafold_url",
    "fetch_alphafold_model",
    # Database integration
    "map_ids_uniprot",
    "fetch_uniprot_fasta",
    "fetch_interpro_domains",
    # Secondary structure
    "simple_helix_coil_propensity",
]

if _structure_analysis_available:
    __all__.extend([
        "predict_secondary_structure",
        "identify_domains",
        "analyze_protein_stability",
        "predict_protein_family",
        "analyze_post_translational_modifications",
    ])

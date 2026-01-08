"""DNA sequence analysis and genomics module for METAINFORMANT.

This module provides comprehensive tools for DNA sequence analysis, including
FASTA/FASTQ processing, sequence alignment, phylogenetics, population genetics,
motif discovery, and genomic data retrieval from NCBI databases.
"""

from __future__ import annotations

# Import subpackages
from . import (
    sequence,
    alignment,
    expression,
    variation,
    population,
    phylogeny,
    external,
    io,
    integration,
)

# Re-export modules to maintain backward compatibility
# Sequence
from .sequence import core as sequences
from .sequence import composition
from .sequence import motifs
from .sequence import restriction
from .sequence import consensus

# Alignment
from .alignment import pairwise as alignment
from .alignment import msa
from .alignment import distances

# Expression
from .expression import codon
from .expression import transcription
from .expression import translation

# Variation
from .variation import mutations
from .variation import variants

# Population
from .population import core as population
from .population import analysis as population_analysis
from .population import visualization as population_viz

# Phylogeny
from .phylogeny import tree as phylogeny

# External
from .external import ncbi
from .external import entrez
from .external import genomes

# I/O
from .io import fastq

# Integration
from .integration import rna as rna_integration

# Direct imports of commonly used classes
from .io.fastq import FastqRecord

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "sequence",
    "alignment",
    "expression",
    "variation",
    "population",
    "phylogeny",
    "external",
    "io",
    "integration",

    # Module aliases (Backward Compatibility)
    "sequences",
    "composition",
    "motifs",
    "restriction",
    "consensus",
    
    "alignment",
    "msa",
    "distances",
    
    "codon",
    "transcription",
    "translation",
    
    "mutations",
    "variants",
    
    "population",
    "population_analysis",
    "population_viz",
    
    "phylogeny",
    
    "ncbi",
    "entrez",
    "genomes",
    
    "fastq",
    "FastqRecord",
    
    "rna_integration",
]

"""DNA sequence analysis and genomics module for METAINFORMANT.

This module provides comprehensive tools for DNA sequence analysis, including
FASTA/FASTQ processing, sequence alignment, phylogenetics, population genetics,
motif discovery, and genomic data retrieval from NCBI databases.
"""

from __future__ import annotations

# Type checking imports
from typing import TYPE_CHECKING

# Import subpackages
from . import (
    alignment,
    annotation,
    expression,
    external,
    integration,
    io,
    phylogeny,
    population,
    sequence,
    variation,
)

# Alignment
from .alignment import distances, msa
from .alignment import pairwise as alignment

# Annotation
from .annotation import functional, gene_prediction

# Expression
from .expression import codon, transcription, translation

# External
from .external import entrez, genomes, ncbi

# Integration
from .integration import rna as rna_integration

# I/O
from .io import fastq

# Direct imports of commonly used classes
from .io.fastq import FastqRecord

# Phylogeny
from .phylogeny import tree as phylogeny

# Population
from .population import analysis as population_analysis
from .population import core as population
from .population import visualization as population_viz

# Re-export modules to maintain backward compatibility
# Sequence
from .sequence import composition, consensus
from .sequence import core as sequences
from .sequence import kmer, motifs, restriction

# Variation
from .variation import calling, mutations, variants

if TYPE_CHECKING:
    pass

__all__ = [
    # Subpackages
    "annotation",
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
    "kmer",
    "alignment",
    "msa",
    "distances",
    "codon",
    "transcription",
    "translation",
    "mutations",
    "variants",
    "calling",
    "gene_prediction",
    "functional",
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

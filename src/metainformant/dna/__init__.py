"""DNA sequence analysis and genomics module for METAINFORMANT.

This module provides comprehensive tools for DNA sequence analysis, including
FASTA/FASTQ processing, sequence alignment, phylogenetics, population genetics,
motif discovery, and genomic data retrieval from NCBI databases.
"""

from __future__ import annotations

# Import all DNA analysis submodules
from . import (
    alignment,
    codon,
    composition,
    consensus,
    distances,
    entrez,
    fastq,
    genomes,
    motifs,
    msa,
    mutations,
    ncbi,
    phylogeny,
    population,
    population_analysis,
    population_viz,
    restriction,
    rna_integration,
    sequences,
    transcription,
    translation,
    variants,
)

# Optional imports with graceful fallbacks
try:
    from . import population_analysis
except ImportError:
    population_analysis = None

try:
    from . import population_viz
except ImportError:
    population_viz = None

try:
    from . import rna_integration
except ImportError:
    rna_integration = None

# Type checking imports
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    # Core sequence analysis
    "sequences",
    "alignment",
    "phylogeny",
    "population",
    "composition",
    "distances",

    # Specialized analysis
    "codon",
    "motifs",
    "msa",
    "mutations",
    "restriction",
    "consensus",

    # Data formats and I/O
    "fastq",
    "transcription",
    "translation",
    "variants",

    # External data sources
    "genomes",
    "ncbi",
    "entrez",

    # Optional/advanced modules
    "population_analysis",
    "population_viz",
    "rna_integration",
]

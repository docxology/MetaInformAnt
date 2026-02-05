"""METAINFORMANT: Comprehensive Bioinformatics Toolkit for Multi-Omic Analysis.

METAINFORMANT is a Python package providing integrated tools for multi-omic
biological data analysis, including DNA/RNA sequencing, GWAS, proteomics,
epigenomics, and systems biology approaches.

Key Features:
- Multi-omic data integration and analysis
- Statistical and machine learning methods
- Biological network analysis
- Information theory applications
- Publication-quality visualization
- Workflow orchestration and automation

Modules:
- core: Infrastructure utilities (I/O, logging, configuration)
- dna: DNA sequence analysis, alignment, phylogenetics
- rna: RNA-seq workflows and amalgkit integration
- protein: Protein sequence and structure analysis
- gwas: Genome-wide association studies
- math: Mathematical biology and theoretical modeling
- information: Information theory for biological data
- life_events: Life course and temporal analysis
- visualization: Plotting and visualization tools
- networks: Biological network analysis
- multiomics: Cross-omics data integration
- singlecell: Single-cell RNA-seq analysis
- simulation: Synthetic data generation
- quality: Data quality control
- ml: Machine learning for biological data
- ontology: Gene ontology and functional annotation
- phenotype: Phenotypic trait analysis
- ecology: Ecological and community analysis
- epigenome: Epigenomic data analysis
- longread: Long-read sequencing analysis (ONT, PacBio)
- structural_variants: Structural variant detection and analysis
- spatial: Spatial transcriptomics analysis
"""

from __future__ import annotations

__version__ = "0.2.0"
__author__ = "MetaInformAnt Development Team"
__license__ = "Apache License 2.0"

# Type checking imports
from typing import TYPE_CHECKING

# Import all domain modules for easy access
from . import (
    core,
    dna,
    ecology,
    epigenome,
    gwas,
    information,
    life_events,
    longread,
    math,
    menu,
    ml,
    multiomics,
    networks,
    ontology,
    phenotype,
    protein,
    quality,
    rna,
    simulation,
    singlecell,
    spatial,
    structural_variants,
    visualization,
)

if TYPE_CHECKING:
    pass

__all__ = [
    # Package metadata
    "__version__",
    "__author__",
    "__license__",
    # Core infrastructure (always available)
    "core",
    # Molecular analysis modules
    "dna",
    "rna",
    "protein",
    "epigenome",
    "longread",
    # Statistical & ML modules
    "gwas",
    "math",
    "ml",
    "information",
    # Systems biology modules
    "networks",
    "multiomics",
    "singlecell",
    "spatial",
    "structural_variants",
    "simulation",
    # Annotation & metadata modules
    "ontology",
    "phenotype",
    "ecology",
    "life_events",
    # Utilities
    "quality",
    "visualization",
    "menu",
]

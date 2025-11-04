"""DNA domain functionality.

This module provides comprehensive DNA sequence analysis tools including:
- Sequence I/O (FASTA/FASTQ)
- Alignment and phylogenetics
- Population genetics statistics
- Genomic data retrieval (NCBI/Entrez)
- Variant calling and analysis
- Motif discovery and restriction enzyme analysis

Import Pattern:
This module uses lazy imports to avoid optional dependency failures during unrelated commands.
Import specific submodules directly:

    from metainformant.dna import sequences
    sequences.read_fasta(...)

    # Or direct import:
    from metainformant.dna.sequences import read_fasta

Available submodules:
- sequences: FASTA/FASTQ reading and basic sequence operations
- alignment: Pairwise and multiple sequence alignment
- phylogeny: Phylogenetic tree construction (Neighbor-Joining)
- population: Population genetics statistics (π, Tajima's D, Fst)
- population_analysis: Advanced population genetics analysis
- population_viz: Population genetics visualization
- genomes, ncbi, entrez: Genomic data retrieval from NCBI
- fastq: FASTQ-specific processing
- variants: Variant calling and VCF support
- motifs: Pattern discovery and PWM analysis
- restriction: Restriction enzyme analysis
- transcription, translation: Sequence conversion
- mutations: Mutation analysis
- codon: Codon usage analysis
- composition: Sequence composition analysis
- consensus: Consensus sequence generation
- distances: Evolutionary distance calculations
"""

__all__ = [
    # Submodules available for import (lazy loading pattern)
    "genomes",
    "ncbi",
    "entrez",
    "sequences",
    "alignment",
    "phylogeny",
    "population",
    "population_analysis",
    "population_viz",
    "consensus",
    "distances",
    "fastq",
    "motifs",
    "restriction",
    "variants",
    "transcription",
    "translation",
    "mutations",
    "codon",
    "composition",
]

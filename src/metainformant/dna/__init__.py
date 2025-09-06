"""DNA domain functionality."""

# Avoid importing all submodules eagerly to prevent optional dependency failures
# during unrelated commands (e.g., when Biopython or ncbi datasets libs are missing).

# Re-export frequently used modules lazily via simple import wrappers
# Callers should import the specific submodules directly when needed.

__all__ = [
    # Keep list for discoverability; do not import to avoid side effects
    "genomes",
    "ncbi",
    "entrez",
    "sequences",
    "alignment",
    "phylogeny",
    "population",
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

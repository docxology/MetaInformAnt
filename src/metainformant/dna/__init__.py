"""DNA domain functionality."""

# Avoid importing all submodules eagerly to prevent optional dependency failures
# during unrelated commands (e.g., when Biopython or ncbi datasets libs are missing).

# Re-export frequently used modules lazily via simple import wrappers
# Callers should import the specific submodules directly when needed.

# Export commonly used functions for convenience while maintaining lazy import pattern
def _lazy_import(module_name: str, attr_name: str):
    """Lazy import wrapper for DNA module functions."""
    def _get():
        module = __import__(f"metainformant.dna.{module_name}", fromlist=[attr_name])
        return getattr(module, attr_name)
    return property(lambda self: _get())

# Commonly used functions can be accessed via module import pattern:
# from metainformant.dna import sequences
# sequences.read_fasta(...)
# 
# Or via direct import:
# from metainformant.dna.sequences import read_fasta

__all__ = [
    # Keep list for discoverability; do not import to avoid side effects
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

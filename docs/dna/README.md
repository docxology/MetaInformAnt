# DNA Domain Documentation

This directory contains comprehensive documentation for METAINFORMANT's DNA analysis capabilities.

## Overview

The DNA domain provides tools for sequence analysis, alignment, phylogenetics, population genetics, and genomic data integration.

## Documentation Files

### Core DNA Analysis
- **`index.md`**: DNA domain overview and module index
- **`sequences.md`**: DNA sequence manipulation and I/O
- **`alignment.md`**: Pairwise sequence alignment algorithms
- **`msa.md`**: Multiple sequence alignment methods
- **`phylogeny.md`**: Phylogenetic tree construction and analysis

### Population Genetics
- **`population.md`**: Population genetic analysis and diversity metrics
- **`distances.md`**: Evolutionary distance calculations
- **`composition.md`**: Sequence composition and GC content analysis
- **`codon.md`**: Codon usage analysis and genetic code

### Molecular Biology
- **`transcription.md`**: DNA to RNA transcription
- **`translation.md`**: Genetic code translation and ORF finding
- **`motifs.md`**: Motif discovery and pattern analysis
- **`restriction.md`**: Restriction enzyme site mapping
- **`mutations.md`**: Mutation analysis and variant calling

### Genomic Data Integration
- **`ncbi.md`**: NCBI database integration
- **`accessions.md`**: Genome accession validation and retrieval
- **`entrez.md`**: Entrez query construction and execution
- **`variants.md`**: VCF variant parsing and analysis

### Sequence Quality
- **`fastq.md`**: FASTQ format processing and quality analysis

## Related Source Code

- See `src/metainformant/dna/` for implementation details
- See `tests/test_dna_*.py` for comprehensive test coverage
- See `src/metainformant/dna/README.md` for module-specific documentation

## Usage Examples

The DNA domain supports a wide range of biological analyses:

```python
from metainformant.dna import sequences, alignment, phylogeny

# Load and analyze DNA sequences
seqs = sequences.read_fasta("genomes.fasta")
aln = alignment.global_pairwise(seqs.values())
tree = phylogeny.neighbor_joining_tree(aln)

# Population genetic analysis
from metainformant.dna import population
pi = population.nucleotide_diversity(seqs.values())
```

## Integration

DNA analysis integrates with:
- **RNA workflows** for transcriptomic context
- **Protein analysis** for functional annotation
- **Visualization tools** for phylogenetic trees and alignments
- **Statistical methods** for evolutionary analysis

## Testing

Comprehensive tests ensure algorithm correctness:
- Sequence format I/O validation
- Algorithm implementation verification
- Performance benchmarking with real data
- Integration testing across modules

## Contributing

When adding new DNA analysis functionality:
1. Update relevant documentation files
2. Add comprehensive tests in `tests/test_dna_*.py`
3. Update module README and API references
4. Ensure compatibility with existing workflows

This documentation provides complete coverage of METAINFORMANT's DNA analysis capabilities.

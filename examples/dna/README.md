# DNA Analysis Examples

This directory contains educational examples demonstrating METAINFORMANT's DNA sequence analysis capabilities, from basic sequence processing to advanced population genetics.

## Overview

These examples showcase the comprehensive DNA analysis toolkit, covering sequence manipulation, alignment, phylogenetic analysis, and population genetics statistics.

## Examples

### Sequence Processing Basics (`example_sequences.py`)

Learn fundamental DNA sequence operations and analysis.

**Demonstrates:**
- FASTA file reading and parsing
- Reverse complement calculation
- GC content analysis
- Sequence validation and basic statistics
- K-mer counting and frequency analysis

```bash
# Learn basic sequence processing
python examples/dna/example_sequences.py
```

**Output:** `output/examples/dna/sequences_analysis.json`

### Sequence Alignment (`example_alignment.py`)

Master pairwise sequence alignment techniques.

**Demonstrates:**
- Global and local sequence alignment
- Alignment scoring matrices
- Gap penalties and scoring parameters
- Alignment visualization and statistics
- Conserved region identification

```bash
# Learn sequence alignment
python examples/dna/example_alignment.py
```

**Output:** `output/examples/dna/alignment_results.json`

### Phylogenetic Analysis (`example_phylogeny.py`)

Build and analyze phylogenetic trees from DNA sequences.

**Demonstrates:**
- Neighbor-joining tree construction
- UPGMA tree building
- Newick format output and parsing
- Tree visualization and statistics
- Bootstrap support calculation

```bash
# Learn phylogenetic tree building
python examples/dna/example_phylogeny.py
```

**Output:** `output/examples/dna/phylogeny_tree.newick`

### Population Genetics (`example_population.py`)

Calculate population genetics statistics from DNA sequence data.

**Demonstrates:**
- Nucleotide diversity (π) calculation
- Tajima's D statistic
- Fu and Li's D* and F* tests
- Segregating sites identification
- Waterson's theta estimation
- Allele frequency calculations

```bash
# Learn population genetics statistics
python examples/dna/example_population.py
```

**Output:** `output/examples/dna/population_stats.json`

## Learning Progression

1. **Start Here**: `example_sequences.py` - Master sequence I/O and basic analysis
2. **Alignment**: `example_alignment.py` - Understand sequence comparison
3. **Evolution**: `example_phylogeny.py` - Learn tree building and phylogenetics
4. **Populations**: `example_population.py` - Explore population genetics

## Related Documentation

- **DNA Module Docs**: [`docs/dna/`](../../docs/dna/) - Complete DNA analysis documentation
- **Core Examples**: [`examples/core/`](../core/) - Foundational METAINFORMANT concepts
- **Scripts**: [`scripts/dna/`](../../scripts/dna/) - Production DNA analysis workflows

## Data Sources

Examples use simulated DNA sequence data for demonstration. For real data analysis, see:
- **NCBI Integration**: `metainformant.dna.ncbi` for genome downloads
- **SRA Data**: `metainformant.dna.genomes` for sequence data retrieval
- **Production Scripts**: `scripts/dna/run_dna_analysis.py` for full workflows

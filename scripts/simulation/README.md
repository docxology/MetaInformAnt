# Simulation Scripts Documentation

This directory contains module-specific simulation scripts for generating synthetic data across all METAINFORMANT domain modules. Each script provides comprehensive simulation capabilities tailored to its respective domain.

## Overview

Each module in METAINFORMANT has a corresponding simulation script that generates realistic synthetic data for testing, validation, and method development. All scripts follow consistent patterns and use core utilities for I/O, logging, path management, and parameter validation.

## Recent Improvements

All simulation scripts have been enhanced with:

- **Comprehensive Parameter Validation**: All parameters are validated using `core.validation` utilities with type checking and range validation
- **Enhanced Documentation**: Complete docstrings with parameter descriptions, return types, and exception documentation
- **Improved Error Handling**: Clear, contextual error messages with proper exception handling
- **Type Safety**: Full type hints throughout all functions
- **Consistent Patterns**: Standardized structure across all scripts for maintainability
- **Domain-Specific Constraints**: Validation of domain-specific requirements (e.g., minimum sequences for phylogeny)

## Available Scripts

### Core Utilities (`simulate_core.py`)

Generates test data for core utilities including configuration files, workflow test data, and I/O pattern examples.

**Simulation Types:**
- `config`: Generate configuration files (YAML, JSON, TOML)
- `workflow`: Generate workflow test data with dependencies
- `io`: Generate test files in various formats (JSON, CSV, TXT, JSONL)

**Example:**
```bash
python3 scripts/simulation/simulate_core.py --type config --n-configs 5 --config-type yaml
python3 scripts/simulation/simulate_core.py --type workflow --n-steps 10
python3 scripts/simulation/simulate_core.py --type io --n-files 20 --file-types json csv
```

### DNA (`simulate_dna.py`)

Generates synthetic DNA sequences for testing and validation.

**Simulation Types:**
- `sequences`: Basic DNA sequence generation with mutations
- `population`: Population genetics data with diversity
- `alignment`: Sequences for alignment testing
- `phylogeny`: Sequences for phylogenetic tree construction

**Example:**
```bash
python3 scripts/simulation/simulate_dna.py --type sequences --n 100 --length 1000
python3 scripts/simulation/simulate_dna.py --type population --n 50 --length 2000 --diversity 0.01
python3 scripts/simulation/simulate_dna.py --type alignment --n 10 --length 500
```

### RNA (`simulate_rna.py`)

Generates synthetic RNA expression data including count matrices, differential expression, and batch effects.

**Simulation Types:**
- `counts`: Basic RNA-seq count matrix
- `differential`: Differential expression between groups
- `batch`: Batch effect simulation

**Example:**
```bash
python3 scripts/simulation/simulate_rna.py --type counts --num-genes 1000 --num-samples 20
python3 scripts/simulation/simulate_rna.py --type differential --num-genes 500 --n-de-genes 50
python3 scripts/simulation/simulate_rna.py --type batch --num-genes 2000 --n-batches 3
```

### Protein (`simulate_protein.py`)

Generates synthetic protein sequences, structure coordinates, domain annotations, and PPI networks.

**Simulation Types:**
- `sequences`: Basic protein sequence generation
- `structure`: Protein structure coordinates (CA atoms)
- `domains`: Domain annotation simulation
- `ppi`: Protein-protein interaction network

**Example:**
```bash
python3 scripts/simulation/simulate_protein.py --type sequences --n 100 --length 200
python3 scripts/simulation/simulate_protein.py --type structure --n 10 --length 150
python3 scripts/simulation/simulate_protein.py --type ppi --n 100 --interactions 200
```

### Epigenome (`simulate_epigenome.py`)

Generates synthetic epigenetic data including DNA methylation, chromatin accessibility, and ChIP-seq peaks.

**Simulation Types:**
- `methylation`: DNA methylation patterns (CpG sites)
- `chromatin`: Chromatin accessibility tracks (bedGraph)
- `chipseq`: ChIP-seq peak data

**Example:**
```bash
python3 scripts/simulation/simulate_epigenome.py --type methylation --n-sites 10000
python3 scripts/simulation/simulate_epigenome.py --type chromatin --chromosome chr1 --length 1000000
python3 scripts/simulation/simulate_epigenome.py --type chipseq --n-peaks 500
```

### Ontology (`simulate_ontology.py`)

Generates synthetic GO term annotations and enrichment test data.

**Simulation Types:**
- `annotations`: Gene-to-GO-term annotations
- `enrichment`: Enrichment test results

**Example:**
```bash
python3 scripts/simulation/simulate_ontology.py --type annotations --n-genes 1000 --n-terms 100
python3 scripts/simulation/simulate_ontology.py --type enrichment --n-genes 500 --n-significant 10
```

### Phenotype (`simulate_phenotype.py`)

Generates synthetic phenotypic trait data including morphological and behavioral traits.

**Simulation Types:**
- `morphological`: Morphological trait simulation
- `behavioral`: Behavioral trait simulation
- `correlation`: Correlated trait simulation

**Example:**
```bash
python3 scripts/simulation/simulate_phenotype.py --type morphological --n-samples 100
python3 scripts/simulation/simulate_phenotype.py --type behavioral --n-samples 200
python3 scripts/simulation/simulate_phenotype.py --type correlation --n-samples 150 --n-traits 10
```

### Ecology (`simulate_ecology.py`)

Generates synthetic ecological data including species abundance, community composition, and environmental metadata.

**Simulation Types:**
- `abundance`: Species abundance matrix
- `community`: Community composition with diversity
- `environmental`: Environmental metadata

**Example:**
```bash
python3 scripts/simulation/simulate_ecology.py --type abundance --n-species 50 --n-samples 20
python3 scripts/simulation/simulate_ecology.py --type community --n-species 100 --diversity 0.7
python3 scripts/simulation/simulate_ecology.py --type environmental --n-samples 25
```

### Math (`simulate_math.py`)

Generates synthetic data for mathematical biology models.

**Simulation Types:**
- `coalescent`: Coalescent simulation data
- `selection`: Selection experiment simulation
- `popgen`: Population genetics data

**Example:**
```bash
python3 scripts/simulation/simulate_math.py --type coalescent --n-samples 10 --n-loci 1000
python3 scripts/simulation/simulate_math.py --type selection --n-generations 50 --population-size 100
python3 scripts/simulation/simulate_math.py --type popgen --n-sequences 20 --diversity 0.01
```

### Visualization (`simulate_visualization.py`)

Generates synthetic data for various plot types and visualizations.

**Simulation Types:**
- `timeseries`: Time-series data
- `multidim`: Multi-dimensional data with clusters
- `statistical`: Statistical test data

**Example:**
```bash
python3 scripts/simulation/simulate_visualization.py --type timeseries --n-points 100 --n-series 5
python3 scripts/simulation/simulate_visualization.py --type multidim --n-samples 50 --n-features 10
python3 scripts/simulation/simulate_visualization.py --type statistical --n-samples 100 --n-groups 4
```

### Single-Cell (`simulate_singlecell.py`)

Generates synthetic single-cell RNA-seq data including count matrices, cell types, and trajectories.

**Simulation Types:**
- `counts`: Basic scRNA-seq count matrix
- `celltypes`: Count matrix with cell type annotations
- `trajectory`: Trajectory/pseudotime data

**Example:**
```bash
python3 scripts/simulation/simulate_singlecell.py --type counts --n-cells 1000 --n-genes 2000
python3 scripts/simulation/simulate_singlecell.py --type celltypes --n-cells 500 --n-types 5
python3 scripts/simulation/simulate_singlecell.py --type trajectory --n-cells 800 --n-states 10
```

### Quality (`simulate_quality.py`)

Generates synthetic quality control data including FASTQ quality scores, contamination patterns, and quality metrics.

**Simulation Types:**
- `fastq`: FASTQ quality score simulation
- `contamination`: Contamination pattern simulation
- `metrics`: Quality control metrics

**Example:**
```bash
python3 scripts/simulation/simulate_quality.py --type fastq --n-reads 10000 --read-length 100
python3 scripts/simulation/simulate_quality.py --type contamination --n-reads 5000 --contamination-rate 0.1
python3 scripts/simulation/simulate_quality.py --type metrics --n-samples 20
```

### Networks (`simulate_networks.py`)

Generates synthetic biological networks including PPI, regulatory, and pathway networks.

**Simulation Types:**
- `ppi`: Protein-protein interaction network
- `regulatory`: Gene regulatory network
- `pathway`: Pathway network

**Example:**
```bash
python3 scripts/simulation/simulate_networks.py --type ppi --n-nodes 100 --n-edges 200
python3 scripts/simulation/simulate_networks.py --type regulatory --n-nodes 50 --n-edges 150
python3 scripts/simulation/simulate_networks.py --type pathway --n-nodes 30 --n-edges 50
```

### Machine Learning (`simulate_ml.py`)

Generates synthetic feature matrices, labels, and test data for classification and regression.

**Simulation Types:**
- `classification`: Classification dataset
- `regression`: Regression dataset
- `features`: Feature matrix with different types

**Example:**
```bash
python3 scripts/simulation/simulate_ml.py --type classification --n-samples 100 --n-features 50 --n-classes 3
python3 scripts/simulation/simulate_ml.py --type regression --n-samples 200 --n-features 30
python3 scripts/simulation/simulate_ml.py --type features --n-samples 150 --n-features 100
```

### Multi-Omics (`simulate_multiomics.py`)

Generates synthetic multi-omics data including cross-platform and integrated datasets.

**Simulation Types:**
- `cross-platform`: Cross-platform multi-omics data (genomics, transcriptomics, proteomics)
- `integrated`: Integrated multi-omics dataset

**Example:**
```bash
python3 scripts/simulation/simulate_multiomics.py --type cross-platform --n-samples 20 --n-genes 1000
python3 scripts/simulation/simulate_multiomics.py --type integrated --n-samples 30 --n-features 200 --n-omics 3
```

### GWAS (`simulate_gwas.py`)

Generates synthetic GWAS data including variants, genotype matrices, phenotypes, and population structure.

**Simulation Types:**
- `variants`: Variant VCF file
- `genotypes`: Genotype matrix
- `phenotype`: Phenotype with genetic effects
- `population`: Population structure data

**Example:**
```bash
python3 scripts/simulation/simulate_gwas.py --type variants --n-variants 1000 --n-samples 100
python3 scripts/simulation/simulate_gwas.py --type genotypes --n-samples 200 --n-sites 500
python3 scripts/simulation/simulate_gwas.py --type phenotype --n-samples 150 --n-causal 10
python3 scripts/simulation/simulate_gwas.py --type population --n-samples 100 --n-populations 2
```

### Life Events (`simulate_life_events.py`)

Generates synthetic life event sequences and life course patterns.

**Simulation Types:**
- `sequences`: Life event sequences
- `patterns`: Life course patterns

**Example:**
```bash
python3 scripts/simulation/simulate_life_events.py --type sequences --n-sequences 100 --min-events 5
python3 scripts/simulation/simulate_life_events.py --type patterns --n-sequences 200 --n-patterns 5
```

### Information Theory (`simulate_information.py`)

Generates synthetic data for information-theoretic analysis.

**Simulation Types:**
- `sequences`: Sequences with varying information content
- `mutual`: Mutual information test data

**Example:**
```bash
python3 scripts/simulation/simulate_information.py --type sequences --n-sequences 100 --length 1000 --complexity high
python3 scripts/simulation/simulate_information.py --type mutual --n-samples 200 --correlation 0.7
```

## Common Features

All simulation scripts share common features:

1. **Consistent Interface**: All scripts use argparse with similar argument patterns
2. **Output Organization**: Outputs are organized in `output/simulation/<module>/` by default
3. **Reproducibility**: Seed-based random number generation for reproducible results
4. **Metadata**: Summary JSON files track simulation parameters and results
5. **Logging**: Comprehensive logging using core.logging utilities
6. **Parameter Validation**: Comprehensive parameter validation using core.validation utilities
   - Type checking for all parameters
   - Range validation (e.g., GC content in [0,1], positive integers)
   - Domain-specific constraints (e.g., minimum sequences for phylogeny)
7. **Error Handling**: Clear error messages and exception handling with context
8. **Documentation**: Comprehensive docstrings with parameter descriptions and examples
9. **Type Safety**: Full type hints for all functions

## Integration with Analysis Modules

Simulated data can be directly used with corresponding analysis modules:

```python
# Simulate DNA sequences
# python3 scripts/simulation/simulate_dna.py --type sequences --n 100

# Analyze with DNA module
from metainformant.dna import sequences, alignment
seqs = sequences.read_fasta("output/simulation/dna/sequences.fasta")
aln = alignment.global_pairwise(list(seqs.values())[:2])
```

## Getting Help

Each script provides comprehensive help:

```bash
python3 scripts/simulation/simulate_<module>.py --help
```

## Parameter Validation

All scripts include comprehensive parameter validation:

- **Type Checking**: Ensures parameters are of correct types
- **Range Validation**: Validates parameter ranges (e.g., probabilities in [0,1])
- **Domain Constraints**: Enforces domain-specific requirements (e.g., minimum 3 sequences for phylogeny)
- **Clear Error Messages**: Descriptive error messages when validation fails

Example validation errors:
```bash
# Invalid GC content
python3 scripts/simulation/simulate_dna.py --type sequences --n 100 --gc-content 1.5
# Error: gc_content must be in range [0.0, 1.0], got 1.5

# Insufficient sequences for phylogeny
python3 scripts/simulation/simulate_dna.py --type phylogeny --n 2
# Error: n_sequences must be >= 3, got 2
```

## Best Practices

1. **Use Seeds**: Always specify `--seed` for reproducible results
2. **Check Outputs**: Review `simulation_summary.json` for simulation metadata
3. **Start Small**: Test with small parameters first, then scale up
4. **Validate Inputs**: Scripts validate parameters, but verify your use case
5. **Review Logs**: Use `--verbose` flag for detailed logging during development

## See Also

- [`src/metainformant/simulation/README.md`](../../src/metainformant/simulation/README.md) - Core simulation module documentation
- [`src/metainformant/README.md`](../../src/metainformant/README.md) - Main package documentation


# Source Code

Python source code for METAINFORMANT bioinformatics toolkit.

## Structure

```
src/
└── metainformant/          # Main package
    ├── core/               # Shared infrastructure
    ├── dna/                # DNA sequence analysis
    ├── rna/                # RNA-seq workflows
    ├── protein/            # Protein analysis
    ├── gwas/               # GWAS pipeline
    ├── visualization/      # Plotting (57+ types)
    ├── networks/           # Biological networks
    ├── singlecell/         # scRNA-seq
    ├── multiomics/         # Multi-omic integration
    ├── ml/                 # Machine learning
    ├── math/               # Population genetics theory
    ├── information/        # Information theory
    ├── life_events/        # Event sequence analysis
    ├── ontology/           # GO analysis
    ├── phenotype/          # Trait analysis
    ├── ecology/            # Community ecology
    ├── epigenome/          # Methylation, ChIP-seq
    ├── simulation/         # Synthetic data
    ├── quality/            # QC metrics
    └── menu/               # Interactive CLI
```

## Usage

```python
# Import modules
from metainformant.core import io, config, paths
from metainformant.dna import sequences, alignment
from metainformant.rna import workflow
from metainformant.visualization import plots

# All operations use real implementations (NO MOCKING policy)
```

## Key Patterns

- **Core utilities**: Use `metainformant.core` for I/O, config, logging
- **Output**: All results go to `output/` directory
- **Configuration**: YAML configs with environment variable overrides
- **Type hints**: Python 3.11+ with comprehensive typing

## Related

- [Module Documentation](../docs/)
- [API Reference](../docs/cli.md)
- [Tests](../tests/)

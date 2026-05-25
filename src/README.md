# Source Code

Python source code for METAINFORMANT bioinformatics toolkit.

## Structure

```
src/
 metainformant/ # Main package
 core/ # Shared infrastructure
 dna/ # DNA sequence analysis
 rna/ # RNA-seq workflows
 protein/ # Protein analysis
 gwas/ # GWAS pipeline
 visualization/ # Plotting (70+ types)
 networks/ # Biological networks
 singlecell/ # scRNA-seq
 multiomics/ # Multi-omic integration
 ml/ # Machine learning
 math/ # Population genetics theory
 information/ # Information theory
 life_events/ # Event sequence analysis
 ontology/ # GO analysis
 phenotype/ # Trait analysis
 ecology/ # Community ecology
 epigenome/ # Methylation, ChIP-seq
 simulation/ # Synthetic data
 longread/ # Long-read sequencing
 metagenomics/ # Metagenomic analysis
 pharmacogenomics/ # Pharmacogenomic analysis
 spatial/ # Spatial transcriptomics
 structural_variants/ # Structural variant analysis
 quality/ # QC metrics
 menu/ # Interactive CLI
```

## Usage

```python
# Import modules using the current package layout
from metainformant.core import io
from metainformant.core.io import paths
from metainformant.core.utils import config
from metainformant.dna.sequence import core as dna_core
from metainformant.rna.engine import workflow
from metainformant.visualization import plots

# All operations use real implementations (REAL IMPLEMENTATION policy)
```

## Key Patterns

- **Core utilities**: Use `metainformant.core.io`, `metainformant.core.utils.config`, and `metainformant.core.utils.logging`
- **Output**: All results go to `output/` directory
- **Configuration**: YAML configs with environment variable overrides
- **Type hints**: Python 3.11+ with comprehensive typing

## Related

- [Module Documentation](../docs/)
- [API Reference](../docs/cli.md)
- [Tests](../tests/)

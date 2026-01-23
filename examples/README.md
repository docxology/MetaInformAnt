# Examples

Runnable examples demonstrating METAINFORMANT functionality across all domains.

## Directory Structure

```
examples/
├── core/           # Core utilities (I/O, config, logging)
├── dna/            # DNA sequence analysis, phylogeny
├── rna/            # RNA-seq workflows
├── protein/        # Protein analysis
├── gwas/           # GWAS pipeline examples
├── visualization/  # Plotting and figures
├── networks/       # Biological networks
├── singlecell/     # Single-cell analysis
├── multiomics/     # Multi-omic integration
├── ml/             # Machine learning
├── math/           # Population genetics, coalescent
├── information/    # Information theory
├── life_events/    # Life event analysis
├── ontology/       # GO analysis
├── phenotype/      # Trait analysis
├── ecology/        # Community ecology
├── epigenome/      # Methylation, ChIP-seq
├── simulation/     # Synthetic data
├── quality/        # QC examples
├── integration/    # Cross-module examples
└── templates/      # Example templates
```

## Running Examples

```bash
# Run a single example
uv run python examples/core/io_example.py

# Run all examples in a domain
uv run python -m pytest examples/dna/ -v

# Validate all examples work
bash scripts/test_examples/validate_examples.sh
```

## Example Format

Each example follows this pattern:

```python
#!/usr/bin/env python3
"""Example: Brief description.

Demonstrates: Feature being demonstrated.
Output: output/examples/<domain>/<name>/
"""
from pathlib import Path
from metainformant.<domain> import <module>

def main():
    # Setup output directory
    output_dir = Path("output/examples/<domain>/<name>")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Example code
    result = <module>.function(...)

    # Save results
    print(f"Results saved to {output_dir}")

if __name__ == "__main__":
    main()
```

## Key Examples

### Core
- `io_example.py` - JSON/CSV file operations
- `config_example.py` - Configuration loading

### DNA
- `sequence_analysis.py` - GC content, k-mers
- `phylogeny_example.py` - Tree construction
- `population_genetics.py` - Diversity statistics

### RNA
- `workflow_example.py` - Amalgkit workflow
- `quantification.py` - Expression quantification

### GWAS
- `association_test.py` - Association testing
- `visualization.py` - Manhattan/QQ plots

### Visualization
- `plots_gallery.py` - All plot types
- `genomics_plots.py` - Genomic visualizations

## Output

All examples write to `output/examples/<domain>/`:

```
output/examples/
├── core/
├── dna/
├── gwas/
└── ...
```

## Dependencies

Examples may require optional dependencies:

```bash
# Install all optional dependencies
uv pip install -e ".[dev,scientific,ml,networks]"
```

## Related

- [Module Documentation](../docs/)
- [Test Suite](../tests/)
- [Source Code](../src/metainformant/)

# Agent Directives: examples

## Role
Working example scripts demonstrating METAINFORMANT functionality.

## Directory Structure
Examples organized by module:
- `core/` - Core infrastructure examples
- `dna/` - DNA sequence analysis examples
- `ecology/` - Ecology analysis examples
- `epigenome/` - Epigenome analysis examples
- `gwas/` - GWAS pipeline examples
- `information/` - Information theory examples
- `integration/` - Cross-module integration examples
- `life_events/` - Life events analysis examples
- `math/` - Mathematical biology examples
- `ml/` - Machine learning examples
- `multiomics/` - Multi-omics integration examples
- `networks/` - Network analysis examples
- `ontology/` - Ontology analysis examples
- `phenotype/` - Phenotype analysis examples
- `protein/` - Protein analysis examples
- `quality/` - Quality control examples
- `rna/` - RNA-seq and amalgkit examples
- `simulation/` - Simulation examples
- `singlecell/` - Single-cell analysis examples
- `templates/` - Example templates for creating new examples
- `visualization/` - Visualization examples

## Rules and Constraints

### Real Implementations Only
All examples must:
- Use REAL functions and algorithms
- Produce REAL outputs
- Make REAL API calls (or skip gracefully)
- NEVER use mocks, stubs, or placeholder data

### Example Naming
Follow the pattern `example_{feature}.py`:
- `example_alignment.py` - Sequence alignment example
- `example_association.py` - GWAS association example

### Example Structure
```python
"""
Example: {Feature Name}

Demonstrates {feature} functionality in METAINFORMANT.
"""
from metainformant.{module} import {functions}

def main():
    # Example implementation with real operations
    result = actual_function(real_data)
    print(f"Result: {result}")

if __name__ == "__main__":
    main()
```

### Testing Examples
Run examples via:
```bash
uv run python examples/{module}/example_{feature}.py
```

Validate all examples:
```bash
uv run python scripts/validate_examples.py
```

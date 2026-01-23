# Specification: examples

## Scope

Working example scripts demonstrating METAINFORMANT functionality. Contains real, executable code examples for all 19 domain modules. Examples use actual implementations, produce real outputs, and make genuine API calls when required.

## Architecture

- **Dependency Level**: Examples
- **Component Type**: Demonstration Scripts
- **Design Pattern**: Real Implementation Examples (NO MOCKING)

### Directory Structure
```
examples/
├── core/               # Core infrastructure examples
├── dna/                # DNA sequence analysis examples
├── rna/                # RNA-seq and amalgkit examples
├── gwas/               # GWAS pipeline examples
├── {module}/           # Domain-specific examples
├── integration/        # Cross-module integration examples
└── templates/          # Templates for creating new examples
```

## Data Structures

### Example File Types
- **example_{feature}.py**: Single-feature demonstration scripts
- **example_{workflow}.py**: Multi-step workflow examples
- **templates/*.py**: Boilerplate templates for new examples

### Example Naming Convention
```
example_{feature}.py       # Feature demonstration
example_{module}_{use_case}.py  # Module-specific use case
```

### Module Subdirectories
Each domain has a corresponding examples subdirectory:
- core/, dna/, rna/, gwas/, protein/, epigenome/
- networks/, multiomics/, singlecell/, visualization/
- quality/, ml/, math/, information/, ontology/
- phenotype/, ecology/, simulation/, life_events/

## Interface

### Running Examples
```bash
# Run individual example
uv run python examples/dna/example_alignment.py
uv run python examples/gwas/example_association.py

# Validate all examples
uv run python scripts/validate_examples.py

# Benchmark examples
uv run python scripts/benchmark_examples.py
```

### Example Structure
```python
"""
Example: Feature Name

Demonstrates feature functionality in METAINFORMANT.
"""
from metainformant.{module} import {functions}

def main():
    # Real implementation with actual operations
    result = actual_function(real_data)
    print(f"Result: {result}")

if __name__ == "__main__":
    main()
```

### Example Conventions
- **Real Implementations Only**: Use actual functions and algorithms
- **Real Outputs**: Produce genuine results, not placeholders
- **Real API Calls**: Make actual API calls or skip gracefully
- **No Mocks**: Never use mocks, stubs, or fake data
- **Self-Contained**: Each example should run independently
- **Documented**: Include docstrings explaining purpose and usage
- **Output to output/**: Write results to output/ directory

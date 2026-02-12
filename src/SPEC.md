# Specification: src

## Scope

Root source directory containing the metainformant Python package. This is the package distribution root with the main `metainformant/` module containing all 25 domain-specific bioinformatics modules.

## Architecture

- **Dependency Level**: Package Root
- **Component Type**: Source Distribution
- **Package Layout**: src-layout (src/metainformant/)

### Directory Structure
```
src/
├── metainformant/          # Main package directory
│   ├── __init__.py         # Package initialization
│   ├── __main__.py         # CLI entry point
│   ├── core/               # Core infrastructure
│   ├── dna/                # DNA analysis
│   ├── rna/                # RNA analysis
│   └── {module}/           # 16 additional domain modules
└── metainformant.egg-info/ # Package metadata (generated)
```

## Data Structures

### Package Components
- **metainformant/**: Main Python package with 25 domain modules
- **metainformant.egg-info/**: Generated package metadata for pip/uv installation

### Module Organization
The metainformant package contains domain-driven modules:
- Core: core/
- Genomics: dna/, rna/, gwas/, epigenome/
- Proteomics: protein/
- Systems Biology: networks/, multiomics/, singlecell/
- Analysis: quality/, ml/, math/, information/
- Annotation: ontology/, phenotype/
- Specialized: ecology/, simulation/, life_events/, visualization/, menu/

## Interface

### Package Installation
```bash
# Editable install for development
uv pip install -e .

# Install with extras
uv pip install -e ".[dev,scientific,ml,networks]"
```

### Importing
```python
# Import from metainformant package
from metainformant.core import io, config, logging
from metainformant.dna import sequence, alignment
from metainformant.rna import workflow, analysis

# CLI access
python -m metainformant --help
```

### Development Conventions
- All source code resides in src/metainformant/
- Use absolute imports from metainformant.*
- Each module has its own SPEC.md documenting internals
- Follow NO MOCKING policy in all implementations

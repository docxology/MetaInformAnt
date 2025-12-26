# METAINFORMANT Examples

This directory contains educational examples demonstrating METAINFORMANT usage patterns across all biological analysis domains. These examples complement the production scripts in `scripts/` by focusing on learning concepts rather than production workflows.

## Overview

METAINFORMANT examples are organized by domain to match the module structure, providing focused, runnable demonstrations of key functionality. Examples focus on education and concept demonstration, while production workflows belong in `scripts/`.

## Example Organization Framework

```mermaid
graph TD
    A[Examples Directory] --> B[Core Utilities]
    A --> C[Domain-Specific]
    A --> D[Integration Patterns]

    B --> E[Configuration Management]
    B --> F[File I/O Patterns]
    B --> G[Logging Systems]
    B --> H[Path Operations]

    C --> I[DNA Sequence Analysis]
    C --> J[RNA Transcriptomics]
    C --> K[GWAS Association]
    C --> L[Protein Structure]
    C --> M[15+ Domain Modules]

    D --> N[Multi-Omics Integration]
    D --> O[End-to-End Workflows]
    D --> P[Cross-Domain Analysis]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style C fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style D fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Learning Progression"
        Q[Foundation] -.-> B
        R[Domain Knowledge] -.-> C
        S[Integration Skills] -.-> D
    end

    subgraph "Example Types"
        T[Educational] -.-> Q
        U[Demonstrative] -.-> R
        V[Integrative] -.-> S
    end

    subgraph "Output Convention"
        W[output/examples/] -.-> T
        X[Informative Names] -.-> T
        Y[Publication Ready] -.-> T
    end
```

### Learning Pathway Architecture

```mermaid
graph LR
    A[Beginner] --> B[Core Concepts]
    B --> C[Domain Basics]
    C --> D[Advanced Methods]
    D --> E[Integration]

    B --> F[example_config.py]
    B --> G[example_io.py]
    B --> H[example_logging.py]

    C --> I[example_sequences.py]
    C --> J[example_alignment.py]
    C --> K[example_phylogeny.py]

    D --> L[example_gwas.py]
    D --> M[example_multiomics.py]
    D --> N[example_networks.py]

    E --> O[example_complete_workflow.py]
    E --> P[example_cross_domain.py]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style E fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Skill Levels"
        Q[Configuration] -.-> B
        R[File I/O] -.-> B
        S[Logging] -.-> B
        T[Path Management] -.-> B
    end

    subgraph "Domain Skills"
        U[Sequence Analysis] -.-> C
        V[Alignment] -.-> C
        W[Phylogenetics] -.-> C
        X[Population Genetics] -.-> C
    end

    subgraph "Advanced Skills"
        Y[Association Testing] -.-> D
        Z[Systems Biology] -.-> D
        AA[Multi-Omics] -.-> D
        BB[Network Analysis] -.-> D
    end

    subgraph "Integration Skills"
        CC[End-to-End Pipelines] -.-> E
        DD[Cross-Domain Analysis] -.-> E
        EE[Workflow Orchestration] -.-> E
        FF[Production Deployment] -.-> E
    end
```

### Example Execution Framework

```mermaid
graph TD
    A[Example Script] --> B[Import Dependencies]
    B --> C[Setup Configuration]

    C --> D[Load Sample Data]
    D --> E[Execute Analysis]

    E --> F[Generate Results]
    F --> G[Create Visualizations]

    G --> H[Save Outputs]
    H --> I[Display Summary]

    I --> J[Example Complete]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style E fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style J fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Execution Pattern"
        K[Self-Contained] -.-> A
        L[Minimal Setup] -.-> A
        M[Educational Focus] -.-> A
        N[Demonstrates Concept] -.-> A
    end

    subgraph "Output Convention"
        O[output/examples/] -.-> H
        P[Informative Names] -.-> H
        Q[Publication Quality] -.-> H
        R[Structured Format] -.-> H
    end

    subgraph "Quality Standards"
        S[Working Code] -.-> I
        T[Clear Documentation] -.-> I
        U[Educational Value] -.-> I
        V[Best Practices] -.-> I
    end
```

### Domain Coverage Matrix

```mermaid
graph TD
    A[Domain Examples] --> B[Core Infrastructure]
    A --> C[Molecular Biology]
    A --> D[Statistical Genetics]
    A --> E[Systems Biology]
    A --> F[Applied Biology]

    B --> G[Configuration & I/O]
    B --> H[Logging & Paths]
    B --> I[Workflow Orchestration]

    C --> J[DNA Sequences]
    C --> K[RNA Expression]
    C --> L[Protein Structure]
    C --> M[Epigenetics]

    D --> N[GWAS Analysis]
    D --> O[Population Genetics]
    D --> P[Statistical Modeling]

    E --> Q[Networks & Pathways]
    E --> R[Multi-Omics Integration]
    E --> S[Single-Cell Analysis]

    F --> T[Ontology & Annotation]
    F --> U[Phenotype Analysis]
    F --> V[Ecology & Diversity]
    F --> W[Life Events]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style E fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style F fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Coverage Areas"
        X[20+ Example Categories] -.-> A
        Y[50+ Individual Examples] -.-> A
        Z[Educational Progression] -.-> A
    end

    subgraph "Quality Assurance"
        AA[Runnable Code] -.-> A
        BB[Real Data] -.-> A
        CC[Comprehensive Docs] -.-> A
        DD[Testing Validation] -.-> A
    end
```

## Directory Structure

```
examples/
├── README.md                    # This file
├── AGENTS.md                    # AI development documentation
├── core/                        # Core utilities examples
│   ├── README.md
│   ├── example_config.py        # Configuration management
│   ├── example_io.py            # File I/O patterns
│   ├── example_logging.py       # Logging patterns
│   ├── example_paths.py         # Path management
│   └── example_workflow.py      # Basic workflow orchestration
├── dna/                         # DNA analysis examples
│   ├── README.md
│   ├── example_sequences.py     # Sequence processing basics
│   ├── example_alignment.py     # Sequence alignment
│   ├── example_phylogeny.py     # Phylogenetic tree building
│   └── example_population.py    # Population genetics statistics
├── rna/                         # RNA analysis examples
│   ├── README.md
│   ├── example_amalgkit.py      # Amalgkit workflow basics
│   └── example_quantification.py # Expression quantification
├── gwas/                        # GWAS examples
│   ├── README.md
│   ├── example_association.py   # Association testing
│   └── example_visualization.py # GWAS visualization
├── [other domains]/             # Examples for all modules
│   ├── README.md
│   └── example_*.py             # Domain-specific examples
└── integration/                 # Cross-domain examples
    ├── README.md
    ├── example_dna_rna.py       # Multi-omic integration
    ├── example_multiomics.py    # Cross-omics analysis
    └── example_complete_workflow.py # End-to-end pipeline
```

## Examples vs Scripts

METAINFORMANT provides both **examples** (for learning) and **scripts** (for production):

| Type | Purpose | Location | Focus |
|------|---------|----------|--------|
| **Examples** | Learning concepts | `examples/` | Educational, minimal demos |
| **Scripts** | Production workflows | `scripts/` | Full orchestrators, batch processing |

### When to Use Examples
- Learning new concepts and APIs
- Understanding module functionality
- Testing ideas before production
- Educational purposes and tutorials

### When to Use Scripts
- Production data analysis
- Batch processing workflows
- Automated pipelines
- End-to-end data processing

## Core Examples (`examples/core/`)

Foundational examples demonstrating METAINFORMANT core utilities:

### Configuration Management
```bash
# Learn how to load configs with environment overrides
python examples/core/example_config.py
```

### File I/O Patterns
```bash
# See JSON, CSV, and JSONL file operations
python examples/core/example_io.py
```

### Logging Setup
```bash
# Learn structured logging patterns
python examples/core/example_logging.py
```

### Path Management
```bash
# Understand path validation and containment
python examples/core/example_paths.py
```

## Domain Examples

### DNA Analysis (`examples/dna/`)

Learn fundamental DNA sequence analysis concepts:

#### Sequence Processing Basics
```bash
# Learn FASTA reading, reverse complement, GC content
python examples/dna/example_sequences.py
```

#### Sequence Alignment
```bash
# Understand pairwise sequence alignment
python examples/dna/example_alignment.py
```

#### Phylogenetics
```bash
# Build neighbor-joining phylogenetic trees
python examples/dna/example_phylogeny.py
```

#### Population Genetics
```bash
# Calculate nucleotide diversity and Tajima's D
python examples/dna/example_population.py
```

### RNA Analysis (`examples/rna/`)

Explore transcriptomic analysis workflows:

#### Amalgkit Workflow Basics
```bash
# Learn amalgkit workflow setup and planning
python examples/rna/example_amalgkit.py
```

#### Expression Quantification
```bash
# Understand gene expression quantification
python examples/rna/example_quantification.py
```

### GWAS (`examples/gwas/`)

Genome-wide association study examples:

#### Association Testing
```bash
# Learn basic association test implementation
python examples/gwas/example_association.py
```

#### GWAS Visualization
```bash
# Create Manhattan and QQ plots
python examples/gwas/example_visualization.py
```

## Integration Examples (`examples/integration/`)

Cross-domain workflows demonstrating multi-omic analysis:

#### DNA-RNA Integration
```bash
# Use DNA coordinates for RNA analysis
python examples/integration/example_dna_rna.py
```

#### Multi-Omics Analysis
```bash
# Cross-omics correlation and integration
python examples/integration/example_multiomics.py
```

#### Complete Workflow
```bash
# End-to-end bioinformatics pipeline
python examples/integration/example_complete_workflow.py
```

## Running Examples

### Prerequisites

```bash
# Install METAINFORMANT and dependencies
uv sync

# Activate environment
source .venv/bin/activate  # or uv run for commands
```

### Execution

```bash
# Run any example
python examples/core/example_config.py
python examples/dna/example_sequences.py
python examples/gwas/example_association.py

# Examples create output in output/examples/
# Results are saved with informative names
```

### Example Output Structure

All examples follow consistent output patterns:

```
output/examples/
├── core/
│   ├── config_example.json
│   └── io_example.csv
├── dna/
│   ├── sequences_analysis.json
│   ├── alignment_results.json
│   ├── phylogeny_tree.newick
│   └── population_stats.json
├── gwas/
│   ├── association_results.json
│   └── manhattan_plot.png
└── integration/
    └── multiomics_correlations.json
```

## Testing Examples

Examples are tested via `tests/test_examples.py`:

```bash
# Run all example tests
python -m pytest tests/test_examples.py -v

# Run specific example test
python -m pytest tests/test_examples.py::TestExamples::test_dna_sequences -v
```

## Contributing Examples

### Example Guidelines

1. **Educational Value**: Each example should teach a specific METAINFORMANT concept
2. **Self-Contained**: Examples should run independently with minimal setup
3. **Well-Documented**: Include comprehensive comments and docstrings
4. **Follow Patterns**: Use core utilities (io, paths, logging) consistently
5. **Output to `output/examples/`**: Save results with informative names

### Adding New Examples

1. **Choose Topic**: Identify a useful pattern from METAINFORMANT modules
2. **Create Example**: Write complete, runnable code following the template
3. **Add Documentation**: Include README and inline comments
4. **Test Example**: Ensure example works and is well-tested
5. **Update Index**: Add to examples/README.md and domain README

### Example Template

```python
#!/usr/bin/env python3
"""Brief description of what this example demonstrates.

This example shows:
- Key concept 1
- Key concept 2
- Key concept 3

Usage:
    python examples/domain/example_name.py [--options]
"""

from __future__ import annotations

from pathlib import Path
from metainformant.core import io, paths, logging
from metainformant.domain import module_function

def main():
    """Main example function."""
    # Setup
    logger = logging.get_logger(__name__)
    output_dir = Path("output/examples/domain")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Example code with clear comments
    result = module_function(...)

    # Save results
    io.dump_json(result, output_dir / "result.json")

    print(f"Example complete. Results saved to {output_dir}")

if __name__ == "__main__":
    main()
```

## Related Documentation

- **Scripts Documentation**: [`scripts/README.md`](../scripts/README.md) - Production workflow scripts
- **Core Documentation**: [`docs/core/`](../docs/core/) - Core utilities documentation
- **Domain Documentation**: [`docs/<domain>/`](../docs/) - Module-specific documentation
- **Testing**: [`docs/testing.md`](../docs/testing.md) - Testing patterns and examples

These examples provide practical guidance for learning METAINFORMANT across all biological analysis domains.


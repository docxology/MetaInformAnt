# AI Agents in Simulation Scripts Development

This document outlines AI assistance in developing METAINFORMANT's simulation utility scripts and automation tools.

## AI Contributions

### Script Architecture
**Code Assistant Agent** designed:
- Modular simulation script organization
- Consistent command-line interfaces for simulation tools
- Cross-platform compatibility approaches
- Integration with domain-specific simulation modules

### Automation Features
**Code Assistant Agent** implemented:
- Comprehensive simulation orchestration scripts (`simulate_*.py`)
- Unified simulation runner (`run_simulation.py`)
- Domain-specific simulation utilities
- Automated testing and validation of simulated data

### User Experience
**Documentation Agent** contributed to:
- Simulation script documentation
- Usage examples and best practices
- Integration guides for synthetic data generation
- Troubleshooting and validation instructions

## Script Categories

### Core Simulation Scripts
**Code Assistant Agent** created:
- `run_simulation.py`: Unified simulation workflow orchestrator
- `_common.py`: Shared utilities for simulation scripts
- Domain-specific simulation modules for all METAINFORMANT modules

### Domain-Specific Simulation
**Code Assistant Agent** implemented simulation scripts for:
- **DNA**: Sequence and variant simulation (`simulate_dna.py`)
- **RNA**: Expression count simulation (`simulate_rna.py`)
- **Protein**: Protein sequence simulation (`simulate_protein.py`)
- **Epigenome**: Epigenetic modification simulation (`simulate_epigenome.py`)
- **GWAS**: Genotype-phenotype simulation (`simulate_gwas.py`)
- **Math**: Mathematical model simulation (`simulate_math.py`)
- **ML**: Machine learning data simulation (`simulate_ml.py`)
- **Multi-omics**: Integrated multi-omics simulation (`simulate_multiomics.py`)
- **Networks**: Biological network simulation (`simulate_networks.py`)
- **Ontology**: Gene ontology simulation (`simulate_ontology.py`)
- **Phenotype**: Trait simulation (`simulate_phenotype.py`)
- **Ecology**: Community simulation (`simulate_ecology.py`)
- **Quality**: Quality metric simulation (`simulate_quality.py`)
- **Single-cell**: scRNA-seq simulation (`simulate_singlecell.py`)
- **Information**: Information theory simulation (`simulate_information.py`)
- **Life Events**: Event sequence simulation (`simulate_life_events.py`)
- **Visualization**: Plot simulation (`simulate_visualization.py`)
- **Core**: Infrastructure simulation (`simulate_core.py`)

## Orchestrator Design Principles

All simulation scripts follow consistent patterns:
- Standard argument parsing with --model, --output, --seed, --verbose
- Integration with `metainformant.core` utilities (logging, paths, I/O)
- Output directory defaulting to `output/simulation/`
- Parameter validation and error handling
- Reproducible random seed management
- Comprehensive documentation and examples

## Development Standards

### Script Quality
- Executable permissions and proper shebang lines
- Comprehensive error handling and exit codes
- Consistent option parsing and help messages
- Cross-platform shell compatibility
- No hardcoded absolute paths

### Repository Organization
- Simulation scripts organized by biological domain
- Clear separation from output directories
- Integration with testing infrastructure
- Reproducible simulation results

### Maintenance Practices
- Regular updates with simulation method improvements
- Testing on multiple environments
- Clear documentation of simulation parameters
- Version tracking for significant changes

## Integration

Simulation scripts integrate with:
- **Core utilities**: Using `metainformant.core` for I/O and path management
- **Domain modules**: Calling simulation functions from respective modules
- **Testing framework**: Providing synthetic data for validation
- **Output directories**: Writing results to `output/simulation/`

## Simulation Script Centralization Philosophy

### Methods vs Data Separation
- **Scripts** (simulation methods) centralized in `scripts/simulation/`
- **Data** (simulated outputs) distributed in `output/simulation/`
- Single source of truth for all simulation scripts
- Update once, applies everywhere
- Proper version control and tracking

### Benefits
- ✅ No simulation script duplication across domains
- ✅ Clean output directories (data only)
- ✅ Easy maintenance and updates
- ✅ Consistent behavior across all simulations

## Related Documentation

For detailed AI contributions to specific simulation domains:

### Core Simulation
- [`scripts/simulation/README.md`](../simulation/README.md) - Simulation scripts overview

### Domain-Specific Simulation
- Additional AGENTS.md files available in `src/metainformant/*/AGENTS.md` for domain-specific simulation implementations

This simulation script collection provides essential tooling for METAINFORMANT development and synthetic data generation workflows.


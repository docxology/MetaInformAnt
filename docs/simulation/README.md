# Simulation Documentation

This directory contains comprehensive documentation for METAINFORMANT's synthetic data generation and agent-based modeling capabilities.

## Overview

The simulation domain provides tools for generating synthetic biological data and modeling complex biological systems with interacting agents.

## Documentation Files

### Core Simulation Tools
- **`index.md`**: Simulation domain overview and module index
- **`sequences.md`**: Synthetic DNA/RNA/protein sequence generation
- **`rna_counts.md`**: Gene expression count simulation
- **`agents.md`**: Agent-based modeling and GridWorld framework

## Related Source Code

- See `src/metainformant/simulation/` for implementation details
- See `tests/test_simulation_*.py` for comprehensive test coverage
- See `src/metainformant/simulation/README.md` for module-specific documentation

## Usage Examples

The simulation domain supports synthetic data generation:

```python
from metainformant.simulation import sequences, rna

# Generate synthetic sequences
dna_seq = sequences.generate_random_dna(1000)
protein_seq = sequences.generate_random_protein(200)

# Simulate expression data
counts = rna.simulate_counts_negative_binomial(100, 6)

# Agent-based modeling
from metainformant.simulation import GridWorld
world = GridWorld(10, 10, num_agents=5)
world.step()
```

## Integration

Simulation tools integrate with:
- **DNA/RNA analysis** for synthetic data validation
- **Machine learning** for training data generation
- **Mathematical models** for theoretical validation
- **Visualization** for simulation result plotting

## Testing

Comprehensive tests ensure simulation reliability:
- Sequence generation algorithm validation
- Expression distribution correctness
- Agent behavior modeling accuracy
- Integration testing with analysis modules

## Contributing

When adding new simulation functionality:
1. Update algorithm documentation
2. Add comprehensive validation tests
3. Ensure reproducible results with proper seeding
4. Update integration examples

This documentation provides complete coverage of METAINFORMANT's simulation capabilities.

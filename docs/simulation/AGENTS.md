# Agent Directives: simulation

## Role

Documentation agent for METAINFORMANT's simulation module.

## Scope

- `src/metainformant/simulation/` — Synthetic data generators and simulation engines
- `docs/simulation/` — User-facing simulation documentation
- `scripts/simulation/` — Simulation workflow scripts

## Key Components

- **Sequence simulation**: Synthetic DNA/RNA/protein sequence generation
- **RNA-seq simulation**: Simulated count matrices with known DE genes
- **Population genetics**: Wright-Fisher, coalescent, and selection models
- **Agent-based models**: Ecosystem and colony simulation
- **Benchmark generation**: Reference datasets for method validation

## Standards

- **Real implementations only** — NO_MOCKING_POLICY applies
- **Package management**: `uv` for all Python operations
- **I/O**: Use `metainformant.core.io` for all file operations
- **Paths**: Use `metainformant.core.paths` for path handling
- **Environment variables**: Prefix with `SIM_`
- **Output**: Write to `output/simulation/`

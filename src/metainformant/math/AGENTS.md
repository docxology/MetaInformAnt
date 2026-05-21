# Agent Directives: math

**Context**: Mathematical biology and theoretical modeling module for METAINFORMANT.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `bayesian/` — exports: `inference`
- `core/` — exports: `utilities`, `visualization`
- `decision_theory/` — exports: `ddm`
- `epidemiology/` — exports: `models`
- `evolutionary_dynamics/` — exports: `core`, `egt`
- `perception/` — exports: `psychophysics`, `signal_detection`
- `population_genetics/` — exports: `coalescent`, `core`, `demography`, `effective_size`, `fst`
- `quantitative_genetics/` — exports: `core`, `price`

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/math/](../../../docs/math/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Core module**: [../core/AGENTS.md](../core/AGENTS.md)
- **Networks module**: [../networks/AGENTS.md](../networks/AGENTS.md)
- **Simulation module**: [../simulation/AGENTS.md](../simulation/AGENTS.md)

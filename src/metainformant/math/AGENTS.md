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
- Use `metainformant.core.io` for file operations — never `import json` directly
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management

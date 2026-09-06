# Specification: math

## Scope
Mathematical biology and theoretical modeling module for METAINFORMANT.

## Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## Data Structures
- **Sub-packages**: bayesian, core, decision_theory, epidemiology, evolutionary_dynamics, perception, population_genetics, quantitative_genetics
- **Top-level module**: `popgen.py` (compatibility exports for population genetics helpers)
- **Key Concepts**: Refer to Pydantic models in source.

## API Definition
### Exports
- `__init__.py` (re-exports all sub-packages listed above)

# Mathematical Biology Documentation

This directory contains comprehensive documentation for METAINFORMANT's mathematical and theoretical biology tools.

## Overview

The math domain provides implementations of key mathematical frameworks in theoretical biology, population genetics, and decision-making processes.

## Documentation Files

### Core Mathematical Biology
- **`index.md`**: Mathematical biology overview and module index
- **`price.md`**: Price equation and selection analysis
- **`selection.md`**: Kin selection and multilevel selection theory
- **`ddm.md`**: Drift-diffusion models for decision making

### Population Genetics
- **`popgen.md`**: Population genetics theory and models
- **`coalescent.md`**: Coalescent theory and genealogy
- **`ld.md`**: Linkage disequilibrium analysis and genetic maps
- **`epidemiology.md`**: Disease modeling and epidemiology
- **`dynamics.md`**: Population dynamics and ecological modeling

## Related Source Code

- See `src/metainformant/math/` for implementation details
- See `tests/test_math_*.py` for comprehensive test coverage
- See `src/metainformant/math/README.md` for module-specific documentation

## Usage Examples

The math domain supports theoretical and quantitative biological analysis:

```python
from metainformant.math import price_equation, kin_selection_response

# Price equation analysis
cov, trans, total = price_equation.decompose_change(
    trait_values=[1.0, 1.2, 0.9],
    fitness_values=[0.2, 0.4, 0.1],
    trait_changes=[0.25, 0.35, 0.15]
)

# Kin selection analysis
response = kin_selection_response(
    relatedness=0.5,
    benefit=0.4,
    cost=0.1
)
```

## Integration

Mathematical biology integrates with:
- **DNA analysis** for evolutionary modeling
- **Simulation** for theoretical validation
- **Statistical methods** for hypothesis testing
- **Visualization** for mathematical model plotting

## Testing

Comprehensive tests ensure mathematical correctness:
- Verification against known analytical solutions
- Numerical stability testing
- Performance benchmarking
- Edge case handling validation

## Contributing

When adding new mathematical models:
1. Update theoretical background documentation
2. Add comprehensive algorithm tests
3. Include references to original literature
4. Ensure numerical stability and accuracy

This documentation provides complete coverage of METAINFORMANT's mathematical biology capabilities.

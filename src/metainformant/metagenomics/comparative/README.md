# Comparative

Comparative metagenomics analysis providing differential abundance testing, indicator species analysis, effect size ranking, and ML-based biomarker discovery for microbiome studies.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports differential_abundance submodule |
| `differential_abundance.py` | ALDEx2, ANCOM, IndVal, LEfSe, and random forest methods |

## Key Functions

| Function | Description |
|----------|-------------|
| `differential_abundance.differential_abundance()` | Run differential abundance testing between groups |
| `differential_abundance.clr_transform()` | Centered log-ratio transformation for compositional data |
| `differential_abundance.indicator_species()` | IndVal analysis for group-associated taxa |
| `differential_abundance.effect_size_analysis()` | LEfSe-style linear discriminant analysis effect sizes |
| `differential_abundance.biomarker_discovery()` | ML-based biomarker discovery using feature importance |

## Usage

```python
from metainformant.metagenomics.comparative import differential_abundance

results = differential_abundance.differential_abundance(
    counts, groups, taxa_names, method="aldex2_like"
)
indicators = differential_abundance.indicator_species(counts, groups, taxa_names, seed=42)
biomarkers = differential_abundance.biomarker_discovery(counts, groups, taxa_names, seed=42)
```

# Integration

Cross-domain utilities for DNA sequence analysis: DNA-RNA integration (ORF finding, splice site prediction, regulatory element detection, expression correlation) and the DNA-simulation adapter that supplies simulated null distributions for the neutrality tests.

Cross-domain imports are allowed only from adapter surfaces: `integration` and `workflow` modules, CLI entry points and the MCP tool layer (`scripts/quality/check_module_boundaries.py`). This package is where `dna` meets the `simulation` package.

## Contents

| File | Purpose |
|------|---------|
| `rna.py` | DNA-RNA integration: gene structure, regulatory elements, and expression correlation |
| `neutrality_null.py` | DNA-simulation adapter: Monte Carlo null distributions for the neutrality tests of `dna.population.analysis`. `analysis` imports it lazily, when a caller asks for `p_value_method="simulation"`; `__init__.py` does not import it, so `import metainformant.dna` does not load the simulation package |

## Key Functions

| Function | Description |
|----------|-------------|
| `find_open_reading_frames()` | Locate ORFs in a DNA sequence |
| `predict_transcription_start_sites()` | Score potential TSS locations using promoter motifs |
| `find_transcription_factor_binding_sites()` | Scan DNA for known TF motif matches |
| `calculate_codon_usage_bias()` | Compute codon bias statistics for a coding sequence |
| `analyze_gene_structure()` | Identify exons, introns, and UTR regions |
| `correlate_dna_with_rna_expression()` | Correlate sequence features with expression levels |
| `predict_splice_sites()` | Detect canonical donor/acceptor splice signals |
| `analyze_regulatory_elements()` | Find promoters, enhancers, and CpG islands |
| `integrate_dna_rna_data()` | Merge DNA and RNA datasets into a unified result |
| `simulated_neutrality_tests()` | Observed statistic and two-sided Monte Carlo p-value of each requested test (`tajima_d`, `fu_li_d_star`, `fu_li_f_star`, `fay_wu_h`) against a simulated neutral coalescent null |
| `monte_carlo_two_sided_p()` | Twice the smaller one-sided tail of the simulated statistics, with the +1 correction, capped at 1 |

## Usage

```python
from metainformant.dna.integration.rna import (
    find_open_reading_frames,
    correlate_dna_with_rna_expression,
    integrate_dna_rna_data,
)

orfs = find_open_reading_frames("ATGGCCATTGTAATGGGCCGCTGAAAGGGT")
merged = integrate_dna_rna_data(dna_data, rna_data)
```

## Simulated neutrality nulls

`simulated_neutrality_tests(sequences, tests, *, n_simulations, random_state=None, outgroup=None)` returns `{test: (statistic, p_value)}`:

```python
from metainformant.dna.integration.neutrality_null import simulated_neutrality_tests

results = simulated_neutrality_tests(
    sequences, ["tajima_d", "fu_li_f_star"], n_simulations=2000, random_state=1
)
tajima_d, p_value = results["tajima_d"]
```

- **Null**: the standard neutral coalescent for the same number of sequences (constant size, infinite sites, no recombination), drawn by `simulation.models.popgen.simulate_coalescent_spectra` with exactly as many mutations in every replicate as the observed spectrum has (the `ms -s` approach).
- **Groups**: `FOLDED_TESTS` (`tajima_d`, `fu_li_d_star`, `fu_li_f_star`) score the folded spectrum and need no outgroup; `POLARIZED_TESTS` (`fay_wu_h`) need `outgroup` and use the polarized spectrum. Each group is scored against its own simulation.
- **p-value**: twice the smaller one-sided tail, `(1 + #{replicates at least as extreme}) / (n_simulations + 1)`, capped at 1. A sample without mutations has statistic 0 and p-value 1.
- **Errors**: an unknown test name, or a polarized test without an outgroup, raises `ValueError`.

Callers normally use the public options of `dna.population.analysis` (`p_value_method="simulation"`, `n_simulations`, `random_state`) rather than this module; see [Neutrality p-values](../../../../docs/dna/population.md#neutrality-p-values-analysis).

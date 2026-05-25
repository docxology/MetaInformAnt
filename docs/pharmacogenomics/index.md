# Pharmacogenomics Index

High‑level map of the module's documentation.

## Getting started

- **[GETTING_STARTED.md](GETTING_STARTED.md)** — install, load data, one‑liner demo
- **[io.md](../io.md)** — platform‑specific loaders (refer to general spatial docs)

## Core documentation

| File | Purpose |
|------|---------|
| [ARCHITECTURE.md](ARCHITECTURE.md) | Layer diagram, data flow, design decisions |
| [API_REFERENCE.md](API_REFERENCE.md) | Complete function and class catalogue |
| [SPEC.md](SPEC.md) | Technical specification: data model, algorithms, schemas |
| [CONFIGURATION.md](CONFIGURATION.md) | All `spatial.*` keys, env vars, YAML examples |
| [PERFORMANCE.md](PERFORMANCE.md) | Benchmarks, scaling, tuning parameters |
| [TROUBLESHOOTING.md](TROUBLESHOOTING.md) | Symptom → fix table, diagnostics |

## Topic guides

| Guide | Content |
|-------|---------|
| [acmg.md](acmg.md) | ACMG variant classification for PGx |
| [cpic.md](cpic.md) | CPIC guideline loading & lookup |
| [drug_interactions.md](drug_interactions.md) | DDI analysis, polypharmacy risk |
| [reporting.md](reporting.md) | Clinical report generation (text/JSON/HTML) |
| [star_alleles.md](star_alleles.md) | Matching engine internals |

## Development

- **[AGENTS.md](AGENTS.md)** — contributor guide: adding genes, tests, style
- **[PAI.md](PAI.md)** — pre‑implementation: problem, prior art, API surface

## Examples

- **[EXAMPLES.md](EXAMPLES.md)** — 25 runnable code snippets covering single‑sample,
  batch, integration, custom tables, reporting.

## Related integration

- **[../INTEGRATION.md](../../INTEGRATION.md)** — DNA→GWAS→Viz multi‑module pipelines
  (PHGx can plug in as downstream phenotype layer).

## Module API summary

```python-snippet
from metainformant.pharmacogenomics import (
    call_star_alleles,
    determine_diplotype,
    classify_phenotype,
    predict_metabolizer,
    lookup_drug_gene,
    get_dosing_recommendation,
)

from metainformant.pharmacogenomics.annotations import cpic, pharmgkb
from metainformant.pharmacogenomics.clinical import generate_clinical_report, classify_variant_acmg
from metainformant.pharmacogenomics.interaction import analyze_drug_gene_interactions
```

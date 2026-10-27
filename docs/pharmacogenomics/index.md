# Pharmacogenomics Index

High‑level map of the module's documentation.

## Getting started

- **[GETTING_STARTED.md](GETTING_STARTED.md)** — install, load data, one‑liner demo
- **[CAPABILITIES.md](CAPABILITIES.md)** — function and class catalogue by submodule

## Core documentation

| File | Purpose |
|------|---------|
| [ARCHITECTURE.md](ARCHITECTURE.md) | Layer diagram, data flow, design decisions |
| [CAPABILITIES.md](CAPABILITIES.md) | Complete function and class catalogue |
| [SPEC.md](SPEC.md) | Technical specification: data model, algorithms, schemas |
| [CONFIGURATION.md](CONFIGURATION.md) | Function arguments that customise tables, guidelines and reports; logging (no env vars or config keys) |
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
- **[PAI.md](PAI.md)** — pre‑implementation: problem, prior art, data scale, user journeys

## Examples

- **[EXAMPLES.md](EXAMPLES.md)** — 22 code snippets covering single‑sample,
  batch, integration, custom tables, reporting.

## Related integration

- **[INTEGRATION.md](../INTEGRATION.md)** — DNA→GWAS→Viz multi‑module pipelines
  (it does not cover pharmacogenomics; pass variant identifiers from such a pipeline to
  `call_star_alleles`).

## Module API summary

```python-snippet
from metainformant.pharmacogenomics import predict_metabolizer
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles
from metainformant.pharmacogenomics.alleles.diplotype import determine_diplotype
from metainformant.pharmacogenomics.alleles.phenotype import classify_phenotype
from metainformant.pharmacogenomics.annotations.cpic import lookup_drug_gene, get_dosing_recommendation

from metainformant.pharmacogenomics.annotations import cpic, pharmgkb
from metainformant.pharmacogenomics.clinical.reporting import generate_clinical_report
from metainformant.pharmacogenomics.clinical.pathogenicity import classify_variant_acmg
from metainformant.pharmacogenomics.clinical.drug_interaction import analyze_drug_gene_interactions
```

# Getting Started: Pharmacogenomics

Pharmacogenomics (PGx) predicts drug response from genetic variation using star allele
conventions and CPIC clinical guidelines. This module handles CYP450 genes (CYP2D6,
CYP2C19, CYP2C9, etc.), TPMT, NUDT15, SLCO1B1, and many more.

## Quick example

```python
from metainformant.pharmacogenomics import call_star_alleles, determine_diplotype, classify_phenotype
from metainformant.pharmacogenomics.annotations.cpic import get_dosing_recommendation

variants = {'rs1065852', 'rs3892097', 'rs5030655'}  # CYP2D6 variants
alleles = call_star_alleles(variants, gene='CYP2D6')
# → [StarAllele(name='*1', ...), StarAllele(name='*4', ...)]

diplotype = determine_diplotype(alleles[0].name, alleles[1].name, gene='CYP2D6')
# → '*1/*4'  (string)

phenotype = classify_phenotype(diplotype, gene='CYP2D6')
# → {'phenotype': MetabolizerPhenotype.INTERMEDIATE, 'activity_score': 1.0}

rec = get_dosing_recommendation('codeine', phenotype['phenotype'].value, gene='CYP2D6')
print(rec['recommendation'])
# → "Use alternative to codeine" for IM/PM
```

## Installation

The pharmacogenomics module ships with built-in allele definitions and CPIC guidelines.
For full functionality (reports, plots), install optional extras:

```bash
uv pip install metainformant[pharmacogenomics,visualization]
```

Optional extra `deconvolution` enables spatial cell-type deconvolution (not needed for PGx).

## Data sources

| Resource | What it provides | Built-in version | Update method |
|----------|------------------|------------------|---------------|
| CPIC guidelines | Drug–gene dosing recommendations per phenotype | v3.0 (2024-01) | `load_cpic_guidelines(path=...)` or `PG_CPIC_URL` env |
| PharmGKB annotations | Evidence-backed variant annotations | Monthly cache (~9 MB) | Auto-refresh or `refresh_pharmgkb_cache()` |
| Star allele tables | Gene→star-allele definitions with activity scores | 28 genes | Custom JSON via `load_allele_definitions()` |

## Core concepts

**Star alleles** – Standardized haplotype labels (e.g., `*1`, `*4`, `*2A`) that represent
specific variant combinations. The matcher maps observed rsIDs to the set of alleles
that could explain them.

**Diplotype** – The unordered pair of star alleles present on the two chromosomes
(e.g., `*1/*4`). Activity score is the sum of each allele's numeric activity.

**Phenotype** – Clinical interpretation: Poor (PM), Intermediate (IM), Normal (NM),
Rapid (RM), or Ultrarapid (UM) metabolizer. Thresholds are gene-specific.

**CPIC level** – Evidence strength: A (strong), B (moderate), C/D (optional/negative).

## Batch processing

```python-snippet
import pandas as pd
from metainformant import parallel, pharmacogenomics

df = pd.read_csv('cohort.csv')  # requires 'rsids' column (semicolon-separated)
def predict(row):
    return pharmacogenomics.predict_metabolizer(
        set(str(row['rsids']).split(';')), gene='CYP2D6'
    )

results = parallel.map(predict, df.to_dict('records'), workers=8, chunk_size=10000)
df['CYP2D6_phenotype'] = [r['phenotype'] for r in results]
```

For >100k samples increase `chunk_size` to reduce overhead; for massive cohorts switch
to Dask:

```python-snippet
import dask.dataframe as dd
from dask.distributed import Client
client = Client(n_workers=8)
ddf = dd.read_csv('cohort_*.csv')
ddf['phenotype'] = ddf.map_partitions(
    lambda part: [pharmacogenomics.predict_metabolizer(set(vs), gene='CYP2D6')['phenotype']
                  for vs in part['rsids']],
    meta=('phenotype','object')
)
```

## Clinical reporting

```python
from metainformant.pharmacogenomics.clinical.reporting import generate_clinical_report, export_report

patient = {
    'patient_id': 'PT001',
    'name': 'Sample',
    'dob': '1970-01-01',
    'medications': ['codeine', 'clopidogrel'],
    'genotype': {
        'CYP2D6': {'diplotype': '*1/*4', 'phenotype': 'Intermediate Metabolizer'},
        'CYP2C19': {'diplotype': '*1/*1', 'phenotype': 'Normal Metabolizer'},
    },
}
report = generate_clinical_report(patient)
print(export_report(report, format='text'))   # also 'json' and 'html'
```

Output includes header, genotype table per gene, drug recommendations with CPIC
levels, DDI warnings, and disclaimer.

## Configuration

Key settings are under `pharmacogenomics.*` and environment variables starting with
`PG_`:

| Key | Env var | Default | Purpose |
|-----|---------|---------|---------|
| `algorithm` | `PG_ALGORITHM` | `auto` | Matching: `fast` \| `balanced` \| `accurate` |
| `chunk_size` | `PG_CHUNK_SIZE` | `5000` | Batch size for `parallel.map` |
| `parallel` | `PG_PARALLEL` | `false` | Enable process-based parallelism |
| `allele_definitions_dir` | `PG_ALLELE_DEFINITIONS_DIR` | built-in | Directory with custom JSON allele tables |
| `cpic_guidelines_url` | `PG_CPIC_GUIDELINES_URL` | built-in | Override CPIC JSON source URL |
| `pharmgkb.tsv_path` | `PG_PHARMGKB_TSV_PATH` | cache dir | Local PharmGKB TSV for offline use |
| `logging` | `PG_LOGGING` | `INFO` | Module log level (`DEBUG` verbose) |

Set via Python:

```python
from metainformant import config
config.set('pharmacogenomics.algorithm', 'balanced')
config.set('pharmacogenomics.chunk_size', 20000)
```

Or YAML (`~/.hermes/config.yaml`):

```yaml
pharmacogenomics:
  algorithm: balanced
  chunk_size: 5000
  parallel: false
```

See [CONFIGURATION.md](CONFIGURATION.md) for the full reference.

## What genes are supported?

The built-in allele definitions cover ~30 pharmacogenes:

- CYP2D6 (140 alleles), CYP2C19 (30), CYP2C9 (19), CYP3A5 (9)
- TPMT (13), NUDT15 (9), UGT1A1 (10)
- SLCO1B1 (12), VKORC1 (8), DPYD (22)
- ... and more

List available genes:

```python
from metainformant.pharmacogenomics.alleles.star_allele import _BUILTIN_ALLELE_DEFINITIONS
print(sorted(_BUILTIN_ALLELE_DEFINITIONS.keys()))
```

## Interactive pipeline (Hermes CLI)

If you run `hermes` in interactive mode, you can:

1. `/tool pharmacogenomics call_star_alleles --gene CYP2D6 --variants rs1065852,rs3892097`
2. `/tool pharmacogenomics classify_phenotype --diplotype *1/*4 --gene CYP2D6`
3. `/tool pharmacogenomics cpic get_dosing --drug codeine --phenotype "Intermediate Metabolizer"`

Type `/help` inside Hermes to explore tool invocations.

## Next steps

- **[API_REFERENCE.md](API_REFERENCE.md)** — Full class and function catalogue
- **[EXAMPLES.md](EXAMPLES.md)** — 20 detailed usage patterns
- **[PERFORMANCE.md](PERFORMANCE.md)** — Benchmarks and scaling tips
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** — Known failure modes and fixes
- **[CONFIGURATION.md](CONFIGURATION.md)** — All config keys and environment variables

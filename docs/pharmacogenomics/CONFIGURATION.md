# Configuration: Pharmacogenomics

## Environment variables (top precedence)

| Variable | Default | Description |
|----------|---------|-------------|
| `PG_ALGORITHM` | `auto` | Algorithm: `fast` | `balanced` | `accurate` | `auto` |
| `PG_CHUNK_SIZE` | `5000` | Batch size for streaming (>20k records recommended 20k) |
| `PG_PARALLEL` | `false` | Enable multiprocessing — `true` or `false` |
| `PG_ALLELE_DEFINITIONS_DIR` | built-in | Path to directory with per-gene `*.json` custom tables |
| `PG_CPIC_GUIDELINES_URL` | built-in | Override URL to CPIC JSON (rare) |
| `PG_PHARMGKB_API_TOKEN` | none | Optional token for authenticated PharmGKB REST queries |
| `PG_LOGGING` | `INFO` | Module log-level — `DEBUG` shows allele-matching trace |

## Python API

```python
from metainformant import config
config.set('pharmacogenomics.algorithm', 'fast')
config.set('pharmacogenomics.chunk_size', 20000)
config.set('pharmacogenomics.parallel', True)
config.set('pharmacogenomics.allele_definitions_dir', '/data/alleles/')
```

## YAML (~/.hermes/config.yaml)

```yaml
pharmacogenomics:
  algorithm: balanced
  chunk_size: 5000
  parallel: false
  allele_definitions_dir: null
  cpic_guidelines_url: null
```

## Algorithm detail

| value | behavior | per-sample runtime (typical) |
|-------|----------|-----------------------------|
| `fast` | Greedy first subset match; no fragmentation check | 5–8 ms |
| `balanced` | Greedy + optimal count when 2–3 candidate alleles tie | 9–12 ms |
| `accurate` | Exhaustive enumeration of all minimal allele subsets (≤6 options) | 25–40 ms |
| `auto` | Heuristic: ≤2 alleles → fast; ≥3 → balanced | — |

## Chunk size & parallel streaming

When processing large DataFrames, use `parallel.map` or `parallel.batch_map`.

| cohort N | recommended `chunk_size` | live memory impact |
|----------|----------------------|-------------------|
| < 5 000  | 1 000                | negligible        |
| 5 K–100 K | 5 000               | ~60 MB            |
| 100 K–1 M | 20 000              | ~120 MB           |
| > 1 M    | 50 000              | ~200 MB           |

Example workflow:

```python
from metainformant import parallel, pharmacogenomics
def allele_predict(row):
    return pharmacogenomics.predict_metabolizer(set(row['rsids'].split(';')), gene='CYP2D6')

results = parallel.map(allele_predict, df.to_dict('records'),
                        workers=8, chunk_size=20000)
```

Workers inherit config; ensure `PG_PARALLEL=true` for best throughput.

## Custom star allele tables

JSON schema per gene file:

```json
{
  "gene": "CYP2E1",
  "star_alleles": [
    {
      "name": "*1",
      "variants": ["rs2031920"],
      "function": "normal",
      "activity": 1.0
    },
    {
      "name": "*5B",
      "variants": ["rs3813867"],
      "function": "decreased",
      "activity": 0.5
    }
  ]
}
```

Required allele keys: `name` (star allele designation), `variants` (list of dbSNP rsIDs used for matching), `function` (free-text), `activity` (float; usually 0.0 .5 .5 1.0 1.5+).

Load from a directory that holds one or more json files:

```python
from metainformant.pharmacogenomics.alleles.star_allele import load_all_allele_definitions
load_all_allele_definitions('/my/custom/path/')
```

Files are indexed by top-level `"gene"` key; you may store many genes in one directory. After load, every subsequent call uses those definitions.

## CPIC guidelines override

The built-in snapshot is CPIC v3.0 (2024-01). To replace or upgrade:

```python
from metainformant.pharmacogenomics.annotations.cpic import load_cpic_guidelines
load_cpic_guidelines(path='/data/cpic_guidelines_2025.json')
```

The JSON file must be a list of guideline dictionaries with these required keys:

- `drug` (str)
- `gene` (str)
- `cpic_level` (`'A'`|`'B'`|`'C'`|`'D'`)
- `phenotypes` — map phenotype string → `{recommendation, classification, evidence}`
- `guideline_url` (optional)

Alternatively, set `PG_CPIC_GUIDELINES_URL` to a downloadable URL; the module will fetch and cache it once.

## PharmGKB cache management

First network call downloads the latest PharmGKB annotation TSV (~9 MB) into `~/.cache/metainformant/pharmgkb/annotations.tsv`. Future calls are offline.

To force re-download (e.g., monthly update):

```python
from metainformant.pharmacogenomics.annotations.pharmgkb import refresh_pharmgkb_cache
refresh_pharmgkb_cache(force=True)
```

For air-gapped systems, pre-populate the cache directory or set `PG_PHARMGKB_TSV_PATH` to a local file. The file must be a TSV with the official PharmGKB columns.

## Validation & sanity

```python
from metainformant import config
errors = config.validate_schema(module='pharmacogenomics')
if errors:
    for e in errors:
        print(f"{e.path}: {e.message}")
```

Common validation errors: unknown algorithm name; chunk_size non-integer; `parallel` non-boolean; `allele_definitions_dir` not a readable directory.

To revert all keys to built-in defaults:

```python
config.reset(module='pharmacogenomics')
```

## Global logging

Set `PG_LOGGING=DEBUG` to see detailed allele-matching traces: every candidate allele evaluated, scores, and tie-breaking decisions. This is invaluable when troubleshooting a new rare variant.

```bash
export PG_LOGGING=DEBUG
python predict.py
```

## Precedence hierarchy

1. Function-call kwargs (win everything)
2. Per-process `config.set()`
3. `config.yaml` file
4. Environment variables
5. Built-in module defaults

## Scope table — what's configurable

| config key | env var | purpose |
|------------|---------|---------|
| `pharmacogenomics.algorithm` | `PG_ALGORITHM` | matching strategy |
| `pharmacogenomics.chunk_size` | `PG_CHUNK_SIZE` | batch size for streaming |
| `pharmacogenomics.parallel` | `PG_PARALLEL` | enable multiprocessing |
| `pharmacogenomics.allele_definitions_dir` | `PG_ALLELE_DEFINITIONS_DIR` | custom JSON directory |
| `pharmacogenomics.cpic_guidelines_url` | `PG_CPIC_GUIDELINES_URL` | CPIC JSON source |
| `pharmacogenomics.pharmgkb.tsv_path` | `PG_PHARMGKB_TSV_PATH` | local PharmGKB TSV |
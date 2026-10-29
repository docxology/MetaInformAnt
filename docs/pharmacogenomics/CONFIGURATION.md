# Configuration: Pharmacogenomics

The pharmacogenomics module has no environment variables, no `pharmacogenomics.*`
configuration keys and no YAML section: nothing under `src/metainformant/pharmacogenomics`
reads the environment or `metainformant.core.utils.config`. Behaviour is configured with
function arguments, listed in the [scope table](#scope-table--whats-configurable) at the
end of this page. (The `PG_*` variables that
`core.utils.config.load_postgres_config_from_env()` reads — `PG_HOST`, `PG_PORT`,
`PG_DATABASE`, `PG_USER`, `PG_PASSWORD` — configure PostgreSQL and are unrelated to this
module.)

To drive those arguments from a file, load the file yourself (for example with
`metainformant.core.utils.config.load_mapping_from_file`) and pass the values explicitly.

## Not implemented

None of the following exists; nothing in the module reads or provides it:

- Environment variables `PG_ALGORITHM`, `PG_CHUNK_SIZE`, `PG_PARALLEL`,
  `PG_ALLELE_DEFINITIONS_DIR`, `PG_CPIC_GUIDELINES_URL`, `PG_PHARMGKB_API_TOKEN`,
  `PG_PHARMGKB_TSV_PATH` and `PG_LOGGING`.
- `pharmacogenomics.*` config keys, a `pharmacogenomics:` YAML section, and the
  `config.set()`, `config.validate_schema()` and `config.reset()` calls
  (`metainformant.core.utils.config` has none of them) — so there is also no precedence
  hierarchy between layers; a function argument is the only way to configure a call.
- Matching algorithms (`fast`, `balanced`, `accurate`, `auto`); `call_star_alleles` has
  the single algorithm described below.
- Chunk-size, memory and process-parallelism settings.
- Downloads and caches: the module has no network code, so no CPIC or PharmGKB data is
  fetched and there is nothing to refresh.

## Algorithm detail

`call_star_alleles(variants, gene, allele_definitions=None)` has one greedy algorithm:

1. Sort the allele definitions by number of defining variants, most specific first. The
   definition named `*1` is the reference and is never matched on its own variants.
2. Keep every definition whose defining variants are all observed and not already
   consumed by an earlier match, so overlapping definitions resolve to the most specific
   one (observing `rs4149056` and `rs2306283` in SLCO1B1 gives `*15`, not `*5` plus `*1B`).
3. Return the matches sorted by name. If nothing matched, return `[*1]` — or `[]` when
   the definitions contain no `*1`.

An empty or unrecognised variant set therefore reads as the reference allele, not as an
error. The only debug output is one line per matched allele, or one line when nothing
matched (see Global logging).

## Batch processing

There is no chunk-size setting. Map a per-sample function over a cohort with
`metainformant.core.execution.parallel`:

- `parallel.thread_map(func, items, max_workers=8)` — thread pool, one call per item,
  results in input order.
- `parallel.process_map(func, items, max_workers=None)` — process pool for CPU-bound
  work; `func` must be a top-level (picklable) function.
- `parallel.parallel_batch(func, items, batch_size=10, max_workers=4)` — splits `items`
  into batches; `func` takes a list and returns a list.

```python-snippet
from metainformant.core.execution import parallel
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles

def called_alleles(row):
    return [a.name for a in call_star_alleles(set(row['rsids'].split(';')), 'CYP2D6')]

results = parallel.thread_map(called_alleles, df.to_dict('records'), max_workers=8)
```

## Custom star allele tables

JSON file per gene (allele name → definition; the table may also be nested under the gene symbol):

```json
{
  "*1": {
    "defining_variants": [],
    "function": "Normal function",
    "activity_value": 1.0
  },
  "*5B": {
    "defining_variants": ["rs3813867"],
    "function": "Decreased function",
    "activity_value": 0.5
  }
}
```

Allele keys: `defining_variants` (list of dbSNP rsIDs or `chr:pos:ref:alt` identifiers used for matching), `function` (free-text), `activity_value` (float; usually 0.0, 0.5, 1.0, 1.5+; required in every file entry; a missing or empty value raises `ValueError`). `clinical_significance` and `evidence_level` are optional. A TSV with the columns `allele`, `defining_variants` (semicolon-separated), `function` and `activity_value` is also accepted.

Load one gene's table from a file and pass it to the allele caller:

```python-snippet
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles, load_allele_definitions
from metainformant.pharmacogenomics.alleles.diplotype import determine_diplotype

definitions = load_allele_definitions('CYP2E1', filepath='/my/custom/path/CYP2E1.json')
alleles = call_star_alleles({'rs3813867'}, 'CYP2E1', allele_definitions=definitions)

# determine_diplotype leaves alleles absent from its scoring table unscored,
# so give it the custom activity values explicitly
scores = {a.name: a.activity_value for a in definitions}
diplotype = determine_diplotype(alleles[0], '*1', 'CYP2E1', scoring_table=scores)
```

### Choosing the source: `source` and `filepath`

`load_allele_definitions(gene, source=None, filepath=None)` picks the source from its arguments
and refuses contradictory ones, so a file is never ignored and the built-in table is never
substituted for one:

| Arguments | Result |
|-----------|--------|
| neither (or `source='cpic'`) | the built-in table for the gene |
| `filepath` only (or `source='file'` with `filepath`) | the definitions in that file (`FileNotFoundError` if it is missing); a built-in gene's table is replaced for that call |
| `source='file'` without `filepath` | `ValueError` |
| `source='cpic'` with `filepath` | `ValueError` |
| any other `source` | `ValueError` |

`source` is case-insensitive.

### Phenotypes for a custom gene

Definitions are loaded per gene and used only when passed as `allele_definitions`; otherwise `call_star_alleles` uses its built-in definitions. Custom definitions and activity values do not add a phenotype definition. `classify_phenotype`, `get_phenotype_thresholds` and the metabolizer front end raise `ValueError` for a gene without its own supported phenotype bands. Returned thresholds always have `gene_specific=True`. A supported gene's diplotype carrying a custom `scoring_table` score is classified from that score; an allele missing from the applicable table instead yields an indeterminate phenotype, `activity_score=None` and its name in `unknown_alleles`.

### Genes with activity values but no calling definitions

UGT1A1 and CYP1A2 have activity values but no bundled calling definitions. `load_allele_definitions`, `call_star_alleles` and `phased_diplotype` without custom definitions raise a `ValueError` naming that gap. UGT1A1 alleles depend on promoter TA-repeat length, which presence-only variant identifiers cannot express; build its diplotype from allele names called elsewhere or supply a definition file. CYP1A2 can be scored, but its retained activity values are not sourced to a supported CPIC phenotype definition: `classify_phenotype` and `predict_metabolizer` raise for CYP1A2.

## CPIC guidelines override

The built-in guidelines are 9 level-A entries. `load_cpic_guidelines` returns the list parsed from your file; it does not replace the built-in data, so pass the list to the functions that take `guidelines`:

```python-snippet
from metainformant.pharmacogenomics.annotations.cpic import get_dosing_recommendation, load_cpic_guidelines

guidelines = load_cpic_guidelines(filepath='cpic_guidelines_2025.json')
rec = get_dosing_recommendation('codeine', 'Poor Metabolizer', guidelines=guidelines)
```

The JSON file must be a list of guideline dictionaries (an object with a `guidelines` list, or a single guideline object, is also accepted) with these keys:

- `drug` (str)
- `gene` (str)
- `cpic_level` (`'A'`|`'A/B'`|`'B'`|`'C'`|`'D'`) — `lookup_drug_gene` and `get_dosing_recommendation` raise `KeyError` when a matching entry lacks it (`list_actionable_genes` skips such an entry)
- `recommendations` — map phenotype string → `{recommendation, classification, implications}`
- `guideline_url` (optional)

A TSV with the columns `drug`, `gene`, `cpic_level` (or `level`) and `guideline_url` (or `url`) also loads, but its entries have empty `recommendations`, so it only serves `lookup_drug_gene` and `list_actionable_genes`.

`generate_clinical_report` and `analyze_drug_gene_interactions` always use the built-in guidelines; they take no `guidelines` argument.

## PharmGKB annotations

The package ships a curated set of PharmGKB clinical annotations; it does not download or cache anything, so there is nothing to refresh. `query_pharmgkb_annotations(gene=..., drug=..., variant=...)` searches the built-in annotations.

To use your own, more recent export, parse a local PharmGKB clinical-annotations TSV:

```python-snippet
from metainformant.pharmacogenomics.annotations.pharmgkb import parse_clinical_annotations
annotations = parse_clinical_annotations('/data/pharmgkb/clinical_annotations.tsv')
```

The parser reads these official PharmGKB TSV columns: `Clinical Annotation ID`, `Gene`, `Level of Evidence`, `Drug(s)`, `Phenotype Category` and `Annotation Text`. A column it does not find gives empty strings, not an error. The parsed list is returned to you; `query_pharmgkb_annotations` keeps searching only the built-in annotations.

## Global logging

The module logs through the standard `logging` logger `metainformant.pharmacogenomics` (child loggers per module), whose level is inherited from the root logger (WARNING by default). Raise it to see what the allele matcher does: at `DEBUG` it logs one line per matched allele, or a line saying that no variant allele matched and `*1` was returned. This helps when a variant set unexpectedly calls `*1`.

```python
import logging
logging.getLogger('metainformant.pharmacogenomics').setLevel(logging.DEBUG)
```

Alternatively set `CORE_LOG_LEVEL=DEBUG` in the environment and call `metainformant.core.utils.logging.configure_logging_from_env()` once at start-up; that sets the root logger level for every module, not only this one.

## Scope table — what's configurable

| What | Function and argument | Default |
|------|-----------------------|---------|
| Star-allele definitions | `call_star_alleles(variants, gene, allele_definitions=...)`; `load_allele_definitions(gene, source=None, filepath=None)` | built-in tables (`source=None` without `filepath`; `filepath` alone loads that file) |
| Allele activity scores | `determine_diplotype(allele1, allele2, gene, scoring_table=...)`; `calculate_activity_score(diplotype, gene, scoring_table=...)` | built-in tables; an unknown allele leaves score `None` |
| Allele function table for `predict_metabolizer` | `predict_metabolizer(genotype, gene, activity_scores=...)`, one gene's `{allele: value}` dict | built-in table; an unknown allele yields `indeterminate`, a `None` score and `unknown_alleles`; unsupported phenotype genes raise |
| CPIC guidelines | `lookup_drug_gene`, `get_dosing_recommendation`, `list_actionable_genes(guidelines=..., min_level='B')` | built-in guidelines; `min_level='B'` |
| Dose guidelines | `dose_adjustment(metabolizer_status, drug, dose_guidelines=...)`; the status may be `"poor"`, `"Poor Metabolizer"`, `"PM"` or a `MetabolizerPhenotype`, and a custom table may key its phenotypes in any of those forms | built-in dose table |
| Which drugs a report covers | `generate_clinical_report(patient_data, genotypes, drugs=...)` | `drugs=None`: every drug with a built-in CPIC guideline for the genotyped genes; `drugs=[]`: none |
| Drug–drug interactions | `predict_drug_interaction(drug_a, drug_b, interaction_db=...)`; `polypharmacy_risk(medications, interaction_db=...)` | `default_interaction_database()` |
| CYP inhibitor/inducer profiles | `cyp_inhibition_prediction(drug, cyp_profiles=...)` | built-in tables |
| Alternative drugs | `suggest_alternatives(drug, phenotype, alternatives_db=...)` | built-in table |
| Report format and disclaimer | `export_report(report, format='text', output_path=None)`; `add_disclaimer(report, custom_disclaimer=None)` | text; standard disclaimer |
| Plot output | every `plot_*` function: `output_path`, `title`, `figsize` | per function |

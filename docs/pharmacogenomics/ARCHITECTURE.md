# Architecture: Pharmacogenomics

![Layer diagram: Public API → Core → Data → External resources](diagram.png)

## Layer cake

1. **Public API** — `predict_metabolizer()`, `call_star_alleles()`, CPIC helpers
2. **Core logic** — diplotype determination, activity scoring, phenotype classification
3. **Allele data** — built-in JSON star-allele tables + custom loader
4. **Annotations** — CPIC guidelines, PharmGKB, drug labeling
5. **Clinical layer** — report generation, DDI analysis, ACMG classification

## Key design decisions

- **Synchronous by default** — single-sample calls return immediately (~7 ms); batch
  parallelism via `parallel.map` or manual `ThreadPoolExecutor`.
- **Deterministic matching** — given identical variant sets the same top-1 allele always wins
  (subsets scored by size, then by allele activity, then lexical).
- **Immutable allele tables** — once loaded per process they're read-only; custom tables
  replace built-ins globally (process-wide).
- **CPIC snapshot** — v3.0 bundled; upgrade via `load_cpic_guidelines(path)` without restart.
- **Cache-once** — PharmGKB TSV downloaded once to `~/.cache/metainformant/pharmgkb/`.

## Call flow

```
user code
  ↓
predict_metabolizer(variants, gene)
  ↓
call_star_alleles(variants, gene)           ← subset-matching engine
  → StarAllele objects (name, activity)
  ↓
determine_diplotype(a1.name, a2.name, gene) ← CPIC phase 1 rules
  → diplotype string
  ↓
classify_phenotype(diplotype, gene)         ← activity sum + gene-specific thresholds
  → MetabolizerPhenotype enum
  ↓
(optional) get_dosing_recommendation(...)  ← CPIC lookup
  → recommendation dict
```

## File layout

```
src/metainformant/pharmacogenomics/
  alleles/
    star_allele.py      ← subset-matching engine (~400 LOC)
    diplotype.py        ← phase-1/2 CPIC diplotype rules per gene
    phenotype.py        ← activity→phenotype threshold tables
    activity.py         ← _ACTIVITY_SCORE_TABLES (float per star allele)
  annotations/
    cpic.py             ← guideline loader + drug–gene lookup
    pharmgkb.py         ← TSV downloader + annotation query
  clinical/
    pathogenicity.py   ← ACMG classification of PGx variants
    reporting.py       ← text/JSON/HTML report builder
  interaction/
    drug_interactions.py  ← DDI analyzer based on CYP inhibition
  io/ (reserved)
```

## Thread / process safety

- All global caches (`_BUILTIN_ALLELE_DEFINITIONS`, `_ACTIVITY_SCORE_TABLES`) are
  populated on first use (lazy). In a `multiprocessing.Pool` without pre-warming,
  each child pays the ~4 ms JSON parse cost once.
- Best practice: call `load_all_allele_definitions()` in the parent *before* creating
  the pool; copy-on-write shares the tables (Linux/Unix). On Windows or when using
  `spawn` start method, do it inside the child with an `initializer`.

## Extensibility hooks

| Hook | Purpose | How to use |
|------|---------|------------|
| Custom allele JSON | Add/override star allele definitions | `load_allele_definitions(gene='XYZ', filepath='~/xyz.json')` |
| CPIC override | Use newer guidelines | `load_cpic_guidelines(path='cpic_2025.json')` |
| Phenotype thresholds | Research overrides | `classify_phenotype(..., override_thresholds={'UM': 2.5})` |
| Activity scoring | Non-standard gene rules | Patch `alleles.activity.compute_activity_score()` |

## Testing strategy

All public functions are covered by unit tests in `tests/pharmacogenomics/`:
- Allele matching edge cases (empty set, conflicting alleles, copy-number markers)
- Diplotypes for all CYP2D6 stars (*1–*137)
- Phenotype classification sanity across 28 genes
- CPIC guideline lookup integrity
- Report generation has deterministic UUID and timestamp tests

No network calls during tests; PharmGKB fixture is a 200-row TSV sample.

## Performance profile

Single-sample end-to-end latency is ~7 ms on modern x86. Memory footprint is ~4–6 MB
for static tables; per-sample transient dicts add ~2 KB. Parallel workers share base
tables via fork (Linux) so memory scales with number of workers, not dataset size.

See [PERFORMANCE.md](PERFORMANCE.md) for benchmarks, scaling curves, and tuning tips.

## Known limitations

- CYP2D6 CNV detection requires duplication-marker rsIDs (e.g., rs28371685); if
  absent the gene is treated as diploid and copy number defaults to 2.
- Very rare novel variants outside dbSNP or not encoded in a star allele are silently
  ignored — they contribute neither match score nor penalty.
- PharmGKB requires internet on first run unless `PG_PHARMGKB_TSV_PATH` points to a local file.

## Future roadmap

- v0.9: Cython rewrite of subset-matching inner loop (5× speed target)
- v1.0: HL7 FHIR Pharmacogenomics Implementation Guide output
- v1.1: Graph-based diplotype inference allowing phasing uncertainty

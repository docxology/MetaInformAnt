# Development Guide: Pharmacogenomics

## Module boundaries (5 layers)

| Layer | Package | Responsibility |
|-------|---------|----------------|
| Public API | `__init__.py` | `predict_metabolizer`, `call_star_alleles`, … |
| Core logic | `alleles/` | Matching, diplotype, phenotype, activity tables |
| Annotations | `annotations/` | CPIC & PharmGKB loaders + caching |
| Clinical | `clinical/` | Report generation, ACMG classification |
| Interaction | `interaction/` | DDI analysis, polypharmacy risk |

Never import `annotations` from inside `alleles` to avoid circular dependency.

## Adding a new gene (step-by-step)

1. Create JSON: `src/metainformant/pharmacogenomics/alleles/data/<GENE>.json`
   Schema: `{ "gene": "XYZ", "star_alleles": [{ "name":"*1", "variants":["rs123"], "activity":1.0 }] }`
2. Add phenotype thresholds to `alleles/phenotype.py` → `_PHENOTYPE_THRESHOLDS[gene] = { 'PM':0, 'IM':(0,1.0), 'NM':[1.0,2.25), … }`
3. If gene has special CPIC phase rules (like CYP2D6 duplication), extend
   `alleles/diplotype.py` → `_SPECIAL_DIPLOTYPE_RULES[gene]`.
4. Register gene in `alleles/__init__.py` `_BUILTIN_ALLELE_DEFINITIONS[gene] = load_builtin(gene)`.

5a. Unit tests: `tests/pharmacogenomics/test_<gene>.py` — cover:
   - Matching known variant combos → expected star allele names
   - Diplotype canonical ordering (e.g., `*1/*4` vs `*4/*1`)
   - Phenotype classification at threshold boundaries
   run `pytest tests/pharmacogenomics/`

5b. If gene has CNV rules, add test for duplication markers.

## Code style & benchmarks

- Public functions: type hints + Google/NumPy docstring (Params/Returns/Examples).
- Immutable value objects: `@dataclass(frozen=True)`.
- Closed sets: `enum.Enum`.
- **Performance guardrail**: single-sample `predict_metabolizer()` must stay <10 ms on
  i7-12700K. Profile hot path with `cProfile`; any >15 ms function needs optimisation
  or Cython rewrite.
- Cache expensive pure functions with `@disk_cache(ttl=86400)`.

## Adding an annotation source (example: ClinVar PGx)

1. Create `annotations/clinvar_pgx.py`.
2. `load_clinvar_pgx(path=None)` — download to cache if `path=None`, store under
   `~/.cache/metainformant/clinvar/`.
3. `query_clinvar(variant_rsid, gene=None) -> list[dict]`.
4. Register in `annotations/__init__.py`.
5. Add example to `docs/pharmacogenomics/EXAMPLES.md`.

All network calls must be cache-first and never raise on network failure; return
empty list instead.

## Testing philosophy

- Deterministic unit tests only; no network (use fixtures).
- Property tests: for every gene, confirm
  `activity_score = classify_phenotype(determine_diplotype(a,b,gene),gene)` lies in the
  interval corresponding to the returned phenotype.
- Round-trip report export: `text → JSON → object → text` identical.

## Release checklist

- [ ] Allele tables validated against latest PharmVar (external script run)
- [ ] CPIC snapshot version bumped if newer
- [ ] CHANGELOG updated with gene additions / breaking API changes
- [ ] All `docs/pharmacogenomics/EXAMPLES.md` snippets run without error
- [ ] Performance regression suite (`pytest -k perf`) passes
- [ ] Config schema validation (`config.validate_schema('pharmacogenomics')`) clean

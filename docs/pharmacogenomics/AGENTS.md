# Development Guide: Pharmacogenomics

## Module boundaries

| Layer | Package | Responsibility |
|-------|---------|----------------|
| Public API | `__init__.py` | the submodules plus `predict_metabolizer`; every other function is imported from the submodule that defines it |
| Core logic | `alleles/` | star-allele matching, diplotype scoring, phenotype assignment; owns the allele values and the phenotype definitions |
| Front end | `metabolism/` | string and dict front end over `alleles`, plus the dose-guideline table; holds no allele or phenotype data of its own |
| Annotations | `annotations/` | bundled CPIC, PharmGKB and drug-label data with file loaders (no network code) |
| Clinical | `clinical/` | report generation, drug-gene analysis, ACMG classification |
| Interaction | `interaction/` | drug-drug interaction analysis, polypharmacy risk |
| Visualization | `visualization/` | matplotlib plots |

Never import `annotations` from inside `alleles` to avoid a circular dependency.

## Adding a new gene (step-by-step)

1. **Allele values.** Add `"GENE": {"*1": 1.0, ...}` to `_ACTIVITY_SCORE_TABLES` in
   `alleles/diplotype.py`: normal function 1.0, decreased 0.5, no function 0.0, or the
   activity values CPIC publishes for the gene (CYP2D6 uses 0.25 for `*9`, `*10`, `*41`).
   This one table feeds `determine_diplotype`, `classify_phenotype` and the `metabolism`
   front end (`default_allele_function_table()` is a copy of it), so there is no second
   table to update.
2. **Calling definitions (optional).** Add `_BUILTIN_ALLELE_DEFINITIONS["GENE"]` in
   `alleles/star_allele.py`: for each allele `defining_variants` (a frozenset of
   identifiers), `function` and an `activity_value` identical to step 1 (a test enforces
   it). If the gene's alleles cannot be expressed as the presence of variant identifiers,
   leave it out and register the gene with its reason in
   `_GENES_WITHOUT_BUILTIN_DEFINITIONS`, as UGT1A1 and CYP1A2 are; a test keeps that
   registry in step with the two tables.
3. **Phenotype definition.** Add `_PHENOTYPE_THRESHOLDS["GENE"]` in `alleles/phenotype.py`
   as consecutive `(upper_bound, upper_inclusive, phenotype)` bands, copied from the CPIC
   publication that is cited in a comment: the first band starts at 0, each next band starts
   where the previous one ends, and the last ends at `math.inf`. Do not round the bounds
   (CYP2D6 uses 1.25 and 2.25). If CPIC defines the gene by pairs of allele functions, also
   add `_FUNCTION_PHENOTYPES["GENE"]`, keyed by the `(less active, more active)` function
   classes. A gene with no supported phenotype entry raises `ValueError`; no other gene's bands
   are substituted.
4. **Special rules.** Copy-number handling exists for CYP2D6 only (`handle_cyp2d6_cnv`); there
   is no registry of per-gene diplotype rules to extend.

5a. Unit tests in `tests/pharmacogenomics/`: add the gene's CPIC diplotype and phenotype
    pairs to `CPIC_DIPLOTYPES` in `test_cpic_phenotype_definitions.py`, written from the
    publication rather than computed from the package's tables (that file also checks that
    the bands tile every score without gaps, that both phenotype front ends agree, and that
    calling definitions carry the table's values). Cover:
    - Matching known variant combinations to the expected star allele names
    - Diplotype canonical ordering (`*1/*4` and `*4/*1` are the same diplotype)
    - Phenotype classification at band boundaries

    Run `uv run python -m pytest -q tests/pharmacogenomics`.

5b. If the gene has copy-number rules, add a test for the duplication markers.

## Code style

- Public functions: type hints and a Google-style docstring (Args, Returns, Raises).
- Closed sets: `enum.Enum`.
- Results are deterministic functions of the inputs (the report's `generated_at` timestamp
  is the one exception): no random numbers, and no iteration over a `set` of strings whose
  order could change the result. `test_package_determinism.py` scans the package for
  random-number use and compares complete analyses across `PYTHONHASHSEED` values.

## Adding an annotation source

1. Create `annotations/<source>.py` with the data as Python literals and a loader that
   reads a file you pass in (the module has no network code; do not add any).
2. Export the query functions from the new module and register it in
   `annotations/__init__.py`.
3. Add an example to `docs/pharmacogenomics/EXAMPLES.md` and an entry per function to
   `docs/pharmacogenomics/CAPABILITIES.md`.

## Testing philosophy

- Deterministic unit tests only; no network, no mocks.
- Expected phenotypes are written from the CPIC publication named next to them, so a band
  or allele value that drifts from CPIC fails a test. For CYP2C19, SLCO1B1 and UGT1A1 the
  allele-function rules are the reference, and a diplotype's activity score alone does not
  determine the phenotype (CYP2C19 `*9/*17` and `*1/*1` both score 2.0 and differ).
- `predict_metabolizer` and `classify_phenotype` must return the same phenotype and
  activity score for every diplotype of every gene in the allele table.
- Reports: assert the report dictionary (genotype results, recommendations, actions) and
  that recommendations reach the exported text.

## Release checklist

- [ ] `tests/pharmacogenomics/` passes, including the hash-seed determinism check
- [ ] Allele values and phenotype bands of any gene added or changed are traced to the CPIC
      publication cited in the source comments
- [ ] Built-in CPIC guidelines updated if CPIC has published newer ones
- [ ] `docs/pharmacogenomics/CAPABILITIES.md` entries match the docstrings of changed functions
- [ ] CHANGELOG updated with gene additions and behaviour changes

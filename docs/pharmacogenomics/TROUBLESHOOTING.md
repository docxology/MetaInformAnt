# Troubleshooting: Pharmacogenomics

## Symptom → Fix table

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `ValueError: Gene 'XYZ' not found in built-in definitions. Available genes: ...` | Gene not in the built-in allele tables | Load a table with `load_allele_definitions('XYZ', source='file', filepath='xyz.json')` and pass it as `allele_definitions=` |
| Custom allele is unscored or Indeterminate | Its name is absent from the scoring table; `determine_diplotype` does not read `StarAllele.activity_value` directly | Supply explicit values in the definition file (missing values raise), then pass `scoring_table={a.name: a.activity_value for a in definitions}` to `determine_diplotype` |
| `FileNotFoundError: CPIC guideline file not found: ...` | The path given to `load_cpic_guidelines(filepath=...)` does not exist; the built-in guidelines live in the source and need no data files | Correct the path |
| `IndexError: list index out of range` | Code indexed `alleles[1]` but `call_star_alleles` matched fewer than two alleles (one matched allele, or just `*1`) | Pad with the reference allele: `names = [a.name for a in alleles] + ['*1', '*1']`, then use `names[0]` and `names[1]` (assumes heterozygous; pass the same allele twice for homozygous) |
| `ValueError: Invalid diplotype string format: '*4'. Expected '*X/*Y'.` | `classify_phenotype` received a single allele | Pass `'*1/*4'` or the `Diplotype` returned by `determine_diplotype` |
| `TypeError: get_dosing_recommendation() got an unexpected keyword argument 'gene'` | The function takes `(drug, phenotype, guidelines=None)`; there is no gene filter | Select the gene's guideline with `lookup_drug_gene(drug, gene)` and pass `guidelines=[guideline]` |
| `get_dosing_recommendation` returns `None` | No built-in guideline for that drug, or none for that phenotype string (9 drug–gene guidelines are built in) | Check `lookup_drug_gene(drug, gene)`; load your own guidelines with `load_cpic_guidelines(filepath=...)` |
| `ValueError: Gene 'UGT1A1' has built-in activity scores but no built-in star-allele definitions: ...` | UGT1A1 and CYP1A2 have activity values but no calling definitions (UGT1A1 alleles are defined by promoter TA-repeat length) | Build the diplotype from allele names called elsewhere with `determine_diplotype('*1', '*28', 'UGT1A1')`, or supply a definition file with `load_allele_definitions('UGT1A1', filepath=...)` |
| `ValueError: source='file' requires a filepath ...` or `ValueError: filepath is only used with source='file' ...` | `load_allele_definitions` got contradicting arguments, which it refuses so that a file is never ignored | Pass `filepath` alone (or `source='file'` with it) to load a file; omit `filepath` for the built-in table |
| `ValueError: Activity score must be a non-negative number, got nan` | A negative or NaN score reached `alleles.phenotype.predict_metabolizer_status` or `classify_metabolizer`, which refuse to map it to a phenotype | Compute the score from non-negative allele values; check for a missing value upstream |
| Unknown allele produces Indeterminate and score `None` | The allele is absent from the applicable activity table; both phenotype front ends report its name in `unknown_alleles` | Check the name or supply sourced custom activity values; the strict `compute_activity_score` helper still raises for unknown alleles |
| `ValueError: Unrecognised metabolizer phenotype: 'x'` | `MetabolizerPhenotype.from_string` got a label it does not know | Use a CPIC name (`'Poor Metabolizer'`), an abbreviation (`'PM'`) or a token (`'poor'`) |
| CYP2D6 `*1/*4` is intermediate, not normal | Activity score 1.0 falls in CPIC's intermediate band (0 < AS < 1.25) | Nothing to fix: `classify_phenotype` and `predict_metabolizer` agree |
| A report shows `Indeterminate` (`IND`) | Unknown allele, malformed diplotype or unsupported phenotype gene, with no usable supplied phenotype | Inspect `unknown_alleles` and `clinical_significance`; correct the genotype or provide independently supported phenotype information |
| `ImportError: matplotlib is required for pharmacogenomics visualization` | matplotlib is missing from the environment (it is a core dependency) | `uv pip install matplotlib` |

## Diagnostic steps

1. **Check installation**
   ```bash
   uv pip show metainformant
   python -c "import metainformant; print(metainformant.__version__)"
   ```

2. **Enable debug logs**
   ```python
   import logging
   logging.getLogger('metainformant.pharmacogenomics').setLevel(logging.DEBUG)
   ```

3. **Inspect the built-in allele tables**
   ```python
   from metainformant.pharmacogenomics.alleles.star_allele import load_allele_definitions
   try:
       load_allele_definitions('NOT_A_GENE')
   except ValueError as exc:
       print(exc)   # lists every gene with a built-in table
   print(len(load_allele_definitions('CYP2D6')))   # alleles in one table, counting *1
   ```

4. **CPIC validation**
   ```python
   from metainformant.pharmacogenomics.annotations.cpic import load_cpic_guidelines
   gl = load_cpic_guidelines()
   print(f"Loaded {len(gl)} guideline entries")   # 9 built-in
   ```

5. **Activity table sanity**
   ```python-snippet
   from metainformant.pharmacogenomics.alleles.star_allele import load_allele_definitions
   table = {a.name: a.activity_value for a in load_allele_definitions('CYP2D6')}
   print('*1' in table, table['*1'])
   ```

## Common pitfalls

- **Empty or unmatched variant set** returns `[*1]` without an error, which downstream reads as `*1/*1` (normal function) — handle missing genotype data separately (`if not variants: return None`) instead of treating it as "no variants"
- **Identifier format** — variants match by exact string equality with the identifiers in the allele table (rsIDs such as `rs3892097`); a numeric-only `3892097` does not match and is ignored without an error
- **Unsupported phenotype genes** — phenotype APIs raise `ValueError` instead of borrowing CYP2D6 bands. A custom scoring table can score a new gene but does not add a phenotype definition. CYP1A2 can only be scored.
- **Parallel pickle errors** — apply to `parallel.process_map`: the mapped function must be top-level (not a lambda or nested function); `parallel.thread_map` does not pickle
- **Memory spikes** — don't accumulate results in a list inside the mapped function; let `parallel.thread_map` collect

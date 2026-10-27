# Troubleshooting: Pharmacogenomics

## Symptom → Fix table

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `ValueError: Gene 'XYZ' not found in built-in definitions. Available genes: ...` | Gene not in the built-in allele tables | Load a table with `load_allele_definitions('XYZ', source='file', filepath='xyz.json')` and pass it as `allele_definitions=` |
| Custom allele behaves as normal function (activity 1.0) | The JSON/TSV entry has no `activity_value`, which defaults to 1.0 (an empty TSV cell raises `ValueError: could not convert string to float`); or the allele went through `determine_diplotype` without `scoring_table=`, which ignores `StarAllele.activity_value` and scores alleles missing from its built-in table as 1.0 | Add `"activity_value": <float>` to every allele entry, and pass `scoring_table={a.name: a.activity_value for a in definitions}` to `determine_diplotype` |
| `FileNotFoundError: CPIC guideline file not found: ...` | The path given to `load_cpic_guidelines(filepath=...)` does not exist; the built-in guidelines live in the source and need no data files | Correct the path |
| `IndexError: list index out of range` | Code indexed `alleles[1]` but `call_star_alleles` matched fewer than two alleles (one matched allele, or just `*1`) | Pad with the reference allele: `names = [a.name for a in alleles] + ['*1', '*1']`, then use `names[0]` and `names[1]` (assumes heterozygous; pass the same allele twice for homozygous) |
| `ValueError: Invalid diplotype string format: '*4'. Expected '*X/*Y'.` | `classify_phenotype` received a single allele | Pass `'*1/*4'` or the `Diplotype` returned by `determine_diplotype` |
| `TypeError: get_dosing_recommendation() got an unexpected keyword argument 'gene'` | The function takes `(drug, phenotype, guidelines=None)`; there is no gene filter | Select the gene's guideline with `lookup_drug_gene(drug, gene)` and pass `guidelines=[guideline]` |
| `get_dosing_recommendation` returns `None` | No built-in guideline for that drug, or none for that phenotype string (9 drug–gene guidelines are built in) | Check `lookup_drug_gene(drug, gene)`; load your own guidelines with `load_cpic_guidelines(filepath=...)` |
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
- **Genes without a scoring table** — `determine_diplotype` does not raise for an unknown gene; every allele scores 1.0 (a warning is logged) unless you pass `scoring_table=`
- **Parallel pickle errors** — apply to `parallel.process_map`: the mapped function must be top-level (not a lambda or nested function); `parallel.thread_map` does not pickle
- **Memory spikes** — don't accumulate results in a list inside the mapped function; let `parallel.thread_map` collect

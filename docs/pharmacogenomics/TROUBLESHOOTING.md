# Troubleshooting: Pharmacogenomics

## Symptom → Fix table

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `ValueError: Unknown gene 'XYZ'` | Gene not in built-in allele table | Add custom JSON via `load_allele_definitions(gene='XYZ', filepath='xyz.json')` |
| `KeyError: activity` | Custom allele missing `activity` float | Add `"activity": <float>` to every star entry |
| `FileNotFoundError` for CPIC JSON | Editable install missing data files | `uv pip install -e .` or load via explicit `path=` |
| `ValueError: Cannot resolve diplotype` | Only 1 allele detected (empty variant set) | Pad with `*1`: `determine_diplotype(allele_a, '*1', gene=...)` |
| `HTTPError 429` from PharmGKB | Rate limit | Set `PG_PHARMGKB_API_TOKEN` or rely on cached TSV |
| `ImportError: scvi-tools` | Missing `[deconvolution]` extra | `uv pip install metainformant[deconvolution]` |
| `MPLBACKEND` error on headless | No display server | `export MPLBACKEND=Agg` |

## Diagnostic steps

1. **Check installation**
   ```bash
   uv pip show metainformant
   python -c "import metainformant; print(metainformant.__version__)"
   ```

2. **Enable debug logs**
   ```python
   import os; os.environ['PG_LOGGING'] = 'DEBUG'
   ```

3. **Inspect loaded allele count**
   ```python
   from metainformant.pharmacogenomics.alleles.star_allele import _BUILTIN_ALLELE_DEFINITIONS
   print({k: len(v) for k,v in _BUILTIN_ALLELE_DEFINITIONS.items()})
   ```

4. **CPIC validation**
   ```python
   from metainformant.pharmacogenomics.annotations.cpic import load_cpic_guidelines
   gl = load_cpic_guidelines()
   print(f"Loaded {len(gl)} guideline entries")
   ```

5. **Activity table sanity**
   ```python
   from metainformant.pharmacogenomics.alleles.activity import _ACTIVITY_SCORE_TABLES
   table = _ACTIVITY_SCORE_TABLES['CYP2D6']
   print('*1' in table, table['*1'])
   ```

## Common pitfalls

- **Empty variant set** raises `ValueError` — guard with `if not variants: return fallback`
- **Wrong gene case** — use uppercase: `'CYP2D6'` not `'cyp2d6'`
- **rsID format** — must include `rs` prefix; numeric-only strings are ignored
- **Parallel pickle errors** — mapped function must be top-level (not lambda inside another function)
- **Memory spikes** — don't accumulate results in a list inside the mapped function; let `parallel.map` collect

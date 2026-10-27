# Star Allele Matching Engine

Star-allele names (`*1`, `*4`, `*2A`) are the standard labels for pharmacogene haplotypes.
This page describes how `call_star_alleles` matches observed variants to them, how CYP2D6
copy number is handled, and where the allele values and definitions live.

## Matching algorithm (detailed walk-through)

Given an observed variant set V and a gene with allele definitions A₁…Aₙ, each with a set
of defining variants Dᵢ (rsIDs or `chr:pos:ref:alt` identifiers):

```
Step 1 — Set aside the reference: the definition named *1 is never matched on its own
         variants, and a definition without defining variants is skipped
Step 2 — Sort the remaining definitions by |Dᵢ|, descending (ties keep definition order)
Step 3 — Walk the sorted list; keep Dᵢ when Dᵢ ⊆ V minus the variants consumed so far,
         then mark Dᵢ as consumed
Step 4 — If nothing was kept, return [*1] (or [] when the definitions contain no *1)
Step 5 — Return the kept alleles sorted by name
```

Complexity is O(n × m) for n definitions of m variants on average. Trying the most specific
allele first and consuming its variants lets overlapping definitions resolve to the more
specific one: observing both `rs4149056` and `rs2306283` in SLCO1B1 gives `*15` (defined by
both), not `*5` (`rs4149056` alone) plus `*1B` (`rs2306283` alone). An empty or unrecognised
variant set returns `[*1]`, which downstream code reads as the reference allele rather than an
error.

## Duplication and CNV handling (CYP2D6 specific)

`handle_cyp2d6_cnv(variants, copy_number)` takes the copy number from the caller; the marker
identifiers `CYP2D6_DEL` (the defining variant of `*5`) and `CYP2D6_DUP` are not used to infer
it and are removed before alleles are called.

| Copy number | Result |
|-------------|--------|
| 0 | `*5/*5` (deletion), activity 0.0 |
| 1 | the called allele with `*5` (deletion), e.g. `*4/*5` |
| 2 | ordinary calling, a `*1` padded where fewer than two alleles are called |
| 3 or more | the extra copies go to the more active allele, whose activity is multiplied by `copy_number - 1` and whose name gets the suffix `x<copy_number - 1>`: copy number 3 with a `*1`/`*4` genotype gives `*1x2/*4` with activity 2.0 |

`call_star_alleles` and `determine_diplotype` themselves never infer copy number: without an
explicit `copy_number` the gene is treated as diploid. `determine_diplotype` scores the
symbolic placeholders `*1xN` and `*2xN` as 2.0, which is exact for two copies only.

## Activity score reference

Per-gene activity values come from `_ACTIVITY_SCORE_TABLES` in
`src/metainformant/pharmacogenomics/alleles/diplotype.py`, the one allele-function table of
the package (`metabolism.metabolizer_status.default_allele_function_table()` returns a copy).
Normal function is 1.0, decreased 0.5, no function 0.0; CYP2D6 uses CPIC's own activity
values:

| Gene | No function (0.0) | Decreased | Normal (1.0) | Increased |
|------|-------------------|-----------|--------------|-----------|
| CYP2D6 | `*3`–`*8`, `*11`, `*12`, `*15`, `*36` | 0.5: `*14`, `*17`, `*29`; 0.25: `*9`, `*10`, `*41` | `*1`, `*2`, `*33`, `*35` | `*1xN`, `*2xN` (2.0 placeholders) |
| CYP2C19 | `*2`–`*8` | 0.5: `*9`, `*10` | `*1` | 1.5: `*17` |
| CYP2C9 | `*3`, `*6` | 0.5: `*2`, `*4`, `*5`, `*8`, `*11` | `*1` | |
| CYP3A5 | `*3`, `*6`, `*7` | | `*1` | |
| DPYD | `*2A`, `*13` | 0.5: `c.2846A>T`, `HapB3` | `*1` | |
| TPMT | `*2`, `*3A`, `*3B`, `*3C`, `*4` | | `*1` | |
| NUDT15 | `*2`, `*3` | | `*1` | |
| SLCO1B1 | `*5`, `*15` | | `*1`, `*1B` | |
| UGT1A1 | | 0.5: `*6`, `*28`, `*37` | `*1`, `*36` | |
| CYP1A2 | | 0.5: `*1C`, `*1K` | `*1A` | 1.5: `*1F` |

UGT1A1 and CYP1A2 have values only: no calling definitions are bundled for them (see
[CONFIGURATION.md](CONFIGURATION.md#genes-with-activity-values-but-no-calling-definitions)).
The calling definitions (`_BUILTIN_ALLELE_DEFINITIONS`, 8 genes and 45 alleles) carry the same
values for the alleles they define; a test enforces it.

## Custom gene or allele addition

A definition file maps each allele name to its defining variants, function and activity value
(a JSON object, optionally nested under the gene symbol, or a TSV with the columns `allele`,
`defining_variants`, `function` and `activity_value`; see [SPEC.md](SPEC.md#custom-allele-definition-files)):

```python-snippet
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles, load_allele_definitions

definitions = load_allele_definitions('CYP2E1', filepath='output/pharmacogenomics/cyp2e1.json')
alleles = call_star_alleles(variants, 'CYP2E1', allele_definitions=definitions)
```

`filepath` alone loads the file; the definitions are used only when you pass them as
`allele_definitions`. Example 7 in [EXAMPLES.md](EXAMPLES.md) writes and loads such a file end
to end.

## Debugging mismatches

If an expected allele is not returned, compare your variants with each definition:

```python
from metainformant.pharmacogenomics.alleles.star_allele import load_allele_definitions, match_allele_definition

for allele in load_allele_definitions('CYP2D6'):
    result = match_allele_definition(your_variants, allele)
    if result['match_score'] > 0:
        print(allele.name, result['is_match'], sorted(result['missing_variants']))
```

`detect_novel_alleles(variants, known_alleles)` lists observed identifiers that no known allele
explains. Common pitfall: identifiers match by exact string equality with those in the table
(rsIDs such as `rs3892097`), so a numeric-only `3892097`, or an identifier in another format,
matches nothing and is ignored without an error.

## Performance internals

Matching is plain Python on `frozenset` operations. The allele tables are Python literals, so
there is no file to parse on first use, and results do not depend on set iteration order.

## Known edge cases

- **Allele with no defining variants** — never matched on its own. `*1` is the reference and is
  returned only when nothing else matched. The CYP2D6 deletion `*5` is defined by the
  identifier `CYP2D6_DEL`, not by an empty set.
- **Overlapping definitions** — resolved by specificity and variant consumption, as above;
  a variant belongs to at most one called allele.
- **Equal numbers of defining variants** — tried in definition order, and the result is sorted
  by allele name, so the output order is deterministic.
- **More than two alleles** — `call_star_alleles` returns every allele it matched, which can be
  more than the two a diplotype needs; choose the pair yourself (or phase the variants with
  `phased_diplotype`).

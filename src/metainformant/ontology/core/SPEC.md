# Specification: core

## 🎯 Scope
Core ontology functionality: GO analysis, OBO parsing, and type definitions.

## 🧱 Architecture
- **Dependency Level**: Core
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 6 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `go.py`
- `go_api.py`
- `hpo.py`
- `obo.py`
- `types.py`

### Behavior Contracts

#### `go.py`
- `enrich_genes(genes, background, annotations, alpha=0.05, method="fisher")` — `method` is `"fisher"` or
  `"hypergeometric"`; anything else raises `ValueError`. Query genes outside the background are dropped, and each term's
  gene set is intersected with the background, so both methods test the same one-sided upper tail on the same counts
  and give the same p-value. Only terms that share a gene with the query are tested. The returned DataFrame is sorted by
  `p_value` and has the columns `go_term`, `p_value`, `enrichment`, `observed`, `expected`, `term_size` (the term's genes
  in the background), `query_size`, `bonferroni_p` (`p_value` times the number of tested terms, capped at 1) and
  `significant` (`bonferroni_p <= alpha`). Every tested term is returned; the frame is empty when nothing can be tested.
  Without scipy the p-value is the exact upper tail, summed in integer arithmetic.
- `semantic_similarity(term1, term2, term_ic, hierarchy, method="resnik")` — `"resnik"` is the information content of
  the most informative common ancestor in the units of `term_ic` (bits for `calculate_term_ic`), so it is not bounded by
  1. `"lin"` is `2 * IC(MICA) / (IC(term1) + IC(term2))` and `"jiang-conrath"` is `1 - distance / (2 * max(term_ic))`
  clamped to [0, 1]. A term missing from `term_ic`, or two terms without a common ancestor, give 0.0; an unknown method
  raises `ValueError`.

#### `obo.py`
- `parse_obo_synonym(value)` returns `(text, scope, dbxrefs)` for the value of a `synonym:` tag: the text without its
  quotes and with escapes decoded, the scope (`EXACT`, `BROAD`, `NARROW`, `RELATED`, or `None`) and the dbxref list
  (empty when absent). A `!` inside the quotes is text; an unescaped quote closes the text only before the end of the
  value, the scope word, `[`, `{` or `!`; malformed values are read leniently.
- `parse_obo` adds the decoded text of every quoted `synonym:` line to `Term.synonyms`. The OBO writer in
  `query.serialize` emits `synonym: "<escaped text>" EXACT []`, so quotes, backslashes and line breaks round-trip.

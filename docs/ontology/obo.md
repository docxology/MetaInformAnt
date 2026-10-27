# OBO Format Parsing

Parser for the Open Biological and Biomedical Ontologies (OBO) flat-file
format. Converts OBO files into METAINFORMANT `Ontology` objects containing
terms, relationships, and metadata.

## Key Concepts

**OBO format** is a human-readable text format used by the Gene Ontology,
Human Phenotype Ontology, and many other biomedical ontologies. An OBO file
contains:

- A **header** with `format-version`, `date`, `ontology` metadata.
- **[Term] stanzas** defining ontology terms with `id`, `name`, `namespace`,
  `def`, `synonym`, `xref`, `is_a`, and other relationship tags.
- **[Typedef] stanzas** defining relationship types.
- **[Instance] stanzas** (rare) defining ontology instances.

**Parsed Term fields:** Each `[Term]` stanza produces a `Term` object with
`id`, `name`, `definition`, `namespace`, `synonyms` (list), `xrefs` (list),
`alt_ids` (list), `is_obsolete` flag, and a `metadata["relationships"]` list of
`(relation_type, target_id)` tuples. When the definition cites dbxrefs, they are
kept apart from the text in `metadata["definition_xrefs"]`.

**Definitions:** The value of a `def:` tag is a quoted string, then a bracketed
dbxref list, then optional `{...}` trailing modifiers and a `! comment`, for
example `def: "A process that \"x\" does." [GOC:mah, PMID:123]`. The parser
stores the text in `Term.definition` and the dbxrefs in
`metadata["definition_xrefs"]`, a list of strings that is absent when the
definition has no dbxrefs.

- The text is returned without its surrounding quotes and with backslash
  escapes decoded: `\n`, `\t` and `\W` become a newline, a tab and a space, and
  a backslash before any of `:,"\()[]{}` stands for that character. Other
  backslash sequences are kept as written.
- A `!` inside the quotes is part of the text. Trailing modifiers and the
  comment after the dbxref list are dropped.
- An unescaped quote inside the text is kept as text unless it is followed by
  the end of the value, `[`, `{` or `!`.
- Malformed values are read leniently: an unquoted value is kept as written
  (minus a trailing `! comment`), and a quoted string that is never closed runs
  to the end of the value.
- `[Typedef]` definitions go through the same parser. They are stored in the
  ontology metadata under `typedef:<id>` as `definition` and, when present,
  `definition_xrefs`.

**Writing definitions:** `save_ontology(onto, path, format="obo")` (in
`metainformant.ontology.query.serialize`) writes each definition as
`def: "<escaped text>" [<escaped dbxrefs>]`. The dbxrefs come from
`metadata["definition_xrefs"]`, and the list is written as `[]` when there are
none. Backslashes and double quotes in the text are escaped, and line breaks and
tabs are written as `\n` and `\t`, so the value stays on one line; in the
dbxrefs, backslashes, commas, closing brackets, double quotes and line breaks
are escaped. A definition and its dbxrefs therefore survive `save_ontology`
followed by `load_ontology` (or `parse_obo`); CR and CRLF line breaks come back
as LF.

**Relationship extraction:** After all terms are parsed, `is_a`, `part_of`,
`regulates`, `negatively_regulates`, `positively_regulates`, `has_part`,
`occurs_in`, `happens_during`, and `ends_during` tags are converted into
`Relationship` objects linking source and target term IDs.

## Function Reference

### parse_obo

```python-snippet
def parse_obo(path: str | Path) -> Ontology
```

Parses an OBO file and returns an `Ontology` object containing all terms and
relationships. Reads the file as UTF-8 text, parses the header, iterates
stanzas, builds `Term` objects, and extracts relationships.

Raises `metainformant.core.utils.errors.ValidationError` if the file does not
exist or a `[Term]` stanza has an `id` but no `name`, and `UnicodeDecodeError`
if the file is not valid UTF-8. Unrecognised lines and tags are skipped rather
than rejected, and a file without stanzas parses to an empty ontology.

```pycon
>>> ontology = parse_obo("data/go.obo")
>>> len(ontology)
50000
```

### validate_obo_format

```python-snippet
def validate_obo_format(path: str | Path) -> tuple[bool, List[str]]
```

Validates an OBO file. Checks for:
- Presence of `format-version` in the header (first 10 lines).
- At least one `[Term]` stanza.
- Successful parse via `parse_obo`, which builds the whole ontology in memory;
  an exception it raises is reported as `Parsing failed: ...`.

Returns `(is_valid, error_messages)`; raises `ValidationError` if the path does
not exist.

### get_obo_statistics

```python-snippet
def get_obo_statistics(path: str | Path) -> Dict[str, Any]
```

Collects file-level statistics without building a full ontology graph.
Returns:
- `total_lines` -- line count of the file.
- `term_count`, `typedef_count`, `instance_count` -- stanza counts.
- `relationship_count` -- `is_a`, `part_of`, `regulates`, `has_part`,
  `occurs_in` tag count.
- `obsolete_count` -- terms marked `is_obsolete: true`.

### parse_obo_definition / quote_obo_string / format_obo_xref_list

```python-snippet
def parse_obo_definition(value: str) -> Tuple[str, List[str]]
def quote_obo_string(text: str) -> str
def format_obo_xref_list(xrefs: Iterable[str]) -> str
```

`parse_obo_definition` splits everything after `def:` into
`(definition_text, dbxrefs)` using the rules under Definitions above.
`quote_obo_string` renders text as a quoted OBO string, surrounding quotes
included, and `format_obo_xref_list` renders dbxrefs as `[a, b]`. They produce
the two halves of the `def:` line the OBO writer emits, and
`parse_obo_definition` reads their output back to the original text and items.

```pycon
>>> parse_obo_definition('"A process that \\"x\\" does." [GOC:mah, PMID:123]')
('A process that "x" does.', ['GOC:mah', 'PMID:123'])
>>> quote_obo_string('A "quoted" word')
'"A \\"quoted\\" word"'
>>> format_obo_xref_list(["GOC:mah", "PMID:1,2"])
'[GOC:mah, PMID:1\\,2]'
```

## Data Types

The parser produces objects from `metainformant.ontology.core.types`:

| Type | Fields |
|------|--------|
| `Term` | id, name, definition, namespace, synonyms, xrefs, is_obsolete, metadata |
| `Relationship` | source, target, relation_type |
| `Ontology` | terms (dict), relationships (list), metadata (dict) |

Factory functions `create_term`, `create_relationship`, `create_ontology` are
used to construct these objects.

## Usage Example

```python
from metainformant.ontology.core.obo import parse_obo, validate_obo_format, get_obo_statistics

# Validate before loading
is_valid, errors = validate_obo_format("data/go-basic.obo")
if not is_valid:
    print("Validation errors:", errors)

# Quick statistics
stats = get_obo_statistics("data/go-basic.obo")
print(f"Terms: {stats['term_count']}, Relationships: {stats['relationship_count']}")

# Full parse
ontology = parse_obo("data/go-basic.obo")
term = ontology.terms["GO:0008150"]
print(f"{term.id}: {term.name} ({term.namespace})")
print(term.definition, term.metadata.get("definition_xrefs", []))
```

## Related Modules

- `metainformant.ontology.query.query` -- traversal, ancestors, descendants
- `metainformant.ontology.core.types` -- Term, Relationship, Ontology data types
- `metainformant.ontology.core.go` -- Gene Ontology specific utilities
- `metainformant.ontology.query.serialize` -- ontology serialisation

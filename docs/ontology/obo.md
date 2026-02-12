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
`is_obsolete` flag, and a `metadata["relationships"]` list of
`(relation_type, target_id)` tuples.

**Relationship extraction:** After all terms are parsed, `is_a`, `part_of`,
`regulates`, `negatively_regulates`, `positively_regulates`, `has_part`,
`occurs_in`, `happens_during`, and `ends_during` tags are converted into
`Relationship` objects linking source and target term IDs.

## Function Reference

### parse_obo

```python
def parse_obo(path: str | Path) -> Ontology
```

Parses an OBO file and returns an `Ontology` object containing all terms and
relationships. Reads the file via `metainformant.core.io.read_text`, parses
the header, iterates stanzas, builds `Term` objects, and extracts
relationships.

Raises `FileNotFoundError` if the file does not exist; `ValueError` on
format errors.

```python
>>> ontology = parse_obo("data/go.obo")
>>> len(ontology)
50000
```

### validate_obo_format

```python
def validate_obo_format(path: str | Path) -> tuple[bool, List[str]]
```

Validates an OBO file without fully loading it into memory. Checks for:
- Presence of `format-version` in the header (first 10 lines).
- At least one `[Term]` stanza.
- Successful parse via `parse_obo`.

Returns `(is_valid, error_messages)`.

### get_obo_statistics

```python
def get_obo_statistics(path: str | Path) -> Dict[str, Any]
```

Collects file-level statistics without building a full ontology graph.
Returns:
- `total_lines` -- line count of the file.
- `term_count`, `typedef_count`, `instance_count` -- stanza counts.
- `relationship_count` -- `is_a`, `part_of`, `regulates`, `has_part`,
  `occurs_in` tag count.
- `obsolete_count` -- terms marked `is_obsolete: true`.

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
```

## Related Modules

- `metainformant.ontology.query.query` -- traversal, ancestors, descendants
- `metainformant.ontology.core.types` -- Term, Relationship, Ontology data types
- `metainformant.ontology.core.go` -- Gene Ontology specific utilities
- `metainformant.ontology.query.serialize` -- ontology serialisation

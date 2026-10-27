"""Regression tests for ``load_ontology(format="json")`` on ontologies parsed from OBO.

A term parsed from OBO keeps its raw ``(relation, target)`` pairs in ``Term.metadata["relationships"]``. The JSON loader
splatted ``Term.metadata`` into ``create_term(...)``, which already has a ``relationships`` argument, so every such file
failed with ``TypeError: create_term() got multiple values for keyword argument 'relationships'``.

While rebuilding the loader a second round-trip defect surfaced: ``add_term`` records the ``is_a`` link of each
``is_a_parents`` entry and the saved relationship list repeats it, so every ``is_a`` link of a term-built ontology came
back twice. The oracle is the ontology that was saved: the loaded one must equal it. No test doubles.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import pytest

from metainformant.ontology.core.obo import parse_obo
from metainformant.ontology.core.types import Ontology, Relationship, Term
from metainformant.ontology.query.serialize import load_ontology, save_ontology

OBO_TEXT = """format-version: 1.2
data-version: releases/2026-01-01
ontology: demo

[Term]
id: GO:0000001
name: root process
namespace: biological_process
def: "A root process." [GOC:mah]

[Term]
id: GO:0000002
name: child process
namespace: biological_process
def: "A child \\"quoted\\" process." [GOC:mah, PMID:123]
synonym: "kid process" EXACT []
xref: Reactome:R-HSA-1
alt_id: GO:0000099
is_a: GO:0000001 ! root process
part_of: GO:0000001 ! root process

[Term]
id: GO:0000003
name: grandchild process
namespace: biological_process
is_obsolete: true
is_a: GO:0000002 ! child process
regulates: GO:0000001 ! root process

[Typedef]
id: part_of
name: part of
"""


def _jsonable(value: Any) -> Any:
    """What a value looks like after a trip through JSON (tuples become lists)."""
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _term_view(term: Term) -> dict[str, Any]:
    return {
        "id": term.term_id,
        "name": term.name,
        "namespace": term.namespace,
        "definition": term.definition,
        "synonyms": term.synonyms,
        "xrefs": term.xrefs,
        "alt_ids": term.alt_ids,
        "is_a_parents": term.is_a_parents,
        "relationships": term.relationships,
        "subsets": term.subsets,
        "is_obsolete": term.is_obsolete,
        "metadata": _jsonable(term.metadata),
    }


def _relationship_counts(onto: Ontology) -> Counter[tuple[str, str, str]]:
    return Counter((rel.source, rel.target, rel.relation_type) for rel in onto.relationships)


def _assert_equal_ontologies(loaded: Ontology, original: Ontology) -> None:
    assert {tid: _term_view(term) for tid, term in loaded.terms.items()} == {
        tid: _term_view(term) for tid, term in original.terms.items()
    }
    assert _relationship_counts(loaded) == _relationship_counts(original)
    assert loaded.parents_of == original.parents_of
    assert loaded.children_of == original.children_of
    assert loaded.get_roots() == original.get_roots()
    assert loaded.get_leaves() == original.get_leaves()
    assert _jsonable(loaded.metadata) == _jsonable(original.metadata)


def test_an_obo_parsed_ontology_survives_a_json_round_trip(tmp_path: Path) -> None:
    obo_path = tmp_path / "demo.obo"
    obo_path.write_text(OBO_TEXT)
    original = parse_obo(obo_path)
    assert original.terms["GO:0000002"].metadata["relationships"], "scenario needs the metadata key that used to clash"
    json_path = tmp_path / "demo.json"

    save_ontology(original, json_path, format="json")
    loaded = load_ontology(json_path, format="json")

    _assert_equal_ontologies(loaded, original)
    assert loaded.terms["GO:0000002"].definition == 'A child "quoted" process.'
    assert ancestors_of(loaded, "GO:0000003") == {"GO:0000002", "GO:0000001"}


def ancestors_of(onto: Ontology, term_id: str) -> set[str]:
    return set(onto.get_ancestors(term_id))


def test_the_bundled_mini_go_file_round_trips(tmp_path: Path) -> None:
    mini = Path(__file__).resolve().parent.parent / "data" / "ontology" / "go_mini.obo"
    original = parse_obo(mini)
    json_path = tmp_path / "go_mini.json"

    save_ontology(original, json_path, format="json")

    _assert_equal_ontologies(load_ontology(json_path, format="json"), original)


def test_is_a_links_are_not_duplicated_for_an_ontology_built_with_add_term(tmp_path: Path) -> None:
    original = Ontology(metadata={"id": "demo"})
    original.add_term(Term("A", "root"))
    original.add_term(Term("B", "child", is_a_parents=["A"]))
    original.add_term(Term("C", "grandchild", is_a_parents=["B"]))
    original.add_relationship(Relationship("C", "A", "part_of", metadata={"evidence": "IDA"}))
    assert len(original.relationships) == 3
    path = tmp_path / "demo.json"

    save_ontology(original, path, format="json")
    loaded = load_ontology(path, format="json")

    assert _relationship_counts(loaded) == _relationship_counts(original)
    assert len(loaded.relationships) == 3
    part_of = [rel for rel in loaded.relationships if rel.relation_type == "part_of"]
    assert [rel.metadata for rel in part_of] == [{"evidence": "IDA"}]


def test_repeated_relationships_are_kept_as_they_were_saved(tmp_path: Path) -> None:
    original = Ontology()
    original.add_term(Term("A", "root"))
    original.add_term(Term("B", "child", is_a_parents=["A"]))
    original.add_relationship(Relationship("B", "A", "is_a"))  # a second, explicit copy of the same link
    original.add_relationship(Relationship("B", "A", "part_of"))
    original.add_relationship(Relationship("B", "A", "part_of"))
    path = tmp_path / "repeats.json"

    save_ontology(original, path, format="json")
    loaded = load_ontology(path, format="json")

    assert (
        _relationship_counts(loaded)
        == _relationship_counts(original)
        == Counter({("B", "A", "is_a"): 2, ("B", "A", "part_of"): 2})
    )


@pytest.mark.parametrize("key", ["terms", "relationships", "source", "name"])
def test_metadata_keys_that_look_like_arguments_do_not_break_loading(tmp_path: Path, key: str) -> None:
    original = Ontology(metadata={key: "x", "id": "demo"})
    original.add_term(Term("A", "root", metadata={key: ["y"], "relationships": [["is_a", "B"]]}))
    path = tmp_path / "keys.json"

    save_ontology(original, path, format="json")
    loaded = load_ontology(path, format="json")

    assert loaded.metadata == original.metadata
    assert loaded.terms["A"].metadata == original.terms["A"].metadata

from __future__ import annotations

from pathlib import Path

from metainformant.ontology.go import load_go_obo, write_go_summary
from metainformant.ontology.query import ancestors, descendants, subgraph


def _data_path() -> Path:
    return Path("tests/data/ontology/go_mini.obo").resolve()


def test_load_go_and_traverse(tmp_path: Path) -> None:
    onto = load_go_obo(_data_path())
    assert onto.num_terms() >= 4

    # Choose a known chain in the mini data
    child = "GO:0000002"
    parent = "GO:0000001"

    assert parent in ancestors(onto, child)
    assert child in descendants(onto, parent)

    # Subgraph retains connectivity
    sg = subgraph(onto, roots=[parent])
    assert sg.has_term(parent)
    assert sg.has_term(child)

    # Write a summary under output/
    out = write_go_summary(onto)
    assert out.exists()



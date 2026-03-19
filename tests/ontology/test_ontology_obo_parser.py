from __future__ import annotations

from metainformant.ontology.core.obo import parse_obo


def test_parse_obo_minimal_fields() -> None:
    from pathlib import Path
    path = Path(__file__).parent.parent / "data" / "ontology" / "go_mini.obo"
    onto = parse_obo(path)
    # Minimal assertions on a few IDs in the mini fixture
    for tid in [
        "GO:0000001",
        "GO:0000002",
        "GO:0000003",
        "GO:0000004",
    ]:
        assert onto.has_term(tid)

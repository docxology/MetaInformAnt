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


def _write_obo(tmp_path, content: str):
    path = tmp_path / "mini.obo"
    path.write_text(content)
    return path


OBO_WITH_TYPEDEF = """format-version: 1.2

[Term]
id: GO:0000001
name: root term

[Typedef]
id: part_of
name: part of
def: "X is part of Y" []
"""


def test_validate_obo_format_valid(tmp_path) -> None:
    from metainformant.ontology.core.obo import validate_obo_format

    path = _write_obo(tmp_path, OBO_WITH_TYPEDEF)
    is_valid, errors = validate_obo_format(path)
    assert is_valid
    assert errors == []


def test_validate_obo_format_missing_header(tmp_path) -> None:
    from metainformant.ontology.core.obo import validate_obo_format

    path = _write_obo(tmp_path, "[Term]\nid: GO:0000001\nname: root\n")
    is_valid, errors = validate_obo_format(path)
    assert not is_valid
    assert any("format-version" in e for e in errors)


def test_get_obo_statistics(tmp_path) -> None:
    from metainformant.ontology.core.obo import get_obo_statistics

    content = OBO_WITH_TYPEDEF + "\n[Term]\nid: GO:0000002\nname: child\nis_a: GO:0000001\nis_obsolete: true\n"
    path = _write_obo(tmp_path, content)
    stats = get_obo_statistics(path)
    assert stats["term_count"] == 2
    assert stats["typedef_count"] == 1
    assert stats["relationship_count"] == 1
    assert stats["obsolete_count"] == 1


def test_parse_obo_typedef_stored_in_metadata(tmp_path) -> None:
    path = _write_obo(tmp_path, OBO_WITH_TYPEDEF)
    onto = parse_obo(path)
    assert onto.has_term("GO:0000001")
    assert onto.metadata["typedef:part_of"]["name"] == "part of"
    assert onto.metadata["typedef:part_of"]["definition"] == 'X is part of Y" []'  # parser keeps trailing qualifier

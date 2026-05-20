from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.core.utils.errors import ValidationError
from metainformant.ontology.core.go import load_go_obo, validate_go_ontology, write_go_summary
from metainformant.ontology.query.query import ancestors, descendants, get_subontology


def _data_path() -> Path:
    return Path(__file__).parent.parent / "data" / "ontology" / "go_mini.obo"


def test_load_go_and_traverse(tmp_path: Path) -> None:
    onto = load_go_obo(_data_path())
    assert len(onto) >= 4

    # Choose a known chain in the mini data
    child = "GO:0000002"
    parent = "GO:0000001"

    assert parent in ancestors(onto, child)
    assert child in descendants(onto, parent)

    # Subgraph retains connectivity
    sg = get_subontology(onto, [parent])
    assert sg.has_term(parent)
    assert sg.has_term(child)

    # Write a summary under output/
    out = write_go_summary(onto)
    assert out.exists()

    # Validate GO ontology
    is_valid, errors = validate_go_ontology(onto)


def test_validate_go_ontology(tmp_path: Path) -> None:
    """Test GO ontology validation."""
    onto = load_go_obo(_data_path())
    is_valid, errors = validate_go_ontology(onto)
    # Test data should be valid
    assert isinstance(is_valid, bool)
    assert isinstance(errors, list)


def test_load_go_obo_nonexistent_file() -> None:
    """Test loading nonexistent OBO file raises error."""
    with pytest.raises(ValidationError, match="does not exist"):
        load_go_obo("nonexistent_file.obo")


def test_write_go_summary_custom_path(tmp_path: Path) -> None:
    """Test writing GO summary to custom path."""
    onto = load_go_obo(_data_path())
    custom_path = tmp_path / "custom_summary.json"
    summary_path = write_go_summary(onto, dest=custom_path)
    assert summary_path == custom_path
    assert summary_path.exists()

    # Verify summary contains expected fields
    import json

    summary_data = json.loads(summary_path.read_text())
    assert "term_count" in summary_data
    assert "namespaces" in summary_data
    assert "num_roots" in summary_data
    assert "num_leaves" in summary_data

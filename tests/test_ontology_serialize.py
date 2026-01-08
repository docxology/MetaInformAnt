"""Tests for ontology serialization module.

All tests follow NO_MOCKING_POLICY and use real implementations.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.ontology.serialize import load_ontology, save_ontology
from metainformant.ontology.types import Ontology, Term


class TestSerialization:
    """Test ontology serialization."""

    def test_save_and_load_ontology(self, tmp_path: Path):
        """Test saving and loading ontology."""
        # Create test ontology
        onto = Ontology()
        term1 = Term(
            term_id="GO:001",
            name="term1",
            namespace="test",
            definition="Test term 1",
            alt_ids=["GO:001a"],
            is_a_parents=["GO:002"],
            relationships={"part_of": ["GO:003"]},
            synonyms=["synonym1"],
            xrefs=["xref:001"],
            subsets=["subset1"],
        )
        term2 = Term(term_id="GO:002", name="term2", namespace="test")
        onto.add_term(term1)
        onto.add_term(term2)

        # Save
        save_path = tmp_path / "test_onto.json"
        saved_path = save_ontology(onto, save_path)
        assert saved_path.exists()
        assert saved_path == save_path

        # Load
        loaded = load_ontology(save_path)
        assert loaded.num_terms() == 2
        assert loaded.has_term("GO:001")
        assert loaded.has_term("GO:002")

        # Verify term attributes preserved
        loaded_term1 = loaded.get_term("GO:001")
        assert loaded_term1 is not None
        assert loaded_term1.name == "term1"
        assert loaded_term1.namespace == "test"
        assert loaded_term1.definition == "Test term 1"
        assert loaded_term1.alt_ids == ["GO:001a"]
        assert loaded_term1.is_a_parents == ["GO:002"]
        assert loaded_term1.relationships == {"part_of": ["GO:003"]}
        assert loaded_term1.synonyms == ["synonym1"]
        assert loaded_term1.xrefs == ["xref:001"]
        assert loaded_term1.subsets == ["subset1"]

    def test_save_load_roundtrip(self, tmp_path: Path):
        """Test round-trip serialization preserves structure."""
        onto = Ontology()
        term = Term(
            term_id="GO:001",
            name="test",
            namespace="test",
            is_a_parents=["GO:002"],
        )
        onto.add_term(term)

        # Round trip
        save_path = tmp_path / "roundtrip.json"
        save_ontology(onto, save_path)
        loaded = load_ontology(save_path)

        # Verify structure
        assert loaded.num_terms() == onto.num_terms()
        assert loaded.has_term("GO:001")
        assert "GO:002" in loaded.parents_of.get("GO:001", set())

    def test_load_nonexistent_file(self):
        """Test loading nonexistent file raises error."""
        from metainformant.core.utils.errors import IOError

        with pytest.raises(IOError, match="not found"):
            load_ontology("nonexistent.json")

    def test_save_invalid_ontology(self, tmp_path: Path):
        """Test saving empty ontology."""
        onto = Ontology()
        save_path = tmp_path / "empty.json"
        save_ontology(onto, save_path)
        
        loaded = load_ontology(save_path)
        assert loaded.num_terms() == 0


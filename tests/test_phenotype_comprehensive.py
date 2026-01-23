"""Comprehensive tests for phenotype module.

Tests cover AntWiki data loading and parsing functionality.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metainformant.phenotype.data.antwiki import load_antwiki_json
from metainformant.core.utils.errors import IOError as CoreIOError, ValidationError


class TestAntWikiLoading:
    """Tests for AntWiki JSON loading."""

    def test_load_list_format(self, tmp_path: Path):
        """Test loading JSON file with list format."""
        data = [
            {"species": "Camponotus pennsylvanicus", "measurements": {"length": 10.0}},
            {"species": "Formica rufa", "measurements": {"length": 8.0}},
        ]

        json_file = tmp_path / "antwiki.json"
        json_file.write_text(json.dumps(data))

        loaded = load_antwiki_json(json_file)

        assert len(loaded) == 2
        assert loaded[0]["species"] == "Camponotus pennsylvanicus"

    def test_load_dict_format(self, tmp_path: Path):
        """Test loading JSON file with single dict format."""
        data = {"species": "Camponotus pennsylvanicus", "measurements": {"length": 10.0}}

        json_file = tmp_path / "antwiki.json"
        json_file.write_text(json.dumps(data))

        loaded = load_antwiki_json(json_file)

        assert len(loaded) == 1
        assert isinstance(loaded, list)
        assert loaded[0]["species"] == "Camponotus pennsylvanicus"

    def test_load_empty_list(self, tmp_path: Path):
        """Test loading empty JSON."""
        json_file = tmp_path / "empty.json"
        json_file.write_text(json.dumps([]))

        loaded = load_antwiki_json(json_file)

        assert loaded == []

    def test_invalid_json(self, tmp_path: Path):
        """Test handling invalid JSON."""
        json_file = tmp_path / "invalid.json"
        json_file.write_text("not valid json")

        with pytest.raises((CoreIOError, json.JSONDecodeError, ValueError)):
            load_antwiki_json(json_file)

    def test_complex_antwiki_structure(self, tmp_path: Path):
        """Test loading complex AntWiki data structure."""
        data = [
            {
                "species": "Camponotus pennsylvanicus",
                "measurements": {"worker_length_mm": [6.0, 13.0], "head_width_mm": [1.8, 3.2]},
                "traits": ["arboreal", "carnivorous", "polygynous"],
                "distribution": "North America",
            },
            {
                "species": "Formica rufa",
                "measurements": {"worker_length_mm": [4.5, 9.0]},
                "traits": ["terrestrial", "omnivorous"],
            },
        ]

        json_file = tmp_path / "antwiki.json"
        json_file.write_text(json.dumps(data))

        loaded = load_antwiki_json(json_file)

        assert len(loaded) == 2
        assert "measurements" in loaded[0]
        assert "traits" in loaded[0]
        assert len(loaded[0]["traits"]) == 3
        assert loaded[0]["measurements"]["worker_length_mm"] == [6.0, 13.0]

    def test_validation_errors(self, tmp_path: Path):
        """Test validation errors for invalid data structures."""
        # Test invalid top-level structure
        json_file = tmp_path / "invalid.json"
        json_file.write_text(json.dumps("not a list or dict"))

        with pytest.raises(ValidationError):
            load_antwiki_json(json_file)

        # Test missing species field
        invalid_data = [{"measurements": {"length": 10.0}}]
        json_file.write_text(json.dumps(invalid_data))

        with pytest.raises(ValidationError):
            load_antwiki_json(json_file)

        # Test invalid measurements type
        invalid_data2 = [{"species": "Test", "measurements": "not a dict"}]
        json_file.write_text(json.dumps(invalid_data2))

        with pytest.raises(ValidationError):
            load_antwiki_json(json_file)

        # Test invalid traits type
        invalid_data3 = [{"species": "Test", "traits": "not a list"}]
        json_file.write_text(json.dumps(invalid_data3))

        with pytest.raises(ValidationError):
            load_antwiki_json(json_file)

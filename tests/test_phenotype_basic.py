"""Basic tests for phenotype module functionality."""

import json
from pathlib import Path

import pytest

from metainformant.phenotype.data.antwiki import load_antwiki_json
from metainformant.core.utils.errors import IOError as CoreIOError, ValidationError


class TestAntWiki:
    """Test AntWiki phenotype data loading."""

    def test_load_antwiki_json_with_list(self, tmp_path: Path):
        """Test loading AntWiki JSON with list format."""
        test_data = [
            {"species": "Camponotus pennsylvanicus", "trait": "body_length", "value": 12.5},
            {"species": "Formica rufa", "trait": "head_width", "value": 2.1}
        ]

        test_json_file = tmp_path / "antwiki_test.json"
        test_json_file.write_text(json.dumps(test_data))

        result = load_antwiki_json(test_json_file)

        assert len(result) == 2
        assert result[0]["species"] == "Camponotus pennsylvanicus"
        assert result[1]["species"] == "Formica rufa"

    def test_load_antwiki_json_with_dict(self, tmp_path: Path):
        """Test loading AntWiki JSON with dict format."""
        test_data = {"species": "Lasius niger", "trait": "worker_number", "value": 5000}

        test_json_file = tmp_path / "antwiki_test.json"
        test_json_file.write_text(json.dumps(test_data))

        result = load_antwiki_json(test_json_file)

        assert len(result) == 1
        assert result[0]["species"] == "Lasius niger"

    def test_load_antwiki_json_empty(self, tmp_path: Path):
        """Test loading empty AntWiki JSON."""
        test_data = []

        test_json_file = tmp_path / "antwiki_test.json"
        test_json_file.write_text(json.dumps(test_data))

        result = load_antwiki_json(test_json_file)

        assert len(result) == 0

    def test_load_antwiki_json_invalid_format(self, tmp_path: Path):
        """Test loading AntWiki JSON with invalid format."""
        test_json_file = tmp_path / "invalid.json"
        test_json_file.write_text("invalid json")

        # Should raise IOError for invalid JSON
        with pytest.raises((CoreIOError, json.JSONDecodeError, ValueError)):
            load_antwiki_json(test_json_file)

    def test_load_antwiki_json_nonexistent_file(self):
        """Test loading from nonexistent file."""
        nonexistent_file = Path("nonexistent_file.json")

        # Should raise FileNotFoundError for missing file
        with pytest.raises(FileNotFoundError):
            load_antwiki_json(nonexistent_file)

    def test_antwiki_json_structure(self, tmp_path: Path):
        """Test that AntWiki data has expected structure."""
        # Create valid test data
        test_data = [
            {"species": "Camponotus pennsylvanicus", "measurements": {"length": 10.0}, "traits": ["arboreal"]}
        ]
        test_file = tmp_path / "antwiki_test.json"
        test_file.write_text(json.dumps(test_data))
        
        # Load test data
        result = load_antwiki_json(test_file)

        if len(result) > 0:
            # Check that each entry has expected fields
            entry = result[0]
            assert "species" in entry or "taxon" in entry
            assert isinstance(entry, dict)
            
    def test_antwiki_json_validation(self, tmp_path: Path):
        """Test data validation."""
        # Test missing species field
        invalid_data = [{"measurements": {"length": 10.0}}]
        test_file = tmp_path / "invalid.json"
        test_file.write_text(json.dumps(invalid_data))
        
        with pytest.raises(ValidationError):
            load_antwiki_json(test_file)
            
        # Test with validate=False should work
        result = load_antwiki_json(test_file, validate=False)
        assert len(result) == 1

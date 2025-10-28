"""Basic tests for phenotype module functionality."""

import json
from pathlib import Path

import pytest

from metainformant.phenotype import antwiki


class TestAntWiki:
    """Test AntWiki phenotype data loading."""

    def setup_method(self):
        """Setup test environment."""
        self.test_data_dir = Path("tests/data/phenotype")
        self.test_json_file = self.test_data_dir / "antwiki_dataset_sorted_final_01.json"

    def test_load_antwiki_json_with_list(self):
        """Test loading AntWiki JSON with list format."""
        test_data = [
            {"species": "Camponotus pennsylvanicus", "trait": "body_length", "value": 12.5},
            {"species": "Formica rufa", "trait": "head_width", "value": 2.1}
        ]

        with open(self.test_json_file, 'w') as f:
            json.dump(test_data, f)

        result = antwiki.load_antwiki_json(self.test_json_file)

        assert len(result) == 2
        assert result[0]["species"] == "Camponotus pennsylvanicus"
        assert result[1]["species"] == "Formica rufa"

    def test_load_antwiki_json_with_dict(self):
        """Test loading AntWiki JSON with dict format."""
        test_data = {"species": "Lasius niger", "trait": "worker_number", "value": 5000}

        with open(self.test_json_file, 'w') as f:
            json.dump(test_data, f)

        result = antwiki.load_antwiki_json(self.test_json_file)

        assert len(result) == 1
        assert result[0]["species"] == "Lasius niger"

    def test_load_antwiki_json_empty(self):
        """Test loading empty AntWiki JSON."""
        test_data = []

        with open(self.test_json_file, 'w') as f:
            json.dump(test_data, f)

        result = antwiki.load_antwiki_json(self.test_json_file)

        assert len(result) == 0

    def test_load_antwiki_json_invalid_format(self):
        """Test loading AntWiki JSON with invalid format."""
        test_data = "invalid json"

        with open(self.test_json_file, 'w') as f:
            f.write(test_data)

        result = antwiki.load_antwiki_json(self.test_json_file)

        # Should return empty list for invalid data
        assert len(result) == 0

    def test_load_antwiki_json_nonexistent_file(self):
        """Test loading from nonexistent file."""
        nonexistent_file = Path("nonexistent_file.json")

        result = antwiki.load_antwiki_json(nonexistent_file)

        # Should return empty list for missing file
        assert len(result) == 0

    def test_antwiki_json_structure(self):
        """Test that AntWiki data has expected structure."""
        # Load real test data
        result = antwiki.load_antwiki_json(self.test_json_file)

        if len(result) > 0:
            # Check that each entry has expected fields
            entry = result[0]
            assert "species" in entry or "taxon" in entry
            assert isinstance(entry, dict)

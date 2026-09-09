import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

# Paths
SCRIPT_PATH = Path("scripts/rna/normalize_tissue_metadata.py")
MAPPING_PATH = Path("config/amalgkit/tissue_mapping.yaml")
PATCHES_PATH = Path("config/amalgkit/tissue_patches.yaml")


class TestTissueNormalizationScript:
    def test_script_execution(self, tmp_path):
        """Test that the normalization script runs successfully."""
        if not SCRIPT_PATH.exists():
            pytest.skip("Normalization script not found")

        # Create dummy metadata
        metadata_path = tmp_path / "metadata.tsv"
        output_path = tmp_path / "metadata_normalized.tsv"

        df = pd.DataFrame(
            {
                "run": ["SRR1", "SRR2", "SRR3", "SRR4"],
                "bioproject": ["PRJ1", "PRJ1", "PRJ2", "PRJ3"],
                "tissue": ["brain", "Brain", "whole body", "unknown_tissue"],
                "biosample": ["BS1", "BS2", "BS3", "BS4"],
            }
        )
        df.to_csv(metadata_path, sep="\t", index=False)

        # Run script
        cmd = [
            sys.executable,
            str(SCRIPT_PATH),
            "--input",
            str(metadata_path),
            "--output",
            str(output_path),
            "--mapping",
            str(MAPPING_PATH),
            "--patches",
            str(PATCHES_PATH),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        # Check output
        assert output_path.exists()
        df_out = pd.read_csv(output_path, sep="\t")

        assert "tissue_normalized" in df_out.columns
        # "brain" -> "brain"
        # "Brain" -> "brain" (case insensitive)
        # "whole body" -> "whole_body" (synonym)

        assert df_out.loc[0, "tissue_normalized"] == "brain"
        assert df_out.loc[1, "tissue_normalized"] == "brain"
        assert df_out.loc[2, "tissue_normalized"] == "whole_body"

    def test_script_with_patches(self, tmp_path):
        """Test script interactions with patches."""
        if not SCRIPT_PATH.exists() or not PATCHES_PATH.exists():
            pytest.skip("Script or patches missing")

        # Create metadata that relies on patches
        # e.g. PRJNA339620 -> mushroom_body (from tissue_patches.yaml)
        metadata_path = tmp_path / "metadata_patched.tsv"
        output_path = tmp_path / "metadata_patched_out.tsv"

        df = pd.DataFrame(
            {"run": ["SRR_PATCH_TEST"], "bioproject": ["PRJNA339620"], "tissue": [""]}  # Should map to mushroom_body
        )
        df.to_csv(metadata_path, sep="\t", index=False)

        cmd = [
            sys.executable,
            str(SCRIPT_PATH),
            "--input",
            str(metadata_path),
            "--output",
            str(output_path),
            "--mapping",
            str(MAPPING_PATH),
            "--patches",
            str(PATCHES_PATH),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        df_out = pd.read_csv(output_path, sep="\t")
        assert df_out.loc[0, "tissue_normalized"] == "mushroom_body"


class TestTissueNormalizerModule:
    """Direct module API tests for metainformant.rna.amalgkit.tissue_normalizer."""

    def _write_mapping(self, path: Path, payload: object) -> Path:
        with open(path, "w", encoding="utf-8") as handle:
            yaml.dump(payload, handle)
        return path

    def test_load_tissue_mapping_empty_file_returns_empty_dict(self, tmp_path: Path):
        """An existing-but-empty mapping file degrades to no mapping."""
        from metainformant.rna.amalgkit.tissue_normalizer import load_tissue_mapping

        mapping_path = tmp_path / "tissue_mapping.yaml"
        mapping_path.write_text("", encoding="utf-8")
        assert load_tissue_mapping(mapping_path) == {}

    def test_load_tissue_missing_file_returns_empty_dict(self, tmp_path: Path):
        """A missing mapping file returns {} rather than raising."""
        from metainformant.rna.amalgkit.tissue_normalizer import load_tissue_mapping

        assert load_tissue_mapping(tmp_path / "absent.yaml") == {}

    def test_load_tissue_mapping_ignores_non_list_values(self, tmp_path: Path):
        """Scalar keys (comments-as-keys) are filtered out of the mapping."""
        from metainformant.rna.amalgkit.tissue_normalizer import load_tissue_mapping

        mapping_path = self._write_mapping(
            tmp_path / "tissue_mapping.yaml",
            {"brain": ["Brain", "brain"], "note": "scalar entry"},
        )
        assert load_tissue_mapping(mapping_path) == {"brain": ["Brain", "brain"]}

    def test_normalize_tissue_exact_case_insensitive_and_prefix(self):
        """Exact match, case-insensitive match, and prefix fallback all resolve."""
        from metainformant.rna.amalgkit.tissue_normalizer import build_synonym_lookup, normalize_tissue

        lookup = build_synonym_lookup({"fat_body": ["fat body", "FB"], "brain": ["brain"]})
        assert normalize_tissue("brain", lookup) == "brain"
        assert normalize_tissue("  BRAIN ", lookup) == "brain"
        assert normalize_tissue("fat body of 1 queen; kept in a group", lookup) == "fat_body"

    def test_normalize_tissue_unmapped_returns_default(self):
        """Unmapped values return the caller-requested default (default "" contract)."""
        from metainformant.rna.amalgkit.tissue_normalizer import build_synonym_lookup, normalize_tissue

        lookup = build_synonym_lookup({"brain": ["brain"]})
        assert normalize_tissue("wing disc", lookup) == ""
        assert normalize_tissue("wing disc", lookup, default="unknown") == "unknown"
        assert normalize_tissue("", lookup) == ""

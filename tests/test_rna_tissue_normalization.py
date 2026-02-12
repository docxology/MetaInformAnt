
import pytest
import pandas as pd
import subprocess
import sys
from pathlib import Path

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
        
        df = pd.DataFrame({
            "run": ["SRR1", "SRR2", "SRR3", "SRR4"],
            "bioproject": ["PRJ1", "PRJ1", "PRJ2", "PRJ3"],
            "tissue": ["brain", "Brain", "whole body", "unknown_tissue"],
            "biosample": ["BS1", "BS2", "BS3", "BS4"]
        })
        df.to_csv(metadata_path, sep="\t", index=False)
        
        # Run script
        cmd = [
            sys.executable, str(SCRIPT_PATH),
            "--input", str(metadata_path),
            "--output", str(output_path),
            "--mapping", str(MAPPING_PATH),
            "--patches", str(PATCHES_PATH)
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
        
        df = pd.DataFrame({
            "run": ["SRR_PATCH_TEST"],
            "bioproject": ["PRJNA339620"], # Should map to mushroom_body
            "tissue": [""] 
        })
        df.to_csv(metadata_path, sep="\t", index=False)
        
        cmd = [
            sys.executable, str(SCRIPT_PATH),
            "--input", str(metadata_path),
            "--output", str(output_path),
            "--mapping", str(MAPPING_PATH),
            "--patches", str(PATCHES_PATH)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Script failed: {result.stderr}"
        
        df_out = pd.read_csv(output_path, sep="\t")
        assert df_out.loc[0, "tissue_normalized"] == "mushroom_body"


import pytest
import yaml
from pathlib import Path
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

CONFIG_DIR = Path("config/amalgkit")
TISSUE_PATCHES = CONFIG_DIR / "tissue_patches.yaml"
TISSUE_MAPPING = CONFIG_DIR / "tissue_mapping.yaml"

class TestTissueConfiguration:
    def test_tissue_patches_integrity(self):
        """Verify tissue_patches.yaml structure and content."""
        if not TISSUE_PATCHES.exists():
            pytest.skip(f"{TISSUE_PATCHES} not found")
            
        with open(TISSUE_PATCHES) as f:
            data = yaml.safe_load(f)
            
        assert "samples" in data, "Missing 'samples' section"
        assert "bioprojects" in data, "Missing 'bioprojects' section"
        
        # Verify samples structure
        if data["samples"]:
            for sample_id, tissue in data["samples"].items():
                assert isinstance(sample_id, str), f"Sample ID {sample_id} must be string"
                assert isinstance(tissue, str), f"Tissue {tissue} for {sample_id} must be string"
                
        # Verify bioprojects structure
        if data["bioprojects"]:
            for bp_id, tissue in data["bioprojects"].items():
                assert isinstance(bp_id, str), f"BioProject ID {bp_id} must be string"
                assert bp_id.startswith(("PRJ", "SRP", "ERP", "DRP")), f"Invalid BioProject ID format: {bp_id}"
                assert isinstance(tissue, str), f"Tissue {tissue} for {bp_id} must be string"

    def test_tissue_mapping_integrity(self):
        """Verify tissue_mapping.yaml structure and valid mapping targets."""
        if not TISSUE_MAPPING.exists():
            pytest.skip(f"{TISSUE_MAPPING} not found")
            
        with open(TISSUE_MAPPING) as f:
            mapping = yaml.safe_load(f)
            
        assert isinstance(mapping, dict), "Root must be a dictionary"
        
        # Check for duplicates across all synonym lists
        seen_synonyms = {}
        
        for canonical, synonyms in mapping.items():
            assert isinstance(synonyms, list), f"Synonyms for {canonical} must be a list"
            assert len(synonyms) > 0, f"No synonyms defined for {canonical}"
            
            # Check consistency
            assert canonical in synonyms, f"Canonical name '{canonical}' should be in its own synonym list for completeness"
            
            for syn in synonyms:
                if syn in seen_synonyms:
                    prev_canonical = seen_synonyms[syn]
                    pytest.fail(f"Duplicate synonym '{syn}' found in '{canonical}' and '{prev_canonical}'")
                seen_synonyms[syn] = canonical

    def test_cross_file_consistency(self):
        """Ensure tissues used in patches exist in mapping."""
        if not TISSUE_PATCHES.exists() or not TISSUE_MAPPING.exists():
            pytest.skip("Config files missing")
            
        with open(TISSUE_PATCHES) as f:
            patches = yaml.safe_load(f)
            
        with open(TISSUE_MAPPING) as f:
            mapping = yaml.safe_load(f)
            
        valid_tissues = set(mapping.keys())
        
        # Check bioproject patches
        if patches.get("bioprojects"):
            for bp, tissue in patches["bioprojects"].items():
                assert tissue in valid_tissues, f"BioProject {bp} uses undefined tissue '{tissue}'. Add to tissue_mapping.yaml first."
                
        # Check sample patches
        if patches.get("samples"):
            for sample, tissue in patches["samples"].items():
                assert tissue in valid_tissues, f"Sample {sample} uses undefined tissue '{tissue}'. Add to tissue_mapping.yaml first."

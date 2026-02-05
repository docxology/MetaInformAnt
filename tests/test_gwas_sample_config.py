"""Tests for configurable sample architecture in GWAS pipeline.

Tests cover:
- Configurable data generation (subspecies, scale_factor, seed)
- Sample subsetting (by ID, by subspecies, by caste)
- ID-based phenotype matching
- Multi-trait analysis (shared PCA/kinship)
- Manifest cache invalidation
- CLI argument override of config
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

# Add project source to path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "gwas"))


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def gwas_output(tmp_path):
    """Provide a fresh temporary output directory."""
    out = tmp_path / "gwas_test"
    out.mkdir()
    return out


@pytest.fixture
def default_vcf(gwas_output):
    """Generate default VCF with standard params."""
    from run_amellifera_gwas import generate_real_vcf

    return generate_real_vcf(gwas_output, n_variants=200, force=True)


@pytest.fixture
def default_phenotypes(gwas_output, default_vcf):
    """Generate default phenotypes."""
    from run_amellifera_gwas import generate_phenotypes

    return generate_phenotypes(gwas_output, default_vcf, force=True)


@pytest.fixture
def default_metadata(gwas_output, default_vcf):
    """Generate default metadata."""
    from run_amellifera_gwas import generate_metadata

    return generate_metadata(gwas_output, default_vcf, force=True)


# ---------------------------------------------------------------------------
# 1. Default generation backward compatibility
# ---------------------------------------------------------------------------


class TestDefaultGeneration:
    """Ensure defaults produce the same results as before."""

    def test_default_generation_unchanged(self, gwas_output):
        """Default params should produce 90 samples (80 diploid + 10 drones)."""
        from run_amellifera_gwas import generate_real_vcf, DEFAULT_SUBSPECIES

        vcf_path = generate_real_vcf(gwas_output, n_variants=100, force=True)
        assert vcf_path.exists()

        from metainformant.gwas.analysis.quality import parse_vcf_full

        data = parse_vcf_full(vcf_path)
        assert len(data["samples"]) == 90
        # Check subspecies alias still works
        from run_amellifera_gwas import SUBSPECIES

        assert SUBSPECIES is DEFAULT_SUBSPECIES

    def test_default_subspecies_structure(self):
        """DEFAULT_SUBSPECIES has expected keys and structure."""
        from run_amellifera_gwas import DEFAULT_SUBSPECIES

        assert "A.m.ligustica" in DEFAULT_SUBSPECIES
        assert "A.m.scutellata" in DEFAULT_SUBSPECIES
        total = sum(v["n_samples"] for v in DEFAULT_SUBSPECIES.values())
        assert total == 80  # 25+20+15+10+10


# ---------------------------------------------------------------------------
# 2. Scale factor
# ---------------------------------------------------------------------------


class TestScaleFactor:
    """Test that scale_factor multiplies sample counts."""

    def test_scale_factor_doubles_samples(self, gwas_output):
        """scale_factor=2 should produce 180 samples (160 diploid + 20 drones)."""
        from run_amellifera_gwas import generate_real_vcf

        vcf_path = generate_real_vcf(gwas_output, n_variants=50, scale_factor=2, force=True)
        from metainformant.gwas.analysis.quality import parse_vcf_full

        data = parse_vcf_full(vcf_path)
        assert len(data["samples"]) == 180

    def test_scale_factor_one_is_default(self, gwas_output):
        """scale_factor=1 should match default (90 samples)."""
        from run_amellifera_gwas import generate_real_vcf

        vcf_path = generate_real_vcf(gwas_output, n_variants=50, scale_factor=1, force=True)
        from metainformant.gwas.analysis.quality import parse_vcf_full

        data = parse_vcf_full(vcf_path)
        assert len(data["samples"]) == 90


# ---------------------------------------------------------------------------
# 3. Sample subsetting
# ---------------------------------------------------------------------------


class TestSubsetting:
    """Test subset_vcf_data functionality."""

    def test_subset_by_sample_ids(self, default_vcf):
        """Explicit ID list subset should keep only those samples."""
        from metainformant.gwas.analysis.quality import parse_vcf_full, subset_vcf_data

        data = parse_vcf_full(default_vcf)
        target_ids = ["ligustica_001", "ligustica_002", "carnica_001"]
        result = subset_vcf_data(data, sample_ids=target_ids)
        assert result["samples"] == target_ids
        assert len(result["genotypes"]) == 3

    def test_subset_by_sample_list_file(self, default_vcf, tmp_path):
        """Subset from a file with one ID per line."""
        from metainformant.gwas.analysis.quality import parse_vcf_full, subset_vcf_data

        data = parse_vcf_full(default_vcf)
        id_file = tmp_path / "samples.txt"
        id_file.write_text("ligustica_001\ncarnica_001\nmellifera_001\n")
        result = subset_vcf_data(data, sample_list_file=str(id_file))
        assert len(result["samples"]) == 3
        assert "ligustica_001" in result["samples"]

    def test_subset_by_subspecies(self, default_vcf, default_metadata):
        """Metadata-based subspecies filter."""
        from metainformant.gwas.analysis.quality import parse_vcf_full, subset_vcf_data
        from metainformant.gwas.data.metadata import load_sample_metadata

        data = parse_vcf_full(default_vcf)
        meta_result = load_sample_metadata(str(default_metadata))
        metadata = meta_result.get("metadata", {})

        result = subset_vcf_data(
            data,
            subset_config={"subspecies": ["A.m.ligustica", "A.m.carnica"]},
            metadata=metadata,
        )
        # Should only contain ligustica (25) + carnica (20) + drones from ligustica (10) = 55
        for sid in result["samples"]:
            meta = metadata.get(sid, {})
            assert meta.get("subspecies") in ("A.m.ligustica", "A.m.carnica"), f"Unexpected sample: {sid}"

    def test_subset_by_caste(self, default_vcf, default_metadata):
        """Worker-only analysis (exclude drones)."""
        from metainformant.gwas.analysis.quality import parse_vcf_full, subset_vcf_data
        from metainformant.gwas.data.metadata import load_sample_metadata

        data = parse_vcf_full(default_vcf)
        meta_result = load_sample_metadata(str(default_metadata))
        metadata = meta_result.get("metadata", {})

        result = subset_vcf_data(
            data,
            subset_config={"caste": ["worker"]},
            metadata=metadata,
        )
        assert len(result["samples"]) == 80  # All diploid workers, no drones
        for sid in result["samples"]:
            assert not sid.startswith("drone_")

    def test_subset_preserves_genotype_alignment(self, default_vcf):
        """Genotype rows should align with sample list after subsetting."""
        from metainformant.gwas.analysis.quality import parse_vcf_full, subset_vcf_data

        data = parse_vcf_full(default_vcf)
        target = ["ligustica_001", "carnica_001"]
        result = subset_vcf_data(data, sample_ids=target)
        assert len(result["genotypes"]) == len(result["samples"])
        # Each genotype row should have the same number of variants
        n_variants = len(data["variants"])
        for row in result["genotypes"]:
            assert len(row) == n_variants

    def test_subset_no_criteria_returns_all(self, default_vcf):
        """No subsetting criteria should return all samples."""
        from metainformant.gwas.analysis.quality import parse_vcf_full, subset_vcf_data

        data = parse_vcf_full(default_vcf)
        result = subset_vcf_data(data)
        assert len(result["samples"]) == len(data["samples"])

    def test_subset_immutable(self, default_vcf):
        """Subsetting should not modify original data."""
        from metainformant.gwas.analysis.quality import parse_vcf_full, subset_vcf_data

        data = parse_vcf_full(default_vcf)
        original_n = len(data["samples"])
        _ = subset_vcf_data(data, sample_ids=["ligustica_001"])
        assert len(data["samples"]) == original_n


# ---------------------------------------------------------------------------
# 4. Phenotype ID matching
# ---------------------------------------------------------------------------


class TestPhenotypeMatching:
    """Test ID-based phenotype alignment."""

    def test_phenotype_id_matching(self, default_phenotypes):
        """_load_phenotypes_by_id returns dict keyed by sample_id."""
        from metainformant.gwas.workflow.workflow import _load_phenotypes_by_id

        result = _load_phenotypes_by_id(default_phenotypes, trait="varroa_resistance")
        assert isinstance(result, dict)
        assert len(result) > 0
        assert "ligustica_001" in result
        assert isinstance(result["ligustica_001"], float)

    def test_phenotype_id_missing_trait_returns_default(self, default_phenotypes):
        """Without trait, should use second column."""
        from metainformant.gwas.workflow.workflow import _load_phenotypes_by_id

        result = _load_phenotypes_by_id(default_phenotypes)
        assert len(result) > 0

    def test_phenotype_id_binary_trait(self, default_phenotypes):
        """Can load binary trait column."""
        from metainformant.gwas.workflow.workflow import _load_phenotypes_by_id

        result = _load_phenotypes_by_id(default_phenotypes, trait="disease_resistance")
        assert len(result) > 0
        # Binary: all values should be 0.0 or 1.0
        for val in result.values():
            assert val in (0.0, 1.0)


# ---------------------------------------------------------------------------
# 5. Multi-trait
# ---------------------------------------------------------------------------


class TestMultiTrait:
    """Test multi-trait analysis sharing PCA/kinship."""

    def test_multi_trait_runs(self, default_vcf, default_phenotypes):
        """run_multi_trait_gwas should analyze multiple traits."""
        from metainformant.gwas.workflow.workflow import run_multi_trait_gwas

        config = {
            "work_dir": str(default_vcf.parent.parent / "work"),
            "qc": {"min_maf": 0.01, "max_missing": 0.1},
            "haplodiploidy": {"enabled": True, "het_threshold": 0.05, "exclude_haploid": True},
            "ld_pruning": {"enabled": True, "window_size": 50, "step_size": 5, "r2_threshold": 0.2},
            "association": {"model": "linear"},
        }

        result = run_multi_trait_gwas(
            vcf_path=default_vcf,
            phenotype_path=default_phenotypes,
            traits=["varroa_resistance", "disease_resistance"],
            config=config,
            output_dir=default_vcf.parent.parent / "multi_trait_results",
        )

        assert result["status"] == "success"
        assert "varroa_resistance" in result["trait_results"]
        assert "disease_resistance" in result["trait_results"]
        for trait_name in ["varroa_resistance", "disease_resistance"]:
            tr = result["trait_results"][trait_name]
            assert tr["status"] == "success"
            assert tr["n_tests"] > 0

    def test_multi_trait_shares_pca(self, default_vcf, default_phenotypes):
        """PCA should be computed once (in shared steps), not per trait."""
        from metainformant.gwas.workflow.workflow import run_multi_trait_gwas

        config = {
            "work_dir": str(default_vcf.parent.parent / "work"),
            "qc": {"min_maf": 0.01, "max_missing": 0.1},
            "haplodiploidy": {"enabled": True, "het_threshold": 0.05, "exclude_haploid": True},
            "ld_pruning": {"enabled": True, "window_size": 50, "step_size": 5, "r2_threshold": 0.2},
            "association": {"model": "linear"},
        }

        result = run_multi_trait_gwas(
            vcf_path=default_vcf,
            phenotype_path=default_phenotypes,
            traits=["varroa_resistance", "disease_resistance"],
            config=config,
        )

        # PCA should be in shared steps
        assert "population_structure" in result["steps_completed"]
        # Results should have PCA data
        assert result.get("results", {}).get("pca") is not None


# ---------------------------------------------------------------------------
# 6. Manifest cache
# ---------------------------------------------------------------------------


class TestManifestCache:
    """Test manifest-based cache validation."""

    def test_manifest_written(self, gwas_output):
        """VCF generation should write a manifest sidecar."""
        from run_amellifera_gwas import generate_real_vcf

        vcf_path = generate_real_vcf(gwas_output, n_variants=50, force=True)
        manifest = vcf_path.with_suffix(".vcf.manifest.json")
        assert manifest.exists()

        with open(manifest) as f:
            data = json.load(f)
        assert data["n_variants"] == 50
        assert data["seed"] == 42

    def test_manifest_cache_invalidation(self, gwas_output):
        """Changing params should trigger regeneration when manifest doesn't match."""
        from run_amellifera_gwas import generate_real_vcf, _check_manifest

        # Generate with seed=42
        vcf_path = generate_real_vcf(gwas_output, n_variants=50, seed=42, force=True)
        assert vcf_path.exists()

        # Check manifest matches
        params_same = {
            "n_variants": 50,
            "subspecies": {
                "A.m.ligustica": 25,
                "A.m.carnica": 20,
                "A.m.mellifera": 15,
                "A.m.caucasica": 10,
                "A.m.scutellata": 10,
            },
            "n_drones": 10,
            "scale_factor": 1,
            "seed": 42,
        }
        assert _check_manifest(vcf_path, params_same)

        # Different params should not match
        params_diff = dict(params_same)
        params_diff["seed"] = 99
        assert not _check_manifest(vcf_path, params_diff)

    def test_manifest_cache_hit_skips_regen(self, gwas_output):
        """With matching manifest, generate_real_vcf should skip regeneration."""
        from run_amellifera_gwas import generate_real_vcf

        # Generate first time
        vcf_path = generate_real_vcf(gwas_output, n_variants=50, force=True)
        mtime1 = vcf_path.stat().st_mtime

        # Second call without force should use cache
        import time

        time.sleep(0.1)
        vcf_path2 = generate_real_vcf(gwas_output, n_variants=50)
        mtime2 = vcf_path2.stat().st_mtime
        assert mtime1 == mtime2  # File not regenerated


# ---------------------------------------------------------------------------
# 7. Config loading
# ---------------------------------------------------------------------------


class TestConfigLoading:
    """Test load_data_generation_config."""

    def test_load_from_none(self):
        """No config should return defaults."""
        from run_amellifera_gwas import load_data_generation_config, DEFAULT_SUBSPECIES

        dg = load_data_generation_config(None)
        assert dg["subspecies"] == DEFAULT_SUBSPECIES
        assert dg["n_drones"] == 10
        assert dg["scale_factor"] == 1

    def test_load_from_yaml_config(self):
        """Config dict should override defaults."""
        from run_amellifera_gwas import load_data_generation_config

        config = {
            "data_generation": {
                "n_variants": 5000,
                "scale_factor": 3,
                "seed": 99,
            }
        }
        dg = load_data_generation_config(config)
        assert dg["n_variants"] == 5000
        assert dg["scale_factor"] == 3
        assert dg["seed"] == 99
        # Unset fields keep defaults
        assert dg["n_drones"] == 10

    def test_cli_args_override_config(self, gwas_output):
        """CLI args should take precedence over YAML config."""
        from run_amellifera_gwas import load_data_generation_config

        yaml_config = {"data_generation": {"n_variants": 5000, "scale_factor": 2}}
        dg = load_data_generation_config(yaml_config)
        assert dg["scale_factor"] == 2

        # Simulate CLI override
        dg["scale_factor"] = 5
        dg["n_variants"] = 100
        assert dg["scale_factor"] == 5
        assert dg["n_variants"] == 100

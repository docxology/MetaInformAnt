"""Test GWAS configuration for P. barbatus species."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.gwas import execute_gwas_workflow, load_gwas_config


def test_load_pbarbatus_config() -> None:
    """Test loading P. barbatus GWAS configuration."""
    config_path = Path("config/gwas/gwas_pbarbatus.yaml")

    if not config_path.exists():
        pytest.skip(f"Configuration file not found: {config_path}")

    cfg = load_gwas_config(config_path)

    # Verify basic configuration
    assert cfg.work_dir is not None
    assert "pbarbatus" in str(cfg.work_dir)

    # Verify genome configuration
    assert cfg.genome is not None
    assert cfg.genome.get("accession") == "GCF_000187915.1"
    assert cfg.genome.get("dest_dir") is not None
    assert "genome" in str(cfg.genome.get("dest_dir"))
    assert cfg.genome.get("include") is not None
    assert "genome" in cfg.genome.get("include", [])

    # Verify variants configuration
    assert cfg.variants is not None

    # Verify QC configuration
    assert cfg.qc is not None
    assert "min_maf" in cfg.qc

    # Verify structure configuration
    assert cfg.structure is not None
    assert cfg.structure.get("compute_pca") is True or cfg.structure.get("compute_pca") is None

    # Verify association configuration
    assert cfg.association is not None
    assert "model" in cfg.association


def test_pbarbatus_config_validation() -> None:
    """Test configuration validation for P. barbatus."""
    config_path = Path("config/gwas/gwas_pbarbatus.yaml")

    if not config_path.exists():
        pytest.skip(f"Configuration file not found: {config_path}")

    cfg = load_gwas_config(config_path)

    # Test validation mode
    result = execute_gwas_workflow(cfg, check=True)

    assert result["status"] == "validated"
    assert "config" in result


def test_pbarbatus_config_parameters() -> None:
    """Test that P. barbatus config has all required parameters."""
    config_path = Path("config/gwas/gwas_pbarbatus.yaml")

    if not config_path.exists():
        pytest.skip(f"Configuration file not found: {config_path}")

    cfg = load_gwas_config(config_path)

    # Check genome parameters
    if cfg.genome:
        assert "accession" in cfg.genome
        assert cfg.genome["accession"] == "GCF_000187915.1"
        assert "dest_dir" in cfg.genome
        assert "include" in cfg.genome

    # Check variants parameters
    assert isinstance(cfg.variants, dict)

    # Check QC parameters
    assert "min_maf" in cfg.qc
    assert "max_missing" in cfg.qc
    assert "hwe_pval" in cfg.qc

    # Check structure parameters
    if cfg.structure:
        assert "compute_pca" in cfg.structure or cfg.structure.get("compute_pca") is None
        assert "n_components" in cfg.structure or cfg.structure.get("n_components") is None

    # Check association parameters
    assert "model" in cfg.association
    assert "trait" in cfg.association

    # Check correction parameters
    assert "method" in cfg.correction
    assert "alpha" in cfg.correction

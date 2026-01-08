"""Tests for GWAS configuration management."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.gwas import GWASWorkflowConfig, load_gwas_config
from metainformant.core.io.io import dump_json


def test_load_gwas_config_basic(tmp_path: Path) -> None:
    """Test loading basic GWAS configuration from YAML."""
    config_file = tmp_path / "test_config.yaml"
    config_content = """
work_dir: output/gwas/test
threads: 4
"""
    config_file.write_text(config_content)

    config = load_gwas_config(config_file)
    assert isinstance(config, GWASWorkflowConfig)
    assert config.work_dir == Path("output/gwas/test").expanduser().resolve()
    assert config.threads == 4


def test_load_gwas_config_with_all_sections(tmp_path: Path) -> None:
    """Test loading complete GWAS configuration with all sections."""
    config_file = tmp_path / "full_config.yaml"
    config_content = """
work_dir: output/gwas/full_test
threads: 8
log_dir: output/gwas/logs

genome:
  accession: GCF_000001405.40
  dest_dir: output/gwas/genome

variants:
  vcf_files:
    - data/variants/test.vcf.gz

qc:
  min_maf: 0.01
  max_missing: 0.05

samples:
  phenotype_file: data/phenotypes.tsv

association:
  model: linear
  trait: height
"""
    config_file.write_text(config_content)

    config = load_gwas_config(config_file)
    assert config.genome is not None
    assert config.genome["accession"] == "GCF_000001405.40"
    assert config.variants["vcf_files"] == ["data/variants/test.vcf.gz"]
    assert config.qc["min_maf"] == 0.01
    assert config.association["model"] == "linear"


def test_load_gwas_config_json(tmp_path: Path) -> None:
    """Test loading GWAS configuration from JSON."""
    config_file = tmp_path / "test_config.json"
    config_data = {
        "work_dir": "output/gwas/test",
        "threads": 6,
        "genome": {
            "accession": "GCF_000001405.40"
        }
    }
    dump_json(config_data, config_file)

    config = load_gwas_config(config_file)
    assert config.threads == 6
    assert config.genome is not None
    assert config.genome["accession"] == "GCF_000001405.40"


def test_load_gwas_config_missing_file() -> None:
    """Test error handling for missing configuration file."""
    missing_file = Path("nonexistent_config.yaml")
    with pytest.raises((FileNotFoundError, ValueError)):
        load_gwas_config(missing_file)


def test_gwas_workflow_config_defaults() -> None:
    """Test GWASWorkflowConfig with default values."""
    config = GWASWorkflowConfig(work_dir=Path("output/test"))
    assert config.threads == 8
    assert config.log_dir is None
    assert config.genome is None
    assert config.variants == {}
    assert config.qc == {}
    assert config.samples == {}
    assert config.structure == {}
    assert config.association == {}
    assert config.correction == {}
    assert config.output == {}


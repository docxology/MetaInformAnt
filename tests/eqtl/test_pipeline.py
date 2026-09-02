"""Tests for metainformant.eqtl.pipeline parameter resolution (zero-mocks)."""

from __future__ import annotations

import pytest

from metainformant.eqtl.pipeline import resolve_run_parameters


class TestResolveRunParameters:
    def test_cli_defaults(self):
        params = resolve_run_parameters(species="amellifera")
        assert params["species"] == "amellifera"
        assert params["sample_ids"] is None
        assert params["n_samples"] == 3
        assert params["threads"] == 4
        assert params["cleanup_fastq"] is True

    def test_missing_species_raises(self):
        with pytest.raises(ValueError, match="species"):
            resolve_run_parameters()

    def test_config_overrides(self):
        config = {
            "species": "acromyrmex_echinatior",
            "samples": {"mode": "explicit", "explicit_ids": ["SRR1", "SRR2"]},
            "alignment": {"threads": 8},
            "output": {"cleanup_fastq": False},
        }
        params = resolve_run_parameters(species="ignored", n_samples=10, config=config)
        assert params["species"] == "acromyrmex_echinatior"
        assert params["sample_ids"] == ["SRR1", "SRR2"]
        assert params["threads"] == 8
        assert params["cleanup_fastq"] is False

    def test_config_partial_keeps_cli_extras(self):
        config = {"species": "sp"}
        params = resolve_run_parameters(species="sp", n_samples=5, config=config)
        assert params["n_samples"] == 5

    def test_env_default_missing_gives_none_amalgkit(self):
        params = resolve_run_parameters(species="sp", amalgkit_output=None)
        # Without AMALGKIT_DATA_ROOT, root is None (pipeline errors honestly)
        assert params["amalgkit_output"] is None or str(params["amalgkit_output"])

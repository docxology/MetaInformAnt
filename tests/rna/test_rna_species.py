"""Tests for configuration-derived RNA species discovery."""

from pathlib import Path

from metainformant.rna.engine.species import (
    configured_data_root,
    discover_species_config_names,
    discover_species_names,
    species_name_from_config,
)


def test_species_name_from_config_accepts_paths() -> None:
    assert species_name_from_config(Path("amalgkit_apis_mellifera.yaml")) == "apis_mellifera"


def test_discovery_excludes_non_species_configs(tmp_path: Path) -> None:
    for name in (
        "amalgkit_apis_mellifera.yaml",
        "amalgkit_solenopsis_invicta.yaml",
        "amalgkit_template.yaml",
        "amalgkit_test.yaml",
        "amalgkit_cross_species.yaml",
    ):
        (tmp_path / name).write_text("{}\n", encoding="utf-8")

    assert discover_species_config_names(tmp_path) == [
        "amalgkit_apis_mellifera.yaml",
        "amalgkit_solenopsis_invicta.yaml",
    ]
    assert discover_species_names(tmp_path) == ["apis_mellifera", "solenopsis_invicta"]


def test_configured_data_root_honors_environment(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(tmp_path / "mounted-data"))
    assert configured_data_root() == tmp_path / "mounted-data"

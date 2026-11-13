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


def _species_yaml(directory: Path, slug: str, declared: str, taxon: int | None = 7460) -> None:
    taxon_line = f"taxon_id: {taxon}\n" if taxon is not None else ""
    (directory / f"amalgkit_{slug}.yaml").write_text(f"species_list:\n  - {declared}\n{taxon_line}", encoding="utf-8")


def test_configured_species_reads_declared_identity(tmp_path: Path) -> None:
    from metainformant.rna.engine.species import configured_species

    _species_yaml(tmp_path, "apis_mellifera", "Apis_mellifera", 7460)
    _species_yaml(tmp_path, "solenopsis_invicta", "Solenopsis_invicta", None)
    (tmp_path / "amalgkit_template.yaml").write_text("species_list: [Homo_sapiens]\n", encoding="utf-8")

    records = configured_species(tmp_path)

    assert [(r.slug, r.display_name, r.taxon_id) for r in records] == [
        ("apis_mellifera", "Apis_mellifera", "7460"),
        ("solenopsis_invicta", "Solenopsis_invicta", None),
    ]
    assert records[0].config_path == tmp_path / "amalgkit_apis_mellifera.yaml"
    assert configured_species(tmp_path / "missing") == []


def test_configured_species_refuses_drift_and_duplicates(tmp_path: Path) -> None:
    import pytest

    from metainformant.rna.engine.species import configured_species

    _species_yaml(tmp_path, "apis_mellifera", "Bombus_terrestris")
    with pytest.raises(ValueError, match="differs from its filename slug"):
        configured_species(tmp_path)

    (tmp_path / "amalgkit_apis_mellifera.yaml").write_text("species_list: [A, B]\n", encoding="utf-8")
    with pytest.raises(ValueError, match="exactly one species_list entry"):
        configured_species(tmp_path)

    _species_yaml(tmp_path, "apis_mellifera", "Apis_mellifera", 1)
    _species_yaml(tmp_path, "bombus_terrestris", "Bombus_terrestris", 1)
    with pytest.raises(ValueError, match="duplicate taxon_id"):
        configured_species(tmp_path)


def test_resolve_species_work_dir_prefers_populated_alias(tmp_path: Path) -> None:
    from metainformant.rna.engine.species import resolve_species_work_dir, species_directory_candidates

    assert species_directory_candidates("apis_mellifera") == ("apis_mellifera", "apis_mellifera_all", "amellifera")
    assert species_directory_candidates("nasonia_vitripennis") == ("nasonia_vitripennis",)
    # Nothing exists yet: the canonical path, not created.
    assert resolve_species_work_dir(tmp_path, "apis_mellifera") == tmp_path / "apis_mellifera" / "work"
    assert not (tmp_path / "apis_mellifera").exists()
    # An empty canonical directory loses to an alias that holds work evidence.
    (tmp_path / "apis_mellifera" / "work").mkdir(parents=True)
    (tmp_path / "amellifera" / "work" / "quant").mkdir(parents=True)
    assert resolve_species_work_dir(tmp_path, "apis_mellifera") == tmp_path / "amellifera" / "work"
    # Evidence in the canonical tree wins.
    (tmp_path / "apis_mellifera" / "work" / "metadata").mkdir()
    assert resolve_species_work_dir(tmp_path, "apis_mellifera") == tmp_path / "apis_mellifera" / "work"
    # A bare alias directory is used only when the canonical one is missing.
    (tmp_path / "pbarbatus").mkdir()
    assert resolve_species_work_dir(tmp_path, "pogonomyrmex_barbatus") == tmp_path / "pbarbatus" / "work"


def test_orchestrator_work_dir_delegates_to_species_rule(monkeypatch, tmp_path: Path) -> None:
    from metainformant.rna.engine.streaming_orchestrator import _species_work_dir

    (tmp_path / "pbarbatus" / "work" / "index").mkdir(parents=True)
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(tmp_path))
    assert _species_work_dir("pogonomyrmex_barbatus") == tmp_path / "pbarbatus" / "work"

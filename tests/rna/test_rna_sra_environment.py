"""Tests for campaign-local SRA subprocess environments."""

from __future__ import annotations

from pathlib import Path

from metainformant.rna.amalgkit.sra_environment import build_sra_environment, ensure_ncbi_settings


def test_build_sra_environment_is_campaign_local_and_repeatable(tmp_path: Path, monkeypatch) -> None:
    data_root = tmp_path / "campaign"
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))
    monkeypatch.delenv("AMALGKIT_NCBI_SETTINGS", raising=False)

    first = build_sra_environment(base_environment={"PATH": "/usr/bin", "KEEP": "yes"})
    second = build_sra_environment(base_environment={"PATH": "/usr/bin", "KEEP": "yes"})

    settings = ensure_ncbi_settings(data_root)
    assert first == second
    assert first["KEEP"] == "yes"
    assert Path(first["NCBI_SETTINGS"]) == settings
    assert Path(first["TMPDIR"]).is_relative_to(data_root)
    assert Path(first["VDB_CONFIG"]).is_relative_to(data_root)
    assert str(data_root / ".sra-cache") in settings.read_text(encoding="utf-8")


def test_explicit_ncbi_settings_override_is_respected(tmp_path: Path, monkeypatch) -> None:
    override = tmp_path / "user-settings.mkfg"
    override.write_text('/repository/user/main/public/root = "/custom/cache"\n', encoding="utf-8")
    monkeypatch.setenv("AMALGKIT_NCBI_SETTINGS", str(override))

    environment = build_sra_environment(tmp_path / "campaign", base_environment={})

    assert environment["NCBI_SETTINGS"] == str(override)

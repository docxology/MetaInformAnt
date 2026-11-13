"""Fresh-host preparation for frozen-cohort quantification.

Real files under ``tmp_path``: a cohort installed by the production
``install_cohort`` (with its index bytes), the species YAML whose bytes the
cohort pins, and the production reference-manifest writer. Acceptance is the
same ``verify_species_reference`` that cohort discovery and imports use.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metainformant.rna.engine.cohort import install_cohort, load_frozen_cohort, verify_species_reference
from metainformant.rna.engine.cohort_host import EXIT_BLOCKED, EXIT_USAGE, main
from tests.rna.frozen_cohort_support import (
    COHORT_ID,
    CONFIG_BYTES,
    CONFIG_NAME,
    INDEX_BYTES,
    INDEX_NAME,
    SPECIES,
    build_frozen_source,
)


def _tree(root: Path) -> dict[str, bytes]:
    return {str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()}


class Host:
    def __init__(self, base: Path) -> None:
        build_frozen_source(base / "frozen")
        self.cohort = base / "cohort"
        install_cohort(base / "frozen", self.cohort, cohort_id=COHORT_ID)
        self.config_dir = base / "config"
        self.config_dir.mkdir()
        (self.config_dir / CONFIG_NAME).write_bytes(CONFIG_BYTES)
        self.root = base / "fresh-root"
        self.index = self.root / SPECIES / "work" / "index" / INDEX_NAME

    def args(self, *extra: str) -> list[str]:
        return [
            "--data-root",
            str(self.root),
            "--config-dir",
            str(self.config_dir),
            "--cohort",
            str(self.cohort),
            *extra,
        ]


@pytest.fixture
def host(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Host:
    monkeypatch.delenv("AMALGKIT_DATA_ROOT", raising=False)
    return Host(tmp_path)


def test_plan_on_an_empty_root_writes_nothing(host: Host, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(host.args()) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["mode"] == "plan"
    assert report["species"][0]["actions"] == ["copy_index", "write_reference_manifest"]
    assert report["species"][0]["state"] == "to_prepare"
    assert not host.root.exists()


def test_execute_places_cohort_index_bytes_and_a_verified_manifest_idempotently(
    host: Host, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    assert main(host.args("--execute")) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["states"] == {"ready": 1}
    assert host.index.read_bytes() == INDEX_BYTES
    manifest = json.loads((host.root / SPECIES / "work" / "reference" / "reference_aliases.json").read_text())
    assert manifest["status"] == "complete" and manifest["species"] == SPECIES
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(host.root))
    check = verify_species_reference(load_frozen_cohort(host.cohort), host.root, SPECIES, index_dir=host.index.parent)
    assert check.ok, check.reason
    monkeypatch.delenv("AMALGKIT_DATA_ROOT")

    before = _tree(host.root)
    assert main(host.args("--execute")) == 0
    again = json.loads(capsys.readouterr().out)
    assert again["species"][0]["actions"] == [] and again["states"] == {"ready": 1}
    assert _tree(host.root) == before, "a prepared root is left byte-identical"
    assert main(host.args()) == 0
    assert json.loads(capsys.readouterr().out)["species"][0]["state"] == "ready"


def test_different_index_bytes_are_never_replaced(host: Host, capsys: pytest.CaptureFixture[str]) -> None:
    host.index.parent.mkdir(parents=True)
    host.index.write_bytes(b"a rebuilt index is not the cohort index\n")
    assert main(host.args("--execute")) == EXIT_BLOCKED
    record = json.loads(capsys.readouterr().out)["species"][0]
    assert record["state"] == "blocked" and "never replaced" in record["reason"]
    assert host.index.read_bytes() == b"a rebuilt index is not the cohort index\n"
    assert not (host.root / SPECIES / "work" / "reference").exists()


def test_configuration_drift_and_unknown_species_are_usage_errors(
    host: Host, capsys: pytest.CaptureFixture[str]
) -> None:
    (host.config_dir / CONFIG_NAME).write_bytes(CONFIG_BYTES + b"threads: 2\n")
    assert main(host.args("--execute")) == EXIT_USAGE
    assert "configuration digest" in capsys.readouterr().out
    (host.config_dir / CONFIG_NAME).write_bytes(CONFIG_BYTES)
    assert main(host.args("--species", "no_such_species")) == EXIT_USAGE
    assert not host.root.exists()

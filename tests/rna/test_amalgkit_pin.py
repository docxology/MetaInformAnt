"""The Amalgkit pin is declared once, verified by commit, and never installed at runtime."""

from __future__ import annotations

import importlib.metadata
import logging
import re
import tomllib
from pathlib import Path

import pytest

from metainformant.rna import amalgkit
from metainformant.rna.amalgkit import _amalgkit_impl
from metainformant.rna.amalgkit.amalgkit import InstalledAmalgkit, verify_installed_amalgkit
from metainformant.rna.core.deps import check_step_dependencies, ensure_required_amalgkit
from metainformant.rna.engine.workflow_core import AmalgkitWorkflowConfig, load_workflow_config

REPO_ROOT = Path(__file__).resolve().parents[2]
AMALGKIT_REQUIREMENT = re.compile(r"^amalgkit @ git\+(?P<url>https://\S+?)@(?P<rev>[0-9a-f]{40})$")


def _pinned_build(**changes: object) -> InstalledAmalgkit:
    fields: dict[str, object] = {
        "version": amalgkit.REQUIRED_AMALGKIT_VERSION,
        "url": amalgkit.AMALGKIT_REPOSITORY_URL,
        "commit_id": amalgkit.AMALGKIT_SOURCE_REVISION,
        "requested_revision": amalgkit.AMALGKIT_SOURCE_REVISION,
        "editable": False,
    }
    fields.update(changes)
    return InstalledAmalgkit(**fields)  # type: ignore[arg-type]


def _distribution_fingerprint(name: str) -> tuple[str, str | None, int]:
    """Version, recorded origin and METADATA mtime: any reinstall changes at least one."""

    distribution = importlib.metadata.distribution(name)
    metadata = next(Path(str(distribution.locate_file(f))) for f in distribution.files or [] if f.name == "METADATA")
    return distribution.version, distribution.read_text("direct_url.json"), metadata.stat().st_mtime_ns


def test_pyproject_extras_pin_the_runtime_revision() -> None:
    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    extras = pyproject["project"]["optional-dependencies"]
    pins = {}
    for extra in ("rna", "all"):
        (requirement,) = [entry for entry in extras[extra] if entry.startswith("amalgkit")]
        match = AMALGKIT_REQUIREMENT.match(requirement)
        assert match is not None, requirement
        pins[extra] = (match["url"], match["rev"])
    assert pins["rna"] == pins["all"] == (amalgkit.AMALGKIT_REPOSITORY_URL, amalgkit.AMALGKIT_SOURCE_REVISION)


def test_lock_resolves_the_runtime_revision() -> None:
    lock = tomllib.loads((REPO_ROOT / "uv.lock").read_text(encoding="utf-8"))
    (entry,) = [package for package in lock["package"] if package["name"] == "amalgkit"]
    assert entry["version"] == amalgkit.REQUIRED_AMALGKIT_VERSION
    source = entry["source"]["git"]
    assert source.startswith(amalgkit.AMALGKIT_REPOSITORY_URL + "?rev=")
    assert source.endswith("#" + amalgkit.AMALGKIT_SOURCE_REVISION)


def test_derived_constants_name_the_commit_not_a_tag() -> None:
    revision = amalgkit.AMALGKIT_SOURCE_REVISION
    assert re.fullmatch(r"[0-9a-f]{40}", revision)
    assert amalgkit.AMALGKIT_INSTALL_SPEC == f"git+{amalgkit.AMALGKIT_REPOSITORY_URL}@{revision}"
    assert amalgkit.AMALGKIT_RELEASE_TAG == f"master@{revision[:12]}"
    # Upstream does not tag patch builds; a v<version> tag reference would not resolve.
    assert f"v{amalgkit.REQUIRED_AMALGKIT_VERSION}" not in amalgkit.AMALGKIT_INSTALL_SPEC
    assert amalgkit.AMALGKIT_SYNC_HINT == "uv sync --frozen --extra rna"


def test_verification_accepts_only_the_pinned_commit() -> None:
    ok, message = verify_installed_amalgkit(_pinned_build())
    assert ok, message
    ok, message = verify_installed_amalgkit(_pinned_build(url=amalgkit.AMALGKIT_REPOSITORY_URL + ".git/"))
    assert ok, message

    rejected = {
        "another commit": _pinned_build(commit_id="0" * 40),
        "same version string from an unidentified build": _pinned_build(commit_id=None, url=None),
        "another version": _pinned_build(version="0.16.60"),
        "a fork": _pinned_build(url="https://github.com/docxology/amalgkit"),
        "an editable checkout": _pinned_build(editable=True),
    }
    for label, build in rejected.items():
        ok, message = verify_installed_amalgkit(build)
        assert not ok, label
        assert amalgkit.AMALGKIT_SYNC_HINT in message, label


def test_installed_distribution_records_the_pinned_commit() -> None:
    pytest.importorskip("amalgkit")
    installed = amalgkit.installed_amalgkit_distribution()
    assert installed is not None
    assert (installed.version, installed.commit_id) == (
        amalgkit.REQUIRED_AMALGKIT_VERSION,
        amalgkit.AMALGKIT_SOURCE_REVISION,
    )
    assert verify_installed_amalgkit() == verify_installed_amalgkit(installed)
    assert verify_installed_amalgkit()[0]


def test_version_mismatch_fails_closed_without_touching_the_environment(caplog: pytest.LogCaptureFixture) -> None:
    pytest.importorskip("amalgkit")
    before = {name: _distribution_fingerprint(name) for name in ("amalgkit", "pandas")}

    with caplog.at_level(logging.INFO):
        ok, message, info = amalgkit.ensure_cli_available(auto_install=True, min_version="9999.0.0")

    assert not ok
    assert amalgkit.AMALGKIT_SYNC_HINT in message
    assert info is not None and info["attempted"] is False and info["valid"] is False
    assert {name: _distribution_fingerprint(name) for name in ("amalgkit", "pandas")} == before
    assert not [record for record in caplog.records if "pip install" in record.getMessage()]


def test_runtime_code_contains_no_package_installation() -> None:
    sources = [
        Path(_amalgkit_impl.__file__),
        REPO_ROOT / "src" / "metainformant" / "rna" / "core" / "deps.py",
        REPO_ROOT / "src" / "metainformant" / "rna" / "engine" / "workflow_execution.py",
    ]
    for source in sources:
        text = source.read_text(encoding="utf-8")
        assert '"pip", "install"' not in text, source
        assert "uv pip install" not in text, source


def test_pinned_requirement_is_reported_ready() -> None:
    pytest.importorskip("amalgkit")
    ok, message = ensure_required_amalgkit()
    assert ok, message


def test_auto_install_key_is_accepted_and_ignored_with_one_warning(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    assert AmalgkitWorkflowConfig(work_dir=tmp_path).auto_install_amalgkit is False

    config_path = tmp_path / "species.yaml"
    config_path.write_text(
        f"work_dir: {tmp_path / 'work'}\nthreads: 1\nauto_install_amalgkit: true\nspecies_list:\n  - Test_species\n",
        encoding="utf-8",
    )
    config = load_workflow_config(config_path)
    assert config.auto_install_amalgkit is True

    monkeypatch.setattr(_amalgkit_impl, "_auto_install_warning_issued", False)
    with caplog.at_level(logging.WARNING):
        first = check_step_dependencies("merge", {}, config)
        second = check_step_dependencies("sanity", {}, config)
    warnings = [record for record in caplog.records if "auto_install_amalgkit is ignored" in record.getMessage()]
    assert len(warnings) == 1
    assert amalgkit.AMALGKIT_SYNC_HINT in warnings[0].getMessage()
    # The outcome is the same as with the key absent: availability alone decides.
    absent = AmalgkitWorkflowConfig(work_dir=tmp_path / "work", threads=1, species_list=["Test_species"])
    assert first[0] == second[0] == check_step_dependencies("merge", {}, absent)[0]


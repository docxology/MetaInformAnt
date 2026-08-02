"""Tests for the real-implementation policy scanner."""

from __future__ import annotations

from pathlib import Path

from scripts.quality.verify_real_implementation_policy import (
    iter_policy_files,
    scan_repo,
)


def test_repo_real_implementation_policy_scan_passes() -> None:
    """Repository text should avoid old policy names and test-double APIs."""
    assert scan_repo() == []


def test_policy_file_discovery_prunes_submodules_and_symlinked_directories(
    tmp_path: Path,
) -> None:
    """Independent or external trees must never enter the parent scan."""

    root = tmp_path / "repo"
    root.mkdir()
    (root / "included.py").write_text("print('included')\n", encoding="utf-8")

    submodule = root / "projects" / "submodule"
    submodule.mkdir(parents=True)
    (submodule / ".git").write_text("gitdir: elsewhere\n", encoding="utf-8")
    (submodule / "excluded.py").write_text("Mock()\n", encoding="utf-8")

    external = tmp_path / "external"
    external.mkdir()
    (external / "excluded.py").write_text("Mock()\n", encoding="utf-8")
    (root / "external-link").symlink_to(external, target_is_directory=True)

    discovered = {
        path.relative_to(root) for path in iter_policy_files(root)
    }

    assert discovered == {Path("included.py")}

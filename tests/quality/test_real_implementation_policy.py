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

    discovered = {path.relative_to(root) for path in iter_policy_files(root)}

    assert discovered == {Path("included.py")}


def test_transient_project_state_reports_are_out_of_scan_scope(tmp_path: Path) -> None:
    """Untracked lane/status reports must never gate the policy scan.

    Parallel campaign lanes drop PROJECT_STATE_REPORT_*.md files at the repo
    root that legitimately discuss the policy by name. They are session-local
    artifacts, not tracked source, so the scanner must skip them.
    """

    root = tmp_path / "repo"
    root.mkdir()
    report = root / "PROJECT_STATE_REPORT_2026-09-01_R4_T5.md"
    report.write_text("real-implementation policy reference here\nMock(\n", encoding="utf-8")
    tracked = root / "src.py"
    tracked.write_text("x = 1\n", encoding="utf-8")

    discovered = {path.relative_to(root) for path in iter_policy_files(root)}

    assert discovered == {Path("src.py")}
    assert scan_repo(root) == []


def test_untracked_files_excluded_from_scan_scope(tmp_path: Path) -> None:
    """Untracked worktree files are not repository state for this scan.

    A tracked file with a banned phrase fails; the same content untracked is
    out of scope (it belongs to its authoring lane until committed).
    """
    import subprocess as sp

    root = tmp_path / "repo"
    root.mkdir()
    env_banned = root / "banned.py"
    env_banned.write_text("x = 1  # NO_MOCKING_POLICY\n", encoding="utf-8")
    sp.run(["git", "init", "-q", str(root)], check=True)
    sp.run(["git", "-C", str(root), "add", "banned.py"], check=True)

    assert scan_repo(root) != []

    # Now make the file untracked again - it leaves the scan scope.
    sp.run(["git", "-C", str(root), "rm", "--cached", "-q", "banned.py"], check=True)
    assert env_banned.exists()
    assert scan_repo(root) == []

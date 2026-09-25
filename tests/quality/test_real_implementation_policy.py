"""Tests for the real-implementation policy scanner."""

from __future__ import annotations

import subprocess as sp
from functools import lru_cache
from pathlib import Path

from scripts.quality.verify_real_implementation_policy import (
    MONKEYPATCH_SETATTR_PATTERNS,
    iter_policy_files,
    scan_repo,
)

# Measured 2026-09-24 with the ``monkeypatch-setattr`` rule:
# 110 rebinding sites across 21 tracked files. These are pre-existing
# ``monkeypatch.setattr`` calls that replace functions, methods, lambdas or
# module attributes with test doubles -- a documented violation of
# tests/REAL_IMPLEMENTATION_TESTING_POLICY.md ("pytest.MonkeyPatch for
# function replacement" is prohibited). Mass-migrating them to
# real-implementation tests is deliberate follow-up work, NOT part of the
# scanner change; until then the scanner must keep REPORTING the debt, and
# this assertion pins its size so the debt can neither grow nor silently
# shrink. When a site is migrated to a real implementation, re-measure with
# ``python scripts/quality/verify_real_implementation_policy.py`` and lower
# this number; any other change is a regression.
KNOWN_MONKEYPATCH_SETATTR_SITES = 110


@lru_cache(maxsize=1)
def _repo_violations() -> tuple:
    """Run the I/O-heavy repo-wide scan once per session."""
    return tuple(scan_repo())


def test_repo_real_implementation_policy_scan_passes() -> None:
    """No old policy names and no banned test-double APIs may exist.

    The ``monkeypatch-setattr`` rule intentionally still reports known debt
    (see KNOWN_MONKEYPATCH_SETATTR_SITES); every other rule must be clean.
    """
    remaining = [v for v in _repo_violations() if v.rule != "monkeypatch-setattr"]

    assert remaining == []


def test_monkeypatch_setattr_debt_is_reported_and_bounded() -> None:
    """The scanner must keep reporting the known monkeypatch.setattr debt.

    110 sites across 21 files rebind callables/attributes via
    ``monkeypatch.setattr``. Migration to real implementations is tracked
    follow-up work; the scanner reports them so the repository can see the
    debt, and this test documents the exact allowlist size and its rationale
    (mass migration is out of scope for the scanner change).
    """
    debt = [v for v in _repo_violations() if v.rule == "monkeypatch-setattr"]

    assert debt, "the documented monkeypatch.setattr debt must still be reported"
    assert len(debt) == KNOWN_MONKEYPATCH_SETATTR_SITES, (
        "the monkeypatch.setattr debt changed; re-measure with "
        "scripts/quality/verify_real_implementation_policy.py and update "
        "KNOWN_MONKEYPATCH_SETATTR_SITES with the new count and date"
    )


def test_monkeypatch_setattr_pattern_matches_rebinding_forms() -> None:
    """Call forms of ``setattr`` match; env helpers and non-call refs do not."""
    pattern = MONKEYPATCH_SETATTR_PATTERNS[0]

    assert pattern.search("monkeypatch.setattr(module, 'func', fake)")
    assert pattern.search("monkeypatch.setattr('pkg.mod.attr', fake)")
    assert pattern.search("monkeypatch.setattr(obj, 'method', stub)")
    # A call whose opening paren terminates the line is still matched
    # (line-based scanning); only a reference without a call paren is not.
    assert pattern.search("monkeypatch.setattr(\n    obj, 'attr', value,\n)")
    assert pattern.search("monkeypatch.setattr(os, 'environ', mapping)")
    # Whitelisted environment helpers are never matched.
    assert not pattern.search("monkeypatch.setenv('A', '1')")
    assert not pattern.search("monkeypatch.delenv('A', raising=False)")
    assert not pattern.search("monkeypatch.chdir(tmp_path)")


def test_monkeypatch_setattr_rule_flags_rebinding_but_not_env_config(
    tmp_path: Path,
) -> None:
    """``monkeypatch.setattr`` is a test double; ``setenv``/``delenv`` are not.

    Environment configuration (``monkeypatch.setenv`` / ``delenv``) is allowed
    real configuration under the policy and must never be reported by the
    monkeypatch-setattr rule.
    """

    root = tmp_path / "repo"
    root.mkdir()

    rebinding = root / "rebinding.py"
    rebinding.write_text(
        "def test_rebinds_callable(monkeypatch):\n"
        "    monkeypatch.setattr(mod, 'fetch_data', lambda: [])\n",
        encoding="utf-8",
    )
    env_config = root / "env_config.py"
    env_config.write_text(
        "def test_sets_env(monkeypatch):\n"
        "    monkeypatch.setenv('METAINFORMANT_HOME', '/tmp/x')\n"
        "    monkeypatch.delenv('METAINFORMANT_DEBUG', raising=False)\n",
        encoding="utf-8",
    )

    sp.run(["git", "init", "-q", str(root)], check=True)
    sp.run(["git", "-C", str(root), "add", "rebinding.py", "env_config.py"], check=True)

    violations = scan_repo(root)
    rebinding_hits = [v for v in violations if v.path == Path("rebinding.py")]
    env_hits = [v for v in violations if v.path == Path("env_config.py")]

    assert any(v.rule == "monkeypatch-setattr" for v in rebinding_hits)
    assert not any(v.rule == "monkeypatch-setattr" for v in env_hits)
    assert not any(v.rule == "test-double-api" for v in env_hits)


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
    report.write_text(
        "real-implementation policy reference here\nMock(\n", encoding="utf-8"
    )
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

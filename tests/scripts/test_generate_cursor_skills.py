"""Tests for scripts/package/generate_cursor_skills.py slug assignment.

Zero-mocks: pure functions exercised on real paths; the generated-tree check
inspects the live .cursor/skills tree without stubbing.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "package"))

import generate_cursor_skills as gcs  # noqa: E402


def test_short_path_gets_readable_slug() -> None:
    assert gcs.skill_slug_for_rel(Path("src/metainformant/rna")) == ("metainformant-src-metainformant-rna")


def test_root_slug() -> None:
    assert gcs.skill_slug_for_rel(Path(".")) == "metainformant-root"


def test_long_path_compressed_to_readable_slug() -> None:
    slug = gcs.skill_slug_for_rel(Path("src/metainformant/structural_variants/visualization"))
    assert slug == "metainformant-structural-variants-visualization"
    assert len(slug) <= 64
    assert re.match(r"^[a-z0-9-]+$", slug)


def test_prefix_segment_not_doubled() -> None:
    slug = gcs.skill_slug_for_rel(Path("src/metainformant/visualization/interactive_dashboards"))
    assert "metainformant-metainformant" not in slug
    assert slug == "metainformant-visualization-interactive-dashboards"


def test_extremely_long_single_segment_falls_back_to_digest() -> None:
    long_name = "x" * 80
    slug = gcs.skill_slug_for_rel(Path(f"a/b/{long_name}"))
    assert len(slug) <= 64
    assert slug.startswith("metainformant-")


def test_all_live_assignments_readable_and_unique() -> None:
    files = gcs.iter_agents_files(gcs.REPO_ROOT)
    assert len(files) > 400
    assignments = gcs.assign_slugs(files, gcs.REPO_ROOT)
    slugs = list(assignments.values())
    assert len(slugs) == len(set(slugs)), "slug collision"
    for slug in slugs:
        assert len(slug) <= 64
        assert re.match(r"^[a-z0-9-]+$", slug)


def test_no_hash_named_live_skill_dirs() -> None:
    """Regression: the 8 formerly sha-named wrapper dirs must stay readable."""
    skills = gcs.SKILLS_ROOT
    if not skills.is_dir():
        return
    hashed = [d.name for d in skills.iterdir() if re.search(r"[0-9a-f]{40}", d.name)]
    assert hashed == [], f"hash-named skill dirs reappeared: {hashed}"


def test_check_mode_passes_on_live_tree() -> None:
    assert gcs.run_check(gcs.REPO_ROOT, allow_uninitialized_submodules=True) == 0


def test_docstring_prose_backticks_do_not_fail_validation() -> None:
    """Regression: docstring prose may backtick-quote function names."""
    module_dir = gcs.REPO_ROOT / "src" / "metainformant" / "popgen"
    if not module_dir.is_dir():
        return
    prose_body = (
        "## Module surface (generated, validated)\n"
        "Purpose: exposed via :func:`analyze_dataset`, with scripts as orchestrator.\n"
        "- Public submodules: `workflow`."
    )
    assert gcs.validate_module_skill(module_dir, prose_body) == []


@pytest.fixture()
def scoped_repository(tmp_path: Path) -> Path:
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    for name in ("AGENTS.md", "CLAUDE.md", "docs/REAL_IMPLEMENTATION_POLICY.md"):
        path = repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# Local guidance\n")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(repo),
            "-c",
            "user.name=Test",
            "-c",
            "user.email=test@example.invalid",
            "commit",
            "-qm",
            "fixture",
        ],
        check=True,
    )
    commit = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], check=True, capture_output=True, text=True
    ).stdout.strip()
    subprocess.run(
        ["git", "-C", str(repo), "update-index", "--add", "--cacheinfo", f"160000,{commit},projects/private"],
        check=True,
    )
    for guide in gcs.iter_agents_files(repo):
        gcs.write_skill(guide, repo, gcs.skill_slug_for_rel(guide.parent.relative_to(repo)), False)
    target = repo / "projects/private/AGENTS.md"
    slug = gcs.skill_slug_for_rel(target.parent.relative_to(repo))
    skill = repo / ".cursor/skills" / slug / "SKILL.md"
    skill.parent.mkdir(parents=True)
    skill.write_text(gcs.render_skill_content(target, repo, slug))
    return repo


def test_strict_scope_refuses_missing_submodule_while_explicit_scope_reports_deferral(
    scoped_repository: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert gcs.run_check(scoped_repository) == 1
    assert gcs.run_check(scoped_repository, allow_uninitialized_submodules=True) == 0
    assert "1 wrappers deferred for 1 uninitialized submodules" in capsys.readouterr().out


def test_scoped_check_still_rejects_a_genuine_orphan(scoped_repository: Path) -> None:
    target = scoped_repository / "docs/unregistered/AGENTS.md"
    slug = gcs.skill_slug_for_rel(target.parent.relative_to(scoped_repository))
    skill = scoped_repository / ".cursor/skills" / slug / "SKILL.md"
    skill.parent.mkdir()
    skill.write_text(gcs.render_skill_content(target, scoped_repository, slug))
    assert gcs.run_check(scoped_repository, allow_uninitialized_submodules=True) == 1


def test_deferred_wrapper_identity_is_checked(scoped_repository: Path) -> None:
    skill = scoped_repository / ".cursor/skills/metainformant-projects-private/SKILL.md"
    skill.write_text(skill.read_text().replace("name: metainformant-projects-private", "name: wrong"))
    assert gcs.run_check(scoped_repository, allow_uninitialized_submodules=True) == 1


def test_partial_submodule_is_not_silently_deferred(scoped_repository: Path) -> None:
    project = scoped_repository / "projects/private"
    project.mkdir(parents=True)
    (project / "README.md").write_text("# Populated project\n")
    assert gcs.run_check(scoped_repository, allow_uninitialized_submodules=True) == 1
    guide = project / "AGENTS.md"
    guide.write_text("# Project guidance\n")
    gcs.write_skill(guide, scoped_repository, "metainformant-projects-private", False)
    assert gcs.run_check(scoped_repository) == 0


def test_submodule_symlink_cannot_escape_scope(scoped_repository: Path, tmp_path: Path) -> None:
    project = scoped_repository / "projects/private"
    project.parent.mkdir(exist_ok=True)
    outside = tmp_path / "outside"
    outside.mkdir()
    project.symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match="escapes"):
        gcs.run_check(scoped_repository, allow_uninitialized_submodules=True)

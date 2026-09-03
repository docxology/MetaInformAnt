"""Tests for scripts/package/generate_cursor_skills.py slug assignment.

Zero-mocks: pure functions exercised on real paths; the generated-tree check
inspects the live .cursor/skills tree without stubbing.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "package"))

import generate_cursor_skills as gcs  # noqa: E402


def test_short_path_gets_readable_slug() -> None:
    assert gcs.skill_slug_for_rel(Path("src/metainformant/rna")) == (
        "metainformant-src-metainformant-rna"
    )


def test_root_slug() -> None:
    assert gcs.skill_slug_for_rel(Path(".")) == "metainformant-root"


def test_long_path_compressed_to_readable_slug() -> None:
    slug = gcs.skill_slug_for_rel(
        Path("src/metainformant/structural_variants/visualization")
    )
    assert slug == "metainformant-structural-variants-visualization"
    assert len(slug) <= 64
    assert re.match(r"^[a-z0-9-]+$", slug)


def test_prefix_segment_not_doubled() -> None:
    slug = gcs.skill_slug_for_rel(
        Path("src/metainformant/visualization/interactive_dashboards")
    )
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
    assert gcs.run_check(gcs.REPO_ROOT) == 0


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

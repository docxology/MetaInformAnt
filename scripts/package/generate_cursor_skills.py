#!/usr/bin/env python3
"""Generate Cursor Agent Skills under .cursor/skills/ from AGENTS.md locations.

Usage:
  uv run python scripts/package/generate_cursor_skills.py          # write/update SKILL.md files
  uv run python scripts/package/generate_cursor_skills.py --check  # verify only (exit 1 on drift)
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path
from typing import Iterable

# Repo root = parent of scripts/
REPO_ROOT = Path(__file__).resolve().parents[2]
SKILLS_ROOT = REPO_ROOT / ".cursor" / "skills"

SKIP_DIR_NAMES = frozenset(
    {
        ".git",
        ".venv",
        "node_modules",
        "__pycache__",
        "output",
        "dist",
        "build",
        ".tmp",
        ".uv-cache",
        ".mypy_cache",
        ".pytest_cache",
        ".ruff_cache",
        "htmlcov",
        ".eggs",
    }
)


def should_skip_agents_path(path: Path, repo: Path) -> bool:
    try:
        rel = path.relative_to(repo)
    except ValueError:
        return True
    for part in rel.parts:
        if part in SKIP_DIR_NAMES:
            return True
        if part.endswith(".egg-info"):
            return True
    return False


def iter_agents_files(repo: Path) -> list[Path]:
    found: list[Path] = []
    for p in repo.rglob("AGENTS.md"):
        if p.is_file() and not should_skip_agents_path(p, repo):
            found.append(p)
    return sorted(found, key=lambda x: str(x))


def _sanitize_segment(part: str) -> str:
    s = part.lower().replace("_", "-")
    s = re.sub(r"[^a-z0-9-]+", "-", s)
    s = re.sub(r"-+", "-", s).strip("-")
    return s or "x"


def folder_relative_to_repo(agents_file: Path, repo: Path) -> Path:
    return agents_file.parent.resolve().relative_to(repo.resolve())


def display_folder_label(folder: Path, repo: Path) -> str:
    rel = folder.resolve().relative_to(repo.resolve())
    if rel == Path(".") or not rel.parts:
        return "repository root"
    return rel.as_posix()


def skill_slug_for_rel(rel: Path) -> str:
    if rel == Path("."):
        return "metainformant-root"
    segs = [_sanitize_segment(p) for p in rel.parts]
    base = "metainformant-" + "-".join(segs)
    if len(base) <= 64 and re.match(r"^[a-z0-9-]+$", base):
        return base
    digest = hashlib.sha256(str(rel).encode()).hexdigest()[:50]
    return f"metainformant-{digest}"


def assign_slugs(agents_files: Iterable[Path], repo: Path) -> dict[Path, str]:
    """Map each AGENTS.md path to a unique skill slug."""
    slug_assignments: dict[Path, str] = {}
    used_slugs: dict[str, Path] = {}

    for agents in agents_files:
        rel = folder_relative_to_repo(agents, repo)
        rel_s = str(rel)
        slug = skill_slug_for_rel(rel)
        if slug in used_slugs and used_slugs[slug] != agents:
            digest = hashlib.sha256(rel_s.encode()).hexdigest()[:10]
            alt = f"metainformant-{digest}"
            if len(alt) > 64:
                alt = alt[:64]
            slug = alt
        while slug in used_slugs and used_slugs[slug] != agents:
            digest = hashlib.sha256(f"{rel_s}-{slug}".encode()).hexdigest()[:8]
            slug = (f"metainformant-{digest}")[:64]
        used_slugs[slug] = agents
        slug_assignments[agents] = slug

    return slug_assignments


def ups_to_repo(skill_md: Path, repo: Path) -> str:
    depth = len(skill_md.parent.resolve().relative_to(repo.resolve()).parts)
    return "../" * depth


def relative_link_from_skill_to_repo_file(skill_md: Path, target: Path, repo: Path) -> str:
    """Markdown link path from SKILL.md to a repo file."""
    up = ups_to_repo(skill_md, repo)
    rel_target = target.resolve().relative_to(repo.resolve()).as_posix()
    return f"{up}{rel_target}"


def build_skill_body(agents_file: Path, repo: Path, skill_md: Path) -> str:
    folder = agents_file.parent
    rel_folder = display_folder_label(folder, repo)
    up = ups_to_repo(skill_md, repo)
    agents_link = relative_link_from_skill_to_repo_file(skill_md, agents_file, repo)
    lines = [
        f"# METAINFORMANT — `{rel_folder}`",
        "",
        "Before editing files in this subtree:",
        "",
        f"- Read [`AGENTS.md`]({agents_link}) for this folder (canonical technical context).",
    ]
    readme = folder / "README.md"
    if readme.is_file():
        rl = relative_link_from_skill_to_repo_file(skill_md, readme, repo)
        lines.append(f"- Optional overview: [`README.md`]({rl}).")
    lines.extend(
        [
            f"- Global rules: [`CLAUDE.md`]({up}CLAUDE.md) at repo root (uv, `output/`, `.tmp/`, no mocks).",
            f"- Testing policy: [`docs/NO_MOCKING_POLICY.md`]({up}docs/NO_MOCKING_POLICY.md).",
            "- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.",
            "",
            "Keep changes scoped; match existing patterns in this directory.",
        ]
    )
    return "\n".join(lines) + "\n"


def build_description(rel_folder: str) -> str:
    where = "the repository root" if rel_folder == "repository root" else f"directory {rel_folder}"
    text = (
        f"METAINFORMANT rules for {where}. "
        f"Use when editing, adding tests, or reviewing code under this path. "
        f"Read the linked AGENTS.md first; use uv only, write outputs to output/, no mocks."
    )
    if len(text) > 1024:
        return text[:1021] + "..."
    return text


def write_skill(agents_file: Path, repo: Path, slug: str, dry_run: bool) -> Path:
    skill_dir = SKILLS_ROOT / slug
    skill_md = skill_dir / "SKILL.md"
    rel_folder = display_folder_label(agents_file.parent, repo)
    description = build_description(rel_folder)
    body = build_skill_body(agents_file, repo, skill_md)
    name = slug
    content = f"---\nname: {name}\ndescription: {description}\n---\n\n{body}"
    if not dry_run:
        skill_dir.mkdir(parents=True, exist_ok=True)
        skill_md.write_text(content, encoding="utf-8")
    return skill_md


def parse_frontmatter(skill_md: Path) -> tuple[str | None, str | None]:
    text = skill_md.read_text(encoding="utf-8")
    if not text.startswith("---"):
        return None, None
    end = text.find("\n---\n", 3)
    if end == -1:
        return None, None
    fm = text[3:end]
    name_m = re.search(r"^name:\s*(.+)$", fm, re.MULTILINE)
    desc_m = re.search(r"^description:\s*(.+)$", fm, re.MULTILINE)
    name = name_m.group(1).strip() if name_m else None
    desc = desc_m.group(1).strip() if desc_m else None
    return name, desc


def run_check(repo: Path) -> int:
    agents_files = iter_agents_files(repo)
    assignments = assign_slugs(agents_files, repo)
    errors: list[str] = []

    for agents, slug in assignments.items():
        skill_md = SKILLS_ROOT / slug / "SKILL.md"
        if not skill_md.is_file():
            errors.append(f"Missing skill for {agents}: expected {skill_md}")
            continue
        name, desc = parse_frontmatter(skill_md)
        if not name or not desc:
            errors.append(f"Invalid frontmatter in {skill_md}")
            continue
        if name != slug:
            errors.append(f"name mismatch in {skill_md}: got {name!r}, expected {slug!r}")
        if len(desc) > 1024:
            errors.append(f"description too long in {skill_md} ({len(desc)} chars)")
        if not re.match(r"^[a-z0-9-]+$", name):
            errors.append(f"invalid name field in {skill_md}: {name!r}")

    # Orphan skill directories (has SKILL.md but no AGENTS mapping)
    expected_slugs = set(assignments.values())
    if SKILLS_ROOT.is_dir():
        for child in SKILLS_ROOT.iterdir():
            if child.name == "README.md" or not child.is_dir():
                continue
            sm = child / "SKILL.md"
            if sm.is_file() and child.name not in expected_slugs:
                errors.append(f"Orphan skill directory (no matching AGENTS.md): {child}")

    if errors:
        print("check failed:", file=sys.stderr)
        for e in errors:
            print(f"  {e}", file=sys.stderr)
        return 1
    print(f"check ok: {len(assignments)} skills match AGENTS.md locations")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify skills only")
    parser.add_argument("--dry-run", action="store_true", help="print actions without writing")
    args = parser.parse_args()
    repo = REPO_ROOT

    agents_files = iter_agents_files(repo)
    assignments = assign_slugs(agents_files, repo)

    if args.check:
        return run_check(repo)

    for agents, slug in assignments.items():
        write_skill(agents, repo, slug, args.dry_run)
        if args.dry_run:
            print(f"would write .cursor/skills/{slug}/SKILL.md <- {agents.relative_to(repo)}")

    print(f"wrote {len(assignments)} skills under .cursor/skills/")
    return 0


if __name__ == "__main__":
    sys.exit(main())

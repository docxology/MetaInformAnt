#!/usr/bin/env python3
"""Generate Cursor Agent Skills under .cursor/skills/ from AGENTS.md locations.

Usage:
  uv run python scripts/package/generate_cursor_skills.py          # write/update SKILL.md files
  uv run python scripts/package/generate_cursor_skills.py --check  # verify only (exit 1 on drift)
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import os
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
        "tmp",
        ".uv-cache",
        ".mypy_cache",
        ".pytest_cache",
        ".ruff_cache",
        "htmlcov",
        ".eggs",
    }
)


# Directories skipped inside NESTED git repositories only (not the top-level repo).
# Top-level repo trees keep their canonical skills; nested-repo data/output/results/
# trees (e.g. projects/hymenoptera_amalgkit/data/** per-sample quant dirs from the
# live campaign) carry runtime artifacts, not documentation, and must not be walked.
NESTED_REPO_SKIP_DIR_NAMES = frozenset({"data", "output", "results", "logs", ".downloads"})


def is_nested_git_repo(directory: Path, repo: Path) -> bool:
    """True if ``directory`` is a git repository distinct from the top-level ``repo``."""
    if not (directory / ".git").exists():
        return False
    try:
        return directory.resolve() != repo.resolve()
    except (OSError, ValueError):
        return False


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
    # Prune excluded trees while walking.  A post-filtered rglob still visits
    # every file under preserved runtime artifacts such as ``tmp/``.
    for directory, dirnames, filenames in os.walk(repo, followlinks=False):
        skip = set(SKIP_DIR_NAMES)
        if is_nested_git_repo(Path(directory), repo):
            # Nested repo: also skip its data/output/results/logs trees.
            skip |= NESTED_REPO_SKIP_DIR_NAMES
        dirnames[:] = sorted(name for name in dirnames if name not in skip)
        if "AGENTS.md" in filenames:
            path = Path(directory) / "AGENTS.md"
            if not should_skip_agents_path(path, repo):
                found.append(path)
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


MODULES_SRC = Path("src") / "metainformant"


def module_name_for_rel(rel: Path) -> str | None:
    """Return the metainformant module name if ``rel`` is src/metainformant/<module>."""
    parts = rel.parts
    if len(parts) == 3 and parts[0] == "src" and parts[1] == "metainformant":
        return parts[2]
    return None


def _parse_package_init(init_py: Path) -> tuple[str | None, list[str]]:
    """AST-parse a package __init__.py: (docstring, exported submodule names).

    Static analysis keeps the check gate import-free and deterministic: heavy
    scientific dependencies are never executed just to validate documentation.
    """
    if not init_py.is_file():
        return None, []
    try:
        tree = ast.parse(init_py.read_text(encoding="utf-8"))
    except SyntaxError:
        return None, []
    doc = ast.get_docstring(tree)
    names: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "__all__":
                    if isinstance(node.value, (ast.List, ast.Tuple)):
                        for elt in node.value.elts:
                            if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                                names.append(elt.value)
        elif isinstance(node, ast.ImportFrom) and node.level == 1 and node.module is None:
            for alias in node.names:
                if alias.name != "*":
                    names.append(alias.name)
    seen: set[str] = set()
    ordered = [n for n in names if not (n in seen or seen.add(n))]
    return doc, ordered


def _submodule_exists(module_dir: Path, name: str) -> bool:
    return (module_dir / name / "__init__.py").is_file() or (module_dir / f"{name}.py").is_file()


def build_module_sections(module_dir: Path) -> list[str]:
    """Enriched per-module skill sections validated against the real package."""
    doc, names = _parse_package_init(module_dir / "__init__.py")
    module = module_dir.name
    lines: list[str] = []
    if doc:
        lines.append(f"Purpose: {doc.strip()}")
    submods = [n for n in names if _submodule_exists(module_dir, n)]
    if submods:
        lines.append("Public submodules: " + ", ".join(f"`{n}`" for n in submods) + ".")
    lines.append(
        f"Canonical import: `import metainformant.{module}` "
        f"(submodules: `from metainformant import {module}` then `{module}.<submodule>`)."
    )
    tests_dir = REPO_ROOT / "tests" / module
    if tests_dir.is_dir():
        lines.append(f"Test entry point: `uv run pytest tests/{module} -q` " "(one pytest directory per invocation).")
    return lines


def _validate_package_init(module_dir: Path) -> list[str]:
    """Statically verify a package __init__.py's relative imports resolve.

    ``from . import x`` requires ``x`` to exist on disk (module file or
    package directory); ``from .x import y`` requires ``x`` to exist.
    Re-exported symbols (``y``) are attributes, not files, and are not
    checked here — importing the package is the runtime's job, not the
    documentation gate's.
    """
    errors: list[str] = []
    init_py = module_dir / "__init__.py"
    if not init_py.is_file():
        return errors
    try:
        tree = ast.parse(init_py.read_text(encoding="utf-8"))
    except SyntaxError as exc:
        return [f"{init_py}: syntax error in __init__.py: {exc}"]
    for node in ast.walk(tree):
        if not (isinstance(node, ast.ImportFrom) and node.level == 1):
            continue
        if node.module is None:
            for alias in node.names:
                name = alias.name
                if name != "*" and not _submodule_exists(module_dir, name):
                    errors.append(
                        f"module skill drift: {module_dir.name} __init__ does "
                        f"`from . import {name}` but '{name}' does not exist on disk"
                    )
        else:
            # Dotted path (e.g. from .analysis.association import x): every
            # prefix must resolve on disk, walking packages segment by segment.
            current = module_dir
            resolved = True
            for seg in node.module.split("."):
                pkg_init = current / seg / "__init__.py"
                mod_file = current / f"{seg}.py"
                if pkg_init.is_file():
                    current = current / seg
                elif not mod_file.is_file():
                    resolved = False
                    break
            if not resolved:
                errors.append(
                    f"module skill drift: {module_dir.name} __init__ does "
                    f"`from .{node.module} import ...` but '{node.module}' does not exist on disk"
                )
    return errors


def validate_module_skill(module_dir: Path, skill_body: str) -> list[str]:
    """Fail loudly when a module skill or its package __init__ has drift.

    Two layers:
    1. The generated skill lists only submodules verified to exist on disk,
       so any backtick-quoted submodule name in the body that is missing on
       disk is stale prose.
    2. The package __init__'s own relative imports must resolve statically.
    """
    errors: list[str] = []
    # Only the generated "Public submodules:" line makes submodule claims; the
    # Purpose line embeds the package docstring, whose prose may legitimately
    # backtick-quote function names that are not submodules.
    submodule_line = ""
    for line in skill_body.splitlines():
        if line.startswith("- Public submodules:"):
            submodule_line = line
            break
    doc, names = _parse_package_init(module_dir / "__init__.py")
    for name in names:
        if f"`{name}`" in submodule_line and not _submodule_exists(module_dir, name):
            errors.append(
                f"module skill drift: {module_dir.name} skill advertises submodule "
                f"'{name}' which no longer exists on disk"
            )
    errors.extend(_validate_package_init(module_dir))
    return errors


def _readable_fallback_slug(rel: Path) -> str | None:
    """Readable compressed slug for over-long paths, or None if it cannot fit.

    Strategy: progressively drop leading path segments (keeping the most
    specific tail segments) until the slug fits the 64-character limit.
    Returns None when even the last segment alone is too long, in which case
    the caller falls back to a documented digest slug.
    """
    segs = [_sanitize_segment(p) for p in rel.parts]
    for drop in range(len(segs)):
        kept = segs[drop:]
        # The implicit leading "metainformant" participates in dedupe so a
        # kept "metainformant" segment does not double the package prefix and
        # consecutive duplicate segments collapse (metainformant-x-x -> x).
        collapsed: list[str] = []
        for seg in ["metainformant", *kept]:
            if not collapsed or collapsed[-1] != seg:
                collapsed.append(seg)
        base = "-".join(collapsed)
        if len(base) <= 64 and re.match(r"^[a-z0-9-]+$", base):
            return base
    return None


def skill_slug_for_rel(rel: Path) -> str:
    if rel == Path("."):
        return "metainformant-root"
    segs = [_sanitize_segment(p) for p in rel.parts]
    base = "metainformant-" + "-".join(segs)
    if len(base) <= 64 and re.match(r"^[a-z0-9-]+$", base):
        return base
    # Long paths: compress to a readable slug by dropping leading segments so
    # skill names stay human-greppable (e.g. src/metainformant/<module>/...)
    # instead of opaque digests. Digest fallback only when compression cannot
    # fit; assign_slugs guarantees uniqueness for any collisions.
    compressed = _readable_fallback_slug(rel)
    if compressed is not None:
        return compressed
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
            f"- Global rules: [`CLAUDE.md`]({up}CLAUDE.md) at repo root "
            "(uv, `output/`, `.tmp/`, real implementations).",
            f"- Testing policy: [`docs/REAL_IMPLEMENTATION_POLICY.md`]({up}docs/REAL_IMPLEMENTATION_POLICY.md).",
            "- Use `metainformant.core.io` for file I/O and `metainformant.core.utils.logging` for logs.",
        ]
    )
    module_name = module_name_for_rel(folder_relative_to_repo(agents_file, repo))
    if module_name is not None:
        module_sections = build_module_sections(folder)
        if module_sections:
            lines.append("")
            lines.append("## Module surface (generated, validated)")
            for section in module_sections:
                lines.append(f"- {section}")
    lines.extend(
        [
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
        f"Read the linked AGENTS.md first; use uv only, write outputs to output/, real implementations."
    )
    if len(text) > 1024:
        return text[:1021] + "..."
    return text


def render_skill_content(agents_file: Path, repo: Path, slug: str) -> str:
    """Render the deterministic wrapper for one ``AGENTS.md`` file."""

    skill_dir = SKILLS_ROOT / slug
    skill_md = skill_dir / "SKILL.md"
    rel_folder = display_folder_label(agents_file.parent, repo)
    description = build_description(rel_folder)
    body = build_skill_body(agents_file, repo, skill_md)
    name = slug
    return f"---\nname: {name}\ndescription: {description}\n---\n\n{body}"


def write_skill(agents_file: Path, repo: Path, slug: str, dry_run: bool) -> Path:
    """Write one generated wrapper and return its path."""

    skill_dir = SKILLS_ROOT / slug
    skill_md = skill_dir / "SKILL.md"
    content = render_skill_content(agents_file, repo, slug)
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

        expected = render_skill_content(agents, repo, slug)
        actual = skill_md.read_text(encoding="utf-8")
        if actual != expected:
            errors.append(f"content drift in {skill_md}; regenerate with generate_cursor_skills.py")

        module_name = module_name_for_rel(folder_relative_to_repo(agents, repo))
        if module_name is not None:
            errors.extend(f"{skill_md}: {e}" for e in validate_module_skill(agents.parent, actual))

        required_targets = [repo / "CLAUDE.md", repo / "docs" / "REAL_IMPLEMENTATION_POLICY.md"]
        readme = agents.parent / "README.md"
        if readme.is_file():
            required_targets.append(readme)
        for target in [agents, *required_targets]:
            if not target.is_file():
                errors.append(f"broken canonical link target for {skill_md}: {target}")

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

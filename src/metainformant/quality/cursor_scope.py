"""Explicit checkout scope for generated guidance in optional Git submodules."""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path


def uninitialized_submodules(repo: Path) -> tuple[Path, ...]:
    """Identify only registered, empty Gitlink paths; populated checkouts stay in scope."""
    root = repo.resolve()
    environment = os.environ.copy()
    for name in ("GIT_DIR", "GIT_WORK_TREE", "GIT_COMMON_DIR", "GIT_INDEX_FILE"):
        environment.pop(name, None)
    environment["GIT_OPTIONAL_LOCKS"] = "0"
    result = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--stage", "-z"],
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )
    missing: list[Path] = []
    for record in result.stdout.split("\0"):
        if not record:
            continue
        header, separator, name = record.partition("\t")
        fields = header.split()
        if not separator or len(fields) != 3:
            raise ValueError("malformed Git index entry")
        if fields[0] != "160000":
            continue
        if fields[2] != "0":
            raise ValueError("unmerged submodule index entry")
        path = root / name
        if path.is_symlink() or not path.resolve().is_relative_to(root):
            raise ValueError("submodule path escapes the checked repository")
        if not path.exists() or (path.is_dir() and not any(path.iterdir())):
            missing.append(path.resolve())
    return tuple(missing)


def deferred_agents_target(skill: Path, repo: Path, missing: tuple[Path, ...]) -> Path | None:
    """Resolve the canonical link only when it belongs to an unavailable Gitlink."""
    links = re.findall(r"^- Read \[`AGENTS\.md`\]\(([^)]+)\) for this folder", skill.read_text(), re.MULTILINE)
    if len(links) != 1 or Path(links[0]).is_absolute():
        return None
    target = (skill.parent / links[0]).resolve()
    if target.name != "AGENTS.md" or not target.is_relative_to(repo.resolve()):
        return None
    return target if any(target.is_relative_to(path) for path in missing) else None

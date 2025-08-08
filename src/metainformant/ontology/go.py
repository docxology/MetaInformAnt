from __future__ import annotations

from pathlib import Path


def count_go_scripts(go_dir: Path) -> int:
    return sum(1 for p in go_dir.glob("*.py") if p.is_file())



from __future__ import annotations

import os
from pathlib import Path

from metainformant.core import paths as core_paths


def test_expand_and_resolve_user(tmp_path: Path) -> None:
    home = tmp_path / "home"
    home.mkdir(parents=True, exist_ok=True)
    old_home = os.environ.get("HOME")
    os.environ["HOME"] = str(home)
    try:
        p = core_paths.expand_and_resolve("~/data")
        assert str(p).startswith(str(home))
    finally:
        if old_home is None:
            os.environ.pop("HOME", None)
        else:
            os.environ["HOME"] = old_home


def test_is_within(tmp_path: Path) -> None:
    root = tmp_path / "root"
    sub = root / "a" / "b"
    sub.mkdir(parents=True, exist_ok=True)
    assert core_paths.is_within(sub, root)
    assert not core_paths.is_within(root, sub)

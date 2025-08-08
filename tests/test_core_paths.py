from __future__ import annotations

from pathlib import Path
from metainformant.core import paths as core_paths


def test_expand_and_resolve_user(tmp_path: Path, monkeypatch) -> None:
    home = tmp_path / "home"
    home.mkdir(parents=True, exist_ok=True)
    monkeypatch.setenv("HOME", str(home))
    p = core_paths.expand_and_resolve("~/data")
    assert str(p).startswith(str(home))


def test_is_within(tmp_path: Path) -> None:
    root = tmp_path / "root"
    sub = root / "a" / "b"
    sub.mkdir(parents=True, exist_ok=True)
    assert core_paths.is_within(sub, root)
    assert not core_paths.is_within(root, sub)



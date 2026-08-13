"""Regression tests for bounded external download execution."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from metainformant.core.io import download_robust


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"max_retries": 0}, "max_retries"),
        ({"retry_delay": -1}, "retry_delay"),
        ({"timeout": 0}, "timeout"),
    ],
)
def test_robust_download_rejects_unbounded_controls(kwargs: dict[str, int], message: str, tmp_path: Path) -> None:
    """Invalid retry/timeout controls fail before creating external work."""
    with pytest.raises(ValueError, match=message):
        download_robust.robust_download_url("https://example.invalid/file", tmp_path / "file", **kwargs)


def test_robust_download_passes_wall_clock_timeout_to_subprocess(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The configured timeout bounds the whole curl/wget subprocess."""
    monkeypatch.setattr(download_robust.shutil, "which", lambda name: "/usr/bin/curl" if name == "curl" else None)
    observed: dict[str, object] = {}

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        observed["command"] = command
        observed.update(kwargs)
        return subprocess.CompletedProcess(command, 1, stdout="", stderr="timed out")

    monkeypatch.setattr(download_robust.subprocess, "run", fake_run)
    assert not download_robust.robust_download_url(
        "https://example.invalid/file", tmp_path / "file", max_retries=1, retry_delay=0, timeout=7
    )
    assert observed["timeout"] == 7

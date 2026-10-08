"""Compression tuning preserves real FASTQ bytes and rejects invalid profiles."""

from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from metainformant.rna.engine.fastq_compression import compression_level, pigz_command
from metainformant.rna.retrieval.ena_downloader import verify_gzip_integrity


@pytest.mark.parametrize("level", [1, 6, 9])
def test_native_compression_preserves_fastq(tmp_path: Path, level: int) -> None:
    if shutil.which("pigz") is None:
        pytest.skip("pigz is required for native compression verification")
    source = tmp_path / "reads.fastq"
    payload = b"@read/1\nACGTACGT\n+\nIIIIIIII\n" * 1000
    source.write_bytes(payload)
    subprocess.run(pigz_command(source, threads=2, level=level), check=True, capture_output=True)
    target = source.with_suffix(".fastq.gz")
    assert verify_gzip_integrity(target)
    assert gzip.decompress(target.read_bytes()) == payload
    target.write_bytes(target.read_bytes()[:-8])
    assert not verify_gzip_integrity(target)


@pytest.mark.parametrize("level", [0, 10, -1, True, 1.5])
def test_invalid_level_cannot_build_a_command(level: int) -> None:
    with pytest.raises(ValueError):
        pigz_command(Path("reads.fastq"), threads=2, level=level)


@pytest.mark.parametrize("threads", [0, -1, True, 1.5])
def test_invalid_thread_budget_cannot_build_a_command(threads: int) -> None:
    with pytest.raises(ValueError):
        pigz_command(Path("reads.fastq"), threads=threads, level=1)


def test_environment_level_is_explicit_and_strict(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("AMALGKIT_PIPELINE_COMPRESSION_LEVEL", raising=False)
    assert compression_level() == 6
    monkeypatch.setenv("AMALGKIT_PIPELINE_COMPRESSION_LEVEL", "1")
    assert compression_level() == 1
    for invalid in ("", "fast", "0", "10"):
        monkeypatch.setenv("AMALGKIT_PIPELINE_COMPRESSION_LEVEL", invalid)
        with pytest.raises(ValueError):
            compression_level()


def test_controller_requires_explicit_frozen_configuration(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[2]
    command = [
        sys.executable,
        str(root / "scripts/rna/acquisition.py"),
        "aws",
        "--campaign-root",
        str(tmp_path),
        "--repo",
        str(root),
        "--bucket",
        "unused",
        "--cohort",
        "unused",
        "--budget",
        "750",
        "--historical-gross",
        "0",
        "--ami",
        "unused",
        "--instance-profile",
        "unused",
        "--worker-compression-level",
        "1",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode != 0
    assert "requires an explicit --config-dir" in result.stderr
    assert list(tmp_path.iterdir()) == []

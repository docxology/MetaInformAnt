from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.rna.engine.streaming_orchestrator import (
    _local_sra_path,
    _sra_source_fingerprint,
    _write_sra_validation_result,
)


@pytest.mark.parametrize("layout", ["sample", "legacy", "native_sample", "native_legacy"])
@pytest.mark.parametrize("suffix", [".sra", ".sra.cache"])
def test_exact_prefetched_archive_is_discovered_without_unrelated_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    layout: str,
    suffix: str,
) -> None:
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(tmp_path))
    out_dir = tmp_path / "species/work/getfastq"
    sample = out_dir / "SRR123"
    root = sample / ".prefetch/SRR123" if "sample" in layout else tmp_path / ".sra-cache/prefetch/SRR123"
    if layout.startswith("native"):
        root /= "SRR123"
    root.mkdir(parents=True)
    unrelated = root / "SRR124.sra"
    unrelated.write_bytes(b"unrelated payload")
    assert _local_sra_path(sample, out_dir, "SRR123") is None
    archive = root / f"SRR123{suffix}"
    archive.write_bytes(b"resumable complete-looking payload")
    assert _local_sra_path(sample, out_dir, "SRR123") == archive
    assert unrelated.read_bytes() == b"unrelated payload"


def test_fresh_prefetch_is_preferred_to_known_invalid_historical_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(tmp_path))
    old = tmp_path / ".sra-cache/sra/SRR123.sra"
    old.parent.mkdir(parents=True)
    old.write_bytes(b"bad historical archive")
    _write_sra_validation_result(
        old,
        fingerprint=_sra_source_fingerprint(old),
        valid=False,
        returncode=3,
        detail="corrupt historical archive",
    )
    sample = tmp_path / "species/work/getfastq/SRR123"
    fresh = sample / ".prefetch/SRR123/SRR123.sra"
    fresh.parent.mkdir(parents=True)
    fresh.write_bytes(b"fresh archive")
    assert _local_sra_path(sample, sample.parent, "SRR123") == fresh
    assert old.read_bytes() == b"bad historical archive"

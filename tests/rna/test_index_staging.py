"""The index copies kallisto actually reads: the exact-stem stage and the opt-in local cache.

Real files under ``tmp_path``; the production staging functions copy, hash,
fsync and promote them. A copy is reused only while its source and the copy
itself are unchanged, so a damaged or stale copy is never handed to kallisto.
"""

from __future__ import annotations

import hashlib
import json
import os
import threading
from pathlib import Path

import pytest

from metainformant.rna.engine.streaming_orchestrator import (
    INDEX_STAGE_MANIFEST,
    _cache_index_directory,
    _staged_exact_index_dir,
)

SPECIES = "bombus_terrestris"
EXACT = b"exact-stem kallisto index " * 32
SUBSPECIES = b"subspecies kallisto index " * 32


def _reference(tmp_path: Path) -> Path:
    """A shared reference directory whose subspecies index prefix-collides with the species index."""

    index_dir = tmp_path / "genome" / "index"
    index_dir.mkdir(parents=True)
    (index_dir / "Bombus_terrestris.idx").write_bytes(EXACT)
    (index_dir / "Bombus_terrestris_audax.idx").write_bytes(SUBSPECIES)
    return index_dir


def _no_temporaries(directory: Path) -> bool:
    return not [path.name for path in directory.iterdir() if path.name.endswith(".tmp")]


def test_the_staged_copy_is_hash_verified_and_recorded(tmp_path: Path) -> None:
    index_dir = _reference(tmp_path)
    staged = Path(_staged_exact_index_dir(str(index_dir), SPECIES))
    assert staged == Path(f"{index_dir}_exact_stem") / SPECIES / "index"
    assert sorted(path.name for path in staged.iterdir()) == sorted(["Bombus_terrestris.idx", INDEX_STAGE_MANIFEST])
    record = json.loads((staged / INDEX_STAGE_MANIFEST).read_bytes())["copies"]["Bombus_terrestris.idx"]
    assert record["sha256"] == hashlib.sha256(EXACT).hexdigest()
    assert (staged / "Bombus_terrestris.idx").read_bytes() == EXACT


def test_an_unchanged_stage_is_reused_without_copying(tmp_path: Path) -> None:
    index_dir = _reference(tmp_path)
    staged = Path(_staged_exact_index_dir(str(index_dir), SPECIES))
    identity = (staged / "Bombus_terrestris.idx").stat().st_ino
    assert Path(_staged_exact_index_dir(str(index_dir), SPECIES)) == staged
    assert (staged / "Bombus_terrestris.idx").stat().st_ino == identity


def test_a_damaged_staged_copy_of_the_same_size_is_staged_again(tmp_path: Path) -> None:
    index_dir = _reference(tmp_path)
    staged = Path(_staged_exact_index_dir(str(index_dir), SPECIES))
    copy = staged / "Bombus_terrestris.idx"
    copy.write_bytes(b"X" * len(EXACT))  # same size: the old size-only check reused it forever

    again = Path(_staged_exact_index_dir(str(index_dir), SPECIES))

    assert again == staged and copy.read_bytes() == EXACT
    assert _no_temporaries(staged)


def test_a_rebuilt_source_index_is_staged_again(tmp_path: Path) -> None:
    index_dir = _reference(tmp_path)
    staged = Path(_staged_exact_index_dir(str(index_dir), SPECIES))
    rebuilt = b"rebuilt kallisto index " * 40
    (index_dir / "Bombus_terrestris.idx").write_bytes(rebuilt)

    Path(_staged_exact_index_dir(str(index_dir), SPECIES))

    assert (staged / "Bombus_terrestris.idx").read_bytes() == rebuilt
    record = json.loads((staged / INDEX_STAGE_MANIFEST).read_bytes())["copies"]["Bombus_terrestris.idx"]
    assert record["sha256"] == hashlib.sha256(rebuilt).hexdigest()


def test_concurrent_staging_threads_never_share_a_temporary_file(tmp_path: Path) -> None:
    index_dir = _reference(tmp_path)
    results: list[str] = []
    errors: list[BaseException] = []

    def stage() -> None:
        try:
            results.append(_staged_exact_index_dir(str(index_dir), SPECIES))
        except BaseException as error:  # noqa: BLE001 - reported below
            errors.append(error)

    threads = [threading.Thread(target=stage) for _ in range(8)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    assert errors == [] and len(set(results)) == 1
    staged = Path(results[0])
    assert (staged / "Bombus_terrestris.idx").read_bytes() == EXACT
    assert _no_temporaries(staged)


def test_a_damaged_local_cache_copy_is_copied_again(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = tmp_path / "external" / "index"
    source.mkdir(parents=True)
    (source / "genome.idx").write_bytes(EXACT)
    monkeypatch.setenv("AMALGKIT_LOCAL_INDEX_CACHE_DIR", str(tmp_path / "cache"))
    cached = Path(_cache_index_directory(str(source), "test_species"))
    copy = cached / "genome.idx"
    recorded_mtime = copy.stat().st_mtime_ns
    copy.write_bytes(b"Y" * len(EXACT))
    os.utime(copy, ns=(recorded_mtime, recorded_mtime))  # same size and mtime: only the ctime tells

    assert Path(_cache_index_directory(str(source), "test_species")) == cached
    assert copy.read_bytes() == EXACT
    assert _no_temporaries(cached)

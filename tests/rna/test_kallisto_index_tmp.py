"""Concurrent kallisto index builds use private temporary directories (MI-11)."""

from __future__ import annotations

import os
import random
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from metainformant.rna.amalgkit.genome_prep import build_kallisto_index, run_kallisto_index

pytestmark = pytest.mark.external_tool


@pytest.fixture
def transcripts(tmp_path: Path) -> Path:
    if shutil.which("kallisto") is None:
        pytest.skip("kallisto is not installed")
    rng = random.Random(11)
    fasta = tmp_path / "transcripts.fa"
    fasta.write_text("".join(f">t{i}\n{''.join(rng.choice('ACGT') for _ in range(400))}\n" for i in range(4)))
    return fasta


def test_parallel_builds_from_one_directory_all_succeed_and_leave_no_scratch(
    tmp_path: Path, transcripts: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # All builds share one working directory, which is what made ./tmp collide.
    shared_cwd = tmp_path / "cwd"
    shared_cwd.mkdir()
    monkeypatch.chdir(shared_cwd)
    indexes = [tmp_path / "index" / f"species_{i}.idx" for i in range(8)]

    with ThreadPoolExecutor(max_workers=len(indexes)) as pool:
        results = list(pool.map(lambda index: run_kallisto_index(transcripts, index), indexes))

    assert [result.returncode for result in results] == [0] * len(indexes), [r.stderr[-300:] for r in results]
    assert all(index.is_file() and index.stat().st_size > 0 for index in indexes)
    assert sorted(path.name for path in (tmp_path / "index").iterdir()) == sorted(index.name for index in indexes)
    assert list(shared_cwd.iterdir()) == []
    assert all("-T" in result.args for result in results)


def test_failed_build_still_removes_its_scratch_directory(tmp_path: Path, transcripts: Path) -> None:
    index = tmp_path / "index" / "broken.idx"
    result = run_kallisto_index(tmp_path / "absent.fa", index)
    assert result.returncode != 0
    assert os.listdir(index.parent) == []


def test_build_kallisto_index_keeps_its_kmer_and_reports_the_index(tmp_path: Path, transcripts: Path) -> None:
    index = tmp_path / "index" / "species.idx"
    assert build_kallisto_index(transcripts, index, kmer_size=25) == index
    assert [path.name for path in index.parent.iterdir()] == ["species.idx"]

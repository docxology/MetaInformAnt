from __future__ import annotations

import os
import random
import subprocess
from pathlib import Path

import pytest
import yaml

from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.acquisition_references import (
    _IndexDigests,
    prepare_reference_inputs,
)


def test_aliases_share_only_unchanged_resolved_index(tmp_path: Path) -> None:
    original = tmp_path / "reference.idx"
    original.write_bytes(b"original reference")
    alias = tmp_path / "alias.idx"
    alias.symlink_to(original.name)
    digests = _IndexDigests()
    expected = sha256_file(original)
    assert digests.digest(original) == digests.digest(alias) == expected
    assert len(digests.values) == 1

    before = original.stat()
    original.write_bytes(b"modified reference")
    os.utime(original, ns=(before.st_atime_ns, before.st_mtime_ns))
    assert original.stat().st_size == before.st_size
    assert digests.digest(alias) == sha256_file(original) != expected

    replacement = tmp_path / "replacement.idx"
    replacement.write_bytes(b"different reference")
    alias.unlink()
    alias.symlink_to(replacement.name)
    assert digests.digest(alias) == sha256_file(replacement)


def test_separate_files_and_calls_are_hashed_independently(tmp_path: Path) -> None:
    paths = [tmp_path / name for name in ("first.idx", "second.idx")]
    for path in paths:
        path.write_bytes(b"same content")
    digests = _IndexDigests()
    assert digests.digest(paths[0]) == digests.digest(paths[1])
    assert len(digests.values) == 2
    fresh = _IndexDigests()
    assert fresh.values == {}
    assert fresh.digest(paths[0]) == sha256_file(paths[0])


def test_disappearing_index_is_never_accepted_from_cache(tmp_path: Path) -> None:
    index = tmp_path / "removed.idx"
    index.write_bytes(b"reference")
    digests = _IndexDigests()
    digests.digest(index)
    index.unlink()
    with pytest.raises(FileNotFoundError):
        digests.digest(index)


@pytest.mark.external_tool
def test_native_reference_aliases_and_frozen_corruption_rejection(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = tmp_path / "data"
    work = data / "test_species/work"
    index = work / "index/Test_species.idx"
    index.parent.mkdir(parents=True)
    metadata = work / "metadata/metadata_selected.tsv"
    metadata.parent.mkdir()
    metadata.write_text("scientific_name\tunused\nTest species\tx\nTest species subspecies\ty\n")
    rng = random.Random(151)
    fasta = tmp_path / "transcripts.fa"
    fasta.write_text(">gene\n" + "".join(rng.choice("ACGT") for _ in range(1000)) + "\n")
    subprocess.run(["kallisto", "index", "-i", str(index), str(fasta)], check=True, capture_output=True)
    configs = tmp_path / "configs"
    configs.mkdir()
    config = configs / "test.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                "species_list": ["Test_species"],
                "genome": {"index_dir": str(index.parent)},
                "steps": {"quant": {"index_dir": str(index.parent)}},
                "reference_aliases": {"Test species subspecies": "Test_species"},
            }
        )
    )
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data))
    monkeypatch.setenv("AMALGKIT_LOCAL_INDEX_CACHE_DIR", str(tmp_path / "index-cache"))
    tasks = [{"species": "test_species", "config_name": config.name, "reference_index_sha256": sha256_file(index)}]
    records = prepare_reference_inputs(data_root=data, config_dir=configs, tasks=tasks)
    assert records[0]["targets"] == ["Test species", "Test species subspecies"]
    assert prepare_reference_inputs(data_root=data, config_dir=configs, tasks=tasks) == records
    alias = index.parent / "Test_species_subspecies.idx"
    assert alias.resolve() == index
    index.write_bytes(b"corrupt frozen index")
    with pytest.raises(ValueError, match="reference differs from frozen"):
        prepare_reference_inputs(data_root=data, config_dir=configs, tasks=tasks)

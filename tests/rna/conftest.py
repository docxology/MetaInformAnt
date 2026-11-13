from __future__ import annotations

import csv
import gzip
import json
import os
import random
import shutil
import stat
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from uuid import uuid4

import pytest

from metainformant.rna.amalgkit.genome_prep import run_kallisto_index
from metainformant.rna.engine.acquisition_manifest import sha256_file


@dataclass(frozen=True)
class NativeAcquisitionFiles:
    data: Path
    work: Path
    index: Path
    metadata: Path
    config_dir: Path
    config: Path
    manifest: Path

    def produce(self, accession: str = "SRR123", *, db_path: Path | None = None, **options: Any) -> dict[str, Any]:
        """Run one manifest task through a fresh streaming producer and return its sample result.

        A fresh producer reads the current environment, as each launch of the
        live lane does. The species' reference aliases are prepared first, as
        the producer does before it admits a species' samples; ``options`` are
        passed to ``process_single_sample``.
        """
        import yaml

        from metainformant.rna.engine.streaming_orchestrator import (
            StreamingPipelineOrchestrator,
            _ensure_reference_alias_indexes,
        )

        tasks = [json.loads(line) for line in self.manifest.read_text().splitlines()]
        task = next(task for task in tasks if task["accession"] == accession)
        config = yaml.safe_load((self.config_dir / task["config_name"]).read_text()) or {}
        with self.metadata.open(newline="") as handle:
            targets = sorted({row["scientific_name"] for row in csv.DictReader(handle, delimiter="\t")})
        ready, missing, _ = _ensure_reference_alias_indexes(config, task["species"], targets)
        assert ready, missing
        producer = StreamingPipelineOrchestrator(
            self.config_dir, self.data / "logs", db_path=db_path or self.data / "pipeline_progress.db"
        )
        try:
            producer.db.init_species(task["species"], [str(entry["accession"]) for entry in tasks])
            return producer.process_single_sample(
                accession,
                int(task["batch_index"]),
                self.work / "getfastq",
                self.config_dir / task["config_name"],
                task["species"],
                1,
                task["expected_paired"],
                **options,
            )
        finally:
            producer.db.close()


@pytest.fixture
def native_acquisition_factory(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    def build(
        *,
        accessions: tuple[str, ...] = ("SRR123",),
        reference_target: str = "Test species",
        paired_metadata: bool = False,
        sequence_length: int = 500,
        vary_read_positions: bool = False,
    ) -> NativeAcquisitionFiles:
        if not shutil.which("kallisto") or not shutil.which("amalgkit"):
            pytest.skip("Real Kallisto and Amalgkit executables are required")
        monkeypatch.delenv("AMALGKIT_MAX_RAW_BYTES", raising=False)
        data = tmp_path / "data"
        for name, value in {
            "AMALGKIT_DATA_ROOT": str(data),
            "AMALGKIT_MIN_EXTERNAL_FREE_GB": "0",
            "AMALGKIT_MIN_SYSTEM_FREE_GB": "0",
            "AMALGKIT_RECLAIM_RAW_AFTER_QUANT": "no",
            "AMALGKIT_PIPELINE_FASTQ_THREADS": "1",
            "AMALGKIT_PIPELINE_COMPRESSION_THREADS": "1",
            "AMALGKIT_PIPELINE_VALIDATION_SLOTS": "1",
            "AMALGKIT_PIPELINE_COMPRESSION_LEVEL": "1",
        }.items():
            monkeypatch.setenv(name, value)
        work = data / "test_species/work"
        index, metadata = (
            work / "index/Test_species.idx",
            work / "metadata/metadata_selected.tsv",
        )
        index.parent.mkdir(parents=True)
        metadata.parent.mkdir()
        rng = random.Random(161)
        sequences = ["".join(rng.choice("ACGT") for _ in range(sequence_length)) for _ in range(3)]
        fasta = tmp_path / "transcripts.fa"
        fasta.write_text("".join(f">gene{i}\n{sequence}\n" for i, sequence in enumerate(sequences)))
        built = run_kallisto_index(fasta, index)
        assert built.returncode == 0, built.stderr
        rows = []
        for accession in accessions:
            raw = work / "getfastq" / accession
            raw.mkdir(parents=True)
            with gzip.open(raw / f"{accession}.fastq.gz", "wt") as handle:
                for number in range(90):
                    start = 100 + (number if vary_read_positions else 0)
                    sequence = sequences[number % 3][start : start + 75]
                    handle.write(f"@{accession}.{number + 1}\n{sequence}\n+\n{'I' * 75}\n")
            rows.append(
                {
                    "run": accession,
                    "scientific_name": reference_target,
                    "lib_layout": "paired" if paired_metadata else "single",
                    "total_spots": "90",
                    "total_bases": "6750",
                    "spot_length": "75",
                    "size": "15000",
                    "sample_group": "test",
                    "is_sampled": "yes",
                    "is_qualified": "yes",
                    "exclusion": "no",
                }
            )
        with metadata.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        config_dir = tmp_path / "config"
        config_dir.mkdir()
        config = config_dir / "amalgkit_test_species.yaml"
        config.write_text(
            "reference_aliases:\n  Test species subspecies: Test_species\n"
            "species_list: [Test_species]\nsteps:\n  quant:\n"
            "    index_dir: output/amalgkit/test_species/work/index\n"
        )
        tasks = [
            {
                "schema": "metainformant.rna.acquisition_task.v1",
                "task_id": f"test_species/{accession}",
                "species": "test_species",
                "accession": accession,
                "batch_index": number,
                "config_name": config.name,
                "config_sha256": sha256_file(config),
                "reference_index_sha256": sha256_file(index),
                "expected_paired": paired_metadata,
                "fastq_bytes": 15000,
            }
            for number, accession in enumerate(accessions, 1)
        ]
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text("".join(json.dumps(task) + "\n" for task in tasks))
        (tmp_path / "snapshot.json").write_text(
            json.dumps(
                {
                    "schema": "metainformant.rna.acquisition_snapshot.v1",
                    "cloud_launch_policy": "checkpointed",
                    "task_count": len(tasks),
                    "manifest_sha256": sha256_file(manifest),
                    "input_files": [
                        {
                            "path": str(path.relative_to(tmp_path)),
                            "sha256": sha256_file(path),
                        }
                        for path in (index, metadata)
                    ],
                }
            )
        )
        return NativeAcquisitionFiles(data, work, index, metadata, config_dir, config, manifest)

    return build


STAGING_BULK_ROOT_VARIABLE = "AMALGKIT_STAGING_TEST_BULK_ROOT"


def _remove_test_tree(path: Path) -> None:
    """Remove a test-created tree, including entries a test made read-only."""

    if path.is_symlink() or path.is_file():
        path.unlink()
        return
    if not path.exists():
        return
    for directory, _subdirs, _files in os.walk(path):
        os.chmod(directory, stat.S_IRWXU)
    shutil.rmtree(path)


@contextmanager
def bulk_staging_allocator() -> Iterator[Callable[[str], Path]]:
    """Allocate uniquely named paths under ``$AMALGKIT_STAGING_TEST_BULK_ROOT``.

    Cross-volume cases need a directory on a second (bulk) volume that the
    operator names through the environment; allocation skips when it is unset.
    Every path handed out is removed on exit whether the body passed, failed or
    never created it, and exit fails if anything named by this allocator
    remains: tests must not leave pilot directories on the operator's volume.
    """

    allocated: list[Path] = []

    def allocate(prefix: str) -> Path:
        bulk = os.environ.get(STAGING_BULK_ROOT_VARIABLE)
        if not bulk:
            pytest.skip(f"Set {STAGING_BULK_ROOT_VARIABLE} for an isolated cross-volume pilot")
        path = Path(bulk) / f"{prefix}-{uuid4().hex}"
        allocated.append(path)
        return path

    try:
        yield allocate
    finally:
        for path in allocated:
            _remove_test_tree(path)
        leftovers = [str(path) for path in allocated if path.exists() or path.is_symlink()]
        leftovers += [
            str(sibling)
            for path in allocated
            if path.parent.is_dir()
            for sibling in path.parent.glob(f"*{path.name.rsplit('-', 1)[-1]}*")
        ]
        assert leftovers == [], f"test left directories under the bulk root: {leftovers}"


@pytest.fixture
def staging_bulk_dir() -> Iterator[Callable[[str], Path]]:
    """Fixture form of :func:`bulk_staging_allocator` (removal checked at teardown)."""

    with bulk_staging_allocator() as allocate:
        yield allocate

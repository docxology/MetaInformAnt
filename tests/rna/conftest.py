from __future__ import annotations

import csv
import gzip
import json
import random
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest

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
        for name in (
            "AMALGKIT_DURABLE_BUCKET",
            "AMALGKIT_DURABLE_COHORT",
            "AMALGKIT_CLOUD_MAX_RAW_BYTES",
            "AMALGKIT_WORKER_DEADLINE_EPOCH",
        ):
            monkeypatch.delenv(name, raising=False)
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
        sequences = [
            "".join(rng.choice("ACGT") for _ in range(sequence_length))
            for _ in range(3)
        ]
        fasta = tmp_path / "transcripts.fa"
        fasta.write_text(
            "".join(f">gene{i}\n{sequence}\n" for i, sequence in enumerate(sequences))
        )
        subprocess.run(
            ["kallisto", "index", "-i", str(index), str(fasta)],
            check=True,
            capture_output=True,
        )
        rows = []
        for accession in accessions:
            raw = work / "getfastq" / accession
            raw.mkdir(parents=True)
            with gzip.open(raw / f"{accession}.fastq.gz", "wt") as handle:
                for number in range(90):
                    start = 100 + (number if vary_read_positions else 0)
                    sequence = sequences[number % 3][start : start + 75]
                    handle.write(
                        f"@{accession}.{number + 1}\n{sequence}\n+\n{'I' * 75}\n"
                    )
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
        return NativeAcquisitionFiles(
            data, work, index, metadata, config_dir, config, manifest
        )

    return build

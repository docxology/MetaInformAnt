"""Deterministic frozen-cohort fixtures shared by the cohort, import and discovery tests.

Everything is real files written under a pytest ``tmp_path``: a frozen
acquisition snapshot (``snapshot.json`` + ``manifest.jsonl`` +
``inventory.json`` + ``inputs/``) exactly as the completion-inventory freeze
lays it out, and data roots with configuration, kallisto index bytes, a
reference manifest and quantification outputs written by the production
``write_quant_provenance``.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

from metainformant.rna.engine.provenance import write_quant_provenance

COHORT_ID = "cohort-test-20261112"
SPECIES = "test_species"
CONFIG_NAME = "amalgkit_test_species.yaml"
INDEX_NAME = "Test_species.idx"
CONFIG_BYTES = b"species_list: [Test_species]\nthreads: 1\n"
INDEX_BYTES = b"deterministic kallisto index fixture bytes\n" * 64
COLUMNS = ("scientific_name", "run", "total_bases", "lib_layout", "is_sampled", "exclusion")
# (run, total_bases, lib_layout); SRR104 is a cohort permanent_drop.
ROWS = (
    ("SRR103", 3_000_000_000, "PAIRED"),
    ("SRR101", 1_000_000_000, "SINGLE"),
    ("SRR104", 4_000_000_000, "SINGLE"),
    ("SRR102", 2_000_000_000, "PAIRED"),
    ("SRR105", 5_000_000_000_000, "SINGLE"),
)
EXCLUDED = ("SRR104",)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def table_bytes(rows: tuple[tuple[str, int, str], ...] = ROWS) -> bytes:
    """A frozen selection table: CRLF line endings, no quoting, every row sampled and not excluded."""

    lines = ["\t".join(COLUMNS)]
    for run, bases, layout in rows:
        lines.append("\t".join(("Test species", run, str(bases), layout, "yes", "no")))
    return ("\r\n".join(lines) + "\r\n").encode()


@dataclass
class FrozenSource:
    root: Path
    table: Path
    index: Path
    tasks: list[dict[str, object]] = field(default_factory=list)


def build_frozen_source(
    root: Path,
    *,
    rows: tuple[tuple[str, int, str], ...] = ROWS,
    excluded: tuple[str, ...] = EXCLUDED,
    config_bytes: bytes = CONFIG_BYTES,
    index_bytes: bytes = INDEX_BYTES,
    batch_override: dict[str, int] | None = None,
) -> FrozenSource:
    """Write a one-species frozen acquisition snapshot like ``completion_inventory.freeze_inventory``."""

    table = root / "inputs" / SPECIES / "work" / "metadata" / "metadata_selected.tsv"
    index = root / "inputs" / SPECIES / "work" / "index" / INDEX_NAME
    table.parent.mkdir(parents=True)
    index.parent.mkdir(parents=True)
    table.write_bytes(table_bytes(rows))
    index.write_bytes(index_bytes)
    tasks: list[dict[str, object]] = []
    for position, (run, bases, layout) in enumerate(rows, start=1):
        if run in excluded:
            continue
        tasks.append(
            {
                "accession": run,
                "batch_index": (batch_override or {}).get(run, position),
                "config_name": CONFIG_NAME,
                "config_sha256": sha(config_bytes),
                "existing_state": "pending",
                "expected_paired": layout == "PAIRED",
                "fastq_bytes": bases // 2,
                "reference_index_sha256": sha(index_bytes),
                "schema": "metainformant.hymenoptera.gcp_task_manifest.v1",
                "species": SPECIES,
                "task_id": f"{SPECIES}/{run}",
                "total_bases": float(bases),
            }
        )
    manifest = root / "manifest.jsonl"
    manifest.write_text("".join(json.dumps(task, sort_keys=True) + "\n" for task in tasks))
    inventory = {
        "schema": "metainformant.hymenoptera.completion_inventory.v1",
        "frozen_at": "2026-10-05T16:15:59.599841+00:00",
        "species": [
            {
                "species": SPECIES,
                "config_name": CONFIG_NAME,
                "taxid": 1,
                "config_sha256": sha(config_bytes),
                "metadata_sha256": sha(table.read_bytes()),
                "index_name": INDEX_NAME,
                "index_sha256": sha(index_bytes),
                "tasks": tasks,
                "excluded": [{"accession": run, "reason": "permanent_drop"} for run in excluded],
            }
        ],
    }
    (root / "inventory.json").write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    snapshot = {
        "schema": "metainformant.rna.acquisition_snapshot.v1",
        "cloud_launch_policy": "checkpointed",
        "amalgkit_version": "0.16.60",
        "amalgkit_release_tag": "v0.16.60",
        "amalgkit_source_revision": "c656a52aacdcee6fd3bf7e8031769ca957204ebc",
        "inventory_sha256": sha((root / "inventory.json").read_bytes()),
        "manifest_sha256": sha(manifest.read_bytes()),
        "input_files": [
            {"path": f"inputs/{SPECIES}/work/metadata/metadata_selected.tsv", "sha256": sha(table.read_bytes())},
            {"path": f"inputs/{SPECIES}/work/index/{INDEX_NAME}", "sha256": sha(index_bytes)},
        ],
        "species": [SPECIES],
        "task_count": len(tasks),
    }
    (root / "snapshot.json").write_text(json.dumps(snapshot, indent=2, sort_keys=True) + "\n")
    return FrozenSource(root, table, index, tasks)


class DataRoot:
    """A data root with a species reference (index + manifest) and quantification outputs."""

    def __init__(self, base: Path, *, index_bytes: bytes = INDEX_BYTES, config_bytes: bytes = CONFIG_BYTES) -> None:
        self.base = base
        self.config_dir = base / "repo_config"
        self.config_dir.mkdir(parents=True)
        self.config = self.config_dir / CONFIG_NAME
        self.config.write_bytes(config_bytes)
        self.index_dir = base / SPECIES / "genome" / "index"
        self.index_dir.mkdir(parents=True)
        self.index = self.index_dir / INDEX_NAME
        self.index.write_bytes(index_bytes)
        self.work = base / SPECIES / "work"
        self.manifest = self.work / "reference" / "reference_aliases.json"
        self.manifest.parent.mkdir(parents=True)
        self.manifest.write_text(
            json.dumps(
                {
                    "aliases": [{"index": str(self.index), "metadata_target": "Test species"}],
                    "index_dir": str(self.index_dir),
                    "kallisto_index": str(self.index),
                    "schema": "metainformant.rna.reference_aliases.v1",
                    "species": SPECIES,
                    "status": "complete",
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        (self.work / "quant").mkdir(parents=True)

    def sample(self, accession: str) -> Path:
        return self.work / "quant" / accession

    def write_sample(self, accession: str, *, counts_offset: int = 0, config: Path | None = None) -> Path:
        """Write a kallisto-shaped output and its sidecar (bound to this root's manifest)."""

        sample = self.sample(accession)
        sample.mkdir(parents=True)
        rows = "".join(f"transcript{i}\t1000\t850.5\t{10 + i + counts_offset}.25\t{1000 + i}.5\n" for i in range(8))
        (sample / f"{accession}_abundance.tsv").write_text("target_id\tlength\teff_length\test_counts\ttpm\n" + rows)
        (sample / f"{accession}_run_info.json").write_text(
            json.dumps({"n_targets": 8, "n_processed": 1000, "n_pseudoaligned": 900, "kallisto_version": "0.52.0"})
        )
        (sample / f"{accession}_abundance.h5").write_bytes(b"HDF5 fixture " + accession.encode() * 16)
        write_quant_provenance(
            sample,
            species=SPECIES,
            run_accession=accession,
            config_path=config or self.config,
            command=["amalgkit", "quant", "--out_dir", str(self.work), "--batch", "1"],
            reference_manifest_path=self.manifest,
            quantification_file=sample / f"{accession}_abundance.tsv",
        )
        return sample


def write_species_yaml(config_dir: Path, index_dir: Path, *, config_bytes: bytes | None = None) -> Path:
    """Write the species YAML used by orchestrator tests (its bytes define the cohort config digest)."""

    path = config_dir / CONFIG_NAME
    path.write_bytes(
        config_bytes
        or (f"species_list: [Test_species]\ngenome:\n  index_dir: {index_dir}\nsteps:\n  quant: {{}}\n").encode()
    )
    return path

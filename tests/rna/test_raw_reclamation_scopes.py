"""Raw-input reclamation scopes: getfastq raw, fasterq staging, the shared VDB cache and the audited CLI.

Everything runs on real files under ``tmp_path``. Reusable quantifications are
written by the production ``write_quant_provenance`` through the frozen-cohort
fixtures, and the CLI is the production entry point.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from metainformant.rna.engine.campaign_lock import LOCK_FILENAME, writer_status
from metainformant.rna.engine.cohort import install_cohort
from metainformant.rna.engine.provenance import write_quant_provenance
from metainformant.rna.engine.raw_cleanup import (
    _is_raw_read_name,
    reclaim_fasterq_staging,
    reclaim_sample_raw_inputs,
    reclaim_sample_sra_cache,
)
from metainformant.rna.engine.raw_reclamation import main
from tests.rna.frozen_cohort_support import (
    COHORT_ID,
    SPECIES,
    DataRoot,
    build_frozen_source,
)


def _tree(root: Path) -> dict[str, bytes]:
    # The writer lock's owner record is rewritten by every holder; everything else must be byte-stable.
    return {
        str(path.relative_to(root)): path.read_bytes()
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.name != LOCK_FILENAME
    }


def _write(path: Path, data: bytes = b"payload") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return path


PRODUCER_LOCK_HOLDER = """
import sys, time
from metainformant.rna.engine.campaign_lock import acquire
lock = acquire(sys.argv[1], "producer")
print("held", flush=True)
time.sleep(600)
"""


@pytest.mark.parametrize(
    "name",
    [
        "SRR1_1.fastq.gz",
        "SRR1_1.fastq.gz.part",
        "SRR1_1.fastq.gz.part.invalid",
        "SRR1_2.fastq.gz.invalid",
        "SRR1_2.fastq.gz.invalid-1730000000000",
        "SRR1.fastq",
        "SRR1.sra",
        "SRR1.sra.cache",
        "SRR1.sra.vdbcache",
        "SRR1.sra.tmp",
        "SRR1.sra.prf",
        "SRR1.sra.lock",
        "SRR1.sra.part",
        "SRR1_2.fastq.gz.safely_removed",
    ],
)
def test_raw_names_include_failed_partials_and_vdb_states(name: str) -> None:
    assert _is_raw_read_name(name)


@pytest.mark.parametrize(
    "name",
    [
        ".transfer_integrity_failures.tsv",
        "SRR1.json",
        "SRR1.sra.cache.json",
        "SRR1_abundance.tsv",
        "SRR1_run_info.json",
        ".metainformant_quant_provenance.json",
        "SRR1.fastq.gz.sha256",
    ],
)
def test_evidence_and_outputs_are_never_raw(name: str) -> None:
    assert not _is_raw_read_name(name)


def test_sample_reclaim_takes_part_invalid_and_keeps_evidence(tmp_path: Path) -> None:
    work = tmp_path / "work"
    sample = work / "getfastq" / "SRR123456"
    _write(sample / "SRR123456_1.fastq.gz.part.invalid", b"x" * 300)
    _write(sample / "SRR123456_1.fastq.gz.part", b"x" * 200)
    _write(sample / ".prefetch" / "SRR123456" / "SRR123456.sra.cache", b"x" * 100)
    evidence = _write(sample / ".transfer_integrity_failures.tsv", b"evidence\n")

    planned = reclaim_sample_raw_inputs(work, "SRR123456", dry_run=True)
    assert planned["files_deleted"] == 3 and planned["bytes_freed"] == 600
    assert (sample / "SRR123456_1.fastq.gz.part.invalid").exists()

    executed = reclaim_sample_raw_inputs(work, "SRR123456")
    assert executed["files_deleted"] == 3
    assert _tree(sample) == {".transfer_integrity_failures.tsv": b"evidence\n"}
    assert not (sample / ".prefetch").exists(), "emptied nested directories are removed"
    assert evidence.exists()


def test_staging_reclaim_is_exact_and_keeps_validation_witnesses(tmp_path: Path) -> None:
    work = tmp_path / "work"
    staging = work / "getfastq" / ".fasterq_tmp"
    _write(staging / "SRR123456" / "output" / "SRR123456.fastq", b"y" * 1000)
    _write(staging / "SRR123456" / "temp" / "fasterq.tmp.timelock.local.1" / "chunk", b"y" * 10)
    other = _write(staging / "SRR1234567" / "output" / "SRR1234567.fastq", b"other run")
    witness = _write(staging / "raw_validation" / "SRR123456.json", b"{}")

    planned = reclaim_fasterq_staging(work, "SRR123456", dry_run=True)
    assert planned["files_deleted"] == 2 and planned["bytes_freed"] == 1010
    assert (staging / "SRR123456").exists()

    executed = reclaim_fasterq_staging(work, "SRR123456")
    assert executed["files_deleted"] == 2
    assert not (staging / "SRR123456").exists()
    assert other.exists() and witness.exists()

    _write(work / "getfastq" / "SRR123456" / "SRR123456.fastq.gz", b"z" * 5)
    _write(staging / "SRR123456" / "output" / "SRR123456.fastq", b"y")
    combined = reclaim_sample_raw_inputs(work, "SRR123456", include_staging=True)
    assert combined["files_deleted"] == 2
    assert not (staging / "SRR123456").exists()


def test_staging_with_a_symlink_is_protected(tmp_path: Path) -> None:
    work = tmp_path / "work"
    outside = _write(tmp_path / "outside" / "precious.fastq", b"keep")
    workspace = work / "getfastq" / ".fasterq_tmp" / "SRR123456"
    _write(workspace / "output" / "SRR123456.fastq", b"y")
    (workspace / "output" / "link").symlink_to(outside.parent, target_is_directory=True)

    result = reclaim_fasterq_staging(work, "SRR123456")
    assert result["files_deleted"] == 0
    assert result["protected_paths"] and result["errors"]
    assert (workspace / "output" / "SRR123456.fastq").exists() and outside.exists()


def test_sra_cache_reclaim_is_accession_exact(tmp_path: Path) -> None:
    cache = tmp_path / ".sra-cache"
    _write(cache / "sra" / "SRR123456.sra.cache", b"c" * 64)
    _write(cache / "sra" / "SRR123456.sra", b"s" * 32)
    longer = _write(cache / "sra" / "SRR1234567.sra", b"other")
    marker = _write(cache / "sra" / ".metainformant_sra_validation" / "SRR123456.sra.cache.json", b"{}")
    refseq = _write(cache / "refseq" / "CM000054.5", b"ref")
    tree = cache / "prefetch" / "SRR123456" / "SRR123456"
    for suffix in (".sra", ".sra.vdbcache", ".sra.tmp", ".sra.prf", ".sra.lock"):
        _write(tree / f"SRR123456{suffix}", b"p" * 8)
    other_tree = _write(cache / "prefetch" / "SRR1234567" / "SRR1234567" / "SRR1234567.sra", b"other")

    planned = reclaim_sample_sra_cache(cache, "SRR123456", dry_run=True)
    assert planned["files_deleted"] == 7 and planned["bytes_freed"] == 64 + 32 + 5 * 8
    assert (cache / "sra" / "SRR123456.sra.cache").exists()

    executed = reclaim_sample_sra_cache(cache, "SRR123456")
    assert executed["files_deleted"] == 7
    assert not (cache / "prefetch" / "SRR123456").exists()
    assert not (cache / "sra" / "SRR123456.sra").exists()
    assert longer.exists() and marker.exists() and refseq.exists() and other_tree.exists()


def test_sra_cache_keeps_foreign_files_and_refuses_symlinks(tmp_path: Path) -> None:
    cache = tmp_path / ".sra-cache"
    foreign = _write(cache / "prefetch" / "SRR123456" / "notes.txt", b"operator note")
    _write(cache / "prefetch" / "SRR123456" / "SRR123456" / "SRR123456.sra", b"p")
    outside = _write(tmp_path / "elsewhere" / "SRR123456.sra", b"outside")
    (cache / "sra").mkdir(parents=True)
    (cache / "sra" / "SRR123456.sra").symlink_to(outside)

    result = reclaim_sample_sra_cache(cache, "SRR123456")
    assert result["files_deleted"] == 1
    assert foreign.exists() and outside.exists()
    assert str(cache.resolve() / "sra" / "SRR123456.sra") in result["protected_paths"]

    with pytest.raises(ValueError):
        reclaim_sample_sra_cache(cache, "SRR*")


# ---------------------------------------------------------------------------
# Command line on a cohort-bound data root
# ---------------------------------------------------------------------------


class Root:
    """A frozen cohort and a data root holding outputs, raw inputs, staging and cache entries."""

    def __init__(self, base: Path) -> None:
        build_frozen_source(base / "frozen")
        install_cohort(base / "frozen", base / "cohort", cohort_id=COHORT_ID)
        self.base = base
        self.cohort = base / "cohort"
        self.root = DataRoot(base / "data")
        self.data = self.root.base
        self.root.write_sample("SRR101")  # reusable
        other = base / "other.yaml"
        other.write_bytes(b"species_list: [Test_species]\nthreads: 2\n")
        self.root.write_sample("SRR103", config=other)  # different configuration bytes: not reusable
        getfastq = self.root.work / "getfastq"
        for run in ("SRR101", "SRR102", "SRR103", "SRR104"):
            _write(getfastq / run / f"{run}_1.fastq.gz", b"r" * 100)
        _write(getfastq / "SRR101" / "SRR101_2.fastq.gz.part.invalid", b"r" * 50)
        _write(getfastq / ".fasterq_tmp" / "SRR101" / "output" / "SRR101.fastq", b"f" * 70)
        _write(getfastq / ".fasterq_tmp" / "SRR102" / "output" / "SRR102.fastq", b"f" * 70)
        _write(getfastq / ".fasterq_tmp" / "raw_validation" / "SRR101.json", b"{}")
        cache = self.data / ".sra-cache"
        _write(cache / "sra" / "SRR101.sra.cache", b"c" * 40)
        _write(cache / "sra" / "SRR102.sra", b"c" * 40)
        _write(cache / "sra" / ".metainformant_sra_validation" / "SRR101.sra.cache.json", b"{}")
        _write(cache / "prefetch" / "SRR104" / "SRR104" / "SRR104.sra", b"c" * 30)
        _write(cache / "prefetch" / "SRR103" / "SRR103" / "SRR103.sra", b"c" * 30)

    def args(self, *extra: str) -> list[str]:
        return [
            "--data-root",
            str(self.data),
            "--config-dir",
            str(self.root.config_dir),
            "--cohort",
            str(self.cohort),
            "--workers",
            "2",
            *extra,
        ]


@pytest.fixture
def root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Root:
    monkeypatch.delenv("AMALGKIT_DATA_ROOT", raising=False)
    return Root(tmp_path)


ALL_SCOPES = ("--include-staging", "--include-sra-cache", "--include-excluded")


def test_cli_dry_run_writes_nothing_and_plans_exact_scopes(root: Root, capsys: pytest.CaptureFixture[str]) -> None:
    before = _tree(root.data)
    assert main(root.args(*ALL_SCOPES)) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["mode"] == "dry-run" and summary["manifest"] is None
    assert summary["classification_counts"] == {"current": 1, "invalid": 1}
    # SRR101 (reusable): 2 raw + 1 staging file; SRR104 (permanent_drop): 1 raw file.
    assert summary["raw_files"] == 4 and summary["raw_bytes"] == 100 + 50 + 70 + 100
    # SRR101 cache file and SRR104 prefetch archive; SRR102/SRR103 caches stay.
    assert summary["sra_cache_files"] == 2 and summary["sra_cache_bytes"] == 70
    assert _tree(root.data) == before


def test_cli_execute_is_exact_audited_and_idempotent(root: Root, capsys: pytest.CaptureFixture[str]) -> None:
    manifest = root.base / "audit" / "reclaim.json"
    assert main(root.args(*ALL_SCOPES, "--execute", "--manifest", str(manifest))) == 0
    capsys.readouterr()
    getfastq = root.root.work / "getfastq"
    cache = root.data / ".sra-cache"
    assert not (getfastq / "SRR101").exists() and not (getfastq / "SRR104").exists()
    assert (getfastq / "SRR102" / "SRR102_1.fastq.gz").exists(), "pending task keeps its raw input"
    assert (getfastq / "SRR103" / "SRR103_1.fastq.gz").exists(), "non-reusable output keeps its raw input"
    assert not (getfastq / ".fasterq_tmp" / "SRR101").exists()
    assert (getfastq / ".fasterq_tmp" / "SRR102" / "output" / "SRR102.fastq").exists()
    assert (getfastq / ".fasterq_tmp" / "raw_validation" / "SRR101.json").exists()
    assert not (cache / "sra" / "SRR101.sra.cache").exists() and not (cache / "prefetch" / "SRR104").exists()
    assert (cache / "sra" / "SRR102.sra").exists() and (cache / "prefetch" / "SRR103").exists()
    assert (cache / "sra" / ".metainformant_sra_validation" / "SRR101.sra.cache.json").exists()
    assert (root.root.sample("SRR101") / "SRR101_abundance.tsv").exists(), "outputs are never touched"

    record = json.loads(manifest.read_text())
    assert record["mode"] == "execute" and record["cohort_id"] == COHORT_ID
    assert record["released"][SPECIES] == {"SRR101": "reusable_quantification", "SRR104": "permanent_exclusion"}
    assert record["totals"]["raw_files"] == 4 and record["totals"]["sra_cache_files"] == 2
    assert (root.data / LOCK_FILENAME).exists(), "--execute ran under the campaign writer lock"
    assert not writer_status(root.data)["held"], "the lock is released when the pass ends"

    after = _tree(root.data)
    assert main(root.args(*ALL_SCOPES, "--execute", "--manifest", str(root.base / "audit" / "again.json"))) == 0
    again = json.loads(capsys.readouterr().out)
    assert again["raw_files"] == 0 and again["sra_cache_files"] == 0
    assert _tree(root.data) == after


def test_cli_execute_default_manifest_lives_under_results(root: Root, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(root.args("--execute")) == 0
    summary = json.loads(capsys.readouterr().out)
    manifest = Path(summary["manifest"])
    assert manifest.parent == (root.data / "results").resolve()
    assert json.loads(manifest.read_text())["scopes"] == {
        "raw": True,
        "fasterq_staging": False,
        "sra_cache": False,
        "permanent_exclusions": False,
    }
    # Without --include-excluded the permanent_drop task keeps its raw input.
    assert (root.root.work / "getfastq" / "SRR104" / "SRR104_1.fastq.gz").exists()


def test_cli_refuses_while_a_producer_holds_the_root(root: Root) -> None:
    # A live producer holds the campaign writer lock (a separate process, like the real one).
    holder = subprocess.Popen(
        [sys.executable, "-c", PRODUCER_LOCK_HOLDER, str(root.data)], stdout=subprocess.PIPE, text=True
    )
    try:
        assert holder.stdout is not None and holder.stdout.readline().strip() == "held"
        before = _tree(root.data)
        with pytest.raises(SystemExit, match="campaign lock held by producer"):
            main(root.args("--execute", "--manifest", str(root.base / "b.json")))
        assert _tree(root.data) == before
        assert not (root.base / "b.json").exists()
    finally:
        holder.kill()
        holder.wait(timeout=30)
        if holder.stdout is not None:
            holder.stdout.close()


def test_cli_refuses_cache_scope_without_cohort_and_config_drift(root: Root, tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="need --cohort"):
        main(["--data-root", str(root.data), "--config-dir", str(root.root.config_dir), "--include-sra-cache"])
    edited = tmp_path / "edited"
    shutil.copytree(root.root.config_dir, edited)
    config = next(edited.glob("*.yaml"))
    config.write_bytes(config.read_bytes() + b"# edited\n")
    before = _tree(root.data)
    with pytest.raises(SystemExit, match="does not hash to the cohort configuration"):
        main(["--data-root", str(root.data), "--config-dir", str(edited), "--cohort", str(root.cohort), "--execute"])
    assert _tree(root.data) == before


def test_cli_without_cohort_binds_configuration_by_content(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Without a cohort, the sidecar's config hash must equal the current YAML bytes (paths are irrelevant)."""

    species = "demo_ant"
    data = tmp_path / "data"
    config_dir = tmp_path / "repo" / "config"
    config = _write(config_dir / f"amalgkit_{species}.yaml", b"species_list: [Demo_ant]\nthreads: 1\n")
    index = _write(data / species / "genome" / "index" / "Demo_ant.idx", b"index bytes\n" * 32)
    manifest = data / species / "work" / "reference" / "reference_aliases.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(
        json.dumps(
            {
                "aliases": [{"index": str(index), "metadata_target": "Demo ant"}],
                "index_dir": str(index.parent),
                "kallisto_index": str(index),
                "schema": "metainformant.rna.reference_aliases.v1",
                "species": species,
                "status": "complete",
            }
        )
    )
    moved_config = _write(tmp_path / "old-checkout" / config.name, config.read_bytes())
    drifted_config = _write(tmp_path / "old-checkout" / "drifted.yaml", b"species_list: [Demo_ant]\nthreads: 9\n")
    for run, recorded_config in (("SRR201", moved_config), ("SRR202", drifted_config)):
        sample = data / species / "work" / "quant" / run
        rows = "".join(f"t{i}\t1000\t850.5\t{10 + i}.25\t{1000 + i}.5\n" for i in range(8))
        _write(sample / f"{run}_abundance.tsv", ("target_id\tlength\teff_length\test_counts\ttpm\n" + rows).encode())
        _write(
            sample / f"{run}_run_info.json",
            json.dumps(
                {"n_targets": 8, "n_processed": 1000, "n_pseudoaligned": 900, "kallisto_version": "0.52.0"}
            ).encode(),
        )
        write_quant_provenance(
            sample,
            species=species,
            run_accession=run,
            config_path=recorded_config,
            command=["amalgkit", "quant"],
            reference_manifest_path=manifest,
            quantification_file=sample / f"{run}_abundance.tsv",
        )
        _write(data / species / "work" / "getfastq" / run / f"{run}.fastq.gz", b"r" * 10)
    shutil.rmtree(tmp_path / "old-checkout")  # recorded config paths no longer exist

    previous_root = os.environ.get("AMALGKIT_DATA_ROOT")
    assert main(["--data-root", str(data), "--config-dir", str(config_dir), "--execute", "--workers", "1"]) == 0
    assert os.environ.get("AMALGKIT_DATA_ROOT") == previous_root, "the CLI must not leak its data root"
    summary = json.loads(capsys.readouterr().out)
    assert summary["raw_files"] == 1
    assert not (data / species / "work" / "getfastq" / "SRR201").exists()
    assert (data / species / "work" / "getfastq" / "SRR202" / "SRR202.fastq.gz").exists()


# ---------------------------------------------------------------------------
# Producer: reclaim after each reusable quantification
# ---------------------------------------------------------------------------


def _producer_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[DataRoot, Path, object]:
    from metainformant.rna.engine.cohort import load_frozen_cohort
    from tests.rna.frozen_cohort_support import write_species_yaml

    data = DataRoot(tmp_path / "data")
    for name, value in {
        "AMALGKIT_DATA_ROOT": str(data.base),
        "AMALGKIT_MIN_EXTERNAL_FREE_GB": "0",
        "AMALGKIT_MIN_SYSTEM_FREE_GB": "0",
        "AMALGKIT_RECLAIM_RAW_AFTER_QUANT": "yes",
    }.items():
        monkeypatch.setenv(name, value)
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    yaml = write_species_yaml(config_dir, data.index_dir)
    build_frozen_source(tmp_path / "frozen", config_bytes=yaml.read_bytes())
    install_cohort(tmp_path / "frozen", tmp_path / "cohort", cohort_id=COHORT_ID)
    data.write_sample("SRR101", config=yaml)
    for run in ("SRR101", "SRR102"):
        _write(data.work / "getfastq" / run / f"{run}_1.fastq.gz", b"r" * 10)
        _write(data.work / "getfastq" / ".fasterq_tmp" / run / "output" / f"{run}.fastq", b"f" * 10)
        _write(data.base / ".sra-cache" / "sra" / f"{run}.sra.cache", b"c" * 10)
    return data, config_dir, load_frozen_cohort(tmp_path / "cohort")


def test_producer_clears_staging_and_cache_after_cohort_quantification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from metainformant.rna.engine.streaming_orchestrator import (
        StreamingPipelineOrchestrator,
    )

    data, config_dir, cohort = _producer_root(tmp_path, monkeypatch)
    orchestrator = StreamingPipelineOrchestrator(
        config_dir, tmp_path / "logs", db_path=tmp_path / "progress.db", cohort=cohort
    )
    try:
        orchestrator._reclaim_after_current_quantification(SPECIES, "SRR101")
        orchestrator._reclaim_after_current_quantification(SPECIES, "SRR102")  # not quantified: untouched
    finally:
        orchestrator.db.close()
    getfastq = data.work / "getfastq"
    assert not (getfastq / "SRR101").exists()
    assert not (getfastq / ".fasterq_tmp" / "SRR101").exists()
    assert not (data.base / ".sra-cache" / "sra" / "SRR101.sra.cache").exists()
    assert (getfastq / "SRR102" / "SRR102_1.fastq.gz").exists()
    assert (getfastq / ".fasterq_tmp" / "SRR102" / "output" / "SRR102.fastq").exists()
    assert (data.base / ".sra-cache" / "sra" / "SRR102.sra.cache").exists()
    assert (data.sample("SRR101") / "SRR101_abundance.tsv").exists()


def test_producer_without_cohort_keeps_the_shared_cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from metainformant.rna.engine.streaming_orchestrator import (
        StreamingPipelineOrchestrator,
    )

    data, config_dir, _cohort = _producer_root(tmp_path, monkeypatch)
    orchestrator = StreamingPipelineOrchestrator(config_dir, tmp_path / "logs", db_path=tmp_path / "progress.db")
    try:
        orchestrator._reclaim_after_current_quantification(SPECIES, "SRR101")
    finally:
        orchestrator.db.close()
    assert not (data.work / "getfastq" / "SRR101").exists()
    assert not (data.work / "getfastq" / ".fasterq_tmp" / "SRR101").exists()
    assert (data.base / ".sra-cache" / "sra" / "SRR101.sra.cache").exists()

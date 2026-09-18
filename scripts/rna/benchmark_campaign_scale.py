#!/usr/bin/env python3
"""Campaign-scale performance benchmark harness (disposable fixtures only).

Characterizes the orchestration-side costs of the current RNA-seq campaign
(discovery, queue accounting, retry/backoff scheduling, hash/provenance I/O,
and SQLite storage growth) at parameterized task-row scales WITHOUT touching
the live producer, the live data root, or the live progress database:

- constructs a DISPOSABLE SQLite progress DB via the real ``ProgressDB``
  implementation, populated with synthetic task rows (default scales
  1e3 / 1e4 / 1e5)
- builds a synthetic species-config / metadata / quantification filesystem
  fixture bound to the disposable DB
- measures discovery throughput, queue-depth computation, retry/backoff
  scheduling cost, hash/provenance I/O, and storage growth
- emits a JSON report and a Markdown report, both carrying a
  machine-provenance block (commit, host, python, seed, scale parameters)

Safety contract (enforced in code):
- the live data root is never opened for reading or writing; any target
  directory at or beneath it is refused
- the default target directory is ``output/campaign_scale_benchmarks/`` inside
  the repository; the fixture is deleted after the run unless
  ``--keep-fixture`` is passed
- no network access, no external tools

Determinism: all fixture content (species names, run accessions, state
assignment, exclusions, abundance payloads) is generated from a fixed seed;
timings are inherently machine-dependent and are reported as descriptive
observations, not statistical claims.

Usage:
    python scripts/rna/benchmark_campaign_scale.py --scales 1000
    python scripts/rna/benchmark_campaign_scale.py                      # full sweep 1e3/1e4/1e5
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import random
import socket
import sqlite3
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from metainformant.rna.engine.progress_db import ProgressDB, classify_sample_error  # noqa: E402
from metainformant.rna.engine.provenance import (  # noqa: E402
    QUANT_STATUS_CURRENT,
    classify_quantification,
    digest_file,
    read_quant_provenance,
    write_quant_provenance,
)
from metainformant.rna.engine.species import discover_species_config_names  # noqa: E402

BENCHMARK_SCHEMA = "metainformant.rna.benchmark.campaign_scale.v1"
DEFAULT_SEED = 20260917
DEFAULT_SCALES = (1_000, 10_000, 100_000)
DEFAULT_SPECIES_COUNT = 8
DEFAULT_REPEATS = 3
DEFAULT_PROVENANCE_SAMPLE_CAP = 200

# Never opened: the guard below compares path strings lexically (no stat/scan).
LIVE_DATA_ROOT = Path("/Volumes/external_drive/Data/amalgkit")
LIVE_DB_PATH = LIVE_DATA_ROOT / "pipeline_progress.db"

# Deterministic state distribution for a mid-campaign snapshot.
STATE_WEIGHTS: Tuple[Tuple[str, float], ...] = (
    ("pending", 0.58),
    ("downloading", 0.05),
    ("downloaded", 0.14),
    ("quantifying", 0.05),
    ("quantified", 0.15),
    ("failed", 0.02),
    ("quarantined", 0.01),
)

# Realistic failure texts exercising classify_sample_error() marker order.
FAILURE_TEXTS: Tuple[str, ...] = (
    "Download Failed (all sources): curl (56) OpenSSL SSL_read connection reset",
    "fasterq-dump timeout after 3600s while extracting SRA",
    "Quant timeout: kallisto exceeded wall clock budget",
    "Quantification Failed: exit status 1",
    "Quant exception: Traceback (most recent call last)",
)

EXCLUSION_REASON_DETAIL = "terminal-failure audit M-02: environmental class, retryable elsewhere"


# ---------- Safety guard ----------


def _lexical_abspath(raw: str | os.PathLike[str]) -> Path:
    """Resolve a path lexically (no filesystem access) to an absolute path."""
    return Path(os.path.abspath(os.fspath(os.path.expanduser(str(raw)))))


def guard_target_dir(raw_target: str | os.PathLike[str]) -> Path:
    """Refuse any target directory at or beneath the live data root.

    Comparison is lexical: the live data root itself is never stat'ed, read,
    or created by this harness.
    """

    resolved = _lexical_abspath(raw_target)
    live = _lexical_abspath(LIVE_DATA_ROOT)
    if resolved == live or live in resolved.parents:
        raise SystemExit(f"refusing benchmark target inside the live data root: {resolved} " f"(protected: {live})")
    return resolved


# ---------- Machine provenance ----------


def _git(repo: Path, *args: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(repo), *args],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    output = result.stdout.strip()
    return output or None


def machine_provenance(*, seed: int, scales: Sequence[int], target_dir: Path, repo: Path = REPO_ROOT) -> dict:
    """Machine-provenance block shared by the JSON and Markdown reports."""

    head = _git(repo, "rev-parse", "HEAD")
    status = _git(repo, "status", "--porcelain")
    return {
        "schema": BENCHMARK_SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": head,
        "git_branch": _git(repo, "rev-parse", "--abbrev-ref", "HEAD"),
        # None when git status could not be captured in time (unknown).
        "git_dirty": bool(status) if status is not None else None,
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python_version": sys.version.split()[0],
        "seed": seed,
        "scales": list(scales),
        "target_dir": str(_lexical_abspath(target_dir)),
        "execution_context": "local disposable fixture (not a hosted run)",
    }


# ---------- Synthetic fixture ----------


@dataclass
class CampaignFixture:
    """Disposable synthetic campaign fixture bound to one progress DB."""

    target_dir: Path
    config_dir: Path
    metadata_dir: Path
    quant_root: Path
    db_path: Path
    species: List[str]
    task_rows: int
    seed: int
    cohort: Dict[str, List[str]] = field(default_factory=dict)
    state_map: Dict[str, Dict[str, List[str]]] = field(default_factory=dict)
    sample_dirs: Dict[str, List[Path]] = field(default_factory=dict)


def synthetic_species_names(count: int) -> List[str]:
    """Deterministic species identifiers (no reserved template/test markers)."""

    return [f"formica_bench_{i:04d}" for i in range(1, count + 1)]


def synthetic_run_accessions(index: int) -> str:
    return f"SRR{index:010d}"


def write_synthetic_configs(config_dir: Path, species: Sequence[str]) -> None:
    """Write minimal runnable-shaped species configs for the discovery scan."""

    config_dir.mkdir(parents=True, exist_ok=True)
    for name in species:
        work = f"output/amalgkit/{name}/work"
        config_dir.joinpath(f"amalgkit_{name}.yaml").write_text(
            "\n".join(
                [
                    "# METAINFORMANT benchmark fixture (disposable, synthetic)",
                    f"work_dir: {work}",
                    f"log_dir: output/amalgkit/{name}/logs",
                    "threads: 8",
                    "auto_install_amalgkit: false",
                    "filters:",
                    "  require_tissue: false",
                    "species_list:",
                    f"  - {name.capitalize()}",
                    f"taxon_id: {100000 + (sum(name.encode()) % 89999)}",
                    "",
                ]
            ),
            encoding="utf-8",
        )


def write_synthetic_metadata(
    metadata_dir: Path,
    species: Sequence[str],
    task_rows: int,
    rng: random.Random,
) -> Dict[str, List[str]]:
    """Write per-species metadata TSVs; return the parsed cohort (species -> SRR ids)."""

    metadata_dir.mkdir(parents=True, exist_ok=True)
    cohort: Dict[str, List[str]] = {name: [] for name in species}
    run_index = 0
    remaining = task_rows
    for position, name in enumerate(species):
        is_last = position == len(species) - 1
        count = remaining if is_last else max(1, task_rows // len(species))
        count = min(count, remaining)
        remaining -= count
        rows = ["run_accession\tscientific_name\tlibrary_layout"]
        for _ in range(count):
            run_index += 1
            accession = synthetic_run_accessions(run_index)
            cohort[name].append(accession)
            layout = rng.choice(("PAIRED", "SINGLE"))
            rows.append(f"{accession}\t{name.capitalize()}\t{layout}")
        metadata_dir.joinpath(f"metadata_{name}.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    return cohort


def assign_states(cohort: Dict[str, List[str]], rng: random.Random) -> Dict[str, Dict[str, List[str]]]:
    """Deterministically partition each species' cohort across the state machine."""

    weights = [weight for _, weight in STATE_WEIGHTS]
    state_map: Dict[str, Dict[str, List[str]]] = {}
    for name, srr_ids in cohort.items():
        total = len(srr_ids)
        counts = [int(total * weight) for weight in weights]
        counts[-1] = total - sum(counts[:-1])
        states: List[str] = []
        for (state, _), count in zip(STATE_WEIGHTS, counts):
            states.extend([state] * count)
        rng.shuffle(states)
        bucket: Dict[str, List[str]] = {state: [] for state, _ in STATE_WEIGHTS}
        for srr_id, state in zip(srr_ids, states):
            bucket[state].append(srr_id)
        state_map[name] = bucket
    return state_map


def abundance_tsv_payload(rng: random.Random, gene_count: int = 64) -> str:
    """Deterministic kallisto-shaped abundance table (small, fixed size)."""

    rows = ["target_id\tlength\teff_length\test_counts\ttpm"]
    for gene in range(gene_count):
        length = 400 + rng.randrange(2_600)
        eff = round(length - rng.uniform(20.0, 60.0), 3)
        est = round(rng.uniform(1.0, 5_000.0), 2)
        tpm = round(rng.uniform(0.01, 400.0), 4)
        rows.append(f"gene_{gene:04d}\t{length}\t{eff}\t{est}\t{tpm}")
    return "\n".join(rows) + "\n"


def write_synthetic_quant_fixtures(
    quant_root: Path,
    cohort: Dict[str, List[str]],
    config_dir: Path,
    *,
    sample_cap: int,
    rng: random.Random,
) -> Dict[str, List[Path]]:
    """Create per-sample quant dirs (abundance.tsv + provenance sidecar).

    Only the first ``sample_cap`` runs of each species get a filesystem
    fixture; the DB rows exist for the full cohort.
    """

    quant_root.mkdir(parents=True, exist_ok=True)
    sample_dirs: Dict[str, List[Path]] = {}
    for name, srr_ids in cohort.items():
        config_path = config_dir / f"amalgkit_{name}.yaml"
        species_dirs: List[Path] = []
        for srr_id in srr_ids[:sample_cap]:
            sample_dir = quant_root / name / "work" / "quant" / srr_id
            sample_dir.mkdir(parents=True, exist_ok=True)
            sample_dir.joinpath("abundance.tsv").write_text(abundance_tsv_payload(rng), encoding="utf-8")
            write_quant_provenance(
                sample_dir,
                species=name,
                run_accession=srr_id,
                config_path=config_path,
                command=["kallisto", "quant", "-i", "benchmark_index", f"{srr_id}_1.fastq", f"{srr_id}_2.fastq"],
                quantification_file="abundance.tsv",
            )
            species_dirs.append(sample_dir)
        sample_dirs[name] = species_dirs
    return sample_dirs


def populate_progress_db(
    db: ProgressDB,
    cohort: Dict[str, List[str]],
    state_map: Dict[str, Dict[str, List[str]]],
    rng: random.Random,
) -> Dict[str, int]:
    """Fill the disposable DB with the synthetic campaign snapshot."""

    for name, srr_ids in cohort.items():
        db.init_species(name, srr_ids)

    exclusions: Dict[str, List[Dict[str, str]]] = {}
    for name, bucket in state_map.items():
        for state in ("downloading", "downloaded", "quantifying", "quantified", "quarantined"):
            if bucket[state]:
                db.bulk_set_state(name, bucket[state], state)
        for idx, srr_id in enumerate(bucket["failed"]):
            error_text = FAILURE_TEXTS[(idx + len(srr_id)) % len(FAILURE_TEXTS)]
            db.set_state(name, srr_id, "failed", error=error_text)
        # ~1% of the cohort gets a permanent_drop exclusion; a smaller slice
        # gets re_download markers (never blocking eligibility).
        dropped = bucket["pending"][: max(1, len(bucket["pending"]) // 100)] if bucket["pending"] else []
        if dropped:
            exclusions[name] = [
                {"srr_id": srr_id, "reason_code": "permanent_drop", "reason_detail": EXCLUSION_REASON_DETAIL}
                for srr_id in dropped
            ]
        requeue = bucket["pending"][len(dropped) : len(dropped) + max(1, len(bucket["pending"]) // 200)]
        if requeue:
            exclusions[name] = (exclusions.get(name) or []) + [
                {"srr_id": srr_id, "reason_code": "re_download", "reason_detail": "stale partial transfer"}
                for srr_id in requeue
            ]
        if exclusions.get(name):
            db.record_exclusions(name, exclusions[name])

    # Provenance audit rows for a bounded subset (quantified snapshot).
    audit_rows = 0
    for name, bucket in state_map.items():
        for srr_id in bucket["quantified"][:DEFAULT_PROVENANCE_SAMPLE_CAP]:
            db.record_quantification_audit(
                name,
                srr_id,
                status=QUANT_STATUS_CURRENT,
                reason="exact current runtime",
                contract_id=f"bench_{rng.randrange(16**32):032x}",
                observed_amalgkit_version="0.16.60",
                observed_release_tag="v0.16.60",
                observed_source_revision="bench-fixture",
            )
            audit_rows += 1

    # Backdate part of the downloading cohort so reset_stale_downloading()
    # has real work during the retry/backoff benchmark.
    stale_rows = 0
    for name, bucket in state_map.items():
        downloading_ids = bucket["downloading"]
        stale_ids = downloading_ids[: max(1, len(downloading_ids) // 2)] if downloading_ids else []
        if not stale_ids:
            continue
        with sqlite3.connect(str(db.db_path)) as conn:
            conn.executemany(
                "UPDATE samples SET updated_at = datetime('now', '-2 hours') WHERE species = ? AND srr_id = ?",
                [(name, srr_id) for srr_id in stale_ids],
            )
        stale_rows += len(stale_ids)

    return {
        "exclusions": sum(len(v) for v in exclusions.values()),
        "audit_rows": audit_rows,
        "stale_backdated": stale_rows,
    }


def build_campaign_fixture(
    target_dir: Path,
    task_rows: int,
    *,
    species_count: int = DEFAULT_SPECIES_COUNT,
    seed: int = DEFAULT_SEED,
    provenance_sample_cap: int = DEFAULT_PROVENANCE_SAMPLE_CAP,
) -> CampaignFixture:
    """Build the disposable config/metadata/quant fixture and progress DB."""

    rng = random.Random(seed)
    config_dir = target_dir / "config" / "amalgkit"
    metadata_dir = target_dir / "metadata"
    quant_root = target_dir / "amalgkit_root"
    species = synthetic_species_names(species_count)
    write_synthetic_configs(config_dir, species)
    cohort = write_synthetic_metadata(metadata_dir, species, task_rows, rng)
    state_map = assign_states(cohort, rng)
    sample_dirs = write_synthetic_quant_fixtures(
        quant_root, cohort, config_dir, sample_cap=provenance_sample_cap, rng=rng
    )
    db_path = target_dir / "pipeline_progress.db"
    db = ProgressDB(db_path=db_path)
    try:
        populate_progress_db(db, cohort, state_map, rng)
    finally:
        db.close()
    return CampaignFixture(
        target_dir=target_dir,
        config_dir=config_dir,
        metadata_dir=metadata_dir,
        quant_root=quant_root,
        db_path=db_path,
        species=species,
        task_rows=task_rows,
        seed=seed,
        cohort=cohort,
        state_map=state_map,
        sample_dirs=sample_dirs,
    )


def _db_storage_bytes(fixture: CampaignFixture) -> Dict[str, int]:
    """Checkpoint WAL and report on-disk bytes of the disposable DB files."""

    conn = sqlite3.connect(str(fixture.db_path))
    try:
        conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        conn.commit()
    finally:
        conn.close()

    def _size(path: Path) -> int:
        return path.stat().st_size if path.is_file() else 0

    return {
        "db_bytes": _size(fixture.db_path),
        "wal_bytes": _size(fixture.db_path.with_name(fixture.db_path.name + "-wal")),
        "shm_bytes": _size(fixture.db_path.with_name(fixture.db_path.name + "-shm")),
    }


def _fixture_fs_bytes(fixture: CampaignFixture) -> int:
    total = 0
    for path in fixture.target_dir.rglob("*"):
        if path.is_file() and not path.name.startswith("pipeline_progress.db"):
            total += path.stat().st_size
    return total


# ---------- Timing helpers ----------


def _median_of(repeats: int, fn: Callable[[], Any]) -> tuple[float, Any]:
    """Run ``fn`` repeats times; return (median_seconds, last_result)."""

    samples: List[float] = []
    result = None
    for _ in range(repeats):
        start = time.perf_counter()
        result = fn()
        samples.append(time.perf_counter() - start)
    return statistics.median(samples), result


def _rate(operations: int, seconds: float) -> float | None:
    """Operations per second; None when the measurement was too coarse."""

    return round(operations / seconds, 1) if seconds > 0 else None


# ---------- Benchmarks ----------


def bench_discovery(fixture: CampaignFixture, repeats: int) -> dict:
    """Config discovery plus synthetic cohort enumeration throughput."""

    seconds, names = _median_of(repeats, lambda: discover_species_config_names(fixture.config_dir))
    rows = sum(len(ids) for ids in fixture.cohort.values())

    def enumerate_cohort() -> List[Tuple[str, str]]:
        pairs: List[Tuple[str, str]] = []
        for name, srr_ids in fixture.cohort.items():
            for srr_id in srr_ids:
                pairs.append((name, srr_id))
        return pairs

    cohort_seconds, _ = _median_of(repeats, enumerate_cohort)
    return {
        "config_files_scanned": len(names),
        "config_scan_seconds": round(seconds, 6),
        "config_scans_per_second": _rate(repeats * len(names), seconds * repeats),
        "cohort_rows": rows,
        "cohort_enumeration_seconds": round(cohort_seconds, 6),
        "cohort_rows_per_second": _rate(rows, cohort_seconds),
    }


def bench_queue_depth(fixture: CampaignFixture, db: ProgressDB, repeats: int) -> dict:
    """Dashboard queue-depth computation over the full disposable DB."""

    def compute_queue_depth() -> Dict[str, int]:
        totals = db.get_total_counts()
        depth: Dict[str, int] = {}
        for name in fixture.species:
            excluded = db.get_excluded_srr_ids(name, reason_code="permanent_drop")
            pending = db.get_samples(name, "pending")
            depth[name] = len([srr for srr in pending if srr not in excluded])
        queued = sum(depth.values())
        assert queued <= sum(totals.values())
        return depth

    seconds, depth = _median_of(repeats, compute_queue_depth)
    seconds_counts, counts = _median_of(repeats, db.get_counts)
    return {
        "total_task_rows": sum(sum(states.values()) for states in counts.values()),
        "queue_depth_seconds": round(seconds, 6),
        "queue_depth_per_second": _rate(1, seconds),
        "get_counts_seconds": round(seconds_counts, 6),
        "pending_queue_total": sum(depth.values()),
    }


def bench_retry_backoff(
    fixture: CampaignFixture,
    db: ProgressDB,
    *,
    subset_size: int,
) -> dict:
    """Retry scheduling cost: state transitions with error classification."""

    subset: List[Tuple[str, str]] = []
    for name in fixture.species:
        pending = fixture.state_map[name]["pending"]
        take = min(len(pending), max(0, subset_size - len(subset)))
        subset.extend((name, srr_id) for srr_id in pending[:take])
        if len(subset) >= subset_size:
            break

    retry_classes = {
        "transfer_all_sources_failed",
        "extraction_timeout",
        "quantification_timeout",
        "quantification_failed",
        "quantification_exception",
    }
    start = time.perf_counter()
    for idx, (name, srr_id) in enumerate(subset):
        error_text = FAILURE_TEXTS[idx % len(FAILURE_TEXTS)]
        class_name = classify_sample_error(error_text)
        db.set_state(name, srr_id, "failed", error=error_text)
        db.set_state(name, srr_id, "pending", error=None)
        assert class_name in retry_classes
    transition_seconds = time.perf_counter() - start

    stale_start = time.perf_counter()
    reset_rows = db.reset_stale_downloading(stale_seconds=3600)
    stale_seconds = time.perf_counter() - stale_start

    return {
        "transition_pairs": len(subset),
        "transition_seconds": round(transition_seconds, 6),
        "transitions_per_second": _rate(len(subset) * 2, transition_seconds),
        "stale_downloading_resets": reset_rows,
        "stale_scan_seconds": round(stale_seconds, 6),
        "stale_resets_per_second": _rate(max(reset_rows, 1), stale_seconds),
    }


def bench_provenance_io(fixture: CampaignFixture) -> dict:
    """Hash and provenance-sidecar I/O over the synthetic quant fixtures."""

    all_dirs = [path for paths in fixture.sample_dirs.values() for path in paths]
    if not all_dirs:
        return {"sample_dirs": 0}

    hashed_bytes = 0
    hash_start = time.perf_counter()
    digests: List[str] = []
    for sample_dir in all_dirs:
        digest = digest_file(sample_dir / "abundance.tsv")
        assert digest is not None
        digests.append(digest)
        hashed_bytes += (sample_dir / "abundance.tsv").stat().st_size
    hash_seconds = time.perf_counter() - hash_start

    read_start = time.perf_counter()
    payloads = [read_quant_provenance(sample_dir) for sample_dir in all_dirs]
    read_seconds = time.perf_counter() - read_start

    classify_start = time.perf_counter()
    statuses: Dict[str, int] = {}
    for sample_dir, payload in zip(all_dirs, payloads):
        assert payload is not None
        classification = classify_quantification(sample_dir, payload["run_accession"], verify_content=True)
        statuses[classification["status"]] = statuses.get(classification["status"], 0) + 1
    classify_seconds = time.perf_counter() - classify_start

    return {
        "sample_dirs": len(all_dirs),
        "hashed_bytes": hashed_bytes,
        "hash_seconds": round(hash_seconds, 6),
        "hash_megabytes_per_second": _rate(hashed_bytes, hash_seconds) / 1_048_576,
        "sidecar_reads_per_second": _rate(len(all_dirs), read_seconds),
        "classifications_per_second": _rate(len(all_dirs), classify_seconds),
        "classification_status_counts": statuses,
    }


# ---------- One scale ----------


def run_scale_benchmark(
    task_rows: int,
    *,
    target_dir: Path,
    species_count: int,
    seed: int,
    repeats: int,
    provenance_sample_cap: int,
    subset_size: int,
) -> dict:
    """Build one disposable fixture and measure every benchmark at this scale."""

    build_start = time.perf_counter()
    fixture = build_campaign_fixture(
        target_dir,
        task_rows,
        species_count=species_count,
        seed=seed,
        provenance_sample_cap=provenance_sample_cap,
    )
    build_seconds = time.perf_counter() - build_start
    storage = _db_storage_bytes(fixture)
    fs_bytes = _fixture_fs_bytes(fixture)

    db = ProgressDB(db_path=fixture.db_path)
    try:
        discovery = bench_discovery(fixture, repeats)
        queue_depth = bench_queue_depth(fixture, db, repeats)
        retry = bench_retry_backoff(fixture, db, subset_size=min(subset_size, max(task_rows // 5, 20)))
    finally:
        db.close()
    provenance_io = bench_provenance_io(fixture)

    return {
        "task_rows": task_rows,
        "species_count": species_count,
        "fixture_build_seconds": round(build_seconds, 6),
        "storage": {
            **storage,
            "fixture_fs_bytes": fs_bytes,
            "db_bytes_per_row": round(storage["db_bytes"] / task_rows, 2),
        },
        "discovery": discovery,
        "queue_depth": queue_depth,
        "retry_backoff": retry,
        "provenance_io": provenance_io,
    }


# ---------- Rendering ----------


def render_markdown(report: dict) -> str:
    """Human-readable Markdown report with the machine-provenance block."""

    lines: List[str] = []
    lines.append("# Campaign-scale benchmark report")
    lines.append("")
    lines.append(
        "Descriptive observations from a disposable synthetic fixture; these are "
        "local measurements only (not a hosted run) and carry no statistical or "
        "significance claims."
    )
    lines.append("")
    lines.append(
        "| scale (task rows) | build (s) | DB bytes/row | queue-depth (s) | retry transitions/s | classifications/s |"
    )
    lines.append("|---|---|---|---|---|---|")
    for scale in report["scales"]:
        storage = scale["storage"]
        queue = scale["queue_depth"]
        retry = scale["retry_backoff"]
        provenance_io = scale["provenance_io"]
        lines.append(
            f"| {scale['task_rows']:,} | {scale['fixture_build_seconds']} | {storage['db_bytes_per_row']} "
            f"| {queue['queue_depth_seconds']} | {retry['transitions_per_second']} "
            f"| {provenance_io.get('classifications_per_second', 'n/a')} |"
        )
    lines.append("")
    lines.append("## Metric definitions (telemetry)")
    lines.append("")
    lines.append(
        "- **config scans/s**: sorted `amalgkit_*.yaml` glob + marker filter "
        "(`discover_species_config_names`) per second, median of repeats.\n"
        "- **cohort rows/s**: parse of the synthetic per-species metadata TSVs into "
        "(species, run-accession) pairs per second.\n"
        "- **queue depth (s)**: full dashboard pass — `get_counts()`, "
        "`get_total_counts()`, per-species pending `get_samples()` intersected with "
        "permanent-drop exclusions.\n"
        "- **retry transitions/s**: `set_state(failed, error)` + `set_state(pending)` "
        "pairs per second, error text classified via `classify_sample_error`.\n"
        "- **stale resets**: rows reset by `reset_stale_downloading(3600)` on a "
        "backdated downloading cohort.\n"
        "- **hash MB/s / classifications/s**: `digest_file` over the abundance "
        "payloads and full `classify_quantification(verify_content=True)` passes.\n"
        "- **DB bytes/row**: checkpointed SQLite size (db + WAL) per task row."
    )
    lines.append("")
    lines.append("## Budgets (advisory, descriptive)")
    lines.append("")
    lines.append(
        "- Queue-depth dashboard pass should stay comfortably below the 5 s "
        "orchestrator stall-watchdog granularity at every campaign scale.\n"
        "- Phase-1 discovery plus queue accounting must not re-digest the quant "
        "corpus; the reconciliation path classifies by provenance contract "
        "instead (see `ProgressDB.reconcile` docstring).\n"
        "- Retry scheduling is per-sample-commit bound; keep batches bounded by "
        "the measured transitions/s when sizing backoff sweeps."
    )
    lines.append("")
    lines.append("## Machine provenance")
    lines.append("")
    lines.append("```json")
    lines.append(json.dumps(report["provenance"], indent=2, sort_keys=True))
    lines.append("```")
    lines.append("")
    lines.append("## Full measurements (JSON)")
    lines.append("")
    lines.append("```json")
    lines.append(json.dumps({k: v for k, v in report.items() if k != "provenance"}, indent=2, sort_keys=True))
    lines.append("```")
    lines.append("")
    return "\n".join(lines)


def build_report(scale_results: List[dict], *, seed: int, scales: Sequence[int], target_dir: Path) -> dict:
    return {
        "schema": BENCHMARK_SCHEMA,
        "provenance": machine_provenance(seed=seed, scales=scales, target_dir=target_dir),
        "safety": {
            "live_data_root": str(LIVE_DATA_ROOT),
            "live_data_root_read": False,
            "live_db_opened": False,
            "network_access": False,
        },
        "scales": scale_results,
    }


# ---------- CLI ----------


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Campaign-scale benchmark harness (disposable fixtures)")
    parser.add_argument(
        "--target-dir",
        type=Path,
        default=None,
        help=(
            "Benchmark scratch directory "
            "(default: output/campaign_scale_benchmarks/run_<UTC>; never under the live data root)"
        ),
    )
    parser.add_argument(
        "--scales",
        type=str,
        default=",".join(str(s) for s in DEFAULT_SCALES),
        help="Comma-separated task-row scales (default: 1000,10000,100000)",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Fixture RNG seed (default: 20260917)")
    parser.add_argument("--species-count", type=int, default=DEFAULT_SPECIES_COUNT, help="Synthetic species count")
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS, help="Timing repeats per read benchmark")
    parser.add_argument(
        "--provenance-sample-cap",
        type=int,
        default=DEFAULT_PROVENANCE_SAMPLE_CAP,
        help="Max per-species quant fixture sample dirs (bounds hash I/O)",
    )
    parser.add_argument(
        "--keep-fixture",
        action="store_true",
        help="Retain the disposable fixture directory after the run (default: delete)",
    )
    args = parser.parse_args(argv)

    scales = [int(item.strip()) for item in args.scales.split(",") if item.strip()]
    if not scales or any(scale <= 0 for scale in scales):
        parser.error("--scales must be positive integers")

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    target_dir = args.target_dir or REPO_ROOT / "output" / "campaign_scale_benchmarks" / f"run_{stamp}"
    target_dir = guard_target_dir(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    scale_results = []
    for index, task_rows in enumerate(scales):
        scale_dir = target_dir / f"scale_{task_rows}"
        result = run_scale_benchmark(
            task_rows,
            target_dir=scale_dir,
            species_count=args.species_count,
            seed=args.seed,
            repeats=max(1, args.repeats),
            provenance_sample_cap=args.provenance_sample_cap,
            subset_size=200,
        )
        scale_results.append(result)
        print(
            f"[scale {task_rows}] build={result['fixture_build_seconds']}s "
            f"db_bytes_per_row={result['storage']['db_bytes_per_row']} "
            f"queue_depth={result['queue_depth']['queue_depth_seconds']}s "
            f"transitions/s={result['retry_backoff']['transitions_per_second']}"
        )
        if not args.keep_fixture:
            _remove_tree(scale_dir)

    report = build_report(scale_results, seed=args.seed, scales=scales, target_dir=target_dir)
    json_path = target_dir / "report.json"
    md_path = target_dir / "report.md"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path.write_text(render_markdown(report), encoding="utf-8")
    print(f"report: {json_path}")
    print(f"report: {md_path}")
    return 0


def _remove_tree(path: Path) -> None:
    import shutil

    shutil.rmtree(path, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())

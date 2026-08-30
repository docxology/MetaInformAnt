"""Streaming RNA-seq Orchestrator implementation.

This module provides robust orchestration for running the amalgkit pipeline
with ENA-first download streaming, concurrent processing, and automatic cleanup.
"""

from __future__ import annotations

import concurrent.futures
import csv
import errno
import gzip
import json
import logging
import math
import os
import re
import shutil
import stat
import subprocess
import threading
import time
import urllib.request
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

from metainformant.core.utils import logging as log_utils
from metainformant.rna.amalgkit.sra_environment import (
    build_sra_environment,
    ensure_ncbi_settings,
    resolve_data_root,
)
from metainformant.rna.amalgkit.tissue_normalizer import apply_tissue_normalization
from metainformant.rna.core.sample_utils import find_quantification_file
from metainformant.rna.engine.progress_db import ProgressDB
from metainformant.rna.engine.provenance import (
    QUANT_STATUS_CURRENT,
    QUANT_STATUS_VERSION_DRIFT,
    classify_quantification,
    is_current_metadata,
    is_reusable_quantification,
    write_metadata_provenance,
    write_quant_provenance,
)
from metainformant.rna.engine.raw_cleanup import reclaim_sample_raw_inputs
from metainformant.rna.retrieval.ena_downloader import (
    ENADownloader,
    _run_command_in_process_group,
    verify_gzip_integrity,
)

# Determine project root if possible, or use relative paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
NESTED_PROJECT_ROOT = PROJECT_ROOT / "projects/hymenoptera_amalgkit"
CONFIG_DIR = (
    NESTED_PROJECT_ROOT / "config/amalgkit"
    if (NESTED_PROJECT_ROOT / "config/amalgkit").is_dir()
    else PROJECT_ROOT / "config/amalgkit"
)
# Keep the repository-local default for tests and fresh checkouts, while
# allowing production runs to route the large tree to the mounted volume.
LOG_DIR = Path(os.environ.get("AMALGKIT_LOG_DIR", str(PROJECT_ROOT / "output/amalgkit")))


logger = log_utils.get_logger(__name__)

REQUANTIFICATION_POLICIES = frozenset({"preserve", "version-drift", "all"})
DEFAULT_DISCOVERY_WORKERS = 4

_INDEX_CACHE_LOCK = threading.Lock()
_SRA_VALIDATION_SCHEMA = 1


@dataclass(frozen=True)
class PipelineResourceProfile:
    """Bounded resources assigned to one streaming campaign.

    ``workers`` is the number of sample tasks that may be active at once.
    Quantification has its own semaphore: ``quant_slots`` is the maximum
    number of Kallisto/Amalgkit subprocesses, and
    ``quant_threads_per_worker`` is the thread count passed to each one.  The
    download fallback and FASTQ compression are separate child processes, so
    their per-process and concurrent-slot budgets are recorded explicitly
    instead of being hidden hard-coded multipliers. ``validation_slots``
    limits full local gzip scans because they compete with Kallisto for the
    mounted volume. Keeping task workers separate from quant/fallback slots
    lets network I/O feed the CPU-bound stage without silently increasing the
    quantification or external-disk contention budget.
    """

    requested_workers: int
    requested_quant_threads: int
    workers: int
    quant_slots: int
    quant_threads_per_worker: int
    fasterq_threads: int
    fasterq_slots: int
    compression_threads: int
    validation_slots: int
    max_in_flight: int
    host_cpu_count: int

    @property
    def effective_quant_threads(self) -> int:
        """Return the maximum simultaneous Amalgkit quant thread count."""

        return self.quant_slots * self.quant_threads_per_worker

    @property
    def peak_stage_threads(self) -> int:
        """Return the conservative process-thread ceiling across stages."""

        fallback_threads = self.fasterq_slots * max(
            self.fasterq_threads,
            self.compression_threads,
        )
        return max(self.effective_quant_threads, fallback_threads)


def _resource_int(value: Any, default: int, label: str) -> int:
    """Parse a positive resource setting, falling back with a warning."""

    if value is None or value == "":
        return default
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        logger.warning("Invalid %s=%r; using %d", label, value, default)
        return default
    if parsed <= 0:
        logger.warning("Invalid %s=%r; using %d", label, value, default)
        return default
    return parsed


def _format_duration(seconds: int) -> str:
    """Render a timeout compactly for logs and persisted failure messages."""

    if seconds % 3600 == 0:
        return f"{seconds // 3600}h"
    if seconds % 60 == 0:
        return f"{seconds // 60}m"
    return f"{seconds}s"


def build_pipeline_resource_profile(
    workers: int,
    threads: int,
    *,
    quant_slots: Optional[int] = None,
    fasterq_threads: Optional[int] = None,
    fasterq_slots: Optional[int] = None,
    compression_threads: Optional[int] = None,
    validation_slots: Optional[int] = None,
    max_in_flight: Optional[int] = None,
    cpu_count: Optional[int] = None,
) -> PipelineResourceProfile:
    """Build a reproducible, bounded resource profile.

    ``threads`` is a total quantification budget, not a per-worker value.
    ``workers`` controls all active sample tasks (including local validation
    and acquisition), while a quantification semaphore limits Kallisto to at
    most ``min(workers, threads)`` slots.  Stage-specific settings can be
    supplied explicitly or through
    ``AMALGKIT_PIPELINE_QUANT_SLOTS``,
    ``AMALGKIT_PIPELINE_FASTQ_THREADS``,
    ``AMALGKIT_PIPELINE_FASTQ_SLOTS``,
    ``AMALGKIT_PIPELINE_COMPRESSION_THREADS``, and
    ``AMALGKIT_PIPELINE_VALIDATION_SLOTS``, and
    ``AMALGKIT_PIPELINE_MAX_IN_FLIGHT``.
    """

    if workers <= 0:
        raise ValueError("workers must be positive")
    if threads <= 0:
        raise ValueError("threads must be positive")

    host_cpu_count = max(1, cpu_count or (os.cpu_count() or 1))
    # Keep sample-task concurrency independent from the total quantification
    # budget.  This allows more I/O/validation tasks to feed a bounded set of
    # Kallisto processes without turning a low thread budget into hidden
    # quantification oversubscription.
    configured_quant_slots: Any = quant_slots
    if configured_quant_slots is None:
        configured_quant_slots = os.environ.get("AMALGKIT_PIPELINE_QUANT_SLOTS")
    quant_slots = min(
        workers,
        threads,
        _resource_int(configured_quant_slots, min(workers, threads), "AMALGKIT_PIPELINE_QUANT_SLOTS"),
    )
    quant_threads_per_worker = max(1, threads // quant_slots)

    default_fasterq_threads = max(1, min(4, quant_threads_per_worker))
    default_compression_threads = max(1, min(2, quant_threads_per_worker // 2 or 1))
    configured_fasterq: Any = fasterq_threads
    if configured_fasterq is None:
        configured_fasterq = os.environ.get("AMALGKIT_PIPELINE_FASTQ_THREADS")
    configured_fasterq_slots: Any = fasterq_slots
    if configured_fasterq_slots is None:
        configured_fasterq_slots = os.environ.get("AMALGKIT_PIPELINE_FASTQ_SLOTS")
    configured_compression: Any = compression_threads
    if configured_compression is None:
        configured_compression = os.environ.get("AMALGKIT_PIPELINE_COMPRESSION_THREADS")
    configured_validation: Any = validation_slots
    if configured_validation is None:
        configured_validation = os.environ.get("AMALGKIT_PIPELINE_VALIDATION_SLOTS")
    configured_in_flight: Any = max_in_flight
    if configured_in_flight is None:
        configured_in_flight = os.environ.get("AMALGKIT_PIPELINE_MAX_IN_FLIGHT")

    profile = PipelineResourceProfile(
        requested_workers=workers,
        requested_quant_threads=threads,
        workers=workers,
        quant_slots=quant_slots,
        quant_threads_per_worker=quant_threads_per_worker,
        fasterq_threads=_resource_int(
            configured_fasterq,
            default_fasterq_threads,
            "AMALGKIT_PIPELINE_FASTQ_THREADS",
        ),
        fasterq_slots=min(
            workers,
            _resource_int(
                configured_fasterq_slots,
                min(4, workers),
                "AMALGKIT_PIPELINE_FASTQ_SLOTS",
            ),
        ),
        compression_threads=_resource_int(
            configured_compression,
            default_compression_threads,
            "AMALGKIT_PIPELINE_COMPRESSION_THREADS",
        ),
        validation_slots=max(
            1,
            min(
                workers,
                _resource_int(
                    configured_validation,
                    min(4, workers),
                    "AMALGKIT_PIPELINE_VALIDATION_SLOTS",
                ),
            ),
        ),
        max_in_flight=max(
            workers,
            _resource_int(
                configured_in_flight,
                workers * 2,
                "AMALGKIT_PIPELINE_MAX_IN_FLIGHT",
            ),
        ),
        host_cpu_count=host_cpu_count,
    )
    if profile.workers > profile.quant_slots:
        logger.info(
            "Using %d sample-task workers with %d quantification slots; local validation/acquisition "
            "can overlap with the bounded Kallisto budget",
            profile.workers,
            profile.quant_slots,
        )
    if profile.effective_quant_threads > profile.host_cpu_count:
        logger.warning(
            "Quantification budget is %d threads on a host reporting %d CPUs; "
            "expect contention unless this host has an intentionally larger scheduler allocation",
            profile.effective_quant_threads,
            profile.host_cpu_count,
        )
    if profile.peak_stage_threads > profile.host_cpu_count:
        logger.warning(
            "Peak stage process budget is %d threads on a host reporting %d CPUs; "
            "fallback extraction/compression may contend with quantification",
            profile.peak_stage_threads,
            profile.host_cpu_count,
        )
    return profile


def _species_name_from_config(config_name: str) -> str:
    """Derive the output species key from an amalgkit config filename."""
    return config_name.replace("amalgkit_", "").replace(".yaml", "")


def _data_root() -> Path:
    """Return the streaming data root, honoring the external volume override."""

    return resolve_data_root()


def _ensure_ncbi_settings_for_data_root(data_root: Path) -> Path:
    """Compatibility facade for the centralized campaign SRA environment."""

    return ensure_ncbi_settings(data_root)


def _species_work_dir(species_name: str) -> Path:
    """Return the active streaming work directory for a species.

    The external data root contains two historical aliases (Apis and
    Pogonomyrmex).  Prefer a directory that already contains current workflow
    evidence over an empty canonical directory created by an interrupted
    preparation attempt.  This keeps resume runs attached to the existing
    quantifications without copying or deleting data.
    """

    root = _data_root()
    candidates = [species_name]
    if species_name == "apis_mellifera":
        candidates.extend(["apis_mellifera_all", "amellifera"])
    elif species_name == "pogonomyrmex_barbatus":
        candidates.append("pbarbatus")

    # A populated work tree is stronger evidence than directory existence.
    # The marker set is intentionally limited to cheap existence checks; the
    # orchestrator must remain responsive while discovering thousands of rows.
    for candidate in candidates:
        work_dir = root / candidate / "work"
        if work_dir.is_dir() and any(
            (work_dir / marker).exists() for marker in ("metadata", "index", "quant", "merge")
        ):
            return work_dir

    for candidate in candidates:
        species_dir = root / candidate
        if species_dir.is_dir():
            return species_dir / "work"
    return root / species_name / "work"


def _remap_species_workspace_value(value: Any, species_name: str, active_dir_name: str) -> Any:
    """Remap repository-relative paths to an active external alias.

    Configurations remain canonical and reviewable in Git.  When a mounted
    data root uses an established alias, the runtime YAML must point every
    species-specific path (work, logs, FASTQ, merge, filters, and genome) at
    that alias.  Non-path values and shared/reference paths are left intact.
    """

    if isinstance(value, dict):
        return {key: _remap_species_workspace_value(item, species_name, active_dir_name) for key, item in value.items()}
    if isinstance(value, list):
        return [_remap_species_workspace_value(item, species_name, active_dir_name) for item in value]
    if not isinstance(value, str) or active_dir_name == species_name:
        return value

    prefixes = (f"output/amalgkit/{species_name}", f"./output/amalgkit/{species_name}")
    for prefix in prefixes:
        if value == prefix or value.startswith(prefix + "/"):
            suffix = value[len(prefix) :].lstrip("/")
            replacement = f"output/amalgkit/{active_dir_name}"
            return f"{replacement}/{suffix}" if suffix else replacement
    return value


def _runtime_config_path(config_path: Path, species_name: str) -> Path:
    """Return a config path whose workspace paths target the active alias.

    For canonical data directories the checked-in YAML is used directly.  An
    alias-specific YAML is written only under the active external work tree;
    it is deterministic, ignored by the repository, and allows all current
    workflow stages to share the same paths as streaming quantification.
    """

    active_work_dir = _species_work_dir(species_name)
    active_dir_name = active_work_dir.parent.name
    if active_dir_name == species_name:
        return config_path

    try:
        source = yaml.safe_load(config_path.read_text()) or {}
        remapped = _remap_species_workspace_value(source, species_name, active_dir_name)
        payload = yaml.safe_dump(remapped, sort_keys=False)
    except (OSError, yaml.YAMLError) as exc:
        logger.warning("Unable to construct alias runtime config for %s: %s", species_name, exc)
        return config_path

    runtime_dir = active_work_dir / ".metainformant_runtime"
    runtime_path = runtime_dir / config_path.name
    try:
        runtime_dir.mkdir(parents=True, exist_ok=True)
        if not runtime_path.exists() or runtime_path.read_text() != payload:
            runtime_path.write_text(payload)
        logger.info("Using alias-aware runtime config for %s: %s", species_name, runtime_path)
        return runtime_path
    except OSError as exc:
        logger.warning("Unable to persist alias runtime config for %s: %s", species_name, exc)
        return config_path


def _resolve_metadata_path(work_dir: Path) -> Path:
    """Prefer the selected cohort, then fall back to the discovery table.

    ``metadata.tsv`` is the discovery universe; ``metadata_selected.tsv`` is
    the Amalgkit selection output consumed by acquisition and quantification.
    """
    selected_path = work_dir / "metadata/metadata_selected.tsv"
    if selected_path.exists():
        return selected_path
    return work_dir / "metadata/metadata.tsv"


def _filter_metadata_by_size(df: Any, max_gb: float) -> Any:
    """Return metadata rows whose total_bases fit the per-sample size budget."""
    import pandas as pd

    filtered_df = df.copy()
    filtered_df["total_bases"] = pd.to_numeric(filtered_df["total_bases"], errors="coerce").fillna(0)
    max_bases = max_gb * 1e9
    return filtered_df[filtered_df["total_bases"] <= max_bases].copy().sort_values("total_bases")


def _sample_run_column(df: Any) -> str:
    """Resolve the sample accession column used by amalgkit metadata."""
    return "run" if "run" in df.columns else "run_accession"


_RUN_ACCESSION_PATTERN = re.compile(r"^(?:SRR|ERR|DRR)\d+$")

# These are the final raw-input names that Amalgkit can consume directly.  A
# ``.part`` file is intentionally a separate scheduler tier: it is not a
# valid quantification input, but resuming it is cheaper and safer than
# starting a new transfer.
_RAW_INPUT_SUFFIXES = (".fastq.gz", ".fastq", ".fq.gz", ".fq", ".sra", ".sra.cache")
_RAW_PARTIAL_SUFFIXES = tuple(f"{suffix}.part" for suffix in _RAW_INPUT_SUFFIXES)
_FASTQ_INPUT_SUFFIXES = (".fastq.gz", ".fastq", ".fq.gz", ".fq")
_RAW_VALIDATION_SCHEMA = "metainformant.rna.raw_validation.v1"


def _library_layout_is_paired(value: Any) -> bool | None:
    """Resolve an Amalgkit library-layout value without guessing silently."""

    if not isinstance(value, str):
        return None
    normalized = value.strip().lower().replace("_", "-")
    if normalized in {"paired", "paired-end", "pair", "pe"}:
        return True
    if normalized in {"single", "single-end", "se"}:
        return False
    return None


def _raw_fastq_path(sample_dir: Path, accession: str, mate: int | None = None) -> Path | None:
    """Return one non-empty local FASTQ path for an accession and mate."""

    stem = f"{accession}_{mate}" if mate is not None else accession
    for suffix in _FASTQ_INPUT_SUFFIXES:
        candidate = sample_dir / f"{stem}{suffix}"
        try:
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
        except OSError:
            continue
    return None


def _is_accession_fastq_file(path: Path, accession: str) -> bool:
    """Return whether a final FASTQ name belongs to the exact accession.

    Completed mates are useful independently resumable transfer checkpoints
    even though a paired library is not quantifiable until both are present.
    Diagnostic suffixes and ``.part`` files intentionally do not match.
    """

    stems = (accession, f"{accession}_1", f"{accession}_2")
    return any(path.name == f"{stem}{suffix}" for stem in stems for suffix in _FASTQ_INPUT_SUFFIXES)


def _raw_validation_marker(sample_dir: Path, accession: str) -> Path:
    """Return the durable validation marker for one sample's raw inputs.

    The marker lives beside ``getfastq`` under the species work directory, not
    inside the sample directory.  This keeps the sample directory removable
    after quantification while retaining a small audit trail for any future
    reacquisition.  The path is deterministic for both supported
    ``work/getfastq`` and ``work/fastq/getfastq`` layouts.
    """

    return sample_dir.parent.parent / "raw_validation" / f"{accession}.json"


def _raw_validation_marker_matches(
    marker: Path,
    candidates: list[Path],
    expected_paired: bool | None,
) -> bool:
    """Return whether a cached full validation still describes the files."""

    try:
        payload = json.loads(marker.read_text())
    except (OSError, ValueError, TypeError):
        return False
    if payload.get("schema") != _RAW_VALIDATION_SCHEMA:
        return False
    declared_layout = "paired" if expected_paired is True else "single" if expected_paired is False else "unknown"
    recorded_layout = payload.get("declared_library_layout")
    if recorded_layout is not None and recorded_layout != declared_layout:
        return False
    effective_layout = "paired" if len(candidates) == 2 else "single"
    recorded_effective_layout = payload.get("effective_library_layout")
    if recorded_effective_layout is not None and recorded_effective_layout != effective_layout:
        return False
    expected_reconciliation = (
        # Historical witness name retained for compatibility. This records
        # the observed layout, not an invalidating runtime version.
        "amalgkit_0.16.38_paired_metadata_single_fastq"
        if declared_layout == "paired" and effective_layout == "single"
        else "none"
    )
    recorded_reconciliation = payload.get("layout_reconciliation")
    if recorded_reconciliation is not None and recorded_reconciliation != expected_reconciliation:
        return False
    # Older witnesses proved gzip integrity but did not prove that a lone
    # FASTQ under paired metadata was a complete read-1-only stream rather
    # than malformed or non-adjacent paired data. Force one bounded full
    # record scan before such a witness can be reused.
    if (
        expected_reconciliation == "amalgkit_0.16.38_paired_metadata_single_fastq"
        and payload.get("layout_validation") != "full_fastq_record_scan_read1_only"
    ):
        return False
    recorded = payload.get("files")
    if not isinstance(recorded, list) or len(recorded) != len(candidates):
        return False
    by_name = {entry.get("name"): entry for entry in recorded if isinstance(entry, dict)}
    if len(by_name) != len(candidates):
        return False
    for path in candidates:
        entry = by_name.get(path.name)
        if entry is None:
            return False
        try:
            stat = path.stat()
        except OSError:
            return False
        if entry.get("size") != stat.st_size or entry.get("mtime_ns") != stat.st_mtime_ns:
            return False
    return True


def _write_raw_validation_marker(
    marker: Path,
    candidates: list[Path],
    expected_paired: bool | None,
    *,
    layout_validation: str = "payload_stream_and_filename_cardinality",
) -> None:
    """Atomically record a successful local FASTQ validation."""

    files = []
    for path in candidates:
        stat = path.stat()
        files.append({"name": path.name, "size": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    declared_layout = "paired" if expected_paired is True else "single" if expected_paired is False else "unknown"
    effective_layout = "paired" if len(candidates) == 2 else "single"
    payload = {
        "schema": _RAW_VALIDATION_SCHEMA,
        "validation": "full_gzip_stream" if all(path.name.endswith(".gz") for path in candidates) else "nonempty_fastq",
        "declared_library_layout": declared_layout,
        "effective_library_layout": effective_layout,
        "layout_reconciliation": (
            "amalgkit_0.16.38_paired_metadata_single_fastq"
            if declared_layout == "paired" and effective_layout == "single"
            else "none"
        ),
        "layout_validation": layout_validation,
        "files": files,
    }
    # Content-deterministic payload: no wall-clock fields.  Any restart-varying
    # byte in this marker would make a still-valid cached validation look
    # rewritten and force one bounded full record rescan per sample on every
    # resume.  Recency lives in the orchestrator log and the progress DB.
    marker.parent.mkdir(parents=True, exist_ok=True)
    temporary = marker.with_name(f".{marker.name}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, marker)
    except OSError:
        temporary.unlink(missing_ok=True)
        raise


def _local_fastq_payload_is_valid(path: Path) -> bool:
    """Validate a local FASTQ file without treating plain FASTQ as gzip."""

    if path.name.endswith(".gz"):
        return verify_gzip_integrity(path)
    try:
        with path.open("rb") as handle:
            return bool(handle.read(1))
    except OSError:
        return False


def _paths_share_filesystem(first: Path, second: Path) -> bool:
    """Return whether two existing paths reside on the same filesystem."""

    return first.stat().st_dev == second.stat().st_dev


def _sra_source_fingerprint(source: Path) -> dict[str, int]:
    """Return the source fields that make an SRA validation witness current."""

    source_stat = source.stat()
    return {
        "size": source_stat.st_size,
        "mtime_ns": source_stat.st_mtime_ns,
    }


def _sra_validation_marker_path(source: Path) -> Path:
    """Return the hidden validation-witness path for one cached SRA archive."""

    return source.parent / ".metainformant_sra_validation" / f"{source.name}.json"


def _cached_sra_validation_result(source: Path) -> bool | None:
    """Read a current valid/invalid SRA witness, if one exists."""

    marker = _sra_validation_marker_path(source)
    try:
        payload = json.loads(marker.read_text(encoding="utf-8"))
        if payload.get("schema") != _SRA_VALIDATION_SCHEMA:
            return None
        if payload.get("source") != source.name:
            return None
        if payload.get("fingerprint") != _sra_source_fingerprint(source):
            return None
        result = payload.get("valid")
        return result if isinstance(result, bool) else None
    except (OSError, ValueError, TypeError):
        return None


def _write_sra_validation_result(
    source: Path,
    *,
    fingerprint: dict[str, int],
    valid: bool,
    returncode: int,
    detail: str,
) -> None:
    """Persist an atomic source-bound SRA validation witness."""

    marker = _sra_validation_marker_path(source)
    payload = {
        "schema": _SRA_VALIDATION_SCHEMA,
        "source": source.name,
        "fingerprint": fingerprint,
        "valid": valid,
        "validator": "vdb-validate",
        "returncode": returncode,
        "detail": detail[:2000],
    }
    # Content-deterministic payload: no wall-clock fields.  Any restart-varying
    # byte would make a still-valid vdb-validate witness look rewritten and
    # force an expensive re-validation of the cached SRA archive on every
    # resume.  Recency lives in the orchestrator log and the progress DB.
    marker.parent.mkdir(parents=True, exist_ok=True)
    temporary = marker.with_name(f".{marker.name}.{os.getpid()}.{threading.get_ident()}.tmp")
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, marker)
    finally:
        temporary.unlink(missing_ok=True)


def _campaign_cached_sra_path(accession: str) -> Path | None:
    """Return a complete-looking accession file in the campaign VDB cache."""

    cache_dir = _data_root() / ".sra-cache" / "sra"
    for suffix in (".sra", ".sra.cache"):
        candidate = cache_dir / f"{accession}{suffix}"
        try:
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
        except OSError:
            continue
    return None


def _local_sra_path(sample_dir: Path, out_dir: Path, accession: str) -> Path | None:
    """Find an exact local SRA archive in nested, flat, or campaign-cache layout.

    ``prefetch`` may materialize a complete accession as either ``.sra`` or
    ``.sra.cache`` under the campaign-local VDB cache.  Those files are valid
    local databases even though they are outside Amalgkit's per-species
    ``getfastq`` tree.  Prefer them before ENA so an interrupted or previous
    campaign does not pay for the same transfer twice.  Presence is only a
    fast-path hint; ``fasterq-dump`` remains the authoritative extraction and
    layout check, and a failed extraction falls through to normal recovery.
    """

    candidates = [sample_dir / f"{accession}.sra", out_dir / f"{accession}.sra"]
    cached_sra = _campaign_cached_sra_path(accession)
    if cached_sra is not None:
        candidates.append(cached_sra)
    seen: set[Path] = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        try:
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
        except OSError:
            continue
    return None


def _validated_local_fastq_inputs(
    sample_dir: Path,
    accession: str,
    expected_paired: bool | None,
    *,
    assume_payload_valid: bool = False,
) -> list[Path]:
    """Find validated local FASTQs without making an ENA metadata request.

    Existing files are only adopted when the metadata layout is compatible:
    paired libraries use both numbered mates when present, while one
    unnumbered FASTQ is retained for Amalgkit's explicit
    paired-metadata/single-file reconciliation. A lone numbered mate remains
    incomplete. Single libraries require the unnumbered FASTQ, and unknown
    layouts are accepted only when the on-disk naming is unambiguous. Gzip
    streams are fully checked before bypassing the downloader; incomplete
    ``.part`` files therefore remain resumable rather than being mistaken for
    usable input.
    """

    if not sample_dir.is_dir():
        return []

    mate_1 = _raw_fastq_path(sample_dir, accession, 1)
    mate_2 = _raw_fastq_path(sample_dir, accession, 2)
    unnumbered = _raw_fastq_path(sample_dir, accession)

    if expected_paired is True:
        if mate_1 is not None or mate_2 is not None:
            if mate_1 is None or mate_2 is None:
                return []
            candidates = [mate_1, mate_2]
        elif unnumbered is not None:
            candidates = [unnumbered]
        else:
            return []
    elif expected_paired is False:
        if unnumbered is None:
            return []
        candidates = [unnumbered]
    elif mate_1 is not None or mate_2 is not None:
        if mate_1 is None or mate_2 is None:
            return []
        candidates = [mate_1, mate_2]
    elif unnumbered is not None:
        candidates = [unnumbered]
    else:
        return []

    marker = _raw_validation_marker(sample_dir, accession)
    if _raw_validation_marker_matches(marker, candidates, expected_paired):
        return candidates

    layout_validation = "payload_stream_and_filename_cardinality"
    layout_proved_payload = False
    if expected_paired is True and len(candidates) == 1 and candidates[0].name == f"{accession}.fastq.gz":
        # A completed one-file checkpoint may have been left by an interrupted
        # older run before the raw witness was written. Classify it here,
        # before adoption: true interleaved reads are split atomically, while
        # Amalgkit's read-1-only case is accepted only after a complete
        # FASTQ record scan. This path also centralizes classification for new
        # ENA and fasterq outputs, avoiding a second multi-gigabyte scan.
        try:
            normalized = _normalize_interleaved_paired_fastq_layout(
                candidates,
                accession,
                sample_dir,
                expected_paired,
            )
        except (OSError, ValueError) as exc:
            logger.warning(
                "Rejecting one-file paired input for %s after layout validation: %s",
                accession,
                exc,
            )
            return []
        candidates = normalized
        layout_proved_payload = True
        layout_validation = (
            "full_fastq_record_scan_interleaved_split" if len(candidates) == 2 else "full_fastq_record_scan_read1_only"
        )

    if (
        not assume_payload_valid
        and not layout_proved_payload
        and not all(_local_fastq_payload_is_valid(path) for path in candidates)
    ):
        return []
    try:
        _write_raw_validation_marker(
            marker,
            candidates,
            expected_paired,
            layout_validation=layout_validation,
        )
    except OSError as exc:
        # A validation cache is an optimization, never a reason to reject a
        # usable read.  The next run will conservatively revalidate it.
        logger.debug("Could not persist raw validation marker %s: %s", marker, exc)
    return candidates


def _promote_staged_fastq_file(staged: Path, destination: Path) -> None:
    """Promote one validated FASTQ atomically, including across filesystems.

    ``os.replace`` is atomic on one filesystem. Local extraction scratch can
    live on another volume, so an ``EXDEV`` result falls back to an exact-size
    copy into a hidden destination-side temporary file followed by an atomic
    rename. The validated source is removed only after promotion succeeds.
    """

    try:
        os.replace(staged, destination)
        return
    except OSError as exc:
        if exc.errno != errno.EXDEV:
            raise

    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.{threading.get_ident()}.copying")
    try:
        temporary.unlink(missing_ok=True)
        source_size = staged.stat().st_size
        shutil.copy2(staged, temporary)
        if temporary.stat().st_size != source_size:
            raise OSError(f"staged FASTQ size mismatch for {staged.name}")
        os.replace(temporary, destination)
        staged.unlink()
    finally:
        temporary.unlink(missing_ok=True)


def _fastq_pair_key(header: str) -> str | None:
    """Return a mate-insensitive identifier from a FASTQ header."""

    if not header.startswith("@"):
        return None
    fields = header[1:].strip().split()
    if not fields:
        return None
    # ENA's interleaved paired FASTQs use either ``read/1 read/2`` or the
    # Casava ``read 1:N:... read 2:N:...`` convention. The first token is
    # stable for both forms. Do not strip ``.1``/``.2``: SRA uses those as
    # consecutive spot identifiers, not mate suffixes.
    return re.sub(r"/[12]$", "", fields[0])


class _PairedFastqNotInterleaved(ValueError):
    """The stream is structurally FASTQ but adjacent records are not mates."""


def _fastq_header_explicit_mate(header: str) -> int | None:
    """Return an explicit mate number from common FASTQ header conventions."""

    fields = header[1:].strip().split()
    if not fields:
        return None
    slash_suffix = re.search(r"/([12])$", fields[0])
    if slash_suffix:
        return int(slash_suffix.group(1))
    if len(fields) > 1:
        casava = re.match(r"([12]):[YN]:", fields[1])
        if casava:
            return int(casava.group(1))
        sra_suffix = re.search(r"/([12])$", fields[1])
        if sra_suffix:
            return int(sra_suffix.group(1))
    return None


def _validate_single_fastq_stream(source: Path) -> int:
    """Validate a complete one-read stream before layout reconciliation."""

    record_count = 0
    with gzip.open(source, "rt", encoding="utf-8", newline="") as reader:
        while True:
            record = [reader.readline() for _ in range(4)]
            if not record[0]:
                break
            if any(not line for line in record):
                raise ValueError(f"truncated FASTQ record in {source.name}")
            if not record[0].startswith("@"):
                raise ValueError(f"invalid FASTQ header in {source.name}")
            if not record[2].startswith("+"):
                raise ValueError(f"invalid FASTQ separator in {source.name}")
            sequence = record[1].rstrip("\r\n")
            quality = record[3].rstrip("\r\n")
            if not sequence or len(sequence) != len(quality):
                raise ValueError(f"sequence/quality length mismatch in {source.name}")
            if _fastq_header_explicit_mate(record[0]) == 2:
                raise ValueError(f"non-adjacent mate-2 record in one-file paired input {source.name}")
            record_count += 1
    if record_count == 0:
        raise ValueError(f"FASTQ contains no complete records: {source.name}")
    return record_count


def _split_interleaved_fastq(
    source: Path,
    accession: str,
    sample_dir: Path,
) -> list[Path]:
    """Split one interleaved paired FASTQ into numbered mates.

    ENA and some valid SRA archives can expose a paired run as one interleaved
    ``<accession>.fastq.gz`` object even though the portal reports
    ``library_layout=PAIRED``.  Amalgkit/Kallisto require two mate files, so
    convert that representation only after checking every adjacent record
    pair has matching read IDs.  Temporary outputs are never adopted as
    inputs, which keeps interruption resume-safe.
    """

    if not source.is_file() or source.name != f"{accession}.fastq.gz":
        raise ValueError(f"not an unnumbered ENA FASTQ for {accession}: {source.name}")

    output_1 = sample_dir / f"{accession}_1.fastq.gz"
    output_2 = sample_dir / f"{accession}_2.fastq.gz"
    temporary_1 = output_1.with_name(f".{output_1.name}.part")
    temporary_2 = output_2.with_name(f".{output_2.name}.part")
    if output_1.exists() or output_2.exists():
        raise FileExistsError(f"cannot split interleaved FASTQ over existing mate output(s) for {accession}")

    pair_count = 0
    try:
        with gzip.open(source, "rt", encoding="utf-8", newline="") as reader:
            with gzip.open(temporary_1, "wt", compresslevel=1, encoding="utf-8", newline="") as writer_1:
                with gzip.open(temporary_2, "wt", compresslevel=1, encoding="utf-8", newline="") as writer_2:
                    while True:
                        first = [reader.readline() for _ in range(4)]
                        if not first[0]:
                            break
                        second = [reader.readline() for _ in range(4)]
                        if any(not line for line in (*first, *second)):
                            raise ValueError(f"truncated FASTQ record in {source.name}")
                        if not first[0].startswith("@") or not second[0].startswith("@"):
                            raise ValueError(f"invalid FASTQ header in {source.name}")
                        if not first[2].startswith("+") or not second[2].startswith("+"):
                            raise ValueError(f"invalid FASTQ separator in {source.name}")
                        for record in (first, second):
                            sequence = record[1].rstrip("\r\n")
                            quality = record[3].rstrip("\r\n")
                            if not sequence or len(sequence) != len(quality):
                                raise ValueError(f"sequence/quality length mismatch in {source.name}")
                        if _fastq_pair_key(first[0]) != _fastq_pair_key(second[0]):
                            raise _PairedFastqNotInterleaved(
                                f"adjacent FASTQ records are not a mate pair in {source.name}"
                            )
                        first_mate = _fastq_header_explicit_mate(first[0])
                        second_mate = _fastq_header_explicit_mate(second[0])
                        if first_mate == 2 or second_mate == 1:
                            raise ValueError(f"explicit FASTQ mate order is reversed in {source.name}")
                        if first_mate is not None and second_mate is not None and (first_mate, second_mate) != (1, 2):
                            raise ValueError(f"explicit FASTQ mate labels are inconsistent in {source.name}")
                        writer_1.writelines(first)
                        writer_2.writelines(second)
                        pair_count += 1
        if pair_count == 0:
            raise ValueError(f"interleaved FASTQ contains no complete pairs: {source.name}")
        temporary_1.replace(output_1)
        temporary_2.replace(output_2)
        source.unlink()
    except Exception:
        temporary_1.unlink(missing_ok=True)
        temporary_2.unlink(missing_ok=True)
        raise
    return [output_1, output_2]


def _normalize_interleaved_paired_fastq_layout(
    fastq_files: list[Path],
    accession: str,
    sample_dir: Path,
    expected_paired: bool | None,
) -> list[Path]:
    """Normalize a one-file interleaved paired representation."""

    if expected_paired is not True or len(fastq_files) != 1:
        return fastq_files
    source = fastq_files[0]
    if source.name != f"{accession}.fastq.gz":
        return fastq_files
    try:
        return _split_interleaved_fastq(source, accession, sample_dir)
    except _PairedFastqNotInterleaved:
        record_count = _validate_single_fastq_stream(source)
        logger.warning(
            "Run %s is declared paired but contains one validated read-1-only "
            "FASTQ (%d records); retaining the unnumbered file for Amalgkit "
            "0.16.60 layout reconciliation.",
            accession,
            record_count,
        )
        return fastq_files


def _metadata_requires_fastq_input_stats(metadata_path: Path, accession: str) -> bool:
    """Return whether Amalgkit needs read-derived statistics for one accession.

    Amalgkit 0.16.60 infers ``spot_length`` from ``total_bases / total_spots``
    when the public SRA metadata omits the explicit field.  Re-counting a
    validated multi-gigabyte FASTQ merely because ``spot_length`` is missing
    therefore adds substantial I/O/CPU without changing the quantification
    inputs.  Direct FASTQ statistics remain necessary when either count is
    absent or unusable, because ``get_sra_stat`` cannot construct its
    quantification contract without them.
    """

    try:
        with metadata_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if str(row.get("run", "")).strip() != accession:
                    continue
                for column in ("total_spots", "total_bases"):
                    try:
                        value = float(row.get(column, ""))
                    except (TypeError, ValueError):
                        return True
                    if not math.isfinite(value) or value <= 0:
                        return True
                return False
    except (OSError, csv.Error):
        return False
    return False


def _count_fastq_records_and_bases(path: Path) -> tuple[int, int]:
    """Count records and sequence bases in one validated FASTQ payload."""

    opener = gzip.open if path.name.endswith(".gz") else open
    records = 0
    bases = 0
    with opener(path, "rb") as handle:  # type: ignore[arg-type]
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"truncated FASTQ record: {path}")
            records += 1
            bases += len(sequence.rstrip(b"\r\n"))
    return records, bases


def _write_fastq_input_stats(sample_dir: Path, accession: str) -> bool:
    """Write a minimal Amalgkit quant-only stats row for direct FASTQ inputs.

    Direct ENA acquisition does not run Amalgkit's ``getfastq`` stage, so a
    metadata row with zero or missing SRA statistics otherwise fails before
    Kallisto can inspect the already validated FASTQ.  The row is derived from
    the immutable local payload, written atomically, and consumed only as a
    quantification fallback; source metadata is never modified.
    """

    stats_path = sample_dir / "getfastq_stats.tsv"
    try:
        if stats_path.is_file():
            with stats_path.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if str(row.get("run", "")).strip() != accession:
                        continue
                    try:
                        if float(row.get("num_written", "")) > 0:
                            return True
                    except (TypeError, ValueError):
                        break
    except (OSError, csv.Error):
        logger.debug("Unable to inspect existing FASTQ statistics for %s", accession)

    try:
        fastq_files = sorted(
            path
            for path in sample_dir.iterdir()
            if path.is_file() and path.name.startswith(accession) and path.name.endswith(_FASTQ_INPUT_SUFFIXES)
        )
    except OSError:
        return False
    if not fastq_files:
        return False

    try:
        measured = [_count_fastq_records_and_bases(path) for path in fastq_files]
    except (OSError, EOFError, gzip.BadGzipFile, ValueError) as exc:
        logger.warning("Unable to derive FASTQ statistics for %s: %s", accession, exc)
        return False

    paired = [
        result
        for path, result in zip(fastq_files, measured)
        if re.search(r"_[12](?:\.fastq|\.fq)(?:\.gz)?$", path.name)
    ]
    num_written = max((records for records, _ in paired), default=measured[0][0])
    bp_written = sum(bases for _, bases in measured)
    payload = "run\tnum_written\tbp_written\tbp_fastp_in\n" f"{accession}\t{num_written}\t{bp_written}\t{bp_written}\n"
    temporary = stats_path.with_name(f".{stats_path.name}.{os.getpid()}.tmp")
    try:
        sample_dir.mkdir(parents=True, exist_ok=True)
        temporary.write_text(payload, encoding="utf-8")
        os.replace(temporary, stats_path)
    except OSError as exc:
        temporary.unlink(missing_ok=True)
        logger.warning("Unable to persist FASTQ statistics for %s: %s", accession, exc)
        return False
    logger.info(
        "Derived quant-only FASTQ statistics for %s: %d spots, %d bases",
        accession,
        num_written,
        bp_written,
    )
    return True


def _is_valid_run_accession(value: Any) -> bool:
    """Return whether *value* is a concrete SRA/ENA run accession.

    Amalgkit metadata can contain experiment or SRA-primary records without a
    public read-run accession.  Those rows cannot be downloaded or quantified
    by the streaming runner and must not be allowed to become filesystem path
    components (pandas represents blank TSV fields as floating-point NaN).
    """

    return isinstance(value, str) and bool(_RUN_ACCESSION_PATTERN.fullmatch(value.strip()))


def _raw_file_priority(name: str) -> int | None:
    """Return the priority represented by a raw-input filename."""

    if name.endswith(_RAW_INPUT_SUFFIXES):
        return 0
    if name.endswith(_RAW_PARTIAL_SUFFIXES):
        return 1
    return None


def _is_diagnostic_transfer_artifact(path: Path) -> bool:
    """Return whether a failed-transfer artifact should be retained."""

    name = path.name
    return name.endswith(".part") or ".invalid" in name


def _raw_input_state_index(fastq_dir: Path, accessions: set[str]) -> dict[str, tuple[int, int]]:
    """Index raw state and payload size with bounded directory scans.

    The returned tuple is ``(priority, payload_bytes)``: valid FASTQs use
    priority 0, resumable partials priority 1, and fresh acquisitions priority
    2. The size is used only to improve queue order; it never changes the
    provenance or completion contract.
    """

    states = {accession: (2, 0) for accession in accessions}

    def record(path: Path, accession: str) -> None:
        try:
            file_stat = path.stat()
            if not stat.S_ISREG(file_stat.st_mode) or file_stat.st_size <= 0:
                return
            size = file_stat.st_size
        except OSError:
            return
        priority = _raw_file_priority(path.name)
        if priority is None:
            return
        current_priority, current_bytes = states[accession]
        if priority < current_priority:
            states[accession] = (priority, size)
        elif priority == current_priority:
            states[accession] = (current_priority, current_bytes + size)

    try:
        for accession in accessions:
            sample_dir = fastq_dir / accession
            if sample_dir.is_dir():
                for path in sample_dir.iterdir():
                    record(path, accession)

        # Some historical Amalgkit trees keep files directly under
        # ``getfastq``. Read this directory once and recover the accession
        # from the ENA filename instead of rescanning it for every task.
        if fastq_dir.is_dir():
            for path in fastq_dir.iterdir():
                match = re.match(r"^(SRR|ERR|DRR)\d+", path.name)
                if match and match.group(0) in states:
                    record(path, match.group(0))

        # A campaign-local NCBI cache is shared across species and may hold
        # complete SRA databases after an earlier fallback or a separate
        # acquisition pass.  Index exact accession paths without recursively
        # scanning the potentially very large cache tree.
        for accession in accessions:
            cached_sra = _campaign_cached_sra_path(accession)
            if cached_sra is not None:
                record(cached_sra, accession)
    except OSError:
        # Discovery must not abort a campaign because one pre-contract mount
        # entry became unreadable. The normal downloader will report failure.
        return states
    return states


def _raw_input_priority_index(fastq_dir: Path, accessions: set[str]) -> dict[str, int]:
    """Return only the raw-state tier for compatibility with existing callers."""

    return {accession: state[0] for accession, state in _raw_input_state_index(fastq_dir, accessions).items()}


def _sample_raw_input_priority(fastq_dir: Path, srr_id: Any) -> int:
    """Return the acquisition priority for one sample."""

    accession = str(srr_id).strip()
    if not accession:
        return 2
    return _raw_input_priority_index(Path(fastq_dir), {accession}).get(accession, 2)


def _task_has_local_sra(task: Dict[str, Any]) -> bool:
    """Return whether a task has an SRA source that uses the fallback lane."""

    accession = str(task["srr"]).strip()
    fastq_dir = Path(task["fastq_dir"])
    sample_dir = fastq_dir / accession
    expected_paired = task.get("expected_paired")
    unnumbered = _raw_fastq_path(sample_dir, accession)
    mate_1 = _raw_fastq_path(sample_dir, accession, 1)
    mate_2 = _raw_fastq_path(sample_dir, accession, 2)
    if expected_paired is True:
        has_fastq = (mate_1 is not None and mate_2 is not None) or (
            mate_1 is None and mate_2 is None and unnumbered is not None
        )
    elif expected_paired is False:
        has_fastq = unnumbered is not None
    else:
        has_fastq = unnumbered is not None or (mate_1 is not None and mate_2 is not None)
    return not has_fastq and _local_sra_path(sample_dir, fastq_dir, accession) is not None


def _prioritize_tasks_by_raw_state(
    tasks: List[Dict[str, Any]],
    *,
    workers: int | None = None,
    fasterq_slots: int | None = None,
) -> List[Dict[str, Any]]:
    """Prioritize reusable inputs without blocking network workers on SRA I/O.

    Reusable FASTQs remain first. Local SRA databases are admitted in bounded
    fallback-sized chunks, interleaved with partial and new acquisitions, so
    the extraction semaphore does not leave every sample worker waiting on a
    mounted-disk conversion while ENA work could continue independently.
    """

    accessions_by_dir: dict[Path, set[str]] = {}
    for task in tasks:
        fastq_dir = Path(task["fastq_dir"])
        accessions_by_dir.setdefault(fastq_dir, set()).add(str(task["srr"]).strip())
    states_by_dir = {
        fastq_dir: _raw_input_state_index(fastq_dir, accessions) for fastq_dir, accessions in accessions_by_dir.items()
    }
    for task in tasks:
        # Cache the filesystem observation on the transient task record.  A
        # campaign can contain tens of thousands of rows; the index above
        # keeps this lookup O(1) per task after bounded directory scans.
        fastq_dir = Path(task["fastq_dir"])
        accession = str(task["srr"]).strip()
        priority, payload_bytes = states_by_dir[fastq_dir].get(accession, (2, 0))
        task["_raw_input_priority"] = priority
        task["_raw_input_bytes"] = payload_bytes
        task["_uses_local_sra"] = _task_has_local_sra(task)

    def scheduling_key(task: Dict[str, Any]) -> tuple[int, int]:
        priority = task["_raw_input_priority"]
        payload_bytes = task["_raw_input_bytes"]
        if priority == 0:
            # Shortest reusable reads first improves time-to-first-results and
            # prevents a few multi-gigabyte libraries from monopolizing every
            # quant slot while smaller validated samples wait.
            return priority, payload_bytes
        if priority == 1:
            # A larger partial has more acquired data and is usually closer to
            # becoming reusable; finish the most complete resumptions first.
            return priority, -payload_bytes
        return priority, 0

    ordered = sorted(tasks, key=scheduling_key)
    if not ordered:
        return ordered

    # A small initial FASTQ-reuse tier should drain immediately. After that,
    # interleave at most one fallback-sized cache chunk with the non-cache
    # queue. This preserves the preference for downloaded material while
    # keeping network/partial-transfer workers useful during long extraction.
    reusable_fastq = [task for task in ordered if task["_raw_input_priority"] == 0 and not task["_uses_local_sra"]]
    cached_sra = [task for task in ordered if task["_raw_input_priority"] == 0 and task["_uses_local_sra"]]
    non_cached = [task for task in ordered if task["_raw_input_priority"] != 0]
    if not cached_sra or not non_cached:
        return ordered

    worker_budget = max(1, workers or len(ordered))
    fallback_budget = max(1, min(fasterq_slots or 1, worker_budget))
    non_cached_chunk = max(1, worker_budget - fallback_budget)
    result = list(reusable_fastq)
    cache_index = 0
    non_cached_index = 0
    while cache_index < len(cached_sra) or non_cached_index < len(non_cached):
        result.extend(cached_sra[cache_index : cache_index + fallback_budget])
        cache_index += fallback_budget
        result.extend(non_cached[non_cached_index : non_cached_index + non_cached_chunk])
        non_cached_index += non_cached_chunk
    return result


def _filter_valid_run_rows(df: Any, species_name: str, metadata_path: Path) -> Any:
    """Keep downloadable run rows and write an auditable exclusion ledger.

    Rows without a concrete ``run``/``run_accession`` are retained in the
    source metadata provenance only through the ledger; they are excluded from
    the streaming metadata used for batch quantification.  This makes the
    unavailable-data boundary explicit and prevents a missing accession from
    aborting discovery for an entire species.
    """

    run_column = _sample_run_column(df)
    valid_mask = df[run_column].map(_is_valid_run_accession)
    invalid = df.loc[~valid_mask].copy()
    if invalid.empty:
        return df.loc[valid_mask].copy()

    ledger = metadata_path.parent / "validation" / "invalid_run_accessions.tsv"
    ledger.parent.mkdir(parents=True, exist_ok=True)
    invalid.insert(0, "invalid_run_reason", "missing_or_invalid_run_accession")
    invalid.to_csv(ledger, sep="\t", index=False)
    logger.warning(
        "%s: excluding %d metadata rows without a concrete %s accession; " "see %s",
        species_name,
        len(invalid),
        run_column,
        ledger,
    )
    return df.loc[valid_mask].copy()


def _resolve_quant_metadata_path(work_dir: Path) -> Optional[str]:
    """Find the metadata file consumed by amalgkit quant."""
    for metadata_path in [work_dir / "metadata/metadata_selected.tsv", work_dir / "metadata/metadata.tsv"]:
        if metadata_path.exists():
            return str(metadata_path)
    return None


def _resolve_data_config_path(value: Any) -> str:
    """Resolve a config path under the selected external data root."""

    if not isinstance(value, str) or not value:
        return ""
    expanded = os.path.expanduser(os.path.expandvars(value))
    relative = expanded[2:] if expanded.startswith("./") else expanded
    prefix = "output/amalgkit/"
    if relative == "output/amalgkit":
        return str(_data_root())
    if relative.startswith(prefix):
        return str(_data_root() / relative[len(prefix) :])
    return expanded


def _quant_index_candidates(cfg: Dict[str, Any], species_name: str) -> List[str]:
    """Build candidate Kallisto index directories from the current config."""
    quant_index_dir = _resolve_data_config_path(cfg.get("steps", {}).get("quant", {}).get("index_dir", ""))
    genome_dest = _resolve_data_config_path(cfg.get("genome", {}).get("dest_dir", ""))
    genome_index_dir = f"{genome_dest}/index" if genome_dest else ""
    explicit_index_dir = _resolve_data_config_path(cfg.get("genome", {}).get("index_dir", ""))

    candidates = [
        quant_index_dir,
        genome_index_dir,
        explicit_index_dir,
        str(_species_work_dir(species_name) / "index"),
        str(_data_root() / "shared" / "genome" / species_name / "index"),
    ]
    for species in cfg.get("species_list", []):
        candidates.append(str(_data_root() / "shared" / "genome" / species / "index"))
    return [str(Path(candidate).expanduser()) for candidate in candidates if candidate]


def _resolve_index_dir(cfg: Dict[str, Any], species_name: str) -> str:
    """Resolve the first existing index directory containing index files."""
    explicit_index_dir = _resolve_data_config_path(cfg.get("genome", {}).get("index_dir", ""))
    explicit_path = Path(explicit_index_dir).expanduser() if explicit_index_dir else Path()
    if explicit_index_dir and explicit_path.exists():
        return str(explicit_path)

    for candidate in _quant_index_candidates(cfg, species_name):
        if candidate and Path(candidate).exists():
            usable_indexes = [
                index_path
                for index_path in Path(candidate).glob("*.idx")
                if index_path.is_file() and index_path.stat().st_size > 0
            ]
            if usable_indexes:
                return candidate
    return ""


def _normalized_reference_stem(value: Any) -> str:
    """Normalize a metadata target name to Amalgkit's index-file stem."""

    return re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip()).strip("_")


def _find_reference_index(index_dir: Path, stem: str) -> Optional[Path]:
    """Find a non-empty current or compatibility index for one stem."""

    for candidate in (index_dir / f"{stem}.idx", index_dir / f"{stem}_transcripts.idx"):
        try:
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
        except OSError:
            continue
    return None


KALLISTO_DEFAULT_KMER_SIZE = 31


@lru_cache(maxsize=128)
def _kallisto_index_kmer(index_path: str) -> int:
    """Read the k-mer size from a Kallisto index, failing closed to 31.

    Kallisto cannot pseudoalign reads shorter than the index k-mer.  The
    executable project indexes currently use 31-mers, but reading the index
    keeps the eligibility gate correct if a reference is rebuilt with another
    supported k-mer size.  A missing inspection tool must not make discovery
    non-deterministic, so the documented Kallisto default is used as a strict
    fallback and recorded in the eligibility manifest.
    """

    try:
        result = subprocess.run(
            ["kallisto", "inspect", str(index_path)],
            capture_output=True,
            text=True,
            timeout=120,
        )
        combined = f"{result.stdout}\n{result.stderr}"
        match = re.search(r"k-mer length:\s*(\d+)", combined, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
    except (OSError, subprocess.SubprocessError, ValueError) as exc:
        logger.warning("Unable to inspect Kallisto index %s: %s", index_path, exc)
    logger.warning(
        "Using the documented Kallisto k-mer fallback (%d) for %s",
        KALLISTO_DEFAULT_KMER_SIZE,
        index_path,
    )
    return KALLISTO_DEFAULT_KMER_SIZE


def _reference_index_for_kmer(index_dir: Optional[Path], target_names: list[str]) -> Optional[Path]:
    """Choose one readable index representing the active species reference."""

    if index_dir is None:
        return None
    for target in target_names:
        candidate = _find_reference_index(index_dir, _normalized_reference_stem(target))
        if candidate is not None:
            return candidate
    try:
        for candidate in sorted(index_dir.glob("*.idx")):
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
    except OSError:
        return None
    return None


def _filter_kallisto_ineligible_reads(dataframe: Any, kmer_size: int) -> tuple[Any, Any, Any]:
    """Exclude only known reads shorter than the active Kallisto k-mer.

    ``spot_length`` is preferred, with ``total_bases / total_spots`` as the
    documented Amalgkit fallback when the public metadata row omits it.  Rows
    with no usable length estimate are retained for Amalgkit to resolve from
    the downloaded input; silently treating unknown lengths as eligible would
    be less auditable than recording the unresolved value.
    """

    import pandas as pd

    read_length = pd.Series(float("nan"), index=dataframe.index, dtype="float64")
    if "spot_length" in dataframe.columns:
        read_length = pd.to_numeric(dataframe["spot_length"], errors="coerce")
    if {"total_bases", "total_spots"}.issubset(dataframe.columns):
        total_bases = pd.to_numeric(dataframe["total_bases"], errors="coerce")
        total_spots = pd.to_numeric(dataframe["total_spots"], errors="coerce")
        ratio = total_bases.div(total_spots.replace(0, pd.NA))
        read_length = read_length.where(read_length.gt(0), ratio)

    ineligible = read_length.notna() & read_length.lt(int(kmer_size))
    eligible = dataframe.loc[~ineligible].copy()
    excluded = dataframe.loc[ineligible].copy()
    excluded["quant_eligibility_mean_read_length"] = read_length.loc[ineligible].astype(float)
    excluded["quant_eligibility_kmer_size"] = int(kmer_size)
    excluded["quant_eligibility_reason"] = "mean_read_length_below_kallisto_kmer"
    return eligible, excluded, read_length


def _write_quant_eligibility_audit(
    work_dir: Path,
    excluded: Any,
    *,
    kmer_size: int,
    index_path: Optional[Path],
) -> Optional[Path]:
    """Write the fail-closed short-read exclusion table atomically."""

    if excluded is None or excluded.empty:
        return None
    destination = work_dir / "metadata" / "quant_eligibility_exclusions.tsv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    output = excluded.copy()
    output["quant_eligibility_reference_index"] = str(index_path) if index_path else ""
    output["quant_eligibility_kmer_size"] = int(kmer_size)
    temporary = destination.with_name(f".{destination.name}.tmp")
    try:
        output.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, destination)
    except OSError:
        temporary.unlink(missing_ok=True)
        raise
    return destination


def _write_reference_alias_manifest(
    species_name: str,
    work_dir: Path,
    index_dir: Optional[Path],
    targets: list[str],
    aliases: list[dict[str, Any]],
    missing: list[str],
    *,
    kallisto_index: Optional[Path] = None,
    kallisto_kmer_size: Optional[int] = None,
) -> Path:
    """Persist the reference-target audit used by quantification provenance."""

    destination = work_dir / "reference" / "reference_aliases.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Content-deterministic payload: no wall-clock fields.  This manifest's
    # SHA-256 is recorded in per-sample quantification provenance sidecars and
    # re-verified on every resume, so any restart-varying byte (a timestamp)
    # would invalidate every previously quantified sample of the species and
    # quarantine it for redundant re-quantification (observed 2026-08-30:
    # 1,620 apis_mellifera samples re-queued by one restart).  Recency lives
    # in the orchestrator log and the progress DB, not here.
    payload = {
        "schema": "metainformant.rna.reference_aliases.v1",
        "species": species_name,
        "index_dir": str(index_dir) if index_dir else None,
        "targets": targets,
        "aliases": aliases,
        "missing": missing,
        "status": "complete" if not missing else "incomplete",
        "kallisto_index": str(kallisto_index) if kallisto_index else None,
        "kallisto_kmer_size": kallisto_kmer_size,
    }
    temporary = destination.with_name(f".{destination.name}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, destination)
    except OSError:
        temporary.unlink(missing_ok=True)
        raise
    return destination


def _ensure_reference_alias_indexes(
    cfg: Dict[str, Any], species_name: str, target_names: list[str]
) -> tuple[bool, list[str], Optional[Path]]:
    """Ensure every metadata target resolves to a current index.

    Amalgkit resolves an index by metadata target name.  Public metadata can
    contain subspecies or breed labels even when a project intentionally uses
    a species-level reference.  Such mappings must be explicit: this helper
    creates lightweight symlink aliases only for mappings declared in the
    configuration and records them under the external work tree.  Unmapped
    targets fail closed before any reads are downloaded.
    """

    work_dir = _species_work_dir(species_name)
    index_value = _resolve_index_dir(cfg, species_name)
    index_dir = Path(index_value) if index_value else None
    configured = cfg.get("reference_aliases", {})
    aliases_config = configured if isinstance(configured, dict) else {}
    aliases_by_target = {
        str(key).strip().casefold(): str(value).strip()
        for key, value in aliases_config.items()
        if str(key).strip() and str(value).strip()
    }
    targets = sorted({str(target).strip() for target in target_names if str(target).strip()})
    aliases: list[dict[str, Any]] = []
    missing: list[str] = []

    if index_dir is None:
        missing = [f"{target} (no usable index directory)" for target in targets]
    else:
        index_dir.mkdir(parents=True, exist_ok=True)
        for target in targets:
            target_stem = _normalized_reference_stem(target)
            target_index = _find_reference_index(index_dir, target_stem)
            if target_index is not None:
                aliases.append(
                    {
                        "metadata_target": target,
                        "target_stem": target_stem,
                        "reference_stem": target_stem,
                        "index": str(target_index),
                        "mode": "native_index",
                    }
                )
                continue

            source_value = aliases_by_target.get(target.casefold())
            if source_value is None:
                missing.append(target)
                continue
            source_stem = _normalized_reference_stem(source_value)
            source_index = _find_reference_index(index_dir, source_stem)
            if source_index is None:
                missing.append(f"{target} (alias source unavailable: {source_value})")
                continue

            alias_paths = [
                index_dir / f"{target_stem}.idx",
                index_dir / f"{target_stem}_transcripts.idx",
            ]
            conflict = False
            for alias_path in alias_paths:
                try:
                    if alias_path.is_symlink() and alias_path.resolve(strict=False) == source_index.resolve():
                        continue
                    if alias_path.exists():
                        # Never overwrite an existing index with a potentially
                        # different reference.  Require explicit cleanup or a
                        # corrected config instead of silently changing reads.
                        conflict = True
                        continue
                    alias_path.symlink_to(os.path.relpath(source_index, alias_path.parent))
                except OSError as exc:
                    missing.append(f"{target} (unable to create alias: {exc})")
                    conflict = True
            if conflict:
                continue
            aliases.append(
                {
                    "metadata_target": target,
                    "target_stem": target_stem,
                    "reference_stem": source_stem,
                    "source_index": str(source_index),
                    "index_aliases": [str(path) for path in alias_paths],
                    "mode": "explicit_species_level_alias",
                    "note": (
                        "Metadata target lacks an independent configured index; "
                        "reads use the declared reference source."
                    ),
                }
            )

    reference_index = _reference_index_for_kmer(index_dir, targets)
    reference_kmer = _kallisto_index_kmer(str(reference_index)) if reference_index else None
    manifest = _write_reference_alias_manifest(
        species_name,
        work_dir,
        index_dir,
        targets,
        aliases,
        missing,
        kallisto_index=reference_index,
        kallisto_kmer_size=reference_kmer,
    )
    if missing:
        logger.error(
            "%s: reference preflight failed for metadata targets: %s; see %s",
            species_name,
            "; ".join(missing),
            manifest,
        )
        return False, missing, manifest
    if aliases:
        explicit = [entry for entry in aliases if entry.get("mode") == "explicit_species_level_alias"]
        if explicit:
            logger.warning(
                "%s: using %d explicitly declared species-level reference aliases; see %s",
                species_name,
                len(explicit),
                manifest,
            )
    return True, missing, manifest


def _build_quant_command(
    cfg: Dict[str, Any],
    species_name: str,
    batch_index: int,
    threads: int,
    metadata_path: str,
) -> List[str]:
    """Build the amalgkit quant command for one streaming batch."""
    cmd = [
        "amalgkit",
        "quant",
        "--out_dir",
        str(_species_work_dir(species_name)),
        "--metadata",
        metadata_path,
        "--threads",
        str(threads),
        "--batch",
        str(batch_index),
        # A task is scheduled only when no current-method provenance sidecar
        # exists. Replace only explicitly selected re-quantification candidates
        # before writing the new sidecar.
        "--redo",
        "yes",
    ]

    # Keep raw-input deletion behind the provenance gate in
    # ``process_single_sample``.  Amalgkit's own clean-up can run before the
    # current-method sidecar is written, which would make a failed sidecar
    # write indistinguishable from a safely reclaimed sample.  The explicit
    # cleanup helper remains idempotent and removes both reads and redundant
    # ``.safely_removed`` markers after current quantification is proven.
    cmd.extend(["--clean_fastq", "no"])

    index_dir = _resolve_index_dir(cfg, species_name)
    index_dir = _cache_index_directory(index_dir, species_name)
    if index_dir:
        cmd.extend(["--index_dir", index_dir])

    return cmd


def _cache_index_directory(index_dir: str, species_name: str) -> str:
    """Return a validated local copy of a reference index when configured.

    Kallisto repeatedly memory-maps the reference index while the FASTQ files
    remain on the external data volume.  ``AMALGKIT_LOCAL_INDEX_CACHE_DIR`` is
    an opt-in performance lane for hosts where that volume has poor random-read
    latency.  Only ``*.idx`` files are copied, each copy is promoted atomically,
    and a size/mtime manifest invalidates the cache when the source changes.
    With no cache variable, the original path is returned unchanged.
    """

    if not index_dir:
        return index_dir
    cache_value = os.environ.get("AMALGKIT_LOCAL_INDEX_CACHE_DIR", "").strip()
    if not cache_value:
        return index_dir

    source = Path(index_dir).expanduser()
    if not source.is_dir():
        return index_dir
    try:
        source_files = sorted(path for path in source.glob("*.idx") if path.is_file() and path.stat().st_size > 0)
    except OSError as exc:
        logger.warning("Unable to inspect reference index directory %s: %s", source, exc)
        return index_dir
    if not source_files:
        return index_dir

    target = Path(cache_value).expanduser() / species_name / "index"
    manifest_path = target / ".metainformant_index_cache.json"
    expected = {path.name: {"size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns} for path in source_files}

    with _INDEX_CACHE_LOCK:
        try:
            manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
        except (OSError, json.JSONDecodeError, TypeError):
            manifest = {}
        cached_files = manifest.get("files") if isinstance(manifest, dict) else None
        if (
            isinstance(cached_files, dict)
            and manifest.get("source") == str(source)
            and cached_files == expected
            and all(
                (target / name).is_file() and (target / name).stat().st_size == details["size"]
                for name, details in expected.items()
            )
        ):
            return str(target)

        try:
            target.mkdir(parents=True, exist_ok=True)
            for source_file in source_files:
                destination = target / source_file.name
                temporary = target / f".{source_file.name}.{os.getpid()}.tmp"
                shutil.copy2(source_file, temporary)
                os.replace(temporary, destination)
            payload = {"schema": "metainformant.rna.index_cache.v1", "source": str(source), "files": expected}
            temporary_manifest = target / f".{manifest_path.name}.{os.getpid()}.tmp"
            temporary_manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
            os.replace(temporary_manifest, manifest_path)
            logger.info("Using local reference-index cache for %s: %s", species_name, target)
            return str(target)
        except OSError as exc:
            logger.warning("Unable to cache reference indexes for %s; using %s: %s", species_name, source, exc)
            return index_dir


def _build_sample_tasks(
    filtered: Any, srr_col: str, fastq_dir: Path, config_path: Path, species_name: str
) -> List[Dict[str, Any]]:
    """Build task dictionaries from filtered metadata rows."""
    tasks = []
    for i, (_, row) in enumerate(filtered.iterrows()):
        library_layout = row.get("lib_layout")
        if not isinstance(library_layout, str):
            library_layout = row.get("library_layout")
        tasks.append(
            {
                "srr": row[srr_col],
                "batch_idx": i + 1,
                "fastq_dir": fastq_dir,
                "config_path": config_path,
                "species_name": species_name,
                "expected_paired": _library_layout_is_paired(library_layout),
            }
        )
    return tasks


class StreamingPipelineOrchestrator:
    """Orchestrator for streaming ENA download -> Quant processing."""

    def __init__(self, config_dir: Path = CONFIG_DIR, log_dir: Path = LOG_DIR, db_path: Optional[Path] = None):
        self.config_dir = Path(config_dir)
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        # SQLite progress database
        self.db = ProgressDB(db_path) if db_path else ProgressDB()

        # Raw reads are transient once the exact current-contract
        # quantification sidecar has been written.  Keep an explicit opt-out
        # for users who need raw retention, while making the default match the
        # external-volume storage contract.
        self.reclaim_raw_after_quant = os.environ.get(
            "AMALGKIT_RECLAIM_RAW_AFTER_QUANT", "yes"
        ).strip().lower() not in {"0", "false", "no", "off"}
        self.requantification_policy = os.environ.get("AMALGKIT_REQUANTIFICATION_POLICY", "preserve").strip().lower()
        if self.requantification_policy not in REQUANTIFICATION_POLICIES:
            raise ValueError(
                f"Invalid AMALGKIT_REQUANTIFICATION_POLICY={self.requantification_policy!r}; "
                f"expected one of {sorted(REQUANTIFICATION_POLICIES)}"
            )

        # The launcher checks the host volume before starting, but a long
        # campaign can still consume host-side cache/temp space after launch.
        # Keep the same floor in the worker loop so new downloads pause before
        # macOS reaches an emergency-full state.  On Linux/CI, a nonexistent
        # macOS path disables this host-specific guard without affecting the
        # external data-root throttle below.
        self.system_root = Path(os.environ.get("AMALGKIT_SYSTEM_ROOT", "/System/Volumes/Data")).expanduser()
        try:
            self.min_system_free_gb = max(
                0.0,
                float(os.environ.get("AMALGKIT_MIN_SYSTEM_FREE_GB", "4")),
            )
        except ValueError:
            logger.warning("Invalid AMALGKIT_MIN_SYSTEM_FREE_GB; using 4 GiB")
            self.min_system_free_gb = 4.0

        # Keep the mounted-data safety floor explicit and configurable. A
        # temporary test directory or a deliberately local scratch root may
        # live on the system volume rather than the production data volume;
        # callers can set this to zero for such isolated runs without
        # weakening the production default.
        try:
            self.min_external_free_gb = max(
                0.0,
                float(os.environ.get("AMALGKIT_MIN_EXTERNAL_FREE_GB", "8")),
            )
        except ValueError:
            logger.warning("Invalid AMALGKIT_MIN_EXTERNAL_FREE_GB; using 8 GiB")
            self.min_external_free_gb = 8.0

        def _duration_setting(name: str, default: int) -> int:
            raw = os.environ.get(name, str(default)).strip()
            try:
                value = int(raw)
            except ValueError:
                logger.warning("Invalid %s=%r; using %d seconds", name, raw, default)
                return default
            if value <= 0:
                logger.warning("Invalid %s=%r; using %d seconds", name, raw, default)
                return default
            return value

        def _positive_int_setting(name: str, default: int) -> int:
            raw = os.environ.get(name, str(default)).strip()
            try:
                value = int(raw)
            except ValueError:
                logger.warning("Invalid %s=%r; using %d", name, raw, default)
                return default
            if value <= 0:
                logger.warning("Invalid %s=%r; using %d", name, raw, default)
                return default
            return value

        def _nonnegative_int_setting(name: str, default: int) -> int:
            raw = os.environ.get(name, str(default)).strip()
            try:
                value = int(raw)
            except ValueError:
                logger.warning("Invalid %s=%r; using %d", name, raw, default)
                return default
            if value < 0:
                logger.warning("Invalid %s=%r; using %d", name, raw, default)
                return default
            return value

        def _nonnegative_float_setting(name: str, default: float) -> float:
            raw = os.environ.get(name, str(default)).strip()
            try:
                value = float(raw)
            except ValueError:
                logger.warning("Invalid %s=%r; using %.2f", name, raw, default)
                return default
            if value < 0:
                logger.warning("Invalid %s=%r; using %.2f", name, raw, default)
                return default
            return value

        self.download_timeout_seconds = _duration_setting("AMALGKIT_PIPELINE_DOWNLOAD_TIMEOUT_SECONDS", 7200)
        self.download_speed_limit_bytes = _positive_int_setting("AMALGKIT_PIPELINE_DOWNLOAD_SPEED_LIMIT_BYTES", 1024)
        self.download_speed_time_seconds = _duration_setting("AMALGKIT_PIPELINE_DOWNLOAD_SPEED_TIME_SECONDS", 600)
        self.download_retries = _nonnegative_int_setting("AMALGKIT_PIPELINE_ENA_RETRIES", 5)
        self.download_integrity_retries = _nonnegative_int_setting("AMALGKIT_PIPELINE_ENA_INTEGRITY_RETRIES", 1)
        self.invalid_witness_max_bytes = (
            _nonnegative_int_setting("AMALGKIT_PIPELINE_INVALID_WITNESS_MAX_MB", 16) * 1024 * 1024
        )
        self.download_retry_delay_seconds = _duration_setting("AMALGKIT_PIPELINE_ENA_RETRY_DELAY_SECONDS", 5)
        self.ena_api_retries = _nonnegative_int_setting("AMALGKIT_PIPELINE_ENA_API_RETRIES", 2)
        self.ena_api_retry_delay_seconds = _duration_setting("AMALGKIT_PIPELINE_ENA_API_RETRY_DELAY_SECONDS", 2)
        self.ncbi_prefetch_first = os.environ.get("AMALGKIT_NCBI_PREFETCH_FIRST", "no").strip().lower() in {
            "1",
            "true",
            "yes",
            "on",
        }
        self.ncbi_prefetch_timeout_seconds = _duration_setting(
            "AMALGKIT_PIPELINE_NCBI_PREFETCH_TIMEOUT_SECONDS", self.download_timeout_seconds
        )
        self.fasterq_timeout_seconds = _duration_setting("AMALGKIT_PIPELINE_FASTQ_TIMEOUT_SECONDS", 7200)
        self.sra_validation_timeout_seconds = _duration_setting("AMALGKIT_PIPELINE_SRA_VALIDATE_TIMEOUT_SECONDS", 600)
        self.compression_timeout_seconds = _duration_setting("AMALGKIT_PIPELINE_COMPRESSION_TIMEOUT_SECONDS", 1800)
        self.quant_timeout_seconds = _duration_setting("AMALGKIT_PIPELINE_QUANT_TIMEOUT_SECONDS", 7200)

        local_quant_cache = os.environ.get("AMALGKIT_LOCAL_QUANT_SCRATCH_DIR", "").strip()
        self.local_quant_scratch_dir = Path(local_quant_cache).expanduser() if local_quant_cache else None
        local_fasterq_cache = os.environ.get("AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR", "").strip()
        self.local_fasterq_scratch_dir = Path(local_fasterq_cache).expanduser() if local_fasterq_cache else None
        # NCBI documents an 8x--10x accession-size rule of thumb for combined
        # fasterq output and temporary space. Local extraction reserves that
        # full estimate plus one source-sized ephemeral copy when the
        # authoritative SRA is on another volume.
        self.local_fasterq_scratch_multiplier = _nonnegative_float_setting(
            "AMALGKIT_LOCAL_FASTERQ_SCRATCH_MULTIPLIER", 10.0
        )
        self.local_fasterq_scratch_min_gb = _nonnegative_float_setting("AMALGKIT_LOCAL_FASTERQ_SCRATCH_MIN_GB", 2.0)
        self.local_fasterq_remote_reserve_gb = _nonnegative_float_setting(
            "AMALGKIT_LOCAL_FASTERQ_REMOTE_RESERVE_GB", 16.0
        )
        self.local_fastq_stage = os.environ.get("AMALGKIT_LOCAL_FASTQ_STAGE", "").strip().lower() in {
            "1",
            "true",
            "yes",
            "on",
        }
        try:
            self.local_quant_scratch_reserve_gb = max(
                0.0,
                float(
                    os.environ.get(
                        "AMALGKIT_LOCAL_QUANT_SCRATCH_RESERVE_GB",
                        os.environ.get("AMALGKIT_MIN_SYSTEM_FREE_GB", "4"),
                    )
                ),
            )
        except ValueError:
            logger.warning("Invalid AMALGKIT_LOCAL_QUANT_SCRATCH_RESERVE_GB; using 4 GiB")
            self.local_quant_scratch_reserve_gb = 4.0

        # Lock for filesystem operations (symlink cleanup, directory creation)
        self._fs_lock = threading.Lock()
        # Local FASTQ staging and optional fasterq temporary storage can share
        # the host volume. Reserve bytes through one ledger so concurrent
        # stages cannot all observe the same free-space headroom and overfill
        # that volume during a burst of large libraries.
        self._local_scratch_reservation_lock = threading.Lock()
        self._local_scratch_reserved_bytes = 0
        # FASTQ staging occupies the local volume for the entire quantification
        # workspace lifetime, not just while the copy is in progress. Keep one
        # reservation per sample workspace so concurrent quant slots cannot
        # collectively cross the configured host-volume reserve.
        self._local_fastq_stage_reservations: Dict[Path, int] = {}
        # A campaign may use more sample-task workers than quantification
        # slots so local validation and acquisition can feed Kallisto.  This
        # semaphore is installed by ``run_all`` and remains absent for small
        # direct-method calls and unit tests.
        self._quant_semaphore: threading.BoundedSemaphore | None = None
        # SRA extraction is usually external-volume I/O, not network I/O. Keep
        # its concurrency independent from sample workers so a large cache of
        # local archives cannot starve Kallisto or saturate the mounted disk.
        self._fasterq_semaphore: threading.BoundedSemaphore | None = None
        self._raw_validation_semaphore: threading.BoundedSemaphore | None = None

        # Direct calls to ``download_fastq`` (including recovery utilities and
        # tests) still receive an explicit bounded fallback profile.  A full
        # campaign replaces this with the requested profile at the beginning
        # of ``run_all``.
        self._resource_profile = build_pipeline_resource_profile(1, 1)

        # Configure logging to file as well
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        log_file = self.log_dir / f"streaming_orchestrator_{timestamp}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        logger.addHandler(file_handler)

    def query_ena_fastq_urls(self, srr_id: str) -> List[str]:
        """Query ENA API for direct FASTQ download URLs."""
        url = (
            f"https://www.ebi.ac.uk/ena/portal/api/filereport"
            f"?accession={srr_id}&result=read_run"
            f"&fields=run_accession,fastq_ftp&format=tsv"
        )
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:  # nosec B310
                lines = resp.read().decode().strip().split("\n")
                if len(lines) < 2:
                    return []
                ftp_field = lines[1].split("\t")[1] if "\t" in lines[1] else ""
                if not ftp_field:
                    return []
                # Ensure protocol is https
                return [f"https://{p}" if not p.startswith("http") else p for p in ftp_field.split(";") if p]
        except Exception as e:
            logger.warning(f"ENA query failed for {srr_id}: {e}")
            return []

    def _validate_local_sra_archive(
        self,
        source: Path,
        srr_id: str,
    ) -> bool | None:
        """Validate one unchanged local SRA once and cache the result.

        ``False`` is authoritative evidence that the current source fingerprint
        is invalid. ``None`` means validation could not run or timed out, so
        fasterq-dump remains the fallback authority rather than rejecting a
        potentially usable archive.
        """

        cached_result = _cached_sra_validation_result(source)
        if cached_result is not None:
            logger.info(
                "Reusing current SRA validation witness for %s: %s",
                srr_id,
                "valid" if cached_result else "invalid",
            )
            return cached_result

        validator = shutil.which("vdb-validate")
        if validator is None:
            logger.warning(
                "vdb-validate is unavailable for local SRA preflight %s; " "continuing with fasterq-dump validation",
                srr_id,
            )
            return None

        try:
            before = _sra_source_fingerprint(source)
            result = _run_command_in_process_group(
                [validator, str(source)],
                timeout=self.sra_validation_timeout_seconds,
                env=build_sra_environment(_data_root()),
            )
            after = _sra_source_fingerprint(source)
        except subprocess.TimeoutExpired:
            logger.warning(
                "SRA preflight validation timed out for %s (>%s); " "continuing with fasterq-dump validation",
                srr_id,
                _format_duration(self.sra_validation_timeout_seconds),
            )
            return None
        except OSError as exc:
            logger.warning(
                "SRA preflight validation could not inspect %s; " "continuing with fasterq-dump validation: %s",
                srr_id,
                exc,
            )
            return None

        if before != after:
            logger.warning(
                "SRA source changed during preflight validation for %s; "
                "discarding the witness and continuing with fasterq-dump",
                srr_id,
            )
            return None

        detail = (result.stderr or result.stdout or "").strip()
        valid = result.returncode == 0
        try:
            _write_sra_validation_result(
                source,
                fingerprint=after,
                valid=valid,
                returncode=result.returncode,
                detail=detail,
            )
        except OSError as exc:
            logger.warning(
                "Could not persist SRA validation witness for %s: %s",
                srr_id,
                exc,
            )
        if valid:
            logger.info("Validated local SRA archive for %s before extraction.", srr_id)
        else:
            logger.warning(
                "Skipping invalid local SRA archive for %s after vdb-validate " "exit %d: %s",
                srr_id,
                result.returncode,
                detail[:300] or "no validator detail",
            )
        return valid

    def _extract_sra_to_fastq_bounded(
        self,
        source: Path | str,
        sample_dir: Path,
        out_dir: Path,
        srr_id: str,
        expected_paired: bool | None,
    ) -> bool:
        """Validate and extract one SRA run under the fallback-slot budget."""

        semaphore = self._fasterq_semaphore
        if semaphore is None:
            if isinstance(source, Path) and self._validate_local_sra_archive(source, srr_id) is False:
                return False
            return self._extract_sra_to_fastq(
                source,
                sample_dir,
                out_dir,
                srr_id,
                expected_paired,
            )
        with semaphore:
            if isinstance(source, Path) and self._validate_local_sra_archive(source, srr_id) is False:
                return False
            return self._extract_sra_to_fastq(
                source,
                sample_dir,
                out_dir,
                srr_id,
                expected_paired,
            )

    def _extract_sra_to_fastq(
        self,
        source: Path | str,
        sample_dir: Path,
        out_dir: Path,
        srr_id: str,
        expected_paired: bool | None,
    ) -> bool:
        """Extract one local or remote SRA source into validated FASTQs.

        Extraction happens in a per-run staging directory.  Existing sample
        outputs are therefore never overwritten by a partially completed
        ``fasterq-dump`` invocation.  A local ``.sra`` archive remains in
        place until quantification succeeds, so interruption is resumable and
        the current-contract cleanup gate can reclaim it later.
        """

        temp_root = out_dir / ".fasterq_tmp" / srr_id
        staging_dir = temp_root / "output"
        fasterq_temp_dir = temp_root / "temp"
        fasterq_source: Path | str = source
        local_temp_root: Path | None = None
        local_temp_reservation = 0
        source_bytes: int | None = None
        # An interrupted producer can leave only this accession's ephemeral
        # workspace behind. The extraction semaphore is already held here, so
        # no current invocation can own the same staging path. Remove that
        # exact workspace before retrying; never touch the authoritative SRA
        # cache, resumable ENA partials, or promoted sample directory.
        if temp_root.is_symlink() or (temp_root.exists() and not temp_root.is_dir()):
            temp_root.unlink()
        elif temp_root.exists():
            shutil.rmtree(temp_root)
            logger.info("Removed stale fasterq staging workspace for %s: %s", srr_id, temp_root)
        staging_dir.mkdir(parents=True, exist_ok=True)
        if self.local_fasterq_scratch_dir is not None:
            if isinstance(source, Path):
                try:
                    source_bytes = source.stat().st_size
                except OSError as exc:
                    logger.warning(
                        "Unable to size local SRA source for SSD fasterq scratch %s: %s",
                        srr_id,
                        exc,
                    )
                else:
                    # Reserve the documented output/temp workspace estimate
                    # plus one ephemeral source copy. Copying a cross-volume
                    # SRA archive onto the same SSD prevents fasterq-dump's
                    # repeated/random source reads from contending with ENA,
                    # validation, promotion, and quantification on the mounted
                    # data volume.
                    local_temp_reservation = max(
                        int(self.local_fasterq_scratch_min_gb * (1024**3)),
                        math.ceil(source_bytes * (self.local_fasterq_scratch_multiplier + 1.0)),
                    )
            else:
                # A bare accession has no local source size before the SRA
                # Toolkit resolves it. Use a fixed conservative reservation;
                # fasterq-dump retains its own enabled output/temp size check.
                local_temp_reservation = int(self.local_fasterq_remote_reserve_gb * (1024**3))

            if local_temp_reservation and self._reserve_local_scratch(
                self.local_fasterq_scratch_dir,
                local_temp_reservation,
                purpose=f"Local fasterq scratch {srr_id}",
            ):
                try:
                    local_temp_root = self.local_fasterq_scratch_dir / srr_id
                    if local_temp_root.is_symlink() or (local_temp_root.exists() and not local_temp_root.is_dir()):
                        local_temp_root.unlink()
                    elif local_temp_root.exists():
                        shutil.rmtree(local_temp_root)
                    staging_dir = local_temp_root / "output"
                    fasterq_temp_dir = local_temp_root / "temp"
                    staging_dir.mkdir(parents=True, exist_ok=True)
                    fasterq_temp_dir.mkdir(parents=True, exist_ok=True)
                    if isinstance(source, Path) and source_bytes is not None:
                        copying_source: Path | None = None
                        try:
                            if not _paths_share_filesystem(source, local_temp_root):
                                source_dir = local_temp_root / "source"
                                source_dir.mkdir(parents=True, exist_ok=True)
                                staged_source = source_dir / source.name
                                copying_source = source_dir / (
                                    f".{source.name}.{os.getpid()}." f"{threading.get_ident()}.copying"
                                )
                                shutil.copyfile(source, copying_source)
                                copied_bytes = copying_source.stat().st_size
                                if copied_bytes != source_bytes:
                                    raise OSError(
                                        "local SRA source copy size mismatch: "
                                        f"expected {source_bytes}, observed {copied_bytes}"
                                    )
                                os.replace(copying_source, staged_source)
                                fasterq_source = staged_source
                                logger.info(
                                    "Staged SRA source for %s onto local fasterq "
                                    "scratch: %s (%.2f GiB); authoritative source "
                                    "retained at %s",
                                    srr_id,
                                    staged_source,
                                    source_bytes / (1024**3),
                                    source,
                                )
                        except OSError as exc:
                            if copying_source is not None:
                                copying_source.unlink(missing_ok=True)
                            fasterq_source = source
                            logger.warning(
                                "Could not stage SRA source for %s onto local "
                                "fasterq scratch; extracting from the authoritative "
                                "source instead: %s",
                                srr_id,
                                exc,
                            )
                except OSError as exc:
                    logger.warning(
                        "Unable to prepare local fasterq scratch for %s; " "using mounted-volume temporary storage: %s",
                        srr_id,
                        exc,
                    )
                    if local_temp_root is not None:
                        shutil.rmtree(local_temp_root, ignore_errors=True)
                    self._release_local_scratch(
                        local_temp_reservation,
                        purpose=f"Local fasterq scratch {srr_id}",
                    )
                    local_temp_root = None
                    local_temp_reservation = 0
                    fasterq_temp_dir = temp_root / "temp"
                else:
                    logger.info(
                        "Using local fasterq extraction scratch for %s: %s "
                        "(reserved %.2f GiB; ephemeral output and temp are "
                        "local; authoritative SRA/cache and promoted FASTQs "
                        "remain on %s)",
                        srr_id,
                        fasterq_temp_dir,
                        local_temp_reservation / (1024**3),
                        out_dir,
                    )
            else:
                local_temp_reservation = 0
        try:
            fasterq_temp_dir.mkdir(parents=True, exist_ok=True)
            fqd_cmd = [
                "fasterq-dump",
                str(fasterq_source),
                "--outdir",
                str(staging_dir),
                "--temp",
                str(fasterq_temp_dir),
                "--split-3",
                "--threads",
                str(self._resource_profile.fasterq_threads),
                "--skip-technical",
            ]
            fqd_result = _run_command_in_process_group(
                fqd_cmd,
                timeout=self.fasterq_timeout_seconds,
                env=build_sra_environment(_data_root()),
            )
            if fqd_result.returncode != 0:
                logger.error(
                    "fasterq-dump failed for %s (exit %d): %s",
                    srr_id,
                    fqd_result.returncode,
                    fqd_result.stderr[:200],
                )
                return False

            uncompressed = sorted(staging_dir.glob("*.fastq"))
            if not uncompressed:
                logger.error("fasterq-dump produced no FASTQ files for %s", srr_id)
                return False
            compressed_outputs: list[Path] = []
            for fastq in uncompressed:
                compression = _run_command_in_process_group(
                    [
                        "pigz",
                        "-f",
                        "-p",
                        str(self._resource_profile.compression_threads),
                        str(fastq),
                    ],
                    timeout=self.compression_timeout_seconds,
                )
                gz_path = fastq.with_suffix(".fastq.gz")
                if compression.returncode != 0 or not gz_path.is_file() or gz_path.stat().st_size == 0:
                    logger.error("FASTQ compression failed for %s (%s)", srr_id, fastq.name)
                    return False
                compressed_outputs.append(gz_path)

            staged_inputs = _validated_local_fastq_inputs(staging_dir, srr_id, expected_paired)
            if not staged_inputs:
                logger.error(
                    "fasterq-dump output for %s did not match the declared library layout",
                    srr_id,
                )
                return False

            for staged in staged_inputs:
                destination = sample_dir / staged.name
                if destination.exists():
                    if _local_fastq_payload_is_valid(destination):
                        logger.info(
                            "Retaining existing valid FASTQ while promoting complementary "
                            "fallback output for %s: %s",
                            srr_id,
                            destination,
                        )
                        staged.unlink(missing_ok=True)
                        continue
                    invalid = destination.with_name(f"{destination.name}.invalid")
                    counter = 1
                    while invalid.exists():
                        invalid = destination.with_name(f"{destination.name}.invalid.{counter}")
                        counter += 1
                    destination.replace(invalid)
                _promote_staged_fastq_file(staged, destination)

            # The staged gzip payloads were fully validated above and a
            # cross-filesystem promotion verifies the copied byte count before
            # its atomic destination-side rename. Persist the destination
            # layout witness without a second multi-gigabyte gzip scan.
            return bool(
                _validated_local_fastq_inputs(
                    sample_dir,
                    srr_id,
                    expected_paired,
                    assume_payload_valid=True,
                )
            )
        except subprocess.TimeoutExpired:
            logger.error(
                "fasterq-dump timeout for %s (>%s)",
                srr_id,
                _format_duration(self.fasterq_timeout_seconds),
            )
            return False
        except OSError as exc:
            logger.error("SRA extraction filesystem error for %s: %s", srr_id, exc)
            return False
        finally:
            shutil.rmtree(temp_root, ignore_errors=True)
            if local_temp_root is not None:
                shutil.rmtree(local_temp_root, ignore_errors=True)
            if local_temp_reservation:
                self._release_local_scratch(
                    local_temp_reservation,
                    purpose=f"Local fasterq scratch {srr_id}",
                )

    def download_fastq(
        self,
        srr_id: str,
        out_dir: Path,
        expected_paired: bool | None = None,
        *,
        defer_ncbi_fallback: bool = False,
    ) -> bool | None:
        """Download FASTQ files using ENADownloader.

        Uses the ENA Portal API to discover FASTQ URLs and downloads
        via curl with retries. When ``defer_ncbi_fallback`` is true, a run
        without a local SRA archive returns ``None`` after ENA exhaustion so a
        dedicated fallback executor can process it without occupying an ENA
        worker.

        Args:
            srr_id: SRA accession ID (e.g., SRR12345)
            out_dir: Base directory for downloads (sample subdir created automatically)

        Returns:
            True when reusable FASTQs are ready, False after all enabled
            acquisition routes fail, or None when NCBI acquisition was
            deliberately deferred to the scheduler's dedicated fallback lane.
        """
        sample_dir = out_dir / srr_id
        # Remove broken symlinks in the path chain (amalgkit preprocessing creates
        # symlinks like work/getfastq -> /local/path that are invalid inside Docker).
        # Use a lock because 28+ threads may race on the same parent directory.
        with self._fs_lock:
            for check_path in [out_dir, sample_dir]:
                if check_path.is_symlink() and not check_path.exists():
                    logger.warning(f"Removing broken symlink: {check_path} -> {os.readlink(check_path)}")
                    check_path.unlink()
            sample_dir.mkdir(parents=True, exist_ok=True)

        # A valid local sample is the cheapest and most reproducible input.
        # Do this before disk throttling and ENA URL discovery: large active
        # campaigns otherwise spend every worker slot waiting on the ENA API
        # even when the FASTQs are already mounted locally.
        if self._raw_validation_semaphore is None:
            local_fastqs = _validated_local_fastq_inputs(sample_dir, srr_id, expected_paired)
        else:
            with self._raw_validation_semaphore:
                local_fastqs = _validated_local_fastq_inputs(sample_dir, srr_id, expected_paired)
        if local_fastqs:
            logger.info(
                "Reusing %d validated local FASTQ file(s) for %s; skipping ENA acquisition.",
                len(local_fastqs),
                srr_id,
            )
            return True

        local_sra = _local_sra_path(sample_dir, out_dir, srr_id)
        local_sra_attempted = False
        if local_sra is not None:
            logger.info(
                "Extracting existing local SRA archive for %s before any network acquisition: %s",
                srr_id,
                local_sra,
            )
            local_sra_attempted = True
            if self._extract_sra_to_fastq_bounded(
                local_sra,
                sample_dir,
                out_dir,
                srr_id,
                expected_paired,
            ):
                return True
            logger.warning(
                "Local SRA validation/extraction did not produce usable FASTQ files "
                "for %s; continuing with ENA/NCBI recovery.",
                srr_id,
            )

        # Disk space throttle: Prevent starting new enormous downloads if
        # either the external data volume or the host system volume is nearing
        # capacity.  Quantification workers can reclaim external FASTQ space
        # while a paused downloader waits, so this remains resume-safe.
        throttle_attempts = 0
        while True:
            try:
                external_free_gb = shutil.disk_usage(out_dir).free / (1024**3)
                system_free_gb = None
                if self.system_root.is_dir():
                    system_free_gb = shutil.disk_usage(self.system_root).free / (1024**3)
                system_low = system_free_gb is not None and system_free_gb < self.min_system_free_gb
                if external_free_gb < self.min_external_free_gb or system_low:
                    if throttle_attempts % 6 == 0:  # Log every ~5 mins
                        system_text = f"; system {system_free_gb:.1f} GiB free" if system_free_gb is not None else ""
                        logger.warning(
                            "[Disk Throttle] External volume %.1f GiB free "
                            "(threshold %.1f GiB)%s; "
                            "worker pausing to allow quantification/reclamation to proceed...",
                            external_free_gb,
                            self.min_external_free_gb,
                            system_text,
                        )
                    time.sleep(60)  # Wait 1 minute and re-check
                    throttle_attempts += 1
                else:
                    if throttle_attempts > 0:
                        logger.info(
                            "[Disk Throttle] Space available again (external %.1f GiB%s). " "Resuming download for %s.",
                            external_free_gb,
                            (f"; system {system_free_gb:.1f} GiB" if system_free_gb is not None else ""),
                            srr_id,
                        )
                    break
            except Exception as e:
                logger.error(f"Failed to check disk space: {e}")
                break  # Proceed anyway if check fails

        downloader = ENADownloader(
            timeout=self.download_timeout_seconds,
            retries=self.download_retries,
            integrity_retries=self.download_integrity_retries,
            invalid_witness_max_bytes=self.invalid_witness_max_bytes,
            retry_delay_seconds=self.download_retry_delay_seconds,
            speed_limit_bytes=self.download_speed_limit_bytes,
            speed_time_seconds=self.download_speed_time_seconds,
            api_retries=self.ena_api_retries,
            api_retry_delay_seconds=self.ena_api_retry_delay_seconds,
        )
        success, message, downloaded_files = downloader.download_run(srr_id, sample_dir)

        if success:
            for fq in downloaded_files:
                try:
                    file_size = fq.stat().st_size
                except OSError:
                    file_size = 0
                sz_gb = file_size / (1024**3)
                logger.info(f"Downloaded {fq.name}: {sz_gb:.2f} GB")
                if sz_gb == 0:
                    logger.error(f"Downloaded file {fq.name} is empty (0.00 GB). Marking as failed.")
                    success = False
            if success:
                # ENA's transfer validator proves the downloaded payload is a
                # readable gzip stream, but the current-contract witness also
                # needs to bind the accession to the expected single/paired
                # FASTQ layout.  Direct ENA success used to return before this
                # marker was written, leaving an otherwise valid cloud result
                # non-promotable after raw-input reclamation.
                if self._raw_validation_semaphore is None:
                    validated_fastqs = _validated_local_fastq_inputs(
                        sample_dir,
                        srr_id,
                        expected_paired,
                        assume_payload_valid=True,
                    )
                else:
                    with self._raw_validation_semaphore:
                        validated_fastqs = _validated_local_fastq_inputs(
                            sample_dir,
                            srr_id,
                            expected_paired,
                            assume_payload_valid=True,
                        )
                if not validated_fastqs:
                    logger.error(
                        "Downloaded FASTQ layout for %s does not match the declared library layout",
                        srr_id,
                    )
                    success = False
                else:
                    downloaded_files = validated_fastqs
        if not success:
            logger.warning(f"ENA download failed for {srr_id}: {message}")
            logger.info(f"Falling back to NCBI fasterq-dump for {srr_id}...")
            # Retain completed accession FASTQs as per-mate checkpoints. A
            # paired layout remains non-quantifiable until both mates pass the
            # full layout validator, but discarding mate 1 when mate 2 fails
            # needlessly repeats a potentially multi-gigabyte transfer.
            # Retain ENA ``.part`` and diagnostic ``.invalid`` files as well.
            for f in sample_dir.iterdir():
                if (
                    f.is_file()
                    and not _is_diagnostic_transfer_artifact(f)
                    and not f.name.endswith(".sra")
                    and not _is_accession_fastq_file(f, srr_id)
                ):
                    f.unlink(missing_ok=True)

            if defer_ncbi_fallback and not local_sra_attempted:
                logger.info(
                    "Deferring NCBI fallback for %s to the dedicated extraction queue.",
                    srr_id,
                )
                return None
            if local_sra_attempted:
                logger.warning(
                    "Skipping duplicate NCBI fallback extraction for %s: the campaign-local "
                    "SRA archive already failed in this invocation and remains available "
                    "for an idempotent retry.",
                    srr_id,
                )
            else:
                try:
                    success = self._extract_sra_to_fastq_bounded(
                        srr_id,
                        sample_dir,
                        out_dir,
                        srr_id,
                        expected_paired,
                    )
                    if success:
                        logger.info(f"NCBI fallback succeeded for {srr_id}")
                except Exception as e:
                    logger.error(f"NCBI fasterq-dump exception for {srr_id}: {e}")

        if not success:
            # Final cleanup on total failure. Retain canonical completed mates,
            # resumable ENA ``.part`` files, and diagnostic ``.invalid``
            # payloads. A later ENA attempt can skip each validated complete
            # mate, and a fallback promotion can keep it while supplying the
            # missing mate.
            for f in sample_dir.iterdir():
                if (
                    f.is_file()
                    and not _is_diagnostic_transfer_artifact(f)
                    and not f.name.endswith(".sra")
                    and not _is_accession_fastq_file(f, srr_id)
                ):
                    f.unlink(missing_ok=True)

        return success

    def _prefetch_sra_accession(self, srr_id: str) -> Path | None:
        """Resume and verify one NCBI accession before FASTQ extraction."""

        cached = _campaign_cached_sra_path(srr_id)
        if cached is not None:
            return cached
        prefetch = shutil.which("prefetch")
        if prefetch is None:
            logger.error("NCBI prefetch is unavailable for %s", srr_id)
            return None

        root = _data_root() / ".sra-cache" / "prefetch" / srr_id
        root.mkdir(parents=True, exist_ok=True)
        command = [
            prefetch,
            "--transport",
            "http",
            "--resume",
            "yes",
            "--verify",
            "yes",
            "--output-directory",
            str(root),
            srr_id,
        ]
        logger.info("Prefetching NCBI accession %s with resumable verified transport", srr_id)
        try:
            result = _run_command_in_process_group(
                command,
                timeout=self.ncbi_prefetch_timeout_seconds,
                env=build_sra_environment(_data_root()),
            )
        except subprocess.TimeoutExpired:
            logger.error("NCBI prefetch timed out for %s", srr_id)
            return None
        if result.returncode != 0:
            detail = (result.stderr or result.stdout or "").strip()
            logger.error("NCBI prefetch failed for %s (exit %d): %s", srr_id, result.returncode, detail[:300])
            return None

        candidates = [root / f"{srr_id}.sra", root / f"{srr_id}.sra.cache", root / srr_id]
        candidates.extend(sorted(root.glob("*.sra")))
        candidates.extend(sorted(root.glob("*.sra.cache")))
        for candidate in candidates:
            if candidate.is_file() and candidate.stat().st_size > 0:
                logger.info("Prefetch completed for %s: %s", srr_id, candidate)
                return candidate
        cached = _campaign_cached_sra_path(srr_id)
        if cached is not None:
            return cached
        logger.error("NCBI prefetch reported success but produced no SRA archive for %s", srr_id)
        return None

    def _download_fastq_ncbi_only(
        self,
        srr_id: str,
        out_dir: Path,
        expected_paired: bool | None,
    ) -> bool:
        """Run only the bounded NCBI extraction lane for a deferred sample."""

        sample_dir = out_dir / srr_id
        with self._fs_lock:
            sample_dir.mkdir(parents=True, exist_ok=True)
        if self._raw_validation_semaphore is None:
            local_fastqs = _validated_local_fastq_inputs(sample_dir, srr_id, expected_paired)
        else:
            with self._raw_validation_semaphore:
                local_fastqs = _validated_local_fastq_inputs(sample_dir, srr_id, expected_paired)
        if local_fastqs:
            logger.info(
                "Deferred fallback found %d validated FASTQ file(s) for %s; skipping extraction.",
                len(local_fastqs),
                srr_id,
            )
            return True

        local_sra = _local_sra_path(sample_dir, out_dir, srr_id)
        if local_sra is None and self.ncbi_prefetch_first:
            local_sra = self._prefetch_sra_accession(srr_id)
        source: Path | str = local_sra if local_sra is not None else srr_id
        if self.ncbi_prefetch_first and local_sra is None:
            return False
        logger.info("Running dedicated NCBI fallback for %s from %s", srr_id, source)
        success = self._extract_sra_to_fastq_bounded(
            source,
            sample_dir,
            out_dir,
            srr_id,
            expected_paired,
        )
        if success:
            logger.info("Dedicated NCBI fallback succeeded for %s", srr_id)
        return success

    def quant_sample(
        self,
        config_path: Path,
        batch_index: int,
        species_name: str,
        threads: int,
        srr_id: str,
        expected_paired: bool | None = None,
    ) -> tuple:
        """Run amalgkit quant for a single sample using --batch.

        Returns:
            Tuple of (success: bool, error_msg: str or None).
        """
        with open(config_path) as f:
            cfg = yaml.safe_load(f)

        work_dir = _species_work_dir(species_name)
        meta_path = _resolve_quant_metadata_path(work_dir)

        if not meta_path:
            logger.error(f"No metadata found for {species_name} batch {batch_index}")
            return False, f"No metadata found for {species_name}"

        # Direct ENA acquisition bypasses Amalgkit's getfastq stage.  For
        # public metadata rows with missing/zero spot statistics, derive the
        # quant-only fallback from the validated FASTQ before invoking
        # Amalgkit 0.16.60; otherwise it rejects usable reads before Kallisto.
        metadata_file = Path(meta_path)
        sample_fastq_dir = work_dir / "getfastq" / srr_id
        if _metadata_requires_fastq_input_stats(metadata_file, srr_id):
            _write_fastq_input_stats(sample_fastq_dir, srr_id)

        cmd = _build_quant_command(cfg, species_name, batch_index, threads, meta_path)
        scratch_context = self._prepare_local_quant_workspace(
            work_dir,
            species_name,
            srr_id,
            meta_path,
            cmd,
            expected_paired=expected_paired,
        )
        run_cmd = scratch_context[0] if scratch_context else cmd
        scratch_root = scratch_context[1] if scratch_context else None

        log_path = self.log_dir / f"{species_name}_quant.log"

        try:
            with open(log_path, "a") as log_f:
                log_f.write(
                    f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] Quant batch {batch_index} START: {' '.join(run_cmd)}\n"
                )

            result = _run_command_in_process_group(run_cmd, self.quant_timeout_seconds)

            with open(log_path, "a") as log_f:
                log_f.write(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Quant batch {batch_index} END (Exit {result.returncode})\n"
                )
                if result.stdout:
                    log_f.write(f"-- STDOUT --\n{result.stdout}\n")
                if result.stderr:
                    log_f.write(f"-- STDERR --\n{result.stderr}\n")

            if result.returncode == 0:
                if scratch_root is not None and not self._promote_local_quant_output(work_dir, srr_id, scratch_root):
                    return False, f"Unable to promote local quantification output for {srr_id}"
                sample_dir = work_dir / "quant" / srr_id
                quantification_file = find_quantification_file(sample_dir, srr_id)
                if quantification_file is None:
                    return False, f"No abundance output produced for {srr_id}"
                try:
                    reference_manifest = work_dir / "reference" / "reference_aliases.json"
                    write_quant_provenance(
                        sample_dir,
                        species=species_name,
                        run_accession=srr_id,
                        config_path=config_path,
                        command=run_cmd,
                        reference_manifest_path=reference_manifest if reference_manifest.exists() else None,
                        quantification_file=quantification_file,
                    )
                except OSError as exc:
                    return False, f"Unable to write current-method provenance: {exc}"
                return True, None

            # Extract meaningful error from amalgkit output
            error_msg = "Quantification Failed"
            for line in (result.stdout + result.stderr).splitlines():
                line_lower = line.strip().lower()
                if any(
                    kw in line_lower for kw in ["error", "exiting", "no sample", "not found", "failed", "exception"]
                ):
                    error_msg = line.strip()[:120]
                    break
            return False, error_msg
        except subprocess.TimeoutExpired:
            logger.error(
                "Quant timeout batch %s after %s",
                batch_index,
                _format_duration(self.quant_timeout_seconds),
            )
            with open(log_path, "a") as log_f:
                log_f.write(f"TIMEOUT after {self.quant_timeout_seconds}s\n")
            return False, f"Quant timeout (>{_format_duration(self.quant_timeout_seconds)})"
        except Exception as e:
            logger.error(f"Quant exception batch {batch_index}: {e}")
            with open(log_path, "a") as log_f:
                log_f.write(f"EXCEPTION: {e}\n")
            return False, f"Exception: {e}"
        finally:
            if scratch_root is not None:
                self._cleanup_local_quant_workspace(scratch_root)

    def _reserve_local_scratch(
        self,
        scratch_dir: Path,
        payload_bytes: int,
        *,
        purpose: str,
    ) -> bool:
        """Reserve host scratch capacity across all local-I/O optimizations."""

        if payload_bytes <= 0:
            return False
        reserve_bytes = int(self.local_quant_scratch_reserve_gb * (1024**3))
        with self._local_scratch_reservation_lock:
            try:
                scratch_dir.mkdir(parents=True, exist_ok=True)
                free_bytes = shutil.disk_usage(scratch_dir).free
            except OSError as exc:
                logger.warning(
                    "Unable to inspect %s capacity at %s: %s",
                    purpose,
                    scratch_dir,
                    exc,
                )
                return False
            available_after_reservation = free_bytes - self._local_scratch_reserved_bytes - payload_bytes
            if available_after_reservation <= reserve_bytes:
                logger.info(
                    "%s would leave %.2f GiB free, below the %.2f GiB reserve; "
                    "using mounted-volume storage for this sample.",
                    purpose,
                    available_after_reservation / (1024**3),
                    self.local_quant_scratch_reserve_gb,
                )
                return False
            self._local_scratch_reserved_bytes += payload_bytes
            logger.info(
                "[%s] reserved %.2f GiB for this sample; "
                "in-flight reservation %.2f GiB, free before copy %.2f GiB, reserve %.2f GiB",
                purpose,
                payload_bytes / (1024**3),
                self._local_scratch_reserved_bytes / (1024**3),
                free_bytes / (1024**3),
                self.local_quant_scratch_reserve_gb,
            )
            return True

    def _release_local_scratch(self, payload_bytes: int, *, purpose: str) -> None:
        """Release capacity from the shared host-scratch reservation ledger."""

        with self._local_scratch_reservation_lock:
            self._local_scratch_reserved_bytes = max(
                0,
                self._local_scratch_reserved_bytes - max(payload_bytes, 0),
            )
            logger.debug(
                "[%s] released %.2f GiB; in-flight reservation %.2f GiB",
                purpose,
                max(payload_bytes, 0) / (1024**3),
                self._local_scratch_reserved_bytes / (1024**3),
            )

    def _reserve_local_fastq_stage(self, payload_bytes: int) -> bool:
        """Reserve local scratch capacity for one FASTQ staging copy."""

        if self.local_quant_scratch_dir is None:
            return False
        return self._reserve_local_scratch(
            self.local_quant_scratch_dir,
            payload_bytes,
            purpose="Local FASTQ stage",
        )

    def _release_local_fastq_stage(self, payload_bytes: int) -> None:
        """Release a temporary local FASTQ staging reservation."""

        self._release_local_scratch(
            payload_bytes,
            purpose="Local FASTQ stage",
        )

    def _retain_local_fastq_stage_reservation(self, scratch_root: Path, payload_bytes: int) -> None:
        """Keep a successful FASTQ-stage reservation until workspace cleanup."""

        with self._local_scratch_reservation_lock:
            self._local_fastq_stage_reservations[scratch_root.absolute()] = payload_bytes

    def _release_local_fastq_stage_for_workspace(self, scratch_root: Path) -> None:
        """Release the FASTQ bytes retained by one completed/aborted workspace."""

        with self._local_scratch_reservation_lock:
            payload_bytes = self._local_fastq_stage_reservations.pop(scratch_root.absolute(), 0)
        if payload_bytes:
            self._release_local_fastq_stage(payload_bytes)

    def _stage_local_fastq_inputs(
        self,
        work_dir: Path,
        scratch_root: Path,
        srr_id: str,
        expected_paired: bool | None = None,
    ) -> bool:
        """Copy validated FASTQs to local scratch for one quantification.

        The external raw inputs remain authoritative and are never moved. A
        staged copy is promoted only after all files have been copied with
        matching sizes; any failure falls back to the normal mounted-volume
        symlink and leaves the external inputs intact.
        """

        if not self.local_fastq_stage or self.local_quant_scratch_dir is None:
            return False
        source_dir = work_dir / "getfastq" / srr_id
        candidates = _validated_local_fastq_inputs(
            source_dir,
            srr_id,
            expected_paired,
        )
        if not candidates:
            logger.debug("No validated FASTQs available to stage for %s", srr_id)
            return False
        try:
            payload_bytes = sum(path.stat().st_size for path in candidates)
        except OSError as exc:
            logger.warning("Unable to size FASTQs for local staging %s: %s", srr_id, exc)
            return False
        if not self._reserve_local_fastq_stage(payload_bytes):
            return False

        staged_root = scratch_root / "getfastq"
        staged_sample_dir = staged_root / srr_id
        staged_successfully = False
        try:
            staged_sample_dir.mkdir(parents=True, exist_ok=True)
            reused = 0
            for source in candidates:
                destination = staged_sample_dir / source.name
                source_stat = source.stat()
                try:
                    destination_stat = destination.stat()
                except OSError:
                    destination_stat = None
                if (
                    destination_stat is not None
                    and destination_stat.st_size == source_stat.st_size
                    and destination_stat.st_mtime_ns == source_stat.st_mtime_ns
                ):
                    reused += 1
                    continue
                temporary = destination.with_name(f".{destination.name}.{os.getpid()}.{threading.get_ident()}.copying")
                shutil.copy2(source, temporary)
                if temporary.stat().st_size != source_stat.st_size:
                    raise OSError(f"staged FASTQ size mismatch for {source.name}")
                os.replace(temporary, destination)
            logger.info(
                "Staged %d validated FASTQ file(s) locally for %s (%.2f GB; reused %d)",
                len(candidates),
                srr_id,
                payload_bytes / (1024**3),
                reused,
            )
            self._retain_local_fastq_stage_reservation(scratch_root, payload_bytes)
            staged_successfully = True
            return True
        except OSError as exc:
            logger.warning(
                "Local FASTQ staging failed for %s; using mounted inputs: %s",
                srr_id,
                exc,
            )
            shutil.rmtree(staged_root, ignore_errors=True)
            return False
        finally:
            if not staged_successfully:
                self._release_local_fastq_stage(payload_bytes)

    def _prepare_local_quant_workspace(
        self,
        work_dir: Path,
        species_name: str,
        srr_id: str,
        metadata_path: str,
        command: List[str],
        *,
        expected_paired: bool | None = None,
    ) -> tuple[List[str], Path] | None:
        """Prepare a local output workspace and optionally stage its inputs.

        Amalgkit resolves FASTQ inputs relative to ``--out_dir``.  A temporary
        local workspace therefore uses a directory symlink for ``metadata``
        and either a validated local FASTQ copy or a ``getfastq`` symlink.
        The completed sample is copied atomically back into the external work
        tree; no authoritative raw input is moved or deleted by this
        optimization.
        """

        if self.local_quant_scratch_dir is None:
            return None
        scratch_root = self.local_quant_scratch_dir / species_name / srr_id
        try:
            self.local_quant_scratch_dir.mkdir(parents=True, exist_ok=True)
            available_gb = shutil.disk_usage(self.local_quant_scratch_dir).free / (1024**3)
            if available_gb <= self.local_quant_scratch_reserve_gb:
                logger.warning(
                    "Local quant scratch has %.2f GiB free; reserve is %.2f GiB. "
                    "Using external quant output for %s/%s.",
                    available_gb,
                    self.local_quant_scratch_reserve_gb,
                    species_name,
                    srr_id,
                )
                return None
            if scratch_root.is_symlink():
                scratch_root.unlink()
            elif scratch_root.exists() and not scratch_root.is_dir():
                scratch_root.unlink()
            scratch_root.mkdir(parents=True, exist_ok=True)
            # Keep a completed staged FASTQ directory across an interrupted
            # quantification, but reset only transient metadata/output trees.
            # This makes a stopped run resume without copying the same large
            # external FASTQ back to the local SSD.
            for transient_name in ("metadata", "quant"):
                transient = scratch_root / transient_name
                if transient.is_symlink() or (transient.exists() and not transient.is_dir()):
                    transient.unlink()
                elif transient.is_dir():
                    shutil.rmtree(transient)
            staged_getfastq = scratch_root / "getfastq"
            if staged_getfastq.is_symlink() or (staged_getfastq.exists() and not staged_getfastq.is_dir()):
                staged_getfastq.unlink()
            (scratch_root / "metadata").symlink_to(work_dir / "metadata", target_is_directory=True)
            if not self._stage_local_fastq_inputs(
                work_dir,
                scratch_root,
                srr_id,
                expected_paired,
            ):
                # Amalgkit v0.16.60 rejects a symlinked input/output root. If
                # the local FASTQ copy cannot fit under the host reserve, use
                # the canonical external workspace instead of constructing a
                # local workspace that is guaranteed to fail validation.
                self._cleanup_local_quant_workspace(scratch_root)
                logger.info(
                    "Local quant scratch unavailable for %s/%s; using canonical mounted workspace",
                    species_name,
                    srr_id,
                )
                return None
            (scratch_root / "quant").mkdir()
            local_command = list(command)
            local_command[local_command.index("--out_dir") + 1] = str(scratch_root)
            local_command[local_command.index("--metadata") + 1] = str(
                scratch_root / "metadata" / Path(metadata_path).name
            )
            logger.debug("Using local quant scratch for %s/%s: %s", species_name, srr_id, scratch_root)
            return local_command, scratch_root
        except OSError as exc:
            logger.warning(
                "Unable to prepare local quant scratch for %s/%s; using external output: %s",
                species_name,
                srr_id,
                exc,
            )
            if scratch_root.is_dir() and not scratch_root.is_symlink():
                shutil.rmtree(scratch_root, ignore_errors=True)
            self._release_local_fastq_stage_for_workspace(scratch_root)
            return None

    def _promote_local_quant_output(self, work_dir: Path, srr_id: str, scratch_root: Path) -> bool:
        """Copy one validated local quant result into the canonical data-root path."""

        local_sample_dir = scratch_root / "quant" / srr_id
        if find_quantification_file(local_sample_dir, srr_id) is None:
            return False
        external_quant_root = work_dir / "quant"
        external_quant_root.mkdir(parents=True, exist_ok=True)
        external_sample_dir = external_quant_root / srr_id
        temporary = external_quant_root / f".{srr_id}.promote.{os.getpid()}.{threading.get_ident()}"
        try:
            if temporary.exists() or temporary.is_symlink():
                if temporary.is_dir() and not temporary.is_symlink():
                    shutil.rmtree(temporary)
                else:
                    temporary.unlink()
            shutil.copytree(local_sample_dir, temporary, symlinks=False)
            if find_quantification_file(temporary, srr_id) is None:
                shutil.rmtree(temporary, ignore_errors=True)
                return False
            if external_sample_dir.exists() or external_sample_dir.is_symlink():
                archive_dir = work_dir / "archive" / "quant_replaced" / srr_id / str(time.time_ns())
                archive_dir.parent.mkdir(parents=True, exist_ok=True)
                external_sample_dir.replace(archive_dir)
            temporary.replace(external_sample_dir)
            return True
        except OSError as exc:
            logger.error("Unable to promote local quantification for %s: %s", srr_id, exc)
            if temporary.exists() and temporary.is_dir() and not temporary.is_symlink():
                shutil.rmtree(temporary, ignore_errors=True)
            return False

    def _cleanup_local_quant_workspace(self, scratch_root: Path) -> None:
        """Remove only the exact per-sample local scratch directory."""

        cleanup_succeeded = False
        try:
            if scratch_root.is_symlink():
                scratch_root.unlink()
            elif scratch_root.exists():
                shutil.rmtree(scratch_root)
            cleanup_succeeded = True
        except OSError as exc:
            logger.warning("Unable to clean local quant scratch %s: %s", scratch_root, exc)
        if cleanup_succeeded:
            self._release_local_fastq_stage_for_workspace(scratch_root)

    def is_quantified(self, species_name: str, srr_id: str) -> bool:
        """Check for a reusable quantification, including compatible version drift."""
        quant_dir = _species_work_dir(species_name) / "quant" / srr_id
        quant_file = find_quantification_file(quant_dir, srr_id)
        if not quant_dir.exists() or quant_file is None or not quant_file.is_file() or quant_file.stat().st_size <= 100:
            return False
        return is_reusable_quantification(quant_dir, srr_id)

    def _reclaim_after_current_quantification(self, species_name: str, srr_id: str) -> None:
        """Reclaim one sample's raw inputs only after current evidence exists."""

        if not self.reclaim_raw_after_quant or not self.is_quantified(species_name, srr_id):
            return
        cleanup = reclaim_sample_raw_inputs(_species_work_dir(species_name), srr_id)
        if cleanup["errors"]:
            logger.warning(
                "Raw-input reclamation for %s/%s completed with %d errors: %s",
                species_name,
                srr_id,
                len(cleanup["errors"]),
                "; ".join(cleanup["errors"][:3]),
            )
        elif cleanup["files_deleted"]:
            logger.info(
                "Reclaimed %.2f GB of raw inputs for current quantification %s/%s",
                cleanup["bytes_freed"] / (1024**3),
                species_name,
                srr_id,
            )

    def _reset_non_current_quantified_states(self, species_name: str) -> int:
        """Audit quantified outputs without invalidating runtime-version drift."""

        sample_ids = self.db.get_samples(species_name, "quantified")
        if not sample_ids:
            return 0

        profile = getattr(self, "_resource_profile", None)
        validation_slots = max(
            1,
            min(
                len(sample_ids),
                getattr(profile, "validation_slots", 1),
            ),
        )

        # Re-digesting every quantified payload on every resume made the
        # Phase-1 audit O(campaign bytes): a large species (apis_mellifera,
        # ~6.7 GB across 3190 sample dirs) stalled discovery for 15+ minutes
        # on cold external-volume reads, so no task was ever submitted.
        # The contract audit (sidecar, accession, config, reference manifest)
        # remains fail-closed; full content re-verification is opt-in.
        verify_resume_hashes = os.environ.get(
            "AMALGKIT_RESUME_VERIFY_HASHES", ""
        ).strip().lower() in {"1", "true", "yes", "on"}

        def classify(srr_id: str) -> tuple[str, dict[str, Any]]:
            quant_dir = _species_work_dir(species_name) / "quant" / srr_id
            return (
                srr_id,
                classify_quantification(quant_dir, srr_id, verify_content=verify_resume_hashes),
            )

        if validation_slots == 1:
            classifications = [classify(srr_id) for srr_id in sample_ids]
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=validation_slots) as validation_executor:
                classifications = list(validation_executor.map(classify, sample_ids))

        reset = 0
        requantification_candidates = 0
        quarantined = 0
        for srr_id, audit in classifications:
            status = audit["status"]
            if status in {QUANT_STATUS_CURRENT, QUANT_STATUS_VERSION_DRIFT}:
                should_requantify = self.requantification_policy == "all" or (
                    self.requantification_policy == "version-drift" and status == QUANT_STATUS_VERSION_DRIFT
                )
                self.db.record_quantification_audit(
                    species_name,
                    srr_id,
                    status=status,
                    reason=audit["reason"],
                    contract_id=audit.get("contract_id"),
                    observed_amalgkit_version=audit.get("observed_amalgkit_version"),
                    observed_release_tag=audit.get("observed_release_tag"),
                    observed_source_revision=audit.get("observed_source_revision"),
                    state="pending" if should_requantify else "quantified",
                )
                if should_requantify:
                    reset += 1
                    requantification_candidates += 1
            else:
                self.db.record_quantification_audit(
                    species_name,
                    srr_id,
                    status=status,
                    reason=audit["reason"],
                    contract_id=audit.get("contract_id"),
                    observed_amalgkit_version=audit.get("observed_amalgkit_version"),
                    observed_release_tag=audit.get("observed_release_tag"),
                    observed_source_revision=audit.get("observed_source_revision"),
                    state="quarantined",
                )
                reset += 1
                quarantined += 1
        if len(sample_ids) > 1:
            logger.info(
                "Validated %d recorded quantification states for %s with %d provenance-read slot(s); " "policy=%s",
                len(sample_ids),
                species_name,
                validation_slots,
                self.requantification_policy,
            )
        if requantification_candidates:
            logger.info(
                "Queued %d explicit re-quantification candidates for %s",
                requantification_candidates,
                species_name,
            )
        if quarantined:
            logger.info(
                "Quarantined %d incomplete or unverifiable quantification outputs for %s",
                quarantined,
                species_name,
            )
        return reset

    def process_single_sample(
        self,
        srr_id: str,
        batch_idx: int,
        fastq_dir: Path,
        config_path: Path,
        species_name: str,
        threads: int,
        expected_paired: bool | None = None,
        *,
        defer_ncbi_fallback: bool = False,
    ) -> Dict[str, Any]:
        """Processing unit for a single sample."""
        result = {
            "srr": srr_id,
            "batch": batch_idx,
            "downloaded": False,
            "quantified": False,
            "skipped": False,
            "fallback_deferred": False,
            "error": None,
        }

        if self.is_quantified(species_name, srr_id):
            self._reclaim_after_current_quantification(species_name, srr_id)
            result.update({"downloaded": True, "quantified": True, "skipped": True})
            return result

        # Download
        self.db.set_state(species_name, srr_id, "downloading")
        download_outcome = self.download_fastq(
            srr_id,
            fastq_dir,
            expected_paired,
            defer_ncbi_fallback=defer_ncbi_fallback,
        )
        if download_outcome is None:
            self.db.set_state(
                species_name,
                srr_id,
                "pending",
                error="NCBI fallback queued after ENA acquisition failure",
            )
            result["fallback_deferred"] = True
            return result
        if not download_outcome:
            error_msg = "Download Failed (all sources: ENA FTP/HTTP, NCBI)"
            self.db.set_state(species_name, srr_id, "failed", error=error_msg)
            result["error"] = error_msg
            return result
        self.db.set_state(species_name, srr_id, "downloaded")
        result["downloaded"] = True

        # Quantify
        self.db.set_state(species_name, srr_id, "quantifying")
        if self._quant_semaphore is None:
            success, quant_error = self.quant_sample(
                config_path,
                batch_idx,
                species_name,
                threads,
                srr_id,
                expected_paired,
            )
        else:
            with self._quant_semaphore:
                success, quant_error = self.quant_sample(
                    config_path,
                    batch_idx,
                    species_name,
                    threads,
                    srr_id,
                    expected_paired,
                )
        if success:
            self.db.set_state(species_name, srr_id, "quantified")
            self._reclaim_after_current_quantification(species_name, srr_id)
        else:
            error_msg = quant_error or "Quantification Failed"
            self.db.set_state(species_name, srr_id, "failed", error=error_msg)
        result["quantified"] = success
        if not success:
            result["error"] = quant_error or "Quantification Failed"

        return result

    def verify_genome_index(self, config_path: Path, species_name: str) -> bool:
        """Verify Kallisto index exists."""
        try:
            with open(config_path) as f:
                cfg = yaml.safe_load(f)
        except Exception as e:
            logger.error(f"Failed to load config {config_path}: {e}")
            return False

        logger.info(f"Verifying index for {species_name}...")
        for d in _quant_index_candidates(cfg, species_name):
            if not d:
                continue

            path_obj = Path(d)
            abs_path = path_obj.resolve() if path_obj.exists() else path_obj.absolute()

            logger.info(f"  Checking {d} -> {abs_path}")

            if not path_obj.exists():
                logger.warning(f"  Path does not exist: {d}")
                continue

            # Check permissions by attempting listdir
            try:
                files = [
                    candidate
                    for candidate in path_obj.glob("*.idx")
                    if candidate.is_file() and candidate.stat().st_size > 0
                ]
                if files:
                    logger.info(f"  Genome index found in {d}: {[f.name for f in files]}")
                    return True
                else:
                    logger.warning(f"  Path exists but no .idx files found in {d}")
                    try:
                        contents = os.listdir(d)
                        logger.info(f"  Directory contents: {contents[:5]}...")
                    except Exception as e:
                        logger.error(f"  Failed to list directory {d}: {e}")
            except Exception as e:
                logger.error(f"  Permission/Error accessing {d}: {e}")

        logger.error(f"No genome index found for {species_name}")
        return False

    def run_tissue_normalization(self, metadata_path: Path) -> bool:
        """Run tissue normalization atomically and preserve source labels.

        The orchestrator may revisit the same selected table on every resume.
        ``tissue_original`` therefore remains the first observed value instead
        of being overwritten by an already-normalized value. A failed write or
        mapping operation is reported to the caller so discovery can fail
        closed rather than proceeding with an ambiguous metadata state.
        """
        mapping_path = self.config_dir / "tissue_mapping.yaml"
        patches_path = self.config_dir / "tissue_patches.yaml"

        if not mapping_path.exists():
            logger.warning(f"Tissue mapping not found at {mapping_path}, skipping normalization.")
            return True

        logger.info(f"Normalizing tissues in {metadata_path}")
        temporary: Path | None = None
        try:
            import pandas as pd

            df = pd.read_csv(metadata_path, sep="\t", low_memory=False)
            if "tissue" not in df.columns:
                logger.warning("Metadata has no tissue column; skipping tissue normalization.")
                return True
            original_tissue = df["tissue_original"].copy() if "tissue_original" in df.columns else df["tissue"].copy()

            # Apply normalization
            df_norm = apply_tissue_normalization(
                df,
                mapping_path=mapping_path,
                patches_path=patches_path if patches_path.exists() else None,
                output_column="tissue_normalized",
            )

            # Keep the normalized tissue in the canonical field while retaining
            # the source value for provenance.
            if "tissue_normalized" in df_norm.columns:
                df_norm["tissue_original"] = original_tissue
                df_norm["tissue"] = df_norm["tissue_normalized"]
                # Drop temporary column if desired, or keep it

            temporary = metadata_path.with_name(f".{metadata_path.name}.{os.getpid()}.{threading.get_ident()}.tmp")
            df_norm.to_csv(temporary, sep="\t", index=False)
            # Avoid changing the source table when normalization is already
            # idempotent; otherwise promote the complete file atomically.
            if not metadata_path.exists() or temporary.read_bytes() != metadata_path.read_bytes():
                os.replace(temporary, metadata_path)
            else:
                temporary.unlink(missing_ok=True)
            logger.info("Tissue normalization complete and saved.")
            return True

        except Exception as e:
            logger.error(f"Tissue normalization failed: {e}")
            if temporary is not None:
                temporary.unlink(missing_ok=True)
            return False

    def discover_species_tasks(self, config_name: str, max_gb: float, threads: int) -> List[Dict[str, Any]]:
        """Main processing loop for a species.
        Returns a list of tasks instead of executing them immediately.
        """
        source_config_path = self.config_dir / config_name
        if not source_config_path.exists():
            logger.error(f"Config not found: {source_config_path}")
            return []

        species_name = _species_name_from_config(config_name)
        config_path = _runtime_config_path(source_config_path, species_name)
        logger.info(f"=== Processing {species_name} ===")

        # Autonomous Preprocessing Detection
        work_dir = _species_work_dir(species_name)
        metadata_path = _resolve_metadata_path(work_dir)

        needs_metadata = not metadata_path.exists() or not is_current_metadata(work_dir)
        needs_index = not self.verify_genome_index(config_path, species_name)
        needs_prep = needs_metadata or needs_index
        if needs_metadata:
            logger.info(f"Metadata missing for {species_name}. Queuing autonomous preprocessing...")
        if needs_index:
            logger.info(f"Genome index missing for {species_name}. Queuing autonomous preprocessing...")

        if needs_prep:
            logger.info(
                "Running current Amalgkit preparation stages (reference -> metadata -> select) for %s...",
                species_name,
            )

            # Amalgkit 0.16.60 builds the quantification index through the
            # reference-preparation helper; it no longer exposes the removed
            # `config` or `index` CLI stages.
            if needs_index:
                try:
                    from metainformant.rna.engine.workflow import load_workflow_config
                    from metainformant.rna.engine.workflow_planning import prepare_reference_genome

                    runtime_config = load_workflow_config(config_path)
                    if not prepare_reference_genome(runtime_config):
                        logger.error("Reference preparation failed for %s.", species_name)
                        return []
                except Exception as exc:
                    logger.error("Reference preparation exception for %s: %s", species_name, exc)
                    return []

            # Taxonomy acquisition is owned by Amalgkit 0.16.60.  When the
            # campaign launcher provides AMALGKIT_SHARED_DOWNLOAD_DIR, its
            # native taxdump/database cache is reused across species.  The
            # older per-species seed path is intentionally not recreated.
            # Without the shared setting, Amalgkit's own downloader remains
            # responsible for creating its inferred download directory.
            try:
                from metainformant.rna.amalgkit.amalgkit import build_amalgkit_command
                from metainformant.rna.engine.workflow import load_workflow_config
                from metainformant.rna.engine.workflow_planning import (
                    apply_step_defaults,
                    plan_workflow,
                    sanitize_params_for_cli,
                )

                runtime_config = apply_step_defaults(load_workflow_config(config_path))
                planned_steps = dict(plan_workflow(runtime_config))
                for prep_step in ("metadata", "select"):
                    raw_params = planned_steps.get(prep_step)
                    if raw_params is None:
                        logger.error("Current workflow plan omitted %s for %s", prep_step, species_name)
                        return []
                    prep_params = sanitize_params_for_cli(prep_step, dict(raw_params))
                    prep_params["redo"] = "yes"
                    prep_cmd = build_amalgkit_command(prep_step, prep_params)
                    prep_result = subprocess.run(
                        prep_cmd,
                        capture_output=True,
                        text=True,
                        timeout=7200,
                        check=False,
                    )
                    if prep_result.returncode != 0:
                        logger.error(
                            "%s preparation failed for %s:\n%s\n%s",
                            prep_step,
                            species_name,
                            prep_result.stderr,
                            prep_result.stdout,
                        )
                        return []
                logger.info("Metadata/select preparation complete for %s.", species_name)
            except subprocess.TimeoutExpired:
                logger.error("Metadata/select preparation timed out after 2h for %s.", species_name)
                return []
            except Exception as exc:
                logger.error("Metadata/select preparation exception for %s: %s", species_name, exc)
                return []

            # Post-prep explicit verification
            if needs_index and not self.verify_genome_index(config_path, species_name):
                logger.error(f"Genome index still missing for {species_name} after preparation.")
                return []

            metadata_path = _resolve_metadata_path(work_dir)

            if needs_metadata and not metadata_path.exists():
                logger.error(f"No metadata generated for {species_name} after preparation.")
                return []
        # Apply Normalization Step
        if not self.run_tissue_normalization(metadata_path):
            logger.error("Metadata normalization did not complete for %s", species_name)
            return []

        # Some valid Amalgkit 0.16.60 runs emit only ``metadata.tsv`` when
        # selection does not create a separate table.  Materialize the
        # canonical selected-table boundary once so the current metadata
        # witness can remain strict and every later resume can distinguish
        # the source table from the exact acquisition/quantification input.
        selected_metadata_path = work_dir / "metadata" / "metadata_selected.tsv"
        if not selected_metadata_path.is_file():
            try:
                shutil.copy2(metadata_path, selected_metadata_path)
                metadata_path = selected_metadata_path
                logger.info(
                    "Materialized canonical selected metadata for %s: %s",
                    species_name,
                    selected_metadata_path,
                )
            except OSError as exc:
                logger.error(
                    "Unable to materialize selected metadata for %s: %s",
                    species_name,
                    exc,
                )
                return []

        # Load samples
        import pandas as pd

        df = pd.read_csv(metadata_path, sep="\t", low_memory=False)

        filtered = _filter_metadata_by_size(df, max_gb)

        logger.info(f"Samples: {len(df)} total -> {len(filtered)} filtered (<= {max_gb} GB)")

        if filtered.empty:
            logger.info("No samples to process.")
            return []

        before_accession_filter = len(filtered)
        filtered = _filter_valid_run_rows(filtered, species_name, metadata_path)
        if len(filtered) != before_accession_filter:
            logger.info(
                "%s: %d rows remain after concrete run-accession validation",
                species_name,
                len(filtered),
            )
        if filtered.empty:
            logger.warning("No downloadable run accessions remain for %s.", species_name)
            return []

        try:
            runtime_cfg = yaml.safe_load(config_path.read_text()) or {}
        except (OSError, yaml.YAMLError) as exc:
            logger.error("Unable to load runtime config for reference preflight %s: %s", species_name, exc)
            return []
        target_column = next(
            (column for column in ("scientific_name", "organism", "species") if column in filtered.columns),
            None,
        )
        if target_column:
            target_names = sorted(
                str(value).strip() for value in filtered[target_column].dropna().unique() if str(value).strip()
            )
        else:
            target_names = [str(value).strip() for value in runtime_cfg.get("species_list", []) if str(value).strip()]
        references_ready, _, _ = _ensure_reference_alias_indexes(runtime_cfg, species_name, target_names)
        if not references_ready:
            return []

        # Kallisto cannot pseudoalign reads shorter than the active index
        # k-mer.  Exclude only rows with a known incompatible length before
        # acquisition, preserve the complete Amalgkit selection table for
        # auditability, and leave unknown lengths for Amalgkit's own
        # quant-time recovery from downloaded-read statistics.
        index_dir_value = _resolve_index_dir(runtime_cfg, species_name)
        index_dir = Path(index_dir_value) if index_dir_value else None
        reference_index = _reference_index_for_kmer(index_dir, target_names)
        kmer_size = _kallisto_index_kmer(str(reference_index)) if reference_index else KALLISTO_DEFAULT_KMER_SIZE
        quant_eligible, quant_excluded, _ = _filter_kallisto_ineligible_reads(filtered, kmer_size)
        if not quant_excluded.empty:
            selected_snapshot = metadata_path.with_name("metadata_selected_all.tsv")
            try:
                # Refresh the complete-selection snapshot whenever metadata
                # was regenerated.  On a resume without metadata work,
                # preserve the prior full snapshot while the canonical
                # selected table contains only quant-eligible rows.
                if needs_metadata or not selected_snapshot.exists():
                    shutil.copy2(metadata_path, selected_snapshot)
                exclusion_audit = _write_quant_eligibility_audit(
                    work_dir,
                    quant_excluded,
                    kmer_size=kmer_size,
                    index_path=reference_index,
                )
            except OSError as exc:
                logger.error("Unable to preserve quantification eligibility audit for %s: %s", species_name, exc)
                return []
            logger.warning(
                "%s: excluding %d samples before quantification because their known mean read length is below "
                "the Kallisto k-mer size (%d); see %s",
                species_name,
                len(quant_excluded),
                kmer_size,
                exclusion_audit,
            )
            filtered = quant_eligible

        # Register samples in progress DB and reconcile existing results
        srr_col = _sample_run_column(filtered)
        all_srr_ids = filtered[srr_col].tolist()
        self.db.init_species(species_name, all_srr_ids)
        self.db.prune_species(species_name, all_srr_ids)

        quant_dir = _species_work_dir(species_name) / "quant"
        reconciled = self.db.reconcile(species_name, quant_dir)
        if reconciled:
            logger.info(f"Reconciled {reconciled} already-quantified samples from filesystem")
        self._reset_non_current_quantified_states(species_name)

        # Mark all filtered samples as sampled so amalgkit quant processes them
        filtered["is_sampled"] = "yes"

        # Write sorted metadata for batch processing
        filtered.to_csv(metadata_path, sep="\t", index=False)

        # Write the metadata witness only after the final selected table has
        # been written.  The selected table is the acquisition/quantification
        # input and may be filtered for invalid accessions or incompatible
        # read lengths above; hashing it before this write made every restart
        # appear stale and needlessly repeated metadata discovery.
        try:
            rules_path = self.config_dir / "select_rules.tsv"
            write_metadata_provenance(
                work_dir,
                species=species_name,
                config_path=config_path,
                selection_rules_path=rules_path,
            )
        except OSError as exc:
            logger.error("Unable to write metadata provenance for %s: %s", species_name, exc)
            return []

        fastq_dir = work_dir / "getfastq"
        # ``reconcile`` has already established the current filesystem state
        # in the progress database, and the non-current reset above has
        # removed stale quantified states. Re-read that authoritative set
        # once instead of calling ``is_quantified`` for every accession. On
        # large external-volume cohorts, thousands of repeated directory/stat
        # calls here can delay execution before the first reusable sample is
        # even submitted.
        quantified_accessions = set(self.db.get_samples(species_name, "quantified"))
        permanently_excluded = self.db.get_excluded_srr_ids(species_name, reason_code="permanent_drop")
        if permanently_excluded:
            logger.info(
                "Skipping %d permanently excluded accessions for %s per sample_exclusions",
                len(permanently_excluded),
                species_name,
            )
        tasks = [
            task
            for task in _build_sample_tasks(filtered, srr_col, fastq_dir, config_path, species_name)
            if str(task["srr"]).strip() not in quantified_accessions
            and str(task["srr"]).strip() not in permanently_excluded
        ]

        return tasks

    def run_all(
        self,
        species_list: List[str],
        max_gb: float,
        workers: int,
        threads: int,
        *,
        quant_slots: Optional[int] = None,
        fasterq_threads: Optional[int] = None,
        fasterq_slots: Optional[int] = None,
        compression_threads: Optional[int] = None,
        max_in_flight: Optional[int] = None,
        validation_slots: Optional[int] = None,
        requantification_policy: Optional[str] = None,
        discovery_workers: Optional[int] = None,
    ) -> None:
        """Run the current pipeline with bounded, explicit stage resources.

        Species discovery/reference preparation is parallelized in an isolated
        phase. Sample execution starts only after every species discovery
        future has completed, preserving a single global task-prioritization
        boundary and preventing partial discovery from being mistaken for a
        complete cohort. The sample executor still uses a bounded submission
        window, and its active workers and total quantification budget remain
        separate from the fallback ``fasterq-dump`` and ``pigz`` budgets.
        """
        from collections import defaultdict, deque

        if requantification_policy is not None:
            requantification_policy = requantification_policy.strip().lower()
            if requantification_policy not in REQUANTIFICATION_POLICIES:
                raise ValueError(
                    f"Invalid requantification policy {requantification_policy!r}; "
                    f"expected one of {sorted(REQUANTIFICATION_POLICIES)}"
                )
            self.requantification_policy = requantification_policy

        profile = build_pipeline_resource_profile(
            workers,
            threads,
            quant_slots=quant_slots,
            fasterq_threads=fasterq_threads,
            fasterq_slots=fasterq_slots,
            compression_threads=compression_threads,
            validation_slots=validation_slots,
            max_in_flight=max_in_flight,
        )
        self._resource_profile = profile
        self._quant_semaphore = threading.BoundedSemaphore(profile.quant_slots)
        self._fasterq_semaphore = threading.BoundedSemaphore(profile.fasterq_slots)
        self._raw_validation_semaphore = threading.BoundedSemaphore(profile.validation_slots)
        configured_discovery_workers: Any = discovery_workers
        if configured_discovery_workers is None:
            configured_discovery_workers = os.environ.get("AMALGKIT_PIPELINE_DISCOVERY_WORKERS")
        discovery_worker_count = _resource_int(
            configured_discovery_workers,
            DEFAULT_DISCOVERY_WORKERS,
            "AMALGKIT_PIPELINE_DISCOVERY_WORKERS",
        )
        discovery_worker_count = min(max(1, discovery_worker_count), max(1, len(species_list)))
        logger.info(
            "=== Phase 1 & 2: Parallel Species Discovery & Bounded Sample Execution "
            "(requested workers=%d, active workers=%d, discovery workers=%d) ===",
            profile.requested_workers,
            profile.workers,
            discovery_worker_count,
        )
        quantified_by_species: Dict[str, int] = defaultdict(int)

        stale_seconds = int(os.environ.get("AMALGKIT_PIPELINE_STALE_DOWNLOADING_SECONDS", "0"))
        if stale_seconds > 0:
            reset_count = self.db.reset_stale_downloading(stale_seconds)
            if reset_count:
                logger.info(
                    "Reset %d stale downloading rows older than %d seconds to pending before resume.",
                    reset_count,
                    stale_seconds,
                )

        all_tasks: List[Dict[str, Any]] = []
        total_discovered = 0

        discovery_started = time.monotonic()
        discovery_failures: Dict[str, str] = {}
        discovered_by_species: Dict[str, List[Dict[str, Any]]] = {}
        logger.info(
            "=== Phase 1 start: discovering %d species with %d bounded workers; "
            "metadata/reference writes remain species-isolated ===",
            len(species_list),
            discovery_worker_count,
        )
        with concurrent.futures.ThreadPoolExecutor(max_workers=discovery_worker_count) as discovery_executor:
            futures: Dict[concurrent.futures.Future[Any], str] = {}
            discovery_start_times: Dict[str, float] = {}
            for config_name in species_list:
                discovery_start_times[config_name] = time.monotonic()
                futures[discovery_executor.submit(self.discover_species_tasks, config_name, max_gb, threads)] = (
                    config_name
                )
            for future in concurrent.futures.as_completed(futures):
                config_name = futures[future]
                try:
                    discovered_by_species[config_name] = future.result()
                except Exception as exc:
                    discovery_failures[config_name] = str(exc)
                    logger.exception("Fatal error discovering tasks for %s", config_name)
                    continue
                tasks = discovered_by_species[config_name]
                total_discovered += len(tasks)
                logger.info(
                    "Discovery complete for %s: %d tasks in %.1fs (total discovered: %d)",
                    config_name,
                    len(tasks),
                    time.monotonic() - discovery_start_times[config_name],
                    total_discovered,
                )

        discovery_elapsed = max(time.monotonic() - discovery_started, 1e-9)
        logger.info(
            "=== Phase 1 complete: %d species, %d tasks, %.1fs wall time, %.2f species/hour; "
            "discovery workers=%d ===",
            len(species_list) - len(discovery_failures),
            total_discovered,
            discovery_elapsed,
            len(species_list) / discovery_elapsed * 3600,
            discovery_worker_count,
        )
        if discovery_failures:
            details = "; ".join(f"{name}: {message}" for name, message in sorted(discovery_failures.items()))
            raise RuntimeError(
                "Species discovery failed for "
                f"{len(discovery_failures)} of {len(species_list)} configurations; "
                f"no sample tasks were submitted. {details}"
            )

        # Completion order is nondeterministic under parallel discovery. Keep
        # configured species order for stable task batches and reproducible
        # raw-input prioritization.
        for config_name in species_list:
            all_tasks.extend(discovered_by_species.get(config_name, []))

        with (
            concurrent.futures.ThreadPoolExecutor(max_workers=profile.workers) as executor,
            concurrent.futures.ThreadPoolExecutor(max_workers=profile.fasterq_slots) as fallback_executor,
        ):

            logger.info("=== Phase 2 start: all metadata parsed; execution continuing until queue exhaustion ===")

            # Schedule the global executor by raw-input state so mounted data
            # already acquired by earlier campaigns is quantitated and
            # reclaimed first.
            # The sort is stable, retaining metadata size order inside each
            # tier and preserving the original quant batch index.
            prioritized_tasks = _prioritize_tasks_by_raw_state(
                all_tasks,
                workers=profile.workers,
                fasterq_slots=profile.fasterq_slots,
            )
            reusable_count = sum(task["_raw_input_priority"] == 0 for task in prioritized_tasks)
            resumable_count = sum(task["_raw_input_priority"] == 1 for task in prioritized_tasks)
            logger.info(
                "Scheduling %d tasks: %d with existing raw inputs, %d resumable partials, %d new acquisitions; "
                "workers=%d (requested %d), quant_slots=%d, quant_threads=%d total/%d per worker, "
                "validation_slots=%d, fasterq_threads=%d, fasterq_slots=%d in a dedicated fallback queue, "
                "compression_threads=%d, peak_stage_threads=%d, max_in_flight=%d, host_cpus=%d.",
                len(prioritized_tasks),
                reusable_count,
                resumable_count,
                len(prioritized_tasks) - reusable_count - resumable_count,
                profile.workers,
                profile.requested_workers,
                profile.quant_slots,
                profile.effective_quant_threads,
                profile.quant_threads_per_worker,
                profile.validation_slots,
                profile.fasterq_threads,
                profile.fasterq_slots,
                profile.compression_threads,
                profile.peak_stage_threads,
                profile.max_in_flight,
                profile.host_cpu_count,
            )

            # Keep only a bounded number of futures submitted to the executor.
            # The old implementation queued the entire cohort at once; for a
            # 12,000-sample campaign that increased scheduler memory and made
            # interruption/recovery needlessly noisy without increasing
            # concurrency.
            task_iter = iter(prioritized_tasks)
            ready_primary_tasks: deque[Dict[str, Any]] = deque()
            fallback_queue: deque[Dict[str, Any]] = deque()
            primary_futures: Dict[concurrent.futures.Future[Any], Dict[str, Any]] = {}
            fallback_futures: Dict[concurrent.futures.Future[Any], Dict[str, Any]] = {}
            primary_exhausted = False

            def submit_next_primary() -> bool:
                nonlocal primary_exhausted
                if ready_primary_tasks:
                    task = ready_primary_tasks.popleft()
                elif not primary_exhausted:
                    try:
                        task = next(task_iter)
                    except StopIteration:
                        primary_exhausted = True
                        return False
                else:
                    return False
                future = executor.submit(
                    self.process_single_sample,
                    task["srr"],
                    task["batch_idx"],
                    task["fastq_dir"],
                    task["config_path"],
                    task["species_name"],
                    profile.quant_threads_per_worker,
                    task.get("expected_paired"),
                    # A cached SRA is deliberately handled in this lane first:
                    # it may be reusable immediately and, if corrupt, must get
                    # one ENA recovery attempt without repeating the same
                    # expensive local extraction in this invocation. Network
                    # acquisitions defer NCBI work so they never wait behind
                    # the single mounted-disk extraction slot.
                    defer_ncbi_fallback=not task["_uses_local_sra"],
                )
                primary_futures[future] = task
                return True

            def submit_fallbacks() -> None:
                while fallback_queue and len(fallback_futures) < profile.fasterq_slots:
                    task = fallback_queue.popleft()
                    self.db.set_state(
                        task["species_name"],
                        task["srr"],
                        "downloading",
                        error="running dedicated NCBI fallback after ENA acquisition failure",
                    )
                    future = fallback_executor.submit(
                        self._download_fastq_ncbi_only,
                        task["srr"],
                        task["fastq_dir"],
                        task.get("expected_paired"),
                    )
                    fallback_futures[future] = task

            for _ in range(min(profile.max_in_flight, len(prioritized_tasks))):
                if not submit_next_primary():
                    break

            completed = 0
            successful = 0
            execution_started = time.monotonic()
            fallback_error = "Download Failed (all sources: ENA FTP/HTTP, NCBI)"

            def record_terminal(
                task: Dict[str, Any], res: Dict[str, Any] | None, error: Exception | None = None
            ) -> None:
                nonlocal completed, successful
                completed += 1
                srr = task["srr"]
                sp = task["species_name"]
                if error is not None:
                    error_message = f"Unexpected sample task error: {error}"
                    self.db.set_state(sp, srr, "failed", error=error_message)
                    logger.error(
                        "[%d/%d] %s %s: Unexpected error — %s",
                        completed,
                        total_discovered,
                        sp,
                        srr,
                        error,
                    )
                    return
                if res is None:
                    logger.error(
                        "[%d/%d] %s %s: Unexpected empty task result",
                        completed,
                        total_discovered,
                        sp,
                        srr,
                    )
                    return
                if res["quantified"]:
                    successful += 1
                    quantified_by_species[sp] += 1
                    status = "Skipped (Done)" if res.get("skipped") else "Done"
                    elapsed = max(time.monotonic() - execution_started, 1e-9)
                    rate = successful / elapsed * 3600
                    logger.info(
                        "[%d/%d] %s %s: %s (finished-task rate %.2f/hour)",
                        completed,
                        total_discovered,
                        sp,
                        srr,
                        status,
                        rate,
                    )
                elif res["error"]:
                    logger.warning(
                        "[%d/%d] %s %s: Failed (%s)",
                        completed,
                        total_discovered,
                        sp,
                        srr,
                        res["error"],
                    )

            while primary_futures or fallback_futures or fallback_queue or ready_primary_tasks or not primary_exhausted:
                submit_fallbacks()
                while len(primary_futures) < profile.max_in_flight:
                    if not submit_next_primary():
                        break

                active_futures = set(primary_futures) | set(fallback_futures)
                if not active_futures:
                    if primary_exhausted and not fallback_queue and not ready_primary_tasks:
                        break
                    # This should be unreachable, but protects against a
                    # scheduler spin if future routing changes later.
                    logger.error("Execution queues are non-empty but no futures are active; stopping safely.")
                    break
                done, _ = concurrent.futures.wait(
                    active_futures,
                    return_when=concurrent.futures.FIRST_COMPLETED,
                )
                for future in done:
                    if future in primary_futures:
                        task = primary_futures.pop(future)
                        try:
                            res = future.result()
                        except Exception as exc:
                            record_terminal(task, None, exc)
                            continue
                        if res.get("fallback_deferred"):
                            fallback_queue.append(task)
                        else:
                            record_terminal(task, res)
                        continue

                    task = fallback_futures.pop(future)
                    try:
                        fallback_ready = bool(future.result())
                    except Exception as exc:
                        fallback_ready = False
                        logger.error(
                            "Dedicated NCBI fallback raised for %s %s: %s",
                            task["species_name"],
                            task["srr"],
                            exc,
                        )
                    if fallback_ready:
                        # Return the now-local sample to the ordinary lane.
                        # download_fastq() revalidates and reuses the FASTQ,
                        # then quantification competes under the shared CPU
                        # semaphore while the fallback executor immediately
                        # starts the next extraction.
                        ready_primary_tasks.appendleft(task)
                    else:
                        self.db.set_state(
                            task["species_name"],
                            task["srr"],
                            "failed",
                            error=fallback_error,
                        )
                        record_terminal(
                            task,
                            {
                                "quantified": False,
                                "error": fallback_error,
                            },
                        )

        logger.info(
            "=== Acquisition and quantification complete for this invocation; "
            "run the current downstream checkpoint launcher after the producer stops. ==="
        )

        # Global progress summary mapping
        logger.info("=== Final Results Summary ===")
        try:
            all_counts = self.db.get_counts()
            for sp, c in all_counts.items():
                total = sum(c.values())
                q = c.get("quantified", 0)
                f = c.get("failed", 0)
                p = c.get("pending", 0)
                d = c.get("downloaded", 0) + c.get("downloading", 0)
                pct = (q / total * 100) if total > 0 else 0
                logger.info(f"{sp.ljust(30)} | {q}/{total} ({pct:3.0f}%) | " f"Fail: {f:3} | Pend: {p:3} | DL: {d:3}")
        except Exception:
            pass

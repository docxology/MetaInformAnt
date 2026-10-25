"""Read-only, content-bound observations of a moving local acquisition campaign."""

from __future__ import annotations

import hashlib
import json
import math
import stat
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

from metainformant.rna.engine.acquisition_manifest import (
    load_snapshot,
    load_task_selection,
)
from metainformant.rna.engine.provenance import quantification_contract_id


class ObservationError(ValueError):
    pass


class _ChangedInput(ObservationError):
    pass


@dataclass(frozen=True)
class CampaignObservation:
    payload: dict[str, Any]
    sources: dict[str, bytes]


def _integer(value: Any, name: str) -> int:
    if type(value) is not int or value < 0:
        raise ObservationError(f"{name} must be a nonnegative integer")
    return value


def _number(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ObservationError(f"{name} must be finite")
    return float(value)


def _read(path: Path) -> bytes:
    if not stat.S_ISREG(path.stat().st_mode):
        raise ObservationError(f"observation input is not a regular file: {path}")
    with path.open("rb") as handle:
        data = handle.read(64 * 1024**2 + 1)
    if len(data) > 64 * 1024**2:
        raise ObservationError(f"observation input exceeds 64 MiB: {path}")
    return data


def _object(raw: bytes) -> dict[str, Any]:
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise ObservationError("expected a JSON object")
    return value


def _records(raw: bytes) -> tuple[list[dict[str, Any]], bool]:
    incomplete = bool(raw and not raw.endswith(b"\n"))
    complete = raw.rpartition(b"\n")[0] if incomplete else raw
    return [_object(line) for line in complete.splitlines() if line.strip()], incomplete


def _capture_campaign_observation_once(
    operations_root: Path,
    data_root: Path,
    *,
    selection_path: Path,
    journal_path: Path,
    manifest_path: Path | None = None,
    transfer_observation_path: Path | None = None,
    observed_at: float | None = None,
    max_rate_age_seconds: float = 300,
) -> CampaignObservation:
    """Capture stable controls and append-only prefixes without claiming release readiness.

    Startup counts remain producer-reported. Only unique new successes with a
    matching current restore proof, provenance contract and abundance digest
    increment that baseline. No live output, lock, database or producer is changed.
    """

    capture_started_epoch = time.time()
    capture_started_monotonic = time.monotonic()
    operations = operations_root.resolve(strict=True)
    root = data_root.resolve(strict=True)
    manifest = (manifest_path or operations / "frozen/manifest.jsonl").resolve(strict=True)
    now = _number(time.time() if observed_at is None else observed_at, "observed_at")
    max_age = _number(max_rate_age_seconds, "max_rate_age_seconds")
    if max_age < 0:
        raise ObservationError("max_rate_age_seconds must be nonnegative")
    sources: dict[str, bytes] = {}
    paths: dict[str, Path] = {}
    append_only: set[str] = set()
    external_files: dict[str, dict[str, Any]] = {}
    external_stats: dict[Path, tuple[int, ...]] = {}

    def capture(name: str, path: Path, *, prefix: bool = False) -> bytes:
        sources[name] = _read(path)
        paths[name] = path
        if prefix:
            append_only.add(name)
        return sources[name]

    runtime = _object(capture("runtime.json", operations / "active_runtime.json"))
    started = _number(runtime["started_at_epoch"], "runtime start")
    pid = _integer(runtime["pid"], "runtime pid")
    if pid == 0 or not 0 < started <= now:
        raise ObservationError("invalid runtime identity/time")
    scope_path = operations / "active_scope.json"
    scope = _object(capture("scope.json", scope_path))
    if scope_path.stat().st_mtime < started:
        raise ObservationError("scope predates current runtime")
    live = _object(capture("checkpoint_counter.json", operations / "live_progress.json"))
    manifest_raw = capture("manifest.jsonl", manifest)
    snapshot_raw = capture("snapshot.json", manifest.with_name("snapshot.json"))
    snapshot, tasks = load_snapshot(manifest)
    if hashlib.sha256(manifest_raw).hexdigest() != snapshot["manifest_sha256"]:
        raise ObservationError("manifest changed during observation")
    by_id = {task.get("task_id", f"{task['species']}/{task['accession']}"): task for task in tasks}
    total = _integer(scope["scope_total"], "scope_total")
    eligible = _integer(scope["eligible"], "eligible")
    baseline = _integer(scope["verified_existing"], "verified_existing")
    unsupported = scope["unsupported"]
    if not isinstance(unsupported, dict) or set(unsupported) - by_id.keys():
        raise ObservationError("unsupported tasks do not belong to frozen manifest")
    if total != len(tasks) or eligible + len(unsupported) != total or baseline > eligible:
        raise ObservationError("scope counts disagree with frozen manifest")
    selection_path = selection_path.resolve(strict=True)
    journal = journal_path.resolve(strict=True)
    if selection_path.stat().st_mtime < started or journal.stat().st_mtime < started:
        raise ObservationError("partition or journal predates current runtime")
    selection_raw = capture("selection.json", selection_path)
    selection, selected = load_task_selection(
        selection_path,
        snapshot,
        tasks,
        snapshot_sha256=hashlib.sha256(snapshot_raw).hexdigest(),
    )
    if selection is None or set(selection["task_ids"]) & unsupported.keys():
        raise ObservationError("partition includes unsupported tasks")
    selected_ids = {task.get("task_id", f"{task['species']}/{task['accession']}") for task in selected}
    if not journal.resolve().is_relative_to(root):
        raise ObservationError("task journal escapes data root")
    rows, incomplete = _records(capture("task_results.jsonl", journal, prefix=True))
    identities: set[str] = set()
    successful: list[dict[str, Any]] = []
    deferred = 0
    failures: dict[str, int] = {}
    for row in rows:
        identity = row.get("task_id")
        if identity not in selected_ids or identity in identities:
            raise ObservationError("duplicate or out-of-partition journal task")
        task = by_id[identity]
        if row.get("species") != task["species"] or row.get("srr") != task["accession"]:
            raise ObservationError("journal task identity fields disagree")
        if (
            type(row.get("quantified")) is not bool
            or type(row.get("skipped", False)) is not bool
            or type(row.get("download_deferred", False)) is not bool
        ):
            raise ObservationError("journal success/reuse fields must be booleans")
        # The membership check above admits only selected (string) task identities.
        identities.add(cast(str, identity))
        if row["quantified"]:
            successful.append(row)
        else:
            deferred += int(row.get("download_deferred", False))
            reason = str(row.get("error") or "unrecorded")
            failures[reason] = failures.get(reason, 0) + 1

    confirmed: list[str] = []
    unconfirmed: dict[str, str] = {}
    for row in successful:
        if row.get("skipped"):
            continue
        identity = row["task_id"]
        task = by_id[identity]
        sample = root / task["species"] / "work/quant" / task["accession"]
        try:
            proof_path = sample / ".local_restore_verified.json"
            provenance_path = sample / ".metainformant_quant_provenance.json"
            for path in (proof_path, provenance_path):
                if path.is_symlink() or not path.resolve().is_relative_to(root):
                    raise ObservationError("checkpoint evidence escapes data root")
            proof = _object(capture(f"proofs/{identity}.json", proof_path))
            provenance = _object(capture(f"provenance/{identity}.json", provenance_path))
            verified = _number(proof["verified_at_epoch"], "restore time")
            contract = proof["contract_id"]
            if (
                not isinstance(contract, str)
                or len(contract) != 64
                or contract != provenance["quant_contract_id"]
                or contract != quantification_contract_id(provenance)
                or proof["reference_index_sha256"] != task["reference_index_sha256"]
                or provenance["config_sha256"] != task["config_sha256"]
                or provenance["run_accession"] != task["accession"]
                or provenance["species"] != task["species"]
                or not started <= verified <= now
            ):
                raise ObservationError("restore/provenance binding is not current")
            filename = provenance["quantification_file"]
            if not isinstance(filename, str) or Path(filename).name != filename:
                raise ObservationError("abundance filename is not local to sample")
            abundance = sample / filename
            if abundance.is_symlink() or not abundance.resolve().is_relative_to(root):
                raise ObservationError("abundance path escapes data root")
            before = abundance.stat()
            signature = (
                before.st_dev,
                before.st_ino,
                before.st_size,
                before.st_mtime_ns,
                before.st_ctime_ns,
            )
            if not stat.S_ISREG(before.st_mode) or before.st_size == 0:
                raise ObservationError("abundance is not a nonempty regular file")
            with abundance.open("rb") as handle:
                digest = hashlib.file_digest(handle, "sha256").hexdigest()
            after = abundance.stat()
            if signature != (
                after.st_dev,
                after.st_ino,
                after.st_size,
                after.st_mtime_ns,
                after.st_ctime_ns,
            ):
                raise ObservationError("abundance changed during hashing")
            if digest != provenance["quantification_file_sha256"]:
                raise ObservationError("abundance digest disagrees with checkpoint provenance")
            external_files[f"abundance/{identity}.tsv"] = {
                "path": str(abundance),
                "sha256": digest,
                "bytes": before.st_size,
                "captured": False,
                "stat_signature": [str(value) for value in signature],
            }
            external_stats[abundance] = signature
            confirmed.append(identity)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            unconfirmed[identity] = str(exc)
    current = baseline + len(confirmed)
    if current > eligible:
        raise ObservationError("confirmed lower bound exceeds eligible scope")

    transfer: dict[str, Any] = {"available": False, "reason": "no transfer observation"}
    observation_path = transfer_observation_path
    if observation_path is not None and observation_path.exists():
        observation = _object(capture("transfer.json", observation_path))
        measured = _number(observation["observed_at_epoch"], "transfer observation time")
        seconds = _number(observation["seconds"], "transfer window")
        growth = _integer(observation["retained_transfer_growth_bytes"], "retained bytes")
        if observation["pid"] != pid or not started <= measured <= now or measured - seconds < started:
            transfer["reason"] = "transfer observation belongs to another runtime/time"
        elif now - measured > max_age:
            transfer["reason"] = "transfer observation is stale"
        elif seconds <= 0:
            raise ObservationError("transfer window must be positive")
        else:
            transfer = {
                "available": True,
                "seconds": seconds,
                "observed_at_epoch": measured,
                "MiB_s_lower_bound": growth / seconds / 1024**2,
                "retained_bytes": growth,
                "tracked_transfers": _integer(observation["tracked_curl_transfers"], "tracked transfers"),
            }
    reclamation, raw_tail = _records(
        capture("raw_reclamation.jsonl", operations / "raw_reclamation.jsonl", prefix=True)
    )
    reclaimed = sum(_integer(row["bytes_freed"], "reclaimed bytes") for row in reclamation)
    for name, path in paths.items():
        latest = _read(path)
        if not latest.startswith(sources[name]) if name in append_only else latest != sources[name]:
            raise _ChangedInput(f"source changed during capture: {path}")
    for path, captured_signature in external_stats.items():
        latest_stat = path.stat()
        if captured_signature != (
            latest_stat.st_dev,
            latest_stat.st_ino,
            latest_stat.st_size,
            latest_stat.st_mtime_ns,
            latest_stat.st_ctime_ns,
        ):
            raise _ChangedInput(f"abundance changed before capture finished: {path}")
    payload = {
        "schema": "metainformant.rna.campaign_observation.v1",
        "observed_at_epoch": now,
        "capture_started_at_epoch": capture_started_epoch,
        "capture_finished_at_epoch": time.time(),
        "capture_elapsed_seconds": time.monotonic() - capture_started_monotonic,
        "data_root": str(root),
        "operations_root": str(operations),
        "manifest_sha256": snapshot["manifest_sha256"],
        "snapshot_sha256": hashlib.sha256(snapshot_raw).hexdigest(),
        "selection_sha256": hashlib.sha256(selection_raw).hexdigest(),
        "runtime": runtime,
        "run_attribution": "explicit partition and journal; no cryptographic run-to-root binding",
        "scientific_release_ready": False,
        "scope": {
            "total": total,
            "eligible": eligible,
            "unsupported": len(unsupported),
        },
        "counts": {
            "producer_audited_at_start": baseline,
            "current_lower_bound": current,
            "supported_unconfirmed": eligible - current,
            "journal_successes": len(successful),
            "reused": sum(bool(row.get("skipped")) for row in successful),
            "new_restore_bound": len(confirmed),
            "new_unconfirmed": len(unconfirmed),
            "failed_attempts": len(rows) - len(successful),
            "deferred_attempts": deferred,
            "failed_nondeferred_attempts": len(rows) - len(successful) - deferred,
        },
        "confirmed_new_task_ids": sorted(confirmed),
        "unconfirmed_new_tasks": unconfirmed,
        "failure_reasons": failures,
        "journal_incomplete_tail": incomplete,
        "reclamation_incomplete_tail": raw_tail,
        "reclaimed_bytes": reclaimed,
        "reclamation_errors": sum(bool(row["errors"]) for row in reclamation),
        "checkpoint_counter": live,
        "transfer": transfer,
        "limitations": [
            "Startup counts are producer-reported audited counts, not a new all-output audit.",
            "Control files lack a cryptographic run-to-root binding; membership and stable bytes are checked.",
            "Captured append-only prefixes are an operational observation, not a transaction across all live files.",
            "Run accessions and transfer rates do not establish biological replication or a scientific release.",
        ],
        "source_files": {
            **{
                name: {
                    "path": str(paths[name]),
                    "sha256": hashlib.sha256(raw).hexdigest(),
                    "bytes": len(raw),
                    "captured": True,
                    "append_only_prefix": name in append_only,
                    **(
                        {
                            "complete_prefix_bytes": raw.rfind(b"\n") + 1,
                            "excluded_tail_bytes": len(raw) - (raw.rfind(b"\n") + 1),
                        }
                        if name in append_only
                        else {}
                    ),
                }
                for name, raw in sources.items()
            },
            **external_files,
        },
    }
    return CampaignObservation(payload, sources)


def capture_campaign_observation(
    operations_root: Path,
    data_root: Path,
    *,
    selection_path: Path,
    journal_path: Path,
    manifest_path: Path | None = None,
    transfer_observation_path: Path | None = None,
    observed_at: float | None = None,
    max_rate_age_seconds: float = 300,
    stability_attempts: int = 3,
) -> CampaignObservation:
    if type(stability_attempts) is not int or stability_attempts < 1:
        raise ObservationError("stability_attempts must be a positive integer")
    last_error: _ChangedInput | None = None
    for _ in range(stability_attempts):
        try:
            return _capture_campaign_observation_once(
                operations_root,
                data_root,
                selection_path=selection_path,
                journal_path=journal_path,
                manifest_path=manifest_path,
                transfer_observation_path=transfer_observation_path,
                observed_at=observed_at,
                max_rate_age_seconds=max_rate_age_seconds,
            )
        except _ChangedInput as exc:
            last_error = exc
    raise ObservationError(f"source changed across {stability_attempts} capture attempts") from last_error

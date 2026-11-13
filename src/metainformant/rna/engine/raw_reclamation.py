"""Reclaim raw inputs only for verified, reusable quantification outputs.

This is intentionally separate from broad cleanup commands.  It never uses a
progress-table state as evidence: each sample's quantification directory is
classified with the content-addressed reuse contract
(:func:`metainformant.rna.engine.provenance.build_species_reuse_context`), and
only ``current`` or ``version_drift_compatible`` outputs release their inputs.
With ``--cohort`` the frozen cohort's configuration and kallisto index digests
are the expectation, and permanently excluded tasks can release theirs too.

Scopes (each per accession, see :mod:`metainformant.rna.engine.raw_cleanup`):

* always: FASTQ/SRA files and failed partials under ``<species>/work/getfastq``;
* ``--include-staging``: orphaned ``getfastq/.fasterq_tmp/<run>`` workspaces;
* ``--include-sra-cache`` (needs ``--cohort``): the shared VDB cache entries
  ``.sra-cache/sra/<run>.*`` and ``.sra-cache/prefetch/<run>/``, only when every
  cohort task naming the accession is reusable or permanently excluded.

The default is a read-only dry run that prints a JSON summary; ``--manifest``
also writes the full plan.  ``--execute`` refuses while a producer is alive or
a campaign/maintenance lock is held, holds the data-root maintenance lock for
the whole pass, writes the planned manifest before removing anything, and
replaces it with the executed record at the end.
"""

from __future__ import annotations

import argparse
import json
import os
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from metainformant.rna.core.sample_utils import find_quantification_file
from metainformant.rna.engine.campaign_lock import writer_status
from metainformant.rna.engine.cohort import (
    DataRootBusy,
    FrozenCohort,
    data_root_maintenance_lock,
    load_frozen_cohort,
)
from metainformant.rna.engine.provenance import (
    REUSABLE_QUANT_STATUSES,
    SpeciesReuseContext,
    build_species_reuse_context,
)
from metainformant.rna.engine.raw_cleanup import (
    reclaim_sample_raw_inputs,
    reclaim_sample_sra_cache,
    sra_cache_root,
)
from metainformant.rna.engine.species import (
    configured_data_root,
    discover_species_config_names,
)
from metainformant.rna.engine.streaming_orchestrator import _species_work_dir

SCHEMA = "metainformant.rna.raw_reclamation.v2"


def _default_config_dir() -> Path:
    """The Hymenoptera project configuration directory of a source checkout."""

    return Path(__file__).resolve().parents[4] / "projects" / "hymenoptera_amalgkit" / "config" / "amalgkit"


REASON_REUSABLE = "reusable_quantification"
REASON_EXCLUDED = "permanent_exclusion"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--data-root",
        type=Path,
        help="Selected Amalgkit data root (default: AMALGKIT_DATA_ROOT or output/amalgkit)",
    )
    parser.add_argument(
        "--config-dir",
        type=Path,
        default=_default_config_dir(),
    )
    parser.add_argument("--species", action="append", help="Restrict to one or more species identifiers")
    parser.add_argument("--cohort", type=Path, help="Frozen cohort directory: bind reuse to its config/index digests")
    parser.add_argument("--cohort-definition", type=Path, help="Cohort definition file (default: <cohort>/COHORT.json)")
    parser.add_argument("--cohort-definition-sha256", help="Refuse unless the cohort definition has this sha256")
    parser.add_argument(
        "--include-staging", action="store_true", help="Also remove orphaned getfastq/.fasterq_tmp/<run> workspaces"
    )
    parser.add_argument(
        "--include-sra-cache",
        action="store_true",
        help="Also remove the shared .sra-cache entries of released accessions (needs --cohort)",
    )
    parser.add_argument(
        "--include-excluded",
        action="store_true",
        help="Also release inputs of the cohort's permanently excluded tasks (needs --cohort)",
    )
    parser.add_argument("--workers", type=int, default=8, help="Parallel sample classifications (default 8)")
    parser.add_argument("--execute", action="store_true", help="Actually remove inputs; default is a dry run")
    parser.add_argument(
        "--manifest",
        type=Path,
        help="JSON plan/audit path (dry run: optional, nothing is written without it; "
        "--execute default: <data-root>/results/raw_reclamation_<UTC>.json)",
    )
    return parser


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _species_plan(
    config_dir: Path,
    cohort: FrozenCohort | None,
    requested: set[str],
) -> list[tuple[str, str]]:
    """Return ``(species, config_name)`` pairs to scan, refusing unknown requests."""

    if cohort is not None:
        pairs = [(entry.species, entry.config_name) for entry in cohort.species.values()]
    else:
        pairs = [
            (name.removeprefix("amalgkit_").removesuffix(".yaml"), name)
            for name in discover_species_config_names(config_dir)
        ]
    known = {species for species, _ in pairs}
    unknown = requested - known
    if unknown:
        raise SystemExit(f"unknown species: {', '.join(sorted(unknown))}")
    if requested:
        pairs = [pair for pair in pairs if pair[0] in requested]
    if not pairs:
        raise SystemExit(f"no runnable species configurations found under {config_dir}")
    return pairs


def _reuse_context(
    species: str, config_name: str, config_dir: Path, cohort: FrozenCohort | None
) -> SpeciesReuseContext:
    config_path = config_dir / config_name
    if cohort is None:
        return build_species_reuse_context(species, config_path)
    entry = cohort.entry(species)
    context = build_species_reuse_context(
        species,
        config_path,
        cohort_id=cohort.cohort_id,
        reference_index_sha256=entry.index_sha256,
        reference_index_name=entry.index_name,
    )
    if context.config_sha256 != entry.config_sha256:
        raise SystemExit(
            f"{config_path} does not hash to the cohort configuration digest for {species}; refusing to reclaim"
        )
    return context


def _classify(context: SpeciesReuseContext, sample_dir: Path) -> dict[str, Any]:
    run = sample_dir.name
    quant_file = find_quantification_file(sample_dir, run)
    if quant_file is None or not quant_file.is_file() or quant_file.stat().st_size <= 100:
        return {"status": "missing_quantification", "reason": "no non-empty abundance table"}
    try:
        return context.classify(sample_dir, run)
    except OSError as error:
        return {"status": "unreadable", "reason": str(error)}


def _scan_species(species: str, context: SpeciesReuseContext, workers: int) -> tuple[Path, dict[str, dict[str, Any]]]:
    """Classify every accession-named quantification directory of one species."""

    work_dir = _species_work_dir(species)
    quant_root = work_dir / "quant"
    if not quant_root.is_dir():
        return work_dir, {}
    sample_dirs = sorted(
        path
        for path in quant_root.iterdir()
        if path.is_dir() and not path.is_symlink() and path.name[:3] in {"SRR", "ERR", "DRR"}
    )
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        verdicts = list(pool.map(lambda path: _classify(context, path), sample_dirs))
    return work_dir, {path.name: verdict for path, verdict in zip(sample_dirs, verdicts)}


def _sum_into(totals: dict[str, int], scope: str, result: dict[str, Any]) -> None:
    totals[f"{scope}_files"] += int(result["files_deleted"])
    totals[f"{scope}_bytes"] += int(result["bytes_freed"])
    totals[f"{scope}_allocated_bytes"] += int(result["allocated_bytes_freed"])
    totals["errors"] += len(result["errors"])


def _remaining_cache(cache_root: Path) -> dict[str, int]:
    """Count accession-named cache entries that remain after the pass."""

    remaining = {"sra_files": 0, "sra_allocated_bytes": 0, "prefetch_dirs": 0}
    flat = cache_root / "sra"
    if flat.is_dir():
        for path in flat.iterdir():
            if path.is_file() and path.name[:3] in {"SRR", "ERR", "DRR"}:
                remaining["sra_files"] += 1
                remaining["sra_allocated_bytes"] += path.lstat().st_blocks * 512
    prefetch = cache_root / "prefetch"
    if prefetch.is_dir():
        remaining["prefetch_dirs"] = sum(1 for path in prefetch.iterdir() if path.name[:3] in {"SRR", "ERR", "DRR"})
    return remaining


def _describe_writers(status: dict[str, Any]) -> str:
    owner = status.get("owner") or {}
    parts = []
    if status.get("held") and not status.get("inside_lock"):
        parts.append(f"campaign lock held by {owner.get('kind', 'unknown')} pid {owner.get('pid', '?')}")
    for holder in status.get("other_writers") or []:
        parts.append(f"{holder.get('lock')}: {holder.get('reason')} (pid {holder.get('pid')})")
    return "; ".join(parts) or "unknown writer"


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if (args.include_sra_cache or args.include_excluded) and args.cohort is None:
        raise SystemExit("--include-sra-cache and --include-excluded need --cohort (the full task universe)")
    if args.cohort is None and (args.cohort_definition or args.cohort_definition_sha256):
        raise SystemExit("--cohort-definition and --cohort-definition-sha256 need --cohort")
    if args.data_root:
        os.environ["AMALGKIT_DATA_ROOT"] = str(args.data_root.expanduser().resolve())
    data_root = configured_data_root().expanduser().resolve()
    config_dir = args.config_dir.expanduser().resolve()
    cohort = (
        load_frozen_cohort(
            args.cohort.expanduser(),
            definition_path=args.cohort_definition,
            expected_definition_sha256=args.cohort_definition_sha256,
        )
        if args.cohort is not None
        else None
    )
    pairs = _species_plan(config_dir, cohort, set(args.species or []))
    if args.execute:
        # Refuse before the (long) classification pass; the writer lock below is the binding guard.
        status = writer_status(data_root)
        if not status["quiescent"]:
            raise SystemExit(f"refusing to reclaim: another writer holds the data root ({_describe_writers(status)})")

    contexts: dict[str, SpeciesReuseContext] = {}
    work_dirs: dict[str, Path] = {}
    verdicts: dict[str, dict[str, dict[str, Any]]] = {}
    for species, config_name in pairs:
        contexts[species] = _reuse_context(species, config_name, config_dir, cohort)
        work_dirs[species], verdicts[species] = _scan_species(species, contexts[species], args.workers)

    # Released accessions per species, with the evidence that releases them.
    released: dict[str, dict[str, str]] = defaultdict(dict)
    status_counts: Counter[str] = Counter()
    for species, by_run in verdicts.items():
        for run, verdict in by_run.items():
            status_counts[str(verdict.get("status"))] += 1
            if verdict.get("status") in REUSABLE_QUANT_STATUSES:
                released[species][run] = REASON_REUSABLE
        if cohort is not None and args.include_excluded:
            for run in cohort.exclusions_for(species):
                released[species].setdefault(run, REASON_EXCLUDED)

    # The shared cache is keyed by accession only: an accession is released
    # there only when every cohort task naming it is released.
    cache_release: set[str] = set()
    cache_blocked: dict[str, list[str]] = {}
    if cohort is not None and args.include_sra_cache:
        species_by_run: dict[str, set[str]] = defaultdict(set)
        for task in cohort.tasks:
            species_by_run[task.accession].add(task.species)
        for species, runs in released.items():
            for run in runs:
                owners = species_by_run.get(run, {species})
                blocked = []
                for owner in sorted(owners):
                    if run in released.get(owner, {}):
                        continue
                    if owner not in contexts:
                        # Outside the requested species: classify that one directory now.
                        entry = cohort.entry(owner)
                        contexts[owner] = _reuse_context(owner, entry.config_name, config_dir, cohort)
                        work_dirs[owner] = _species_work_dir(owner)
                    sample_dir = work_dirs[owner] / "quant" / run
                    excluded = args.include_excluded and cohort.is_excluded(owner, run)
                    if excluded or (
                        sample_dir.is_dir()
                        and _classify(contexts[owner], sample_dir).get("status") in REUSABLE_QUANT_STATUSES
                    ):
                        continue
                    blocked.append(owner)
                if blocked:
                    cache_blocked[run] = blocked
                else:
                    cache_release.add(run)

    cache_root = sra_cache_root(data_root)
    dry_run = not args.execute
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest_path = args.manifest.expanduser().resolve() if args.manifest else None
    if args.execute and manifest_path is None:
        manifest_path = data_root / "results" / f"raw_reclamation_{timestamp}.json"
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "data_root": str(data_root),
        "config_dir": str(config_dir),
        "cohort_id": cohort.cohort_id if cohort else None,
        "cohort_definition_sha256": cohort.definition_sha256 if cohort else None,
        "scopes": {
            "raw": True,
            "fasterq_staging": bool(args.include_staging),
            "sra_cache": bool(args.include_sra_cache),
            "permanent_exclusions": bool(args.include_excluded),
        },
        "mode": "execute-planned" if args.execute else "dry-run",
        "classification_counts": dict(sorted(status_counts.items())),
        "released": {species: dict(sorted(runs.items())) for species, runs in sorted(released.items())},
        "sra_cache_blocked": dict(sorted(cache_blocked.items())),
    }

    def run_pass() -> dict[str, Any]:
        totals: Counter[str] = Counter()
        records: list[dict[str, Any]] = []
        for species, runs in sorted(released.items()):
            for run, reason in sorted(runs.items()):
                record: dict[str, Any] = {"species": species, "run_accession": run, "reason": reason}
                raw = reclaim_sample_raw_inputs(
                    work_dirs[species], run, dry_run=dry_run, include_staging=args.include_staging
                )
                _sum_into(totals, "raw", raw)
                record["raw"] = raw
                if run in cache_release:
                    cache = reclaim_sample_sra_cache(cache_root, run, dry_run=dry_run)
                    _sum_into(totals, "sra_cache", cache)
                    record["sra_cache"] = cache
                    cache_release.discard(run)
                parts = [raw, record.get("sra_cache") or {}]
                if any(part.get("files_deleted") or part.get("errors") for part in parts):
                    records.append(record)
        totals["released_accessions"] = sum(len(runs) for runs in released.values())
        totals["accessions_with_files"] = len(records)
        return {"totals": dict(sorted(totals.items())), "records": records}

    if dry_run:
        payload.update(run_pass())
        payload["sra_cache_remaining_before"] = _remaining_cache(cache_root) if args.include_sra_cache else None
        if manifest_path is not None:
            _write_json_atomic(manifest_path, payload)
    else:
        try:
            with data_root_maintenance_lock(data_root, purpose="raw-reclamation"):
                assert manifest_path is not None
                _write_json_atomic(manifest_path, payload)
                payload.update(run_pass())
                payload["mode"] = "execute"
                payload["completed_at"] = datetime.now(timezone.utc).isoformat()
                payload["sra_cache_remaining_after"] = _remaining_cache(cache_root) if args.include_sra_cache else None
                _write_json_atomic(manifest_path, payload)
        except DataRootBusy as busy:
            raise SystemExit(f"refusing to reclaim: {busy}") from None

    summary = {
        "manifest": str(manifest_path) if manifest_path else None,
        "mode": payload["mode"],
        "classification_counts": payload["classification_counts"],
        **payload["totals"],
        "sra_cache_blocked": len(cache_blocked),
    }
    print(json.dumps(summary, sort_keys=True))
    return 1 if payload["totals"].get("errors") else 0

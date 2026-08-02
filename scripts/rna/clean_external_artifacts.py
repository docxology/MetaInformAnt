#!/usr/bin/env python3
"""Clear superseded external artifacts without touching current outputs.

The command is intentionally conservative: it inventories before changing
anything, archives superseded Amalgkit stage directories and logs under a dated
data-root archive, and never recursively targets the data root itself.  The
default is a dry run; ``--execute`` is required for the archive operation.
Legacy per-sample quantification directories can be selected explicitly with
``--archive-legacy-quant``. That mode targets only directories containing a
plain ``abundance.tsv`` without a current provenance sidecar.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

SUPERSEDED_WORK_DIRS = {"config", "csca", "curate", "logs", "validation"}
SUPERSEDED_DOWNLOAD_DIRS = {"ete4", "ete_taxonomy"}
STALE_RESULT_NAMES = {
    "pipeline_status_report.md",
    "pipeline_status_report.tsv",
    "pipeline_status_chart.png",
    "pipeline_status_chart_caption.md",
    "manuscript_status.md",
    "evidence_manifest.json",
}


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--execute", action="store_true", help="Archive the listed artifacts; default is dry-run")
    parser.add_argument(
        "--archive-legacy-quant",
        action="store_true",
        help="Also archive per-sample quant directories with plain abundance.tsv and no current provenance",
    )
    parser.add_argument(
        "--legacy-quant-only",
        action="store_true",
        help="Select only legacy quant directories; requires --archive-legacy-quant",
    )
    parser.add_argument("--clear-stale-results", action="store_true", help="Also archive stale report/cross-species outputs")
    parser.add_argument(
        "--archive-superseded-alias",
        action="append",
        default=[],
        help="Archive an explicitly named top-level superseded alias after safety checks",
    )
    parser.add_argument(
        "--keep-path",
        action="append",
        default=[],
        help="Keep an exact data-root-relative path (repeat for active/evidence paths)",
    )
    parser.add_argument("--manifest", type=Path, help="JSON audit path")
    return parser


def _candidates(data_root: Path, include_stale_results: bool, keep_paths: set[Path]) -> list[Path]:
    archive = data_root / "archive"
    selected: set[Path] = set()

    def add(path: Path) -> None:
        resolved = path.resolve(strict=False)
        if any(resolved == kept or kept in resolved.parents for kept in keep_paths):
            return
        selected.add(path)

    for path in data_root.rglob("*.log"):
        if archive not in path.parents:
            add(path)
    for path in data_root.rglob("*.heartbeat.json"):
        if archive not in path.parents:
            add(path)
    for work_dir in data_root.glob("*/work"):
        for name in SUPERSEDED_WORK_DIRS:
            path = work_dir / name
            if path.exists() or path.is_symlink():
                add(path)
        downloads = work_dir / "downloads"
        for name in SUPERSEDED_DOWNLOAD_DIRS:
            path = downloads / name
            if path.exists() or path.is_symlink():
                add(path)
    if include_stale_results:
        results = data_root / "results"
        for name in STALE_RESULT_NAMES:
            path = results / name
            if path.exists() or path.is_symlink():
                add(path)
        for path in results.iterdir() if results.is_dir() else []:
            if path.name.startswith(".") or "aborted" in path.name:
                add(path)
        cross_species = data_root / "cross_species"
        if cross_species.exists() or cross_species.is_symlink():
            add(cross_species)
    return sorted(selected)


def _legacy_quant_candidates(data_root: Path, keep_paths: set[Path]) -> list[Path]:
    """Return exact legacy quant directories eligible for archival.

    A plain ``abundance.tsv`` is not sufficient evidence for the current
    workflow. Only direct children of ``*/work/quant`` with that legacy
    filename and without the current provenance sidecar are selected. Empty,
    partial, accession-qualified, and current-provenance directories remain
    available for resume and downstream validation.
    """

    selected: list[Path] = []
    for quant_root in sorted(data_root.glob("*/work/quant")):
        if not quant_root.is_dir() or quant_root.is_symlink():
            continue
        for sample_dir in sorted(quant_root.iterdir()):
            if (
                not sample_dir.is_dir()
                or sample_dir.is_symlink()
                or sample_dir.name.startswith(".")
                or not (sample_dir / "abundance.tsv").is_file()
                or (sample_dir / ".metainformant_quant_provenance.json").is_file()
            ):
                continue
            resolved = sample_dir.resolve(strict=False)
            if any(resolved == kept or kept in resolved.parents for kept in keep_paths):
                continue
            selected.append(sample_dir)
    return selected


def _active_writer_lock(data_root: Path) -> Path | None:
    """Return an active campaign/finalization lock, if one exists."""

    for lock_name in (".full_campaign.lock", ".finalization.lock"):
        lock_dir = data_root / "results" / lock_name
        if lock_dir.exists():
            return lock_dir
    return None


def _size(path: Path) -> int:
    if path.is_file() or path.is_symlink():
        try:
            return path.stat().st_size
        except OSError:
            return 0
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def _sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError:
        return None
    return digest.hexdigest()


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    data_root = args.data_root.expanduser().resolve()
    if not data_root.is_dir():
        raise SystemExit(f"data root is not a directory: {data_root}")
    keep_paths: set[Path] = set()
    for raw_path in args.keep_path:
        path = Path(raw_path).expanduser()
        resolved = (path if path.is_absolute() else data_root / path).resolve(strict=False)
        if resolved != data_root and data_root not in resolved.parents:
            raise SystemExit(f"--keep-path must be inside the data root: {raw_path}")
        keep_paths.add(resolved)
    if args.legacy_quant_only and not args.archive_legacy_quant:
        raise SystemExit("--legacy-quant-only requires --archive-legacy-quant")
    if args.legacy_quant_only and (args.clear_stale_results or args.archive_superseded_alias):
        raise SystemExit("--legacy-quant-only cannot be combined with stale-results or alias archival")
    if args.archive_legacy_quant and args.execute:
        active_lock = _active_writer_lock(data_root)
        if active_lock is not None:
            raise SystemExit(
                "refusing legacy quant archival while an active writer lock exists: "
                f"{active_lock}"
            )

    candidates = [] if args.legacy_quant_only else _candidates(data_root, args.clear_stale_results, keep_paths)
    if args.archive_legacy_quant:
        candidates.extend(_legacy_quant_candidates(data_root, keep_paths))
        candidates = sorted(set(candidates))
    alias_roots: list[Path] = []
    for alias in args.archive_superseded_alias:
        alias_root = (data_root / alias).resolve()
        if alias_root.parent != data_root or not alias_root.is_dir():
            raise SystemExit(
                "superseded alias must be an existing direct child of data root: "
                f"{alias_root}"
            )
        for path in alias_root.rglob("*"):
            if path.name in {".metainformant_quant_provenance.json", ".metainformant_downstream_provenance.json"}:
                raise SystemExit(f"refusing to archive alias containing current evidence: {path}")
            if path.is_file() and (
                path.name.endswith((".fastq.gz", ".fastq", ".fq.gz", ".fq", ".sra", ".sra.part"))
                or ".fastq.gz.invalid" in path.name
            ):
                raise SystemExit(f"refusing to archive alias containing raw inputs: {path}")
        alias_roots.append(alias_root)
    if alias_roots:
        candidates = [path for path in candidates if not any(alias in path.parents for alias in alias_roots)]
        candidates.extend(
            alias
            for alias in alias_roots
            if not any(alias == kept or alias in kept.parents for kept in keep_paths)
        )
        candidates = sorted(set(candidates))
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive = data_root / "archive" / f"superseded_cleanup_{timestamp}"
    records: list[dict[str, object]] = []
    for source in candidates:
        record: dict[str, object] = {
            "source": str(source),
            "relative_path": str(source.relative_to(data_root)),
            "bytes": _size(source),
            "sha256": _sha256(source),
        }
        destination = archive / source.relative_to(data_root)
        if args.execute:
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(source), str(destination))
            record["action"] = "archived"
            record["destination"] = str(destination)
        else:
            record["action"] = "would_archive"
            record["destination"] = str(destination)
        records.append(record)

    manifest = args.manifest or data_root / "results" / f"superseded_cleanup_{timestamp}.json"
    manifest = manifest.expanduser()
    manifest.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "metainformant.rna.external_cleanup.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "data_root": str(data_root),
        "mode": "execute" if args.execute else "dry-run",
        "archive_legacy_quant": args.archive_legacy_quant,
        "legacy_quant_only": args.legacy_quant_only,
        "include_stale_results": args.clear_stale_results,
        "keep_paths": sorted(str(path) for path in keep_paths),
        "archive": str(archive),
        "candidate_count": len(candidates),
        "records": records,
    }
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"manifest": str(manifest), "mode": payload["mode"], "candidate_count": len(candidates)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

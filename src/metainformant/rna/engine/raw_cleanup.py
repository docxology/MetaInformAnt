"""Safe per-sample reclamation of transient RNA-seq inputs.

The streaming campaign keeps quantification outputs and provenance under the
selected data root.  Once a *current-contract* quantification has succeeded,
the corresponding FASTQ/SRA inputs are reproducible transient inputs and may
be removed to keep the external volume usable.  This module deliberately
touches only exact sample directories and known raw-read filename patterns;
metadata, indexes, quantification tables, and evidence are never candidates.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Iterator

_RUN_ACCESSION = re.compile(r"(?:SRR|ERR|DRR)\d+$")
_RAW_SUFFIXES = (
    ".fastq.gz",
    ".fastq",
    ".fq.gz",
    ".fq",
    ".sra",
    ".sra.part",
    ".fastq.gz.part",
    ".safely_removed",
)


def _is_raw_read_name(name: str) -> bool:
    """Return whether a filename is a recognized raw or failed raw transfer."""

    return name.endswith(_RAW_SUFFIXES) or ".fastq.gz.invalid" in name


def _sample_input_roots(work_dir: Path, run_accession: str) -> Iterator[Path]:
    """Yield supported per-sample input directories without following links."""

    for base in (work_dir / "getfastq", work_dir / "fastq" / "getfastq", work_dir / "fastq"):
        candidate = base / run_accession
        if candidate.is_dir() or candidate.is_symlink():
            yield candidate


def _flat_input_files(work_dir: Path, run_accession: str) -> Iterator[Path]:
    """Yield supported flat downloader outputs for one accession."""

    for base in (work_dir / "getfastq", work_dir / "fastq" / "getfastq", work_dir / "fastq"):
        if not base.is_dir():
            continue
        for path in base.iterdir():
            if path.is_file() and path.name.startswith(run_accession) and _is_raw_read_name(path.name):
                yield path


def reclaim_sample_raw_inputs(
    work_dir: str | Path,
    run_accession: str,
    *,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Remove known raw inputs for one successfully quantified accession.

    The caller is responsible for proving current-contract quantification
    before invoking this function.  A dry run reports the exact files and
    bytes that would be removed.  Symlinked sample directories are reported
    as protected rather than followed, preventing accidental deletion outside
    the selected work tree.
    """

    if not _RUN_ACCESSION.fullmatch(run_accession):
        raise ValueError(f"invalid run accession: {run_accession!r}")

    root = Path(work_dir).expanduser().resolve()
    result: dict[str, Any] = {
        "work_dir": str(root),
        "run_accession": run_accession,
        "dry_run": dry_run,
        "files_deleted": 0,
        "bytes_freed": 0,
        "paths": [],
        "protected_paths": [],
        "errors": [],
    }
    candidates = list(_sample_input_roots(root, run_accession))
    files = list(_flat_input_files(root, run_accession))
    for candidate in candidates:
        if candidate.is_symlink():
            result["protected_paths"].append(str(candidate))
            result["errors"].append(f"refused to follow symlinked raw directory: {candidate}")
            continue
        for path in candidate.rglob("*"):
            if path.is_file() and _is_raw_read_name(path.name):
                files.append(path)

    seen: set[Path] = set()
    for path in files:
        if path.is_symlink():
            result["protected_paths"].append(str(path))
            result["errors"].append(f"refused to remove symlinked raw file: {path}")
            continue
        resolved = path.resolve(strict=False)
        if resolved in seen:
            continue
        seen.add(resolved)
        try:
            size = path.stat().st_size
            result["bytes_freed"] += size
            result["files_deleted"] += 1
            result["paths"].append(str(path))
            if not dry_run:
                path.unlink()
        except OSError as exc:
            result["errors"].append(f"{path}: {exc}")

    if not dry_run:
        for candidate in candidates:
            if candidate.is_dir() and not candidate.is_symlink():
                try:
                    # Only remove directories that are empty after recognized
                    # raw files have been deleted; retain any evidence marker.
                    candidate.rmdir()
                except OSError:
                    pass

    return result

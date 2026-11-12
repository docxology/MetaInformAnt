"""Safe per-sample reclamation of transient RNA-seq inputs.

The streaming campaign keeps quantification outputs and provenance under the
selected data root.  Once a reusable quantification has succeeded, the
corresponding FASTQ/SRA inputs are reproducible transient inputs and may be
removed to keep the external volume usable.  This module deliberately touches
only exact per-accession locations and known raw-read filename patterns;
metadata, indexes, quantification tables, validation witnesses and evidence are
never candidates.

Three per-accession scopes exist:

* :func:`reclaim_sample_raw_inputs` -- FASTQ/SRA files and failed partials in
  the species ``getfastq`` tree (optionally with the accession's orphaned
  ``fasterq-dump`` staging workspace under ``getfastq/.fasterq_tmp/<run>``);
* :func:`reclaim_fasterq_staging` -- only that staging workspace;
* :func:`reclaim_sample_sra_cache` -- the accession's entries in the shared
  campaign VDB cache (``<data-root>/.sra-cache/sra/<run>.sra[.cache]`` and
  ``<data-root>/.sra-cache/prefetch/<run>/``).

Every function refuses a malformed accession, never follows a symbolic link,
reports dry runs with the exact files and bytes, and counts only files that
were actually removed.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any, Iterator

_RUN_ACCESSION = re.compile(r"(?:SRR|ERR|DRR)\d+$")
_RAW_BASE_SUFFIXES = (
    ".fastq.gz",
    ".fastq",
    ".fq.gz",
    ".fq",
    ".sra",
    ".sra.cache",
    ".sra.vdbcache",
    ".vdbcache",
)
# Transfer/extraction states that wrap a raw name: ENA partials (``.part``),
# prefetch temporaries (``.tmp``, ``.prf``, ``.lock``) and rejected payloads
# (``.invalid`` or ``.invalid-<ns>`` written when a payload fails validation).
_TRANSIENT_TRAILERS = re.compile(r"(?:\.part|\.tmp|\.prf|\.lock|\.invalid(?:-\d+)?)+$")
FASTERQ_STAGING_DIRNAME = ".fasterq_tmp"
SRA_CACHE_DIRNAME = ".sra-cache"


def _is_raw_read_name(name: str) -> bool:
    """Return whether a filename is a recognized raw read, VDB file or failed raw transfer."""

    if name.endswith(".safely_removed") or ".fastq.gz.invalid" in name:
        return True
    return _TRANSIENT_TRAILERS.sub("", name).endswith(_RAW_BASE_SUFFIXES)


def _check_accession(run_accession: str) -> None:
    if not _RUN_ACCESSION.fullmatch(run_accession):
        raise ValueError(f"invalid run accession: {run_accession!r}")


def _has_accession_boundary(name: str, run_accession: str) -> bool:
    """Return whether *name* addresses *run_accession*, not a longer accession.

    A valid accession can be a strict prefix of another (``SRR123456`` and
    ``SRR1234567``); a flat downloader output always carries ``_`` or ``.``
    directly after the accession stem.
    """

    return name.startswith(f"{run_accession}_") or name.startswith(f"{run_accession}.")


def _getfastq_bases(work_dir: Path) -> tuple[Path, ...]:
    return (work_dir / "getfastq", work_dir / "fastq" / "getfastq", work_dir / "fastq")


def _sample_input_roots(work_dir: Path, run_accession: str) -> Iterator[Path]:
    """Yield supported per-sample input directories without following links."""

    for base in _getfastq_bases(work_dir):
        candidate = base / run_accession
        if candidate.is_dir() or candidate.is_symlink():
            yield candidate


def _flat_input_files(work_dir: Path, run_accession: str) -> Iterator[Path]:
    """Yield supported flat downloader outputs for one accession."""

    for base in _getfastq_bases(work_dir):
        if not base.is_dir():
            continue
        for path in base.iterdir():
            if path.is_file() and _has_accession_boundary(path.name, run_accession) and _is_raw_read_name(path.name):
                yield path


def _new_result(root: Path, run_accession: str, dry_run: bool, **extra: Any) -> dict[str, Any]:
    return {
        "work_dir": str(root),
        "run_accession": run_accession,
        "dry_run": dry_run,
        "files_deleted": 0,
        "bytes_freed": 0,
        "allocated_bytes_freed": 0,
        "paths": [],
        "protected_paths": [],
        "errors": [],
        **extra,
    }


def _remove_file(path: Path, result: dict[str, Any], dry_run: bool) -> None:
    """Remove one regular file (or plan it) and account apparent and allocated bytes."""

    try:
        info = path.lstat()
        if not dry_run:
            path.unlink()
    except OSError as exc:
        result["errors"].append(f"{path}: {exc}")
        return
    result["files_deleted"] += 1
    result["bytes_freed"] += info.st_size
    result["allocated_bytes_freed"] += getattr(info, "st_blocks", 0) * 512
    result["paths"].append(str(path))


def _tree_has_symlink(root: Path) -> Path | None:
    """Return the first symbolic link inside *root* (``root`` itself included), else ``None``."""

    if root.is_symlink():
        return root
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        for name in (*dirnames, *filenames):
            candidate = Path(dirpath) / name
            if candidate.is_symlink():
                return candidate
    return None


def _remove_tree(root: Path, result: dict[str, Any], dry_run: bool) -> None:
    """Remove every regular file below *root*, then the emptied directories (deepest first).

    The caller has already proven the tree contains no symbolic link. Special
    files (sockets, FIFOs) are left in place and reported, so their
    directories simply stay.
    """

    directories: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root, topdown=True, followlinks=False):
        directories.append(Path(dirpath))
        dirnames.sort()
        for name in sorted(filenames):
            path = Path(dirpath) / name
            if path.is_file():
                _remove_file(path, result, dry_run)
            else:
                result["errors"].append(f"refused to remove a non-regular file: {path}")
    if dry_run:
        return
    for directory in reversed(directories):
        try:
            directory.rmdir()
        except OSError:
            pass


def reclaim_fasterq_staging(
    work_dir: str | Path,
    run_accession: str,
    *,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Remove the accession's orphaned ``fasterq-dump`` workspace(s) under ``getfastq/.fasterq_tmp``.

    The streaming orchestrator stages each extraction in
    ``<getfastq>/.fasterq_tmp/<run>/{output,temp}`` and removes that exact
    workspace before it retries the accession.  A killed producer leaves it
    behind; once the accession is quantified (or permanently excluded) the
    workspace is never read again.  Sibling validation witnesses
    (``.fasterq_tmp/raw_validation``) are not accession directories and are
    never touched.  A workspace containing any symbolic link is protected.
    """

    _check_accession(run_accession)
    root = Path(work_dir).expanduser().resolve()
    result = _new_result(root, run_accession, dry_run, scope="fasterq_staging")
    for base in _getfastq_bases(root):
        staging = base / FASTERQ_STAGING_DIRNAME
        if staging.is_symlink() or not staging.is_dir():
            continue
        workspace = staging / run_accession
        if not (workspace.is_dir() or workspace.is_symlink()):
            continue
        link = _tree_has_symlink(workspace)
        if link is not None:
            result["protected_paths"].append(str(link))
            result["errors"].append(f"refused to clear a staging workspace containing a symlink: {link}")
            continue
        _remove_tree(workspace, result, dry_run)
    return result


def reclaim_sample_raw_inputs(
    work_dir: str | Path,
    run_accession: str,
    *,
    dry_run: bool = False,
    include_staging: bool = False,
) -> dict[str, Any]:
    """Remove known raw inputs for one successfully quantified accession.

    The caller is responsible for proving a reusable quantification (or a
    permanent exclusion) before invoking this function.  A dry run reports the
    exact files and bytes that would be removed.  Symlinked sample directories
    are reported as protected rather than followed, preventing accidental
    deletion outside the selected work tree.  With ``include_staging`` the
    accession's orphaned ``fasterq-dump`` workspace is removed as well (see
    :func:`reclaim_fasterq_staging`).
    """

    _check_accession(run_accession)
    root = Path(work_dir).expanduser().resolve()
    result = _new_result(root, run_accession, dry_run)
    candidates = list(_sample_input_roots(root, run_accession))
    files = list(_flat_input_files(root, run_accession))
    for candidate in candidates:
        if candidate.is_symlink():
            result["protected_paths"].append(str(candidate))
            result["errors"].append(f"refused to follow symlinked raw directory: {candidate}")
            continue
        for dirpath, _dirnames, filenames in os.walk(candidate, followlinks=False):
            for name in filenames:
                if _is_raw_read_name(name):
                    files.append(Path(dirpath) / name)

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
        _remove_file(path, result, dry_run)

    if not dry_run:
        for candidate in candidates:
            if candidate.is_dir() and not candidate.is_symlink():
                # Only directories that are empty after recognized raw files
                # were deleted go away; evidence markers keep theirs.
                for dirpath, _dirnames, _filenames in os.walk(candidate, topdown=False, followlinks=False):
                    try:
                        Path(dirpath).rmdir()
                    except OSError:
                        pass

    if include_staging:
        staging = reclaim_fasterq_staging(root, run_accession, dry_run=dry_run)
        for key in ("files_deleted", "bytes_freed", "allocated_bytes_freed"):
            result[key] += staging[key]
        for key in ("paths", "protected_paths", "errors"):
            result[key].extend(staging[key])
    return result


def sra_cache_root(data_root: str | Path) -> Path:
    """Return the campaign VDB cache directory used by the streaming orchestrator."""

    return Path(data_root) / SRA_CACHE_DIRNAME


def reclaim_sample_sra_cache(
    cache_root: str | Path,
    run_accession: str,
    *,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Remove one accession's entries from the shared campaign VDB cache.

    Candidates are exactly ``<cache>/sra/<run>.<raw suffix>`` (``.sra``,
    ``.sra.cache``, ``.vdbcache`` and their ``.tmp``/``.prf``/``.lock``/
    ``.invalid`` states) and accession-named files below
    ``<cache>/prefetch/<run>/``.  Validation markers
    (``sra/.metainformant_sra_validation``), the ``refseq`` reference cache
    and anything not named for the accession stay.  The orchestrator treats a
    cached archive as authoritative only until quantification succeeds, so the
    caller must prove a reusable quantification (or permanent exclusion) of
    *every* task that names this accession before calling.
    """

    _check_accession(run_accession)
    root = Path(cache_root).expanduser().resolve()
    result = _new_result(root, run_accession, dry_run, scope="sra_cache")
    result["cache_root"] = result.pop("work_dir")

    flat_dir = root / "sra"
    if flat_dir.is_dir() and not flat_dir.is_symlink():
        for path in sorted(flat_dir.iterdir()):
            if not path.name.startswith(f"{run_accession}.") or not _is_raw_read_name(path.name):
                continue
            if path.is_symlink():
                result["protected_paths"].append(str(path))
                result["errors"].append(f"refused to remove symlinked cache file: {path}")
            elif path.is_file():
                _remove_file(path, result, dry_run)

    prefetch_root = root / "prefetch"
    tree = prefetch_root / run_accession
    if prefetch_root.is_dir() and not prefetch_root.is_symlink() and (tree.is_dir() or tree.is_symlink()):
        link = _tree_has_symlink(tree)
        if link is not None:
            result["protected_paths"].append(str(link))
            result["errors"].append(f"refused to clear a prefetch tree containing a symlink: {link}")
        else:
            directories: list[Path] = []
            for dirpath, dirnames, filenames in os.walk(tree, topdown=True, followlinks=False):
                directories.append(Path(dirpath))
                dirnames.sort()
                for name in sorted(filenames):
                    path = Path(dirpath) / name
                    if _has_accession_boundary(name, run_accession) and _is_raw_read_name(name) and path.is_file():
                        _remove_file(path, result, dry_run)
            if not dry_run:
                for directory in reversed(directories):
                    try:
                        directory.rmdir()
                    except OSError:
                        pass
    return result

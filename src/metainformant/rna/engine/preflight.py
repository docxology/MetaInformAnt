"""Campaign environment preflight for the streaming RNA-seq producer.

A producer that starts into a broken environment poisons the cohort instead of
processing it. Observed 2026-09-03: a producer without external-volume write
access ran for twelve minutes, failed 7,291 sample tasks with
``[Errno 1] Operation not permitted``, and a second producer failed its
quantification batches because the bare ``amalgkit`` CLI was not on ``PATH``.
Both failure classes are detectable in milliseconds before any task is
scheduled, so the preflight is mandatory at producer start.

The preflight performs real, minimal probes against the selected environment:

- **Data-root writability.** Write, rename, inspect, and delete a probe file
  under the data root. This exercises the same syscall classes the producer
  depends on (create, rename, stat, unlink) rather than only directory
  existence.
- **Amalgkit CLI resolution.** Resolve the bare ``amalgkit`` command exactly
  the way the producer's quantification subprocesses do, via
  :func:`shutil.which` against the effective ``PATH``.
- **Producer tools.** ``kallisto``, ``fasterq-dump``, ``prefetch``,
  ``vdb-validate`` and ``pigz`` must resolve on the same ``PATH``; a missing
  tool would fail every sample after its download.
- **Kallisto compatibility.** The installed ``kallisto version`` must equal
  the ``kallisto_version`` recorded in one existing ``run_info.json`` per
  species, and a species' kallisto index header must carry the recorded
  ``index_version`` (otherwise new outputs would not be comparable with the
  quantified ones).
- **Free-space floors.** The data root must have at least
  ``AMALGKIT_MIN_EXTERNAL_FREE_GB`` (default 8 GiB) free and the system volume
  ``AMALGKIT_MIN_SYSTEM_FREE_GB`` (default 2 GiB: scratch lives on the data
  root, so the internal disk only needs headroom for logs and the
  environment).
- **Reuse dry-run (optional).** :func:`reuse_dry_run` classifies a small
  sample of recorded quantified outputs per species and refuses to start when
  more than 1% would be quarantined, before any state is changed.

Every check runs even when an earlier check fails, so one run reports the
complete environment state. The module is also runnable directly::

    python -m metainformant.rna.engine.preflight --data-root /path/to/amalgkit
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import re
import shutil
import struct
import subprocess
import sys
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger
from metainformant.rna.amalgkit.sra_environment import resolve_data_root

logger = get_logger(__name__)

PROBE_FILE_NAME = ".metainformant_preflight_probe"
PRODUCER_TOOLS = ("amalgkit", "kallisto", "fasterq-dump", "prefetch", "vdb-validate", "pigz")
"""Executables the streaming producer runs; each must resolve on the producer's ``PATH``."""
DEFAULT_MIN_DATA_ROOT_FREE_GB = 8.0
DEFAULT_MIN_SYSTEM_FREE_GB = 2.0
DEFAULT_SYSTEM_ROOT = "/System/Volumes/Data"
REUSE_DRY_RUN_SEED = 20261112
_KALLISTO_VERSION = re.compile(r"kallisto,?\s+version\s+([0-9][0-9A-Za-z.+-]*)")


class PreflightError(RuntimeError):
    """Raised when one or more campaign preflight checks fail."""


def probe_data_root_writable(data_root: Path) -> Path:
    """Verify ``data_root`` accepts create/rename/stat/unlink operations.

    Returns the probed data root (resolved) on success and raises
    :class:`PreflightError` on the first failed operation. The probe file is
    always removed again; a successful probe leaves no trace.
    """

    root = Path(data_root)
    probe = root / PROBE_FILE_NAME
    moved = root / f"{PROBE_FILE_NAME}.moved"
    try:
        try:
            probe.write_bytes(b"metainformant preflight write probe\n")
        except OSError as exc:
            raise PreflightError(
                f"Data root is not writable (create failed): {root}\n"
                f"  {exc}\n"
                "  Grant the producer's process Full Disk Access / Removable "
                "Volumes permission, or start it from a shell that already has "
                "volume access, then rerun."
            ) from exc
        try:
            os.replace(probe, moved)
            if not moved.is_file():
                raise PreflightError(f"Rename probe vanished after replace: {moved}")
        except OSError as exc:
            raise PreflightError(
                f"Data root does not support rename operations: {root}\n"
                f"  {exc}\n"
                "  The producer quarantines failed transfers with atomic "
                "renames; a read-only or foreign-provenance volume cannot "
                "host a campaign."
            ) from exc
        return root.resolve()
    finally:
        probe.unlink(missing_ok=True)
        moved.unlink(missing_ok=True)


def resolve_amalgkit_cli(search_path: str | None = None) -> str:
    """Resolve the bare ``amalgkit`` command the quantification stage requires.

    ``search_path`` overrides ``PATH`` for the lookup (a colon-separated
    string), mirroring the subprocess environment under test.
    """

    resolved = shutil.which("amalgkit", path=search_path)
    if resolved is None:
        searched = search_path or os.environ.get("PATH", "")
        raise PreflightError(
            "The 'amalgkit' CLI was not found on the effective PATH.\n"
            "  Quantification shells out to bare 'amalgkit'; without it every "
            "quant batch fails with ENOENT.\n"
            "  Add the repository virtualenv bin directory "
            "(.venv/bin) to PATH, or pass --search-path.\n"
            f"  Searched: {searched}"
        )
    return resolved


def resolve_producer_tools(tools: Sequence[str] = PRODUCER_TOOLS, search_path: str | None = None) -> dict[str, str]:
    """Resolve every producer executable on the effective ``PATH``; one error names all missing tools."""

    resolved: dict[str, str] = {}
    missing: list[str] = []
    for tool in tools:
        if tool == "amalgkit":
            try:
                resolved[tool] = resolve_amalgkit_cli(search_path)
            except PreflightError:
                missing.append(tool)
            continue
        found = shutil.which(tool, path=search_path)
        if found is None:
            missing.append(tool)
        else:
            resolved[tool] = found
    if missing:
        searched = search_path if search_path is not None else os.environ.get("PATH", "")
        raise PreflightError(
            f"Required producer tools are not on the effective PATH: {', '.join(missing)}\n"
            "  Every sample would fail after its download. Install them (Homebrew: kallisto, sratoolkit, "
            "pigz; amalgkit comes from the repository virtualenv) or fix PATH.\n"
            f"  Searched: {searched}"
        )
    return resolved


def _gib(value: int) -> float:
    return value / (1024**3)


def _free_gib(path: Path) -> float:
    usage = shutil.disk_usage(path)
    return _gib(usage.free)


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return default
    try:
        value = float(raw)
    except ValueError:
        raise PreflightError(f"{name} must be a number of GiB: {raw!r}") from None
    if not math.isfinite(value) or value < 0:
        raise PreflightError(f"{name} must be a nonnegative number of GiB: {raw!r}")
    return value


def check_free_space(
    data_root: Path,
    *,
    min_data_root_free_gb: float | None = None,
    min_system_free_gb: float | None = None,
    system_root: Path | str | None = None,
) -> dict[str, str]:
    """Check the data-root and system-volume free-space floors (GiB) in Python; return the measured values.

    Floors default to ``AMALGKIT_MIN_EXTERNAL_FREE_GB`` (8) and
    ``AMALGKIT_MIN_SYSTEM_FREE_GB`` (2). The system root defaults to
    ``AMALGKIT_SYSTEM_ROOT`` or ``/System/Volumes/Data`` and is skipped when it
    does not exist (Linux, CI).
    """

    data_floor = (
        min_data_root_free_gb
        if min_data_root_free_gb is not None
        else _env_float("AMALGKIT_MIN_EXTERNAL_FREE_GB", DEFAULT_MIN_DATA_ROOT_FREE_GB)
    )
    system_floor = (
        min_system_free_gb
        if min_system_free_gb is not None
        else _env_float("AMALGKIT_MIN_SYSTEM_FREE_GB", DEFAULT_MIN_SYSTEM_FREE_GB)
    )
    system = Path(system_root or os.environ.get("AMALGKIT_SYSTEM_ROOT", DEFAULT_SYSTEM_ROOT))
    facts: dict[str, str] = {}
    failures: list[str] = []
    try:
        data_free = _free_gib(data_root)
    except OSError as exc:
        raise PreflightError(f"Cannot measure free space on the data root {data_root}: {exc}") from exc
    facts["data_root_free_gib"] = f"{data_free:.1f}"
    if data_free < data_floor:
        failures.append(
            f"Data root {data_root} has {data_free:.1f} GiB free, below the {data_floor:g} GiB floor "
            "(AMALGKIT_MIN_EXTERNAL_FREE_GB)."
        )
    if system.is_dir():
        system_free = _free_gib(system)
        facts["system_free_gib"] = f"{system_free:.1f}"
        if system_free < system_floor:
            failures.append(
                f"System volume {system} has {system_free:.1f} GiB free, below the {system_floor:g} GiB floor "
                "(AMALGKIT_MIN_SYSTEM_FREE_GB); keep TMPDIR and caches on the data root and free space."
            )
    if failures:
        raise PreflightError("\n  ".join(failures))
    return facts


def installed_kallisto_version(kallisto: str) -> str:
    """The version printed by ``kallisto version`` (for example ``0.52.0``)."""

    try:
        completed = subprocess.run([kallisto, "version"], capture_output=True, text=True, timeout=60, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise PreflightError(f"Cannot run '{kallisto} version': {exc}") from exc
    match = _KALLISTO_VERSION.search(f"{completed.stdout}\n{completed.stderr}")
    if completed.returncode != 0 or match is None:
        raise PreflightError(
            f"'{kallisto} version' did not report a version (exit {completed.returncode}): "
            f"{(completed.stdout + completed.stderr).strip()[:200]}"
        )
    return match.group(1)


def kallisto_index_version(index_path: Path) -> int:
    """The index format version stored in the first eight bytes (little-endian) of a kallisto index."""

    with open(index_path, "rb") as handle:
        header = handle.read(8)
    if len(header) != 8:
        raise PreflightError(f"kallisto index is shorter than its header: {index_path}")
    return int(struct.unpack("<Q", header)[0])


def find_run_info(quant_dir: Path, *, max_entries: int = 64) -> Path | None:
    """One ``<run>/<run>_run_info.json`` (or ``run_info.json``) under ``quant_dir``, scanning few entries."""

    try:
        entries = os.scandir(quant_dir)
    except OSError:
        return None
    with entries:
        for checked, entry in enumerate(entries):
            if checked >= max_entries:
                break
            if entry.name.startswith(".") or entry.name.startswith("amalgkit_quant_stage_"):
                continue
            if not entry.is_dir(follow_symlinks=False):
                continue
            for name in (f"{entry.name}_run_info.json", "run_info.json"):
                candidate = Path(entry.path) / name
                if candidate.is_file():
                    return candidate
    return None


def check_kallisto_compatibility(
    kallisto: str,
    quant_dirs: Mapping[str, Path],
    *,
    index_paths: Mapping[str, Path] | None = None,
) -> dict[str, Any]:
    """Require the installed kallisto to match what produced each species' existing outputs.

    For every species with a ``run_info.json`` under its quant directory, the
    installed version must equal the recorded ``kallisto_version`` and, when
    the species' index is given, the index header must equal the recorded
    ``index_version``. Species without outputs are only listed. Returns
    ``{installed, checked, without_outputs}``.
    """

    installed = installed_kallisto_version(kallisto)
    checked: dict[str, dict[str, Any]] = {}
    without: list[str] = []
    failures: list[str] = []
    for species, quant_dir in sorted(quant_dirs.items()):
        run_info_path = find_run_info(Path(quant_dir))
        if run_info_path is None:
            without.append(species)
            continue
        try:
            run_info = json.loads(run_info_path.read_text())
        except (OSError, ValueError) as exc:
            failures.append(f"{species}: unreadable {run_info_path}: {exc}")
            continue
        recorded_version = str(run_info.get("kallisto_version") or "")
        recorded_index = run_info.get("index_version")
        record: dict[str, Any] = {
            "run_info": str(run_info_path),
            "kallisto_version": recorded_version,
            "index_version": recorded_index,
        }
        if recorded_version != installed:
            failures.append(
                f"{species}: installed kallisto {installed} differs from {recorded_version or 'an unrecorded version'} "
                f"recorded in {run_info_path}"
            )
        index_path = (index_paths or {}).get(species)
        if index_path is not None and recorded_index is not None:
            try:
                header = kallisto_index_version(Path(index_path))
            except (OSError, PreflightError) as exc:
                failures.append(f"{species}: cannot read the kallisto index header of {index_path}: {exc}")
            else:
                record["index_header_version"] = header
                if header != int(recorded_index):
                    failures.append(
                        f"{species}: index {index_path} has format version {header}, outputs recorded "
                        f"index_version {recorded_index}"
                    )
        checked[species] = record
    if failures:
        raise PreflightError(
            "Kallisto does not match the existing quantifications (set AMALGKIT_ALLOW_KALLISTO_DRIFT=1 only for "
            "a deliberate toolchain change):\n  " + "\n  ".join(failures)
        )
    return {"installed": installed, "checked": checked, "without_outputs": without}


@dataclass(frozen=True)
class ReuseProbe:
    """One species for :func:`reuse_dry_run`: recorded-quantified runs and how the producer would judge them."""

    species: str
    quant_dir: Path
    recorded_quantified: Sequence[str]
    classify: Callable[[Path, str], Mapping[str, Any]]


def reuse_dry_run(
    probes: Sequence[ReuseProbe],
    *,
    fraction: float = 0.01,
    max_quarantine_fraction: float = 0.01,
    seed: int = REUSE_DRY_RUN_SEED,
    reusable_statuses: frozenset[str] | None = None,
) -> dict[str, Any]:
    """Classify a deterministic sample of each species' recorded-quantified outputs; refuse a mass quarantine.

    Samples ``ceil(fraction * n)`` (at least one) of the runs each probe
    records as quantified and classifies them exactly as resume would. If the
    share that would be quarantined (not reusable, or no output directory)
    exceeds ``max_quarantine_fraction`` in any species, raises
    :class:`PreflightError` naming examples: a provenance or path change must
    never silently send thousands of finished samples back to the network.
    Reads only.
    """

    if reusable_statuses is None:
        from metainformant.rna.engine.provenance import REUSABLE_QUANT_STATUSES

        reusable_statuses = REUSABLE_QUANT_STATUSES
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    rng = random.Random(seed)
    report: dict[str, Any] = {}
    failures: list[str] = []
    for probe in probes:
        runs = sorted(set(probe.recorded_quantified))
        if not runs:
            report[probe.species] = {"recorded": 0, "sampled": 0, "would_quarantine": 0}
            continue
        sample = rng.sample(runs, min(len(runs), max(1, math.ceil(fraction * len(runs)))))
        bad: list[str] = []
        for run in sorted(sample):
            path = Path(probe.quant_dir) / run
            if not path.is_dir():
                bad.append(f"{run}: no output directory")
                continue
            try:
                verdict = probe.classify(path, run)
            except (OSError, ValueError) as exc:
                bad.append(f"{run}: classification failed: {exc}")
                continue
            if verdict.get("status") not in reusable_statuses:
                bad.append(f"{run}: {verdict.get('status')} ({verdict.get('reason')})")
        share = len(bad) / len(sample)
        report[probe.species] = {
            "recorded": len(runs),
            "sampled": len(sample),
            "would_quarantine": len(bad),
            "examples": bad[:5],
        }
        if share > max_quarantine_fraction:
            failures.append(
                f"{probe.species}: {len(bad)} of {len(sample)} sampled quantified outputs would be quarantined "
                f"({share:.1%} > {max_quarantine_fraction:.1%}); e.g. {'; '.join(bad[:3])}"
            )
    if failures:
        raise PreflightError("Reuse dry-run refused the resume:\n  " + "\n  ".join(failures))
    return report


def run_campaign_preflight(
    data_root: Path | str | None = None,
    *,
    search_path: str | None = None,
    tools: Sequence[str] = PRODUCER_TOOLS,
    kallisto_quant_dirs: Mapping[str, Path] | None = None,
    kallisto_index_paths: Mapping[str, Path] | None = None,
    check_space: bool = True,
    min_data_root_free_gb: float | None = None,
    min_system_free_gb: float | None = None,
    system_root: Path | str | None = None,
) -> dict[str, str]:
    """Run all campaign preflight checks and return resolved environment facts.

    All checks run; failures are collected into one :class:`PreflightError`.
    ``data_root`` defaults to the configured ``AMALGKIT_DATA_ROOT`` contract.
    The kallisto compatibility check runs when ``kallisto_quant_dirs`` is
    given (species -> quant directory); ``AMALGKIT_ALLOW_KALLISTO_DRIFT=1``
    turns its failure into a warning.
    """

    facts: dict[str, str] = {}
    failures: list[str] = []

    root = Path(data_root) if data_root is not None else resolve_data_root()
    try:
        facts["data_root"] = str(probe_data_root_writable(root))
    except PreflightError as exc:
        failures.append(str(exc))

    try:
        facts["amalgkit_cli"] = resolve_amalgkit_cli(search_path)
    except PreflightError as exc:
        failures.append(str(exc))

    resolved: dict[str, str] = {}
    others = tuple(tool for tool in tools if tool != "amalgkit")
    if others:
        try:
            resolved = resolve_producer_tools(others, search_path)
        except PreflightError as exc:
            failures.append(str(exc))
        facts.update({f"tool_{name.replace('-', '_')}": path for name, path in resolved.items()})

    if check_space and root.is_dir():
        try:
            facts.update(
                check_free_space(
                    root,
                    min_data_root_free_gb=min_data_root_free_gb,
                    min_system_free_gb=min_system_free_gb,
                    system_root=system_root,
                )
            )
        except PreflightError as exc:
            failures.append(str(exc))

    if kallisto_quant_dirs is not None and "kallisto" in resolved:
        try:
            kallisto_facts = check_kallisto_compatibility(
                resolved["kallisto"], kallisto_quant_dirs, index_paths=kallisto_index_paths
            )
            facts["kallisto_version"] = kallisto_facts["installed"]
            facts["kallisto_species_checked"] = str(len(kallisto_facts["checked"]))
        except PreflightError as exc:
            if os.environ.get("AMALGKIT_ALLOW_KALLISTO_DRIFT", "").strip().lower() in {"1", "yes", "true", "on"}:
                logger.warning("Kallisto drift allowed by AMALGKIT_ALLOW_KALLISTO_DRIFT: %s", exc)
                facts["kallisto_drift_allowed"] = "yes"
            else:
                failures.append(str(exc))

    if failures:
        raise PreflightError("Campaign preflight failed:\n" + "\n".join(f"- {failure}" for failure in failures))
    return facts


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for ``python -m metainformant.rna.engine.preflight``."""

    parser = argparse.ArgumentParser(
        description=(
            "Verify the campaign environment (data-root write access, producer tools, kallisto "
            "compatibility and free space) before starting a producer."
        )
    )
    parser.add_argument(
        "--data-root",
        default=None,
        help="Amalgkit data root to probe (default: AMALGKIT_DATA_ROOT contract)",
    )
    parser.add_argument(
        "--search-path",
        default=None,
        help="Colon-separated PATH override for the tool lookups",
    )
    parser.add_argument(
        "--kallisto-species",
        action="append",
        default=[],
        help="Also compare the installed kallisto with one run_info.json of this species (repeatable)",
    )
    args = parser.parse_args(argv)

    root = Path(args.data_root) if args.data_root is not None else resolve_data_root()
    quant_dirs = {species: root / species / "work" / "quant" for species in args.kallisto_species} or None
    try:
        facts = run_campaign_preflight(root, search_path=args.search_path, kallisto_quant_dirs=quant_dirs)
    except PreflightError as exc:
        print(f"PREFLIGHT FAILED\n{exc}", file=sys.stderr)
        return 1
    for key in sorted(facts):
        print(f"{key}: {facts[key]}")
    print("preflight: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())

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

Every check runs even when an earlier check fails, so one run reports the
complete environment state. The module is also runnable directly::

    python -m metainformant.rna.engine.preflight --data-root /path/to/amalgkit
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

from metainformant.core.utils.logging import get_logger
from metainformant.rna.amalgkit.sra_environment import resolve_data_root

logger = get_logger(__name__)

PROBE_FILE_NAME = ".metainformant_preflight_probe"


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


def run_campaign_preflight(
    data_root: Path | str | None = None,
    *,
    search_path: str | None = None,
) -> dict[str, str]:
    """Run all campaign preflight checks and return resolved environment facts.

    All checks run; failures are collected into one :class:`PreflightError`.
    ``data_root`` defaults to the configured ``AMALGKIT_DATA_ROOT`` contract.
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

    if failures:
        raise PreflightError(
            "Campaign preflight failed:\n" + "\n".join(f"- {failure}" for failure in failures)
        )
    return facts


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for ``python -m metainformant.rna.engine.preflight``."""

    parser = argparse.ArgumentParser(
        description=(
            "Verify the campaign environment (data-root write access and the "
            "amalgkit CLI) before starting a producer."
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
        help="Colon-separated PATH override for the amalgkit lookup",
    )
    args = parser.parse_args(argv)

    try:
        facts = run_campaign_preflight(args.data_root, search_path=args.search_path)
    except PreflightError as exc:
        print(f"PREFLIGHT FAILED\n{exc}", file=sys.stderr)
        return 1
    for key in sorted(facts):
        print(f"{key}: {facts[key]}")
    print("preflight: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())

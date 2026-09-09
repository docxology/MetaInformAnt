"""Cohort funnel accounting for the Hymenoptera RNA-seq campaign.

Builds a data-derived funnel from the campaign progress DB (read-only) and the
per-species amalgkit configuration directory:

- ``configured``: species declared by ``config/amalgkit/amalgkit_*.yaml``
  (excluding the ``template``/``test``/``cross_species`` files), consistent
  with the project's 27-species set.
- ``with_progress``: species with at least one row in the progress DB.
- ``quantified_runs`` / ``failed_runs`` / ``pending_runs``: samples-table row
  counts per run state.
- ``excluded_runs``: rows in the ``sample_exclusions`` table.
- ``active_runs``: runs currently ``downloading`` or ``quantifying``.

Every stage is derived from the data; nothing is hard-coded. The builder fails
closed (raising :class:`CohortFunnelError`, a ``ValueError``) when the DB or
config directory is missing, the DB has no ``samples`` table, or a config file
does not declare a usable ``species_list``. Failure classes come from
:func:`metainformant.rna.engine.progress_db.classify_sample_error`, so this
funnel, the figure generator, and the ``db_failure_classes`` telemetry always
share one durable class set.
"""

from __future__ import annotations

import sqlite3
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import yaml

from metainformant.rna.engine.progress_db import classify_sample_error

__all__ = [
    "CohortFunnelError",
    "FunnelReport",
    "STAGE_NAMES",
    "build_cohort_funnel",
    "iter_failure_trend_rows",
    "render_funnel_lines",
]

#: Config basenames that are not per-species cohort definitions.
EXCLUDED_CONFIG_BASENAMES = frozenset(
    {
        "amalgkit_template.yaml",
        "amalgkit_test.yaml",
        "amalgkit_cross_species.yaml",
    }
)

#: Funnel stages in reporting order.
STAGE_NAMES: tuple[str, ...] = (
    "configured",
    "with_progress",
    "quantified_runs",
    "failed_runs",
    "excluded_runs",
    "pending_runs",
    "active_runs",
)


class CohortFunnelError(ValueError):
    """Cohort accounting cannot proceed from the given inputs (fail closed)."""


@dataclass(frozen=True)
class FunnelReport:
    """Immutable snapshot of the cohort funnel.

    Stage fields appear in :data:`STAGE_NAMES` order; ``reason_codes`` maps
    each observed durable failure class (from ``classify_sample_error``) to
    the number of failed runs in that class. ``db_path`` and ``config_dir``
    record the inputs the snapshot was derived from.
    """

    configured: int
    with_progress: int
    quantified_runs: int
    failed_runs: int
    excluded_runs: int
    pending_runs: int
    active_runs: int
    reason_codes: dict[str, int]
    db_path: Path
    config_dir: Path

    def stage_counts(self) -> dict[str, int]:
        """Return the funnel stages as an ordered ``{stage: count}`` mapping."""

        return {name: getattr(self, name) for name in STAGE_NAMES}

    def to_tsv(self, path: Path) -> None:
        """Write a two-column TSV (``stage_or_class``, ``count``).

        Stages come first in :data:`STAGE_NAMES` order, then a blank line,
        then the observed failure classes sorted by class name — the whole
        file is byte-deterministic for a given snapshot.
        """

        lines = ["stage_or_class\tcount"]
        for name, count in self.stage_counts().items():
            lines.append(f"{name}\t{count}")
        if self.reason_codes:
            lines.append("")
            for class_name in sorted(self.reason_codes):
                lines.append(f"{class_name}\t{self.reason_codes[class_name]}")
        Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def render_funnel_lines(report: FunnelReport) -> list[str]:
    """Render ``cohort_funnel_<stage>: <count>`` summary lines, in stage order."""

    return [f"cohort_funnel_{name}: {report.stage_counts()[name]}" for name in STAGE_NAMES]


def _count_configured_species(config_dir: Path) -> int:
    """Count distinct species declared by the per-species config YAMLs.

    Fails closed if the directory is missing or any per-species config lacks
    a usable ``species_list``.
    """

    config_dir = Path(config_dir)
    if not config_dir.is_dir():
        raise CohortFunnelError(f"amalgkit config directory not found: {config_dir}")
    species: set[str] = set()
    for yaml_path in sorted(config_dir.glob("amalgkit_*.yaml")):
        if yaml_path.name in EXCLUDED_CONFIG_BASENAMES:
            continue
        data = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
        entries = data.get("species_list") if isinstance(data, dict) else None
        if not isinstance(entries, list) or not entries:
            raise CohortFunnelError(f"config {yaml_path.name} declares no species_list")
        for entry in entries:
            if not isinstance(entry, str) or not entry.strip():
                raise CohortFunnelError(f"config {yaml_path.name} has a non-string species_list entry")
            species.add(entry.strip())
    return len(species)


def _failed_rows(db_path: Path) -> list[str | None]:
    """Read the ``error`` text of every failed run, read-only."""

    connection = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        return [row[0] for row in connection.execute("SELECT error FROM samples WHERE state='failed'").fetchall()]
    finally:
        connection.close()


def build_cohort_funnel(
    db_path: Path,
    config_dir: Path,
    *,
    max_gb: float | None = None,
) -> FunnelReport:
    """Build the cohort funnel from the progress DB and config directory.

    Args:
        db_path: Path to the campaign progress SQLite DB; opened read-only.
        config_dir: Directory holding the per-species ``amalgkit_*.yaml``
            configuration files.
        max_gb: Optional upper bound on the DB file size in GiB; a larger DB
            fails closed instead of being read (saturated-volume guard).

    Returns:
        A :class:`FunnelReport` with every stage derived from the data.

    Raises:
        CohortFunnelError: If the DB or config directory is missing, the DB
            has no ``samples`` table, or the DB exceeds ``max_gb``.
    """

    db_path = Path(db_path)
    if not db_path.is_file():
        raise CohortFunnelError(f"progress DB not found: {db_path}")
    if max_gb is not None:
        size_gib = db_path.stat().st_size / 1024**3
        if size_gib > max_gb:
            raise CohortFunnelError(f"progress DB exceeds the {max_gb} GiB limit: {db_path} is {size_gib:.3f} GiB")

    configured = _count_configured_species(Path(config_dir))

    connection = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if "samples" not in tables:
            raise CohortFunnelError(f"progress DB has no samples table: {db_path}")

        def _scalar(query: str) -> int:
            return int(connection.execute(query).fetchone()[0])

        with_progress = _scalar("SELECT COUNT(DISTINCT species) FROM samples")
        quantified_runs = _scalar("SELECT COUNT(*) FROM samples WHERE state='quantified'")
        failed_runs = _scalar("SELECT COUNT(*) FROM samples WHERE state='failed'")
        excluded_runs = _scalar("SELECT COUNT(*) FROM sample_exclusions") if "sample_exclusions" in tables else 0
        pending_runs = _scalar("SELECT COUNT(*) FROM samples WHERE state='pending'")
        active_runs = _scalar("SELECT COUNT(*) FROM samples WHERE state IN ('downloading', 'quantifying')")
    finally:
        connection.close()

    failures = _failed_rows(db_path)
    reason_codes: dict[str, int] = dict(sorted(Counter(classify_sample_error(error) for error in failures).items()))

    return FunnelReport(
        configured=configured,
        with_progress=with_progress,
        quantified_runs=quantified_runs,
        failed_runs=failed_runs,
        excluded_runs=excluded_runs,
        pending_runs=pending_runs,
        active_runs=active_runs,
        reason_codes=reason_codes,
        db_path=db_path,
        config_dir=Path(config_dir),
    )


def iter_failure_trend_rows(db_path: Path) -> Iterable[tuple[str, str | None]]:
    """Yield ``(date(updated_at), error)`` for every failed run, read-only."""

    connection = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        yield from connection.execute(
            "SELECT date(updated_at), error FROM samples" " WHERE state='failed' ORDER BY date(updated_at)"
        ).fetchall()
    finally:
        connection.close()

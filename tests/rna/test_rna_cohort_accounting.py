"""Tests for cohort funnel accounting against deterministic fixtures."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from metainformant.rna.analysis.cohort_accounting import (
    CohortFunnelError,
    FunnelReport,
    STAGE_NAMES,
    build_cohort_funnel,
    iter_failure_trend_rows,
    render_funnel_lines,
)

SAMPLES_SCHEMA = (
    "CREATE TABLE samples ("
    " species     TEXT    NOT NULL,"
    " srr_id      TEXT    NOT NULL,"
    " state       TEXT    NOT NULL DEFAULT 'pending',"
    " error       TEXT,"
    " updated_at  TEXT    NOT NULL DEFAULT (datetime('now')),"
    " PRIMARY KEY (species, srr_id)"
    ")"
)

SAMPLE_ROWS = [
    ("apis_mellifera", "SRR_A1", "quantified", None),
    ("apis_mellifera", "SRR_A2", "quantified", None),
    ("apis_mellifera", "SRR_A3", "pending", None),
    ("apis_mellifera", "SRR_A4", "failed", "kallisto exited with exit code 1"),
    ("apis_mellifera", "SRR_A5", "downloading", None),
    ("apis_mellifera", "SRR_A6", "quantifying", None),
    ("apis_mellifera", "SRR_A7", "failed", "Download Failed (all sources: ENA FTP/HTTP, NCBI)"),
    ("bombus_terrestris", "SRR_B1", "quantified", None),
    ("bombus_terrestris", "SRR_B2", "failed", None),
]

EXCLUSIONS_SCHEMA = (
    "CREATE TABLE sample_exclusions ("
    " species        TEXT    NOT NULL,"
    " srr_id         TEXT    NOT NULL,"
    " reason_code    TEXT    NOT NULL,"
    " reason_detail  TEXT,"
    " recorded_by    TEXT,"
    " recorded_at    TEXT    NOT NULL DEFAULT (datetime('now')),"
    " PRIMARY KEY (species, srr_id)"
    ")"
)

EXCLUSION_ROWS = [
    ("apis_mellifera", "SRR_A9", "permanent_drop"),
    ("bombus_terrestris", "SRR_B9", "re_download"),
]


def _write_config_dir(root: Path) -> Path:
    """Write a two-species amalgkit config dir plus non-cohort files."""
    config_dir = root / "config" / "amalgkit"
    config_dir.mkdir(parents=True, exist_ok=True)
    for slug in ("apis_mellifera", "bombus_terrestris"):
        (config_dir / f"amalgkit_{slug}.yaml").write_text(
            f"work_dir: output/amalgkit/{slug}/work\nspecies_list:\n  - {slug}\n",
            encoding="utf-8",
        )
    for name in ("amalgkit_template.yaml", "amalgkit_test.yaml", "amalgkit_cross_species.yaml"):
        (config_dir / name).write_text("species_list:\n  - Not_A_Cohort_Species\n", encoding="utf-8")
    return config_dir


def _write_progress_db(root: Path, *, with_exclusions: bool = True) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    db_path = root / "pipeline_progress.db"
    connection = sqlite3.connect(db_path)
    try:
        connection.execute(SAMPLES_SCHEMA)
        connection.executemany(
            "INSERT INTO samples (species, srr_id, state, error) VALUES (?, ?, ?, ?)", SAMPLE_ROWS
        )
        if with_exclusions:
            connection.execute(EXCLUSIONS_SCHEMA)
            connection.executemany(
                "INSERT INTO sample_exclusions (species, srr_id, reason_code) VALUES (?, ?, ?)",
                EXCLUSION_ROWS,
            )
        # Fixed dates so trend bucketing is deterministic.
        connection.execute("UPDATE samples SET updated_at='2026-09-01 10:00:00' WHERE srr_id='SRR_A4'")
        connection.execute("UPDATE samples SET updated_at='2026-09-02 11:00:00' WHERE srr_id='SRR_A7'")
        connection.execute("UPDATE samples SET updated_at='2026-09-02 12:00:00' WHERE srr_id='SRR_B2'")
        connection.commit()
    finally:
        connection.close()
    return db_path


def _expected_stage_counts() -> dict[str, int]:
    return {
        "configured": 2,
        "with_progress": 2,
        "quantified_runs": 3,
        "failed_runs": 3,
        "excluded_runs": 2,
        "pending_runs": 1,
        "active_runs": 2,
    }


def test_build_cohort_funnel_stage_counts(tmp_path: Path) -> None:
    """Every stage is derived from the DB rows and config files."""
    db_path = _write_progress_db(tmp_path)
    config_dir = _write_config_dir(tmp_path)

    report = build_cohort_funnel(db_path, config_dir)

    assert report.stage_counts() == _expected_stage_counts()
    # Durable classes via classify_sample_error, deterministic (sorted) order.
    assert report.reason_codes == {
        "transfer_all_sources_failed": 1,
        "unclassified": 1,
        "unrecorded": 1,
    }
    assert isinstance(report, FunnelReport)
    assert report.db_path == db_path
    assert report.config_dir == config_dir


def test_missing_exclusions_table_counts_zero(tmp_path: Path) -> None:
    """A DB without the sample_exclusions table reports excluded_runs=0."""
    db_path = _write_progress_db(tmp_path, with_exclusions=False)
    config_dir = _write_config_dir(tmp_path)

    report = build_cohort_funnel(db_path, config_dir)

    assert report.excluded_runs == 0


def test_render_funnel_lines_order_and_format(tmp_path: Path) -> None:
    """Lines are cohort_funnel_<stage>: <count> in STAGE_NAMES order."""
    db_path = _write_progress_db(tmp_path)
    config_dir = _write_config_dir(tmp_path)
    report = build_cohort_funnel(db_path, config_dir)

    lines = render_funnel_lines(report)

    assert lines == [
        f"cohort_funnel_{name}: {_expected_stage_counts()[name]}" for name in STAGE_NAMES
    ]
    assert lines[0] == "cohort_funnel_configured: 2"
    assert lines[-1] == "cohort_funnel_active_runs: 2"


def test_to_tsv_golden_order(tmp_path: Path) -> None:
    """TSV is stages in order, then a blank-line separated class section."""
    db_path = _write_progress_db(tmp_path)
    config_dir = _write_config_dir(tmp_path)
    report = build_cohort_funnel(db_path, config_dir)
    tsv_path = tmp_path / "cohort_funnel.tsv"

    report.to_tsv(tsv_path)

    assert tsv_path.read_text(encoding="utf-8") == (
        "stage_or_class\tcount\n"
        "configured\t2\n"
        "with_progress\t2\n"
        "quantified_runs\t3\n"
        "failed_runs\t3\n"
        "excluded_runs\t2\n"
        "pending_runs\t1\n"
        "active_runs\t2\n"
        "\n"
        "transfer_all_sources_failed\t1\n"
        "unclassified\t1\n"
        "unrecorded\t1\n"
    )


def test_trend_rows_bucket_by_day(tmp_path: Path) -> None:
    """Trend rows are (date(updated_at), error) pairs in ascending day order."""
    db_path = _write_progress_db(tmp_path)

    assert list(iter_failure_trend_rows(db_path)) == [
        ("2026-09-01", "kallisto exited with exit code 1"),
        ("2026-09-02", "Download Failed (all sources: ENA FTP/HTTP, NCBI)"),
        ("2026-09-02", None),
    ]


def test_missing_db_fails_closed(tmp_path: Path) -> None:
    """A missing progress DB refuses to build a funnel."""
    config_dir = _write_config_dir(tmp_path)

    with pytest.raises(CohortFunnelError, match="progress DB not found"):
        build_cohort_funnel(tmp_path / "absent.db", config_dir)


def test_db_without_samples_table_fails_closed(tmp_path: Path) -> None:
    """A DB without a samples table refuses to build a funnel."""
    db_path = tmp_path / "pipeline_progress.db"
    connection = sqlite3.connect(db_path)
    try:
        connection.execute("CREATE TABLE unrelated (x INTEGER)")
        connection.commit()
    finally:
        connection.close()
    config_dir = _write_config_dir(tmp_path)

    with pytest.raises(CohortFunnelError, match="no samples table"):
        build_cohort_funnel(db_path, config_dir)


def test_missing_config_dir_fails_closed(tmp_path: Path) -> None:
    """A missing config directory refuses to build a funnel."""
    db_path = _write_progress_db(tmp_path)

    with pytest.raises(CohortFunnelError, match="config directory not found"):
        build_cohort_funnel(db_path, tmp_path / "absent_config")


def test_max_gb_guard_fails_closed(tmp_path: Path) -> None:
    """A DB larger than max_gb refuses to be read (saturated-volume guard)."""
    db_path = _write_progress_db(tmp_path)
    config_dir = _write_config_dir(tmp_path)

    with pytest.raises(CohortFunnelError, match="exceeds the"):
        build_cohort_funnel(db_path, config_dir, max_gb=1e-9)


def test_config_without_species_list_fails_closed(tmp_path: Path) -> None:
    """A per-species config without species_list fails closed."""
    db_path = _write_progress_db(tmp_path)
    config_dir = _write_config_dir(tmp_path)
    (config_dir / "amalgkit_formica_exsecta.yaml").write_text("threads: 8\n", encoding="utf-8")

    with pytest.raises(CohortFunnelError, match="declares no species_list"):
        build_cohort_funnel(db_path, config_dir)

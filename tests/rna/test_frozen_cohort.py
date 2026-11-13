"""Frozen cohort definition, installation, verification and the producer task table."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from metainformant.rna.engine.cohort import (
    COHORT_DEFINITION_FILENAME,
    CohortError,
    build_cohort_definition,
    cohort_task_table_bytes,
    cohort_task_table_name,
    default_max_sample_gb,
    install_cohort,
    load_frozen_cohort,
    main,
    quant_batch_run,
    select_cohort_species,
    verify_cohort_configs,
    verify_species_reference,
    write_cohort_definition,
    write_cohort_task_table,
)
from tests.rna.frozen_cohort_support import (
    COHORT_ID,
    CONFIG_BYTES,
    CONFIG_NAME,
    INDEX_BYTES,
    INDEX_NAME,
    SPECIES,
    DataRoot,
    build_frozen_source,
    sha,
)


@pytest.fixture
def installed(tmp_path: Path) -> Path:
    build_frozen_source(tmp_path / "frozen")
    install_cohort(tmp_path / "frozen", tmp_path / "cohort", cohort_id=COHORT_ID)
    return tmp_path / "cohort"


def _make_writable(path: Path) -> None:
    os.chmod(path, 0o644)


def test_definition_is_deterministic_and_describes_the_frozen_snapshot(tmp_path: Path) -> None:
    source = build_frozen_source(tmp_path / "frozen")
    first = build_cohort_definition(source.root, cohort_id=COHORT_ID)
    second = build_cohort_definition(source.root, cohort_id=COHORT_ID)
    digest_a = write_cohort_definition(first, tmp_path / "a.json")
    digest_b = write_cohort_definition(second, tmp_path / "b.json")
    assert digest_a == digest_b
    assert (tmp_path / "a.json").read_bytes() == (tmp_path / "b.json").read_bytes()
    assert first["task_count"] == 4 and first["exclusion_count"] == 1
    assert first["exclusions"] == [{"species": SPECIES, "accession": "SRR104", "reason_code": "permanent_drop"}]
    assert {record["path"] for record in first["files"]} == {
        "snapshot.json",
        "manifest.jsonl",
        "inventory.json",
        f"inputs/{SPECIES}/work/metadata/metadata_selected.tsv",
        f"inputs/{SPECIES}/work/index/{INDEX_NAME}",
    }
    (entry,) = first["species"]
    assert entry["config_sha256"] == sha(CONFIG_BYTES) and entry["index_sha256"] == sha(INDEX_BYTES)


def test_install_copies_read_only_and_loads_verified(tmp_path: Path) -> None:
    source = build_frozen_source(tmp_path / "frozen")
    before = {path: path.read_bytes() for path in source.root.rglob("*") if path.is_file()}
    summary = install_cohort(source.root, tmp_path / "cohort", cohort_id=COHORT_ID)
    assert summary["copied"] == 6 and summary["identical"] == 0
    assert {path: path.read_bytes() for path in source.root.rglob("*") if path.is_file()} == before
    assert not (source.root / COHORT_DEFINITION_FILENAME).exists()
    cohort = load_frozen_cohort(tmp_path / "cohort")
    assert cohort.cohort_id == COHORT_ID and cohort.definition_sha256 == summary["definition_sha256"]
    assert [task.accession for task in cohort.tasks_for(SPECIES)] == ["SRR103", "SRR101", "SRR102", "SRR105"]
    assert [task.batch_index for task in cohort.tasks_for(SPECIES)] == [1, 2, 4, 5]
    assert cohort.exclusions_for(SPECIES) == {"SRR104": "permanent_drop"}
    assert cohort.task(SPECIES, "SRR104") is None and cohort.is_excluded(SPECIES, "SRR104")
    expectation = cohort.expectation(SPECIES, "SRR101")
    assert expectation.reference_index_sha256 == sha(INDEX_BYTES) and expectation.reference_index_name == INDEX_NAME
    for path in (tmp_path / "cohort").rglob("*"):
        if path.is_file():
            assert path.stat().st_mode & 0o222 == 0, path
    # Re-installing is a no-op.
    again = install_cohort(source.root, tmp_path / "cohort", cohort_id=COHORT_ID)
    assert again["copied"] == 0 and again["identical"] == 6


def test_install_removes_the_temporaries_of_a_killed_install(tmp_path: Path) -> None:
    source = build_frozen_source(tmp_path / "frozen")
    install_cohort(source.root, tmp_path / "cohort", cohort_id=COHORT_ID)
    index_dir = tmp_path / "cohort" / "inputs" / SPECIES / "work" / "index"
    current = index_dir / f".{INDEX_NAME}.a1b2c3d4.cohort-install.tmp"
    foreign = index_dir / ".operator-notes"  # not an install temp: left alone
    for path in (current, foreign):
        path.write_bytes(INDEX_BYTES[:100])

    again = install_cohort(source.root, tmp_path / "cohort", cohort_id=COHORT_ID)

    relative = f"inputs/{SPECIES}/work/index"
    assert again["stale_temporaries_removed"] == [f"{relative}/{current.name}"]
    assert not current.exists() and foreign.exists()
    assert again["copied"] == 0 and again["identical"] == 6


def test_install_never_replaces_a_different_existing_file(tmp_path: Path) -> None:
    source = build_frozen_source(tmp_path / "frozen")
    target = tmp_path / "cohort" / "snapshot.json"
    target.parent.mkdir(parents=True)
    target.write_text("{}\n")
    with pytest.raises(FileExistsError):
        install_cohort(source.root, tmp_path / "cohort", cohort_id=COHORT_ID)
    assert target.read_text() == "{}\n"


@pytest.mark.parametrize(
    "relative",
    [f"inputs/{SPECIES}/work/metadata/metadata_selected.tsv", f"inputs/{SPECIES}/work/index/{INDEX_NAME}"],
)
def test_load_refuses_any_changed_input_byte(installed: Path, relative: str) -> None:
    path = installed / relative
    _make_writable(path)
    data = bytearray(path.read_bytes())
    data[-3] ^= 0x01
    path.write_bytes(bytes(data))
    with pytest.raises(CohortError, match="differs from its definition"):
        load_frozen_cohort(installed)


def test_load_refuses_a_changed_definition_and_a_wrong_pin(installed: Path) -> None:
    definition = installed / COHORT_DEFINITION_FILENAME
    digest = sha(definition.read_bytes())
    assert load_frozen_cohort(installed, expected_definition_sha256=digest).definition_sha256 == digest
    with pytest.raises(CohortError, match="pinned"):
        load_frozen_cohort(installed, expected_definition_sha256="0" * 64)
    _make_writable(definition)
    payload = json.loads(definition.read_text())
    payload["exclusions"] = []
    payload["exclusion_count"] = 0
    definition.write_text(json.dumps(payload))
    with pytest.raises(CohortError):
        load_frozen_cohort(installed)


def test_load_refuses_a_missing_or_extra_table_row(tmp_path: Path) -> None:
    source = build_frozen_source(tmp_path / "frozen")
    definition = build_cohort_definition(source.root, cohort_id=COHORT_ID)
    write_cohort_definition(definition, tmp_path / "COHORT.json")
    # A separate (pinned) definition lets the unmodified frozen source be the cohort directory.
    assert len(load_frozen_cohort(source.root, definition_path=tmp_path / "COHORT.json").tasks) == 4
    source.table.write_bytes(source.table.read_bytes() + b"Test species\tSRR999\t1\tSINGLE\tyes\tno\r\n")
    with pytest.raises(CohortError):
        load_frozen_cohort(source.root, definition_path=tmp_path / "COHORT.json")


def test_build_refuses_a_batch_index_that_is_not_the_table_position(tmp_path: Path) -> None:
    source = build_frozen_source(tmp_path / "frozen", batch_override={"SRR102": 3})
    with pytest.raises(CohortError, match="is not the table position 4"):
        build_cohort_definition(source.root, cohort_id=COHORT_ID)


def test_task_table_keeps_frozen_bytes_and_amalgkit_batch_positions(installed: Path, tmp_path: Path) -> None:
    cohort = load_frozen_cohort(installed)
    table = tmp_path / "work" / "metadata" / cohort_task_table_name(COHORT_ID)
    positions = write_cohort_task_table(cohort, SPECIES, table)
    data = table.read_bytes()
    frozen_lines = cohort.metadata_path(SPECIES).read_bytes().splitlines(keepends=True)
    assert data == b"".join(line for line in frozen_lines if b"SRR104" not in line)
    assert data.count(b"\r\n") == 5
    assert positions == {"SRR103": 1, "SRR101": 2, "SRR102": 3, "SRR105": 4}
    for accession, batch in positions.items():
        assert quant_batch_run(table, batch) == accession
    # Identical bytes are not rewritten.
    mtime = table.stat().st_mtime_ns
    assert write_cohort_task_table(cohort, SPECIES, table) == positions
    assert table.stat().st_mtime_ns == mtime
    # Extra exclusions drop rows and renumber the batches.
    extra = write_cohort_task_table(cohort, SPECIES, table, excluded={"SRR103"})
    assert extra == {"SRR101": 1, "SRR102": 2, "SRR105": 3}
    assert quant_batch_run(table, 1) == "SRR101"
    assert cohort_task_table_bytes(cohort, SPECIES)[0] == data


def test_configs_and_species_selection_are_checked(installed: Path, tmp_path: Path) -> None:
    cohort = load_frozen_cohort(installed)
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    (config_dir / CONFIG_NAME).write_bytes(CONFIG_BYTES)
    assert verify_cohort_configs(cohort, config_dir) == {SPECIES: sha(CONFIG_BYTES)}
    (config_dir / CONFIG_NAME).write_bytes(CONFIG_BYTES + b"# edited\n")
    with pytest.raises(CohortError, match="differs from the frozen cohort"):
        verify_cohort_configs(cohort, config_dir)
    assert select_cohort_species(cohort, None) == (SPECIES,)
    with pytest.raises(CohortError, match="not in cohort"):
        select_cohort_species(cohort, ["camponotus_floridanus"])


def test_species_reference_must_hash_to_the_cohort_index(installed: Path, tmp_path: Path) -> None:
    cohort = load_frozen_cohort(installed)
    good = DataRoot(tmp_path / "good")
    check = verify_species_reference(cohort, good.base, SPECIES, index_dir=good.index_dir)
    assert check.ok, check
    bad = DataRoot(tmp_path / "bad", index_bytes=b"another index build\n")
    check = verify_species_reference(cohort, bad.base, SPECIES)
    assert not check.ok and "differ from the cohort index" in check.reason
    missing = verify_species_reference(cohort, tmp_path / "absent", SPECIES)
    assert not missing.ok and "unreadable" in missing.reason


def test_max_sample_budget_follows_the_launcher() -> None:
    assert default_max_sample_gb({}) == 1000.0
    assert default_max_sample_gb({"PIPELINE_MAX_GB": "50"}) == 50.0
    assert default_max_sample_gb({"AMALGKIT_MAX_SAMPLE_GB": "800", "PIPELINE_MAX_GB": "50"}) == 800.0
    with pytest.raises(ValueError):
        default_max_sample_gb({"AMALGKIT_MAX_SAMPLE_GB": "-1"})


def test_cli_define_install_verify(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    source = build_frozen_source(tmp_path / "frozen")
    assert (
        main(["define", "--source", str(source.root), "--cohort-id", COHORT_ID, "--output", str(tmp_path / "d.json")])
        == 0
    )
    defined = json.loads(capsys.readouterr().out)
    assert defined["sha256"] == sha((tmp_path / "d.json").read_bytes())
    assert (
        main(["install", "--source", str(source.root), "--cohort-dir", str(tmp_path / "c"), "--cohort-id", COHORT_ID])
        == 0
    )
    installed = json.loads(capsys.readouterr().out)
    assert installed["definition_sha256"] == defined["sha256"]
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    (config_dir / CONFIG_NAME).write_bytes(CONFIG_BYTES)
    args = ["verify", "--cohort", str(tmp_path / "c"), "--cohort-definition-sha256", defined["sha256"]]
    assert main([*args, "--config-dir", str(config_dir)]) == 0
    assert json.loads(capsys.readouterr().out)["tasks"] == 4
    assert main(["verify", "--cohort", str(tmp_path / "c"), "--cohort-definition-sha256", "f" * 64]) == 2

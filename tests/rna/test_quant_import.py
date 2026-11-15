"""Importing verified quantifications from other data roots (plan, execute, idempotence, conflicts)."""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import pytest

from metainformant.rna.engine.campaign_lock import (
    LOCK_FILENAME as CAMPAIGN_LOCK_FILENAME,
    CampaignLockBusy,
    acquire_writer,
)
from metainformant.rna.engine.cohort import (
    install_cohort,
    load_frozen_cohort,
)
from metainformant.rna.engine.durable_quant import lock_quantification
from metainformant.rna.engine.progress_db import ProgressDB
from metainformant.rna.engine.provenance import (
    QUANT_PROVENANCE_FILENAME,
    QUANT_STATUS_CURRENT,
    REFERENCE_BINDING_FILENAME,
    classify_cohort_quantification,
    read_reference_binding,
)
from metainformant.rna.engine.quant_import import (
    ACTION_EXISTS_DIFFERENT,
    ACTION_EXISTS_IDENTICAL,
    ACTION_IMPORT,
    ACTION_IMPORTED,
    ACTION_SKIPPED,
    DIFFERENCE_DESTINATION_QUANTIFICATION,
    DIFFERENCE_SOURCE_QUANTIFICATION,
    IDENTITY_BASIS_DEEP,
    IDENTITY_BASIS_SHALLOW,
    STAGING_MARKER,
    STAGING_SCHEMA,
    ImportSource,
    ReceiptStore,
    _quantification_against_sidecar,
    execute_imports,
    executed_imports,
    main,
    parse_receipt_store,
    parse_source,
    plan_imports,
    write_plan_manifest,
)
from metainformant.rna.engine.quant_storage import DirectoryStore
from tests.rna.frozen_cohort_support import (
    COHORT_ID,
    CONFIG_NAME,
    INDEX_BYTES,
    SPECIES,
    DataRoot,
    build_frozen_source,
    sha,
)


def _tree(root: Path) -> dict[str, bytes]:
    """Every file's bytes, except the writer lock whose owner record changes with each acquisition."""

    return {
        str(path.relative_to(root)): path.read_bytes()
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.name != CAMPAIGN_LOCK_FILENAME
    }


class World:
    """Cohort + destination root + October-like source root + a receipt store."""

    def __init__(self, base: Path) -> None:
        build_frozen_source(base / "frozen")
        install_cohort(base / "frozen", base / "cohort", cohort_id=COHORT_ID)
        self.cohort = load_frozen_cohort(base / "cohort")
        self.dest = DataRoot(base / "dest")
        self.october = DataRoot(base / "october")
        self.october.write_sample("SRR101")  # October-only, bound by its own root's manifest
        (self.october.sample("SRR101") / "notes.txt").write_text("stray file\n")
        self.october.write_sample("SRR102")  # also in the destination, same bytes
        shutil.copytree(self.october.sample("SRR102"), self.dest.sample("SRR102"))
        self.october.write_sample("SRR103")  # destination holds a different quantification
        self.dest.write_sample("SRR103", counts_offset=5)
        self.october.write_sample("SRR104")  # permanent_drop in the cohort
        self.october.write_sample("SRR900")  # not a cohort task
        other_config = base / "october" / "other.yaml"
        other_config.write_bytes(b"species_list: [Test_species]\nthreads: 2\n")
        self.october.write_sample("SRR105", config=other_config)  # different configuration bytes
        (self.october.work / "quant" / ".import-SRR101-0000").mkdir()  # foreign dot entries are ignored
        self.source = ImportSource("october", self.october.base, "work")
        self.base = base

    def plan(self, **kwargs):
        return plan_imports(self.cohort, self.dest.base, [self.source], workers=1, **kwargs)


@pytest.fixture
def world(tmp_path: Path) -> World:
    return World(tmp_path)


def _by_accession(plan) -> dict[str, object]:
    return {decision.accession: decision for decision in plan.decisions}


def test_parse_source_and_receipt_store_specs() -> None:
    assert parse_source("october=/data/oct") == ImportSource("october", Path("/data/oct"), "work")
    assert parse_source("cloud=/data/c/completed_quant:flat").layout == "flat"
    assert parse_receipt_store("oct=/data/oct/locked-quant:local-20261019") == ReceiptStore(
        "oct", Path("/data/oct/locked-quant"), "local-20261019"
    )
    with pytest.raises(ValueError):
        parse_source("/no/label")
    with pytest.raises(ValueError):
        parse_receipt_store("oct=/no/cohort")


def test_plan_decides_every_sample_and_writes_nothing(world: World, tmp_path: Path) -> None:
    before_dest, before_source = _tree(world.dest.base), _tree(world.october.base)
    plan = world.plan()
    decisions = _by_accession(plan)
    assert decisions["SRR101"].action == ACTION_IMPORT
    assert decisions["SRR101"].evidence["kind"] == "reference_manifest"
    assert decisions["SRR101"].evidence["index_bytes_verified"] is True
    assert {entry["name"] for entry in decisions["SRR101"].files} == {
        QUANT_PROVENANCE_FILENAME,
        "SRR101_abundance.tsv",
        "SRR101_run_info.json",
        "SRR101_abundance.h5",
    }
    assert decisions["SRR102"].action == ACTION_EXISTS_IDENTICAL
    assert decisions["SRR103"].action == ACTION_EXISTS_DIFFERENT
    assert decisions["SRR103"].details["bound"] == "both"
    assert decisions["SRR103"].details["index_bytes_proven"] == "both"
    assert decisions["SRR103"].details["quantification"] == "different"
    assert decisions["SRR103"].reason == "destination holds a different quantification"
    assert (
        decisions["SRR103"].details["source_abundance_sha256"]
        != decisions["SRR103"].details["destination_abundance_sha256"]
    )
    assert (decisions["SRR104"].action, decisions["SRR104"].reason) == (ACTION_SKIPPED, "permanent_drop")
    assert (decisions["SRR900"].action, decisions["SRR900"].reason) == (ACTION_SKIPPED, "not_in_cohort")
    assert (decisions["SRR105"].action, decisions["SRR105"].reason) == (ACTION_SKIPPED, "config_mismatch")
    assert plan.ignored == {"dot_entries": 1}
    manifest = write_plan_manifest(plan, tmp_path / "plans" / "plan.jsonl", [world.source])
    records = [json.loads(line) for line in manifest.read_text().splitlines()]
    assert records[0]["kind"] == "header" and records[-1]["kind"] == "summary"
    assert len([r for r in records if r["kind"] == "decision"]) == 6
    assert not [r for r in records if r["kind"] == "decision" and r["source_relpath"].startswith("/")]
    assert _tree(world.dest.base) == before_dest and _tree(world.october.base) == before_source
    with pytest.raises(ValueError, match="inside a source root"):
        write_plan_manifest(plan, world.october.base / "plan.jsonl", [world.source])


def test_conflict_with_the_same_abundance_is_reported_as_a_runtime_difference(world: World) -> None:
    sidecar = world.dest.sample("SRR102") / QUANT_PROVENANCE_FILENAME
    payload = json.loads(sidecar.read_text())
    payload["command"] = [*payload["command"], "--threads", "2"]
    sidecar.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    before = _tree(world.dest.base)
    decision = _by_accession(world.plan())["SRR102"]
    assert decision.action == ACTION_EXISTS_DIFFERENT
    assert decision.details["quantification"] == "identical"
    assert decision.details["comparison"]["difference"] == "sidecar bytes differ"
    assert decision.reason == "destination holds the same quantification with different runtime files"
    assert _tree(world.dest.base) == before


def test_reference_mismatch_when_index_bytes_differ(world: World, tmp_path: Path) -> None:
    rebuilt = DataRoot(tmp_path / "rebuilt", index_bytes=b"a different kallisto index build\n")
    rebuilt.write_sample("SRR101")
    plan = plan_imports(world.cohort, world.dest.base, [ImportSource("rebuilt", rebuilt.base)], workers=1)
    (decision,) = plan.decisions
    assert (decision.action, decision.reason) == (ACTION_SKIPPED, "reference_mismatch")
    # The destination's own reference is checked too: a destination on other index bytes takes nothing.
    other_dest = DataRoot(tmp_path / "other_dest", index_bytes=b"another index\n")
    plan = plan_imports(world.cohort, other_dest.base, [world.source], workers=1)
    assert _by_accession(plan)["SRR101"].reason == "reference_mismatch"
    assert not plan.references[SPECIES].ok


def test_execute_imports_october_only_samples_portably_and_idempotently(world: World) -> None:
    stale = world.dest.work / "quant" / ".import-SRR101-deadbeef"
    stale.mkdir()
    (stale / STAGING_MARKER).write_text(json.dumps({"schema": STAGING_SCHEMA}))
    unmarked = world.dest.work / "quant" / ".import-unknown"
    unmarked.mkdir()
    conflict_before = _tree(world.dest.sample("SRR103"))
    source_before = _tree(world.october.base)

    plan, manifest, summary = execute_imports(world.cohort, world.dest.base, [world.source], workers=1)
    assert summary["imported"] == 1 and summary["failed"] == 0
    assert summary["stale_staging_removed"] == [stale.name]
    assert not stale.exists() and unmarked.exists()
    assert _by_accession(plan)["SRR101"].action == ACTION_IMPORTED
    imported = world.dest.sample("SRR101")
    source = world.october.sample("SRR101")
    for name in (QUANT_PROVENANCE_FILENAME, "SRR101_abundance.tsv", "SRR101_run_info.json", "SRR101_abundance.h5"):
        assert (imported / name).read_bytes() == (source / name).read_bytes()
    assert not (imported / "notes.txt").exists()
    assert sorted(path.name for path in imported.iterdir()) == sorted(
        [QUANT_PROVENANCE_FILENAME, REFERENCE_BINDING_FILENAME, "SRR101_abundance.tsv", "SRR101_run_info.json"]
        + ["SRR101_abundance.h5"]
    )
    binding = read_reference_binding(imported)
    assert binding["cohort_id"] == COHORT_ID and binding["evidence"]["kind"] == "reference_manifest"
    assert binding["evidence"]["imported_from"] == "october"
    assert not [value for value in binding["evidence"].values() if str(value).startswith("/")]
    assert binding["evidence"]["source_relpath"] == f"{SPECIES}/work/quant/SRR101"
    assert not [p for p in (world.dest.work / "quant").iterdir() if p.name.startswith(".import-SRR101")]
    # The conflicting destination sample and the source are untouched.
    assert _tree(world.dest.sample("SRR103")) == conflict_before
    assert _tree(world.october.base) == source_before
    # The manifest is under the destination and fully parseable.
    assert manifest.parent == world.dest.base / "results" / "imports"
    records = [json.loads(line) for line in manifest.read_text().splitlines()]
    assert [r["kind"] for r in (records[0], records[-1])] == ["header", "summary"]
    assert {r["accession"]: r["action"] for r in records if r["kind"] == "decision"}["SRR101"] == ACTION_IMPORTED
    # The publication was announced (fsynced) before the rename, ahead of its decision record.
    kinds = [(r["kind"], r.get("accession")) for r in records]
    assert kinds.index(("publishing", "SRR101")) < kinds.index(("decision", "SRR101"))

    # Portable: the import classifies as current from its own witness even after the source is gone.
    shutil.rmtree(world.october.base)
    result = classify_cohort_quantification(imported, world.cohort.expectation(SPECIES, "SRR101"), verify_content=True)
    assert result["status"] == QUANT_STATUS_CURRENT and result["evidence"]["kind"] == "binding_witness"


def test_rerun_imports_nothing_new(world: World) -> None:
    execute_imports(world.cohort, world.dest.base, [world.source], workers=1)
    snapshot = _tree(world.dest.base / SPECIES)
    plan, _, summary = execute_imports(world.cohort, world.dest.base, [world.source], workers=1)
    assert summary["imported"] == 0
    decisions = _by_accession(plan)
    assert decisions["SRR101"].action == ACTION_EXISTS_IDENTICAL
    assert decisions["SRR103"].action == ACTION_EXISTS_DIFFERENT
    assert _tree(world.dest.base / SPECIES) == snapshot


def test_receipt_rescues_missing_manifest_evidence_and_detects_tampering(world: World, tmp_path: Path) -> None:
    store_root = tmp_path / "locked-quant"
    lock_quantification(
        DirectoryStore(store_root),
        "local-test",
        SPECIES,
        "SRR101",
        world.october.sample("SRR101"),
        expected_reference_index_sha256=world.cohort.entry(SPECIES).index_sha256,
    )
    store = ReceiptStore("october-locked", store_root, "local-test")
    world.october.manifest.unlink()  # the source root's reference is gone; only the receipt binds the index
    plan = world.plan(receipt_stores=[store])
    decision = _by_accession(plan)["SRR101"]
    assert decision.action == ACTION_IMPORT
    assert decision.evidence["kind"] == "receipt" and decision.evidence["receipt_cohort"] == "local-test"
    assert _by_accession(world.plan())["SRR101"].reason == "no_reference_evidence"

    abundance = world.october.sample("SRR101") / "SRR101_run_info.json"
    abundance.write_text(abundance.read_text().replace("900", "901"))
    assert _by_accession(world.plan(receipt_stores=[store]))["SRR101"].reason == "receipt_mismatch"


def test_execute_refuses_while_another_writer_holds_the_root(world: World) -> None:
    with acquire_writer(world.dest.base, "finalization"):
        before = _tree(world.dest.base)
        with pytest.raises(CampaignLockBusy):
            execute_imports(world.cohort, world.dest.base, [world.source], workers=1)
        assert _tree(world.dest.base) == before


def test_execute_can_adopt_into_the_destination_db(world: World) -> None:
    db_path = world.dest.base / "pipeline_progress.db"
    with ProgressDB(db_path) as db:
        db.init_species(SPECIES, ["SRR101", "SRR102"])
    _, _, summary = execute_imports(world.cohort, world.dest.base, [world.source], adopt_db=db_path, workers=1)
    assert summary["adoption"] == "adopted" and summary["adopted"] == 1
    with ProgressDB(db_path) as db:
        assert db.get_state(SPECIES, "SRR101") == "quantified"
        assert db.get_state(SPECIES, "SRR102") == "pending"


def test_cli_dry_run_and_cohort_hash_mismatch(world: World, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    args = [
        "--data-root",
        str(world.dest.base),
        "--cohort",
        str(world.base / "cohort"),
        "--config-dir",
        str(world.dest.config_dir),
        "--source",
        f"october={world.october.base}",
    ]
    before = _tree(world.dest.base)
    assert main([*args, "--manifest", str(tmp_path / "out" / "plan.jsonl")]) == 0
    assert "import=1" in capsys.readouterr().out
    assert _tree(world.dest.base) == before
    assert main(args) == 2  # a dry run needs an explicit manifest path
    # A changed cohort input byte refuses everything.
    table = world.base / "cohort" / "inputs" / SPECIES / "work" / "metadata" / "metadata_selected.tsv"
    os.chmod(table, 0o644)
    table.write_bytes(table.read_bytes().replace(b"PAIRED", b"PAIRED "))
    assert main([*args, "--manifest", str(tmp_path / "out" / "plan2.jsonl"), "--execute"]) == 2
    assert "differs from its definition" in capsys.readouterr().err
    assert _tree(world.dest.base) == before
    assert not (tmp_path / "out" / "plan2.jsonl").exists()


def test_config_dir_must_match_the_cohort(world: World, tmp_path: Path) -> None:
    config_dir = tmp_path / "edited"
    config_dir.mkdir()
    (config_dir / CONFIG_NAME).write_bytes(b"species_list: [Test_species]\nthreads: 8\n")
    code = main(
        [
            "--data-root",
            str(world.dest.base),
            "--cohort",
            str(world.base / "cohort"),
            "--config-dir",
            str(config_dir),
            "--source",
            f"october={world.october.base}",
            "--manifest",
            str(tmp_path / "plan.jsonl"),
        ]
    )
    assert code == 2


def test_index_fixture_bytes_are_the_cohort_index(world: World) -> None:
    assert world.cohort.entry(SPECIES).index_sha256 == sha(INDEX_BYTES)


def test_staging_left_by_a_killed_import_is_removed_without_a_new_import_decision(world: World) -> None:
    """The staging dir sits in a quantification directory this run imports nothing into: it is removed anyway."""

    execute_imports(world.cohort, world.dest.base, [world.source], workers=1)
    stale = world.dest.work / "quant" / ".import-SRR103-0badc0de"
    (stale / "staged" / "SRR103").mkdir(parents=True)
    (stale / STAGING_MARKER).write_text(json.dumps({"schema": STAGING_SCHEMA}))
    (stale / "staged" / "SRR103" / "SRR103_abundance.tsv").write_bytes(b"partial copy")

    plan, _, summary = execute_imports(world.cohort, world.dest.base, [world.source], workers=1)

    assert not [decision for decision in plan.decisions if decision.action == ACTION_IMPORTED]
    assert summary["stale_staging_removed"] == [stale.name] and not stale.exists()


def test_executed_imports_lists_finished_runs_per_source_root(world: World, tmp_path: Path) -> None:
    assert executed_imports(world.dest.base) == {"finished": {}, "unfinished": []}
    # A dry run (plan manifest), even one written under the destination, is not an executed import.
    plan = world.plan()
    write_plan_manifest(plan, world.dest.base / "results" / "imports" / "dry-run.jsonl", [world.source])
    assert executed_imports(world.dest.base)["finished"] == {}

    _, manifest, summary = execute_imports(world.cohort, world.dest.base, [world.source], workers=1)
    found = executed_imports(world.dest.base, definition_sha256=world.cohort.definition_sha256)
    record = found["finished"][os.path.normpath(str(world.october.base.absolute()))]
    assert record["manifest"] == str(manifest) and record["label"] == "october"
    assert (record["imported"], record["failed"]) == (summary["imported"], summary["failed"]) == (1, 0)
    assert found["unfinished"] == []
    # Another cohort definition does not count.
    assert executed_imports(world.dest.base, definition_sha256="0" * 64)["finished"] == {}

    # A killed import: its manifest has a header and decisions but no summary (and maybe a torn last line).
    lines = manifest.read_text(encoding="utf-8").splitlines()
    killed = manifest.with_name("killed.jsonl")
    killed.write_text("\n".join(lines[:-1]) + '\n{"kind": "deci', encoding="utf-8")
    again = executed_imports(world.dest.base)
    assert again["unfinished"] == [str(killed)]
    assert again["finished"][os.path.normpath(str(world.october.base.absolute()))]["manifest"] == str(manifest)


# ---------------------------------------------------------------------------
# Size-preserving damage (2026-11-15: 4 KiB blocks overwritten in place, sidecar intact)
# ---------------------------------------------------------------------------

_RSYNC_JUNK = b"\x00" * 48 + b"     1,048,576 100%   12.34MB/s    0:00:00 (xfr#1, to-chk=0/1)\r"


def _overwrite_block_in_place(path: Path, *, offset: int = 64, block: int = 4096) -> None:
    """Overwrite ``[offset, end of its 4 KiB block)`` with NULs and rsync progress text; the size never changes."""

    size = path.stat().st_size
    end = min((offset // block + 1) * block, size)
    junk = (_RSYNC_JUNK * ((end - offset) // len(_RSYNC_JUNK) + 1))[: end - offset]
    with path.open("r+b") as handle:
        handle.seek(offset)
        handle.write(junk)
    assert path.stat().st_size == size and junk and path.read_bytes()[offset:end] == junk


def _rewrite_in_place(path: Path, old: bytes, new: bytes) -> None:
    """Replace one same-length byte run (a plausible, still parseable edit); the size never changes."""

    assert len(old) == len(new)
    data = path.read_bytes()
    assert data.count(old) == 1
    path.write_bytes(data.replace(old, new))
    assert len(path.read_bytes()) == len(data)


def _recorded_quantification_sha(sample: Path) -> str:
    return json.loads((sample / QUANT_PROVENANCE_FILENAME).read_bytes())["quantification_file_sha256"]


def test_intact_pair_is_identical_on_the_sidecar_and_both_quantification_hashes(world: World) -> None:
    decision = _by_accession(world.plan())["SRR102"]
    assert decision.action == ACTION_EXISTS_IDENTICAL
    comparison = decision.details["comparison"]
    recorded = _recorded_quantification_sha(world.dest.sample("SRR102"))
    assert comparison["identity_basis"] == IDENTITY_BASIS_SHALLOW == "sidecar_bytes+sizes+quantification_sha256"
    assert comparison["source_quantification_sha256"] == comparison["destination_quantification_sha256"] == recorded
    assert comparison["source_quantification_matches_sidecar"] is True
    assert comparison["destination_quantification_matches_sidecar"] is True
    assert "difference" not in comparison


def test_destination_block_overwritten_in_place_is_never_identical(world: World) -> None:
    damaged = world.dest.sample("SRR102") / "SRR102_abundance.tsv"
    intact = world.october.sample("SRR102") / "SRR102_abundance.tsv"
    _overwrite_block_in_place(damaged)
    # Everything the former default basis (sidecar bytes + sizes) looked at still agrees.
    for name in sorted(path.name for path in world.october.sample("SRR102").iterdir()):
        assert (world.dest.sample("SRR102") / name).stat().st_size == (
            world.october.sample("SRR102") / name
        ).stat().st_size
    assert (world.dest.sample("SRR102") / QUANT_PROVENANCE_FILENAME).read_bytes() == (
        world.october.sample("SRR102") / QUANT_PROVENANCE_FILENAME
    ).read_bytes()
    before_dest, before_source = _tree(world.dest.base), _tree(world.october.base)

    plan = world.plan()
    decision = _by_accession(plan)["SRR102"]
    assert decision.action == ACTION_EXISTS_DIFFERENT
    assert decision.reason == DIFFERENCE_DESTINATION_QUANTIFICATION
    comparison = decision.details["comparison"]
    assert comparison["identity_basis"] == IDENTITY_BASIS_SHALLOW
    assert comparison["difference"] == "destination quantification does not match its sidecar"
    assert comparison["destination_quantification_matches_sidecar"] is False
    assert comparison["destination_quantification_sha256"] == sha(damaged.read_bytes())
    assert comparison["source_quantification_matches_sidecar"] is True
    assert comparison["source_quantification_sha256"] == sha(intact.read_bytes())
    # The sidecar-recorded digests agree, but they no longer describe the destination's table.
    assert decision.details["source_abundance_sha256"] == decision.details["destination_abundance_sha256"]
    assert decision.details["quantification"] == "unknown"
    assert decision.details["source_bound"] is True and decision.details["destination_bound"] is False
    assert decision.details["destination_verdict"].startswith("invalid_source:")
    assert decision.details["bound"] == "source_only"
    assert plan.counts()[f"october/{SPECIES}"].get(ACTION_EXISTS_IDENTICAL) is None
    assert _tree(world.dest.base) == before_dest and _tree(world.october.base) == before_source

    # Execution reports the same conflict and never replaces (or repairs) the damaged destination.
    plan, manifest, summary = execute_imports(world.cohort, world.dest.base, [world.source], workers=1)
    assert summary["imported"] == 1 and summary["failed"] == 0  # SRR101 only
    executed = _by_accession(plan)["SRR102"]
    assert (executed.action, executed.reason) == (ACTION_EXISTS_DIFFERENT, DIFFERENCE_DESTINATION_QUANTIFICATION)
    assert _tree(world.dest.sample("SRR102")) == {
        key.removeprefix(f"{SPECIES}/work/quant/SRR102/"): value
        for key, value in before_dest.items()
        if key.startswith(f"{SPECIES}/work/quant/SRR102/")
    }
    record = next(
        json.loads(line)
        for line in manifest.read_text().splitlines()
        if json.loads(line).get("kind") == "decision" and json.loads(line)["accession"] == "SRR102"
    )
    assert record["action"] == ACTION_EXISTS_DIFFERENT
    assert record["details"]["comparison"]["destination_quantification_matches_sidecar"] is False


def test_source_block_overwritten_in_place_is_not_importable(world: World) -> None:
    # Against an intact destination copy: a conflict naming the source, never "identical".
    _overwrite_block_in_place(world.october.sample("SRR102") / "SRR102_abundance.tsv")
    decision = _by_accession(world.plan())["SRR102"]
    assert (decision.action, decision.reason) == (ACTION_EXISTS_DIFFERENT, DIFFERENCE_SOURCE_QUANTIFICATION)
    comparison = decision.details["comparison"]
    assert comparison["difference"] == "source quantification does not match its sidecar"
    assert comparison["destination_quantification_matches_sidecar"] is True
    assert comparison["source_quantification_matches_sidecar"] is False
    assert decision.details["quantification"] == "unknown"
    assert decision.details["source_bound"] is False and decision.details["destination_bound"] is True

    # Absent from the destination: a still-parseable same-size edit is caught only by the hash, and skipped.
    table = world.october.sample("SRR101") / "SRR101_abundance.tsv"
    _rewrite_in_place(table, b"\t13.25\t", b"\t31.25\t")
    decision = _by_accession(world.plan())["SRR101"]
    assert (decision.action, decision.reason) == (ACTION_SKIPPED, "invalid_source:abundance checksum mismatch")
    plan, _, summary = execute_imports(world.cohort, world.dest.base, [world.source], workers=1)
    assert summary["imported"] == 0 and summary["failed"] == 0
    assert not [d for d in plan.decisions if d.action in (ACTION_IMPORT, ACTION_IMPORTED)]
    assert not world.dest.sample("SRR101").exists()


def test_damage_on_both_sides_names_the_destination_first(world: World) -> None:
    _overwrite_block_in_place(world.dest.sample("SRR102") / "SRR102_abundance.tsv", offset=40)
    _overwrite_block_in_place(world.october.sample("SRR102") / "SRR102_abundance.tsv", offset=96)
    decision = _by_accession(world.plan())["SRR102"]
    assert (decision.action, decision.reason) == (ACTION_EXISTS_DIFFERENT, DIFFERENCE_DESTINATION_QUANTIFICATION)
    comparison = decision.details["comparison"]
    assert comparison["source_quantification_matches_sidecar"] is False
    assert comparison["destination_quantification_matches_sidecar"] is False
    assert decision.details["bound"] == "neither"


def test_damaged_destination_with_a_different_sidecar_reports_the_damage(world: World) -> None:
    """The quantification check runs before the sidecar comparison, so a runtime-only conflict cannot mask damage."""

    sidecar = world.dest.sample("SRR102") / QUANT_PROVENANCE_FILENAME
    payload = json.loads(sidecar.read_text())
    payload["command"] = [*payload["command"], "--threads", "2"]
    sidecar.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    _overwrite_block_in_place(world.dest.sample("SRR102") / "SRR102_abundance.tsv")
    decision = _by_accession(world.plan())["SRR102"]
    assert (decision.action, decision.reason) == (ACTION_EXISTS_DIFFERENT, DIFFERENCE_DESTINATION_QUANTIFICATION)
    assert (
        decision.details["comparison"]["source_sidecar_sha256"]
        != decision.details["comparison"]["destination_sidecar_sha256"]
    )
    assert decision.details["quantification"] == "unknown"


def test_deep_compare_is_unchanged(world: World) -> None:
    decision = _by_accession(world.plan(deep_compare=True))["SRR102"]
    assert decision.action == ACTION_EXISTS_IDENTICAL
    assert decision.details["comparison"] == {
        "identity_basis": IDENTITY_BASIS_DEEP,
        "source_sidecar_sha256": decision.details["comparison"]["source_sidecar_sha256"],
        "destination_sidecar_sha256": decision.details["comparison"]["destination_sidecar_sha256"],
    }
    assert IDENTITY_BASIS_DEEP == "sha256"
    _overwrite_block_in_place(world.dest.sample("SRR102") / "SRR102_abundance.tsv")
    decision = _by_accession(world.plan(deep_compare=True))["SRR102"]
    assert decision.action == ACTION_EXISTS_DIFFERENT
    assert decision.details["comparison"]["difference"] == "file bytes differ: SRR102_abundance.tsv"
    assert "destination_quantification_matches_sidecar" not in decision.details["comparison"]


def test_quantification_must_be_a_regular_file_named_by_the_sidecar(world: World, tmp_path: Path) -> None:
    sample = world.october.sample("SRR101")
    sidecar_bytes = (sample / QUANT_PROVENANCE_FILENAME).read_bytes()
    table = sample / "SRR101_abundance.tsv"
    assert _quantification_against_sidecar(sample, sidecar_bytes) == (sha(table.read_bytes()), True)
    assert _quantification_against_sidecar(sample, b"not json") == (None, False)
    assert _quantification_against_sidecar(sample, b"[]") == (None, False)
    payload = json.loads(sidecar_bytes)
    for name in ("../SRR101/SRR101_abundance.tsv", "", "..", None):
        unsafe = json.dumps({**payload, "quantification_file": name}).encode()
        assert _quantification_against_sidecar(sample, unsafe) == (None, False)
    no_digest = json.dumps({**payload, "quantification_file_sha256": None}).encode()
    assert _quantification_against_sidecar(sample, no_digest) == (None, False)
    # A symlink to an intact copy is not the sample's own table; a missing table matches nothing.
    copy = tmp_path / "elsewhere.tsv"
    shutil.copyfile(table, copy)
    table.unlink()
    assert _quantification_against_sidecar(sample, sidecar_bytes) == (None, False)
    table.symlink_to(copy)
    assert _quantification_against_sidecar(sample, sidecar_bytes) == (None, False)

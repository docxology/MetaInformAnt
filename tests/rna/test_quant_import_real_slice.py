"""Import and cohort checks on a real two-species slice of the Hymenoptera campaign.

The slice (``$RNA_E2_FIXTURE_DIR/real_samples``, built read-only from the campaign volume) holds the real
frozen selection tables and kallisto indexes of atta_cephalotes and dinoponera_quadriceps, real October
(c19) sample directories with their reference-bound receipts, the canonical-root sample directories they
collide with, and ``candb/pipeline_progress.db`` (a copy of the canonical progress DB). Each test copies
what it changes into ``tmp_path``; nothing reads or writes outside the fixture directory and ``tmp_path``
(an audit hook fails the test on any access under ``/Volumes``).
"""

from __future__ import annotations

import json
import os
import shutil
import sqlite3
import sys
from pathlib import Path
from typing import Any

import pytest

from metainformant.rna.engine.cohort import CohortError, load_frozen_cohort
from metainformant.rna.engine.progress_db import ProgressDB
from metainformant.rna.engine.provenance import (
    QUANT_PROVENANCE_FILENAME,
    REFERENCE_BINDING_FILENAME,
    REUSABLE_QUANT_STATUSES,
    classify_cohort_quantification,
)
from metainformant.rna.engine.quant_import import (
    ACTION_EXISTS_DIFFERENT,
    ACTION_EXISTS_IDENTICAL,
    ACTION_IMPORT,
    ACTION_IMPORTED,
    ImportSource,
    ReceiptStore,
    execute_imports,
    plan_imports,
)

FIXTURE_DIR = os.environ.get("RNA_E2_FIXTURE_DIR")
SLICE = Path(FIXTURE_DIR) / "real_samples" if FIXTURE_DIR else None

pytestmark = pytest.mark.skipif(
    SLICE is None or not (SLICE / "SAMPLES.json").is_file(),
    reason="set RNA_E2_FIXTURE_DIR to a directory holding real_samples/ (see build_real_fixture.py)",
)

_VOLUME_ACCESS: list[tuple[str, str]] = []
_WATCHING = [False]


def _audit(event: str, args: tuple[Any, ...]) -> None:
    if _WATCHING[0] and event in {"open", "os.listdir", "os.scandir", "os.rename", "os.remove", "shutil.copyfile"}:
        for value in args[:2]:
            if isinstance(value, (str, os.PathLike)) and os.fspath(value).startswith("/Volumes"):
                _VOLUME_ACCESS.append((event, os.fspath(value)))


sys.addaudithook(_audit)


@pytest.fixture(autouse=True)
def _no_volume_access() -> Any:
    _VOLUME_ACCESS.clear()
    _WATCHING[0] = True
    try:
        yield
    finally:
        _WATCHING[0] = False
    assert not _VOLUME_ACCESS, f"code under test touched the campaign volume: {_VOLUME_ACCESS[:3]}"


def _tree(root: Path) -> dict[str, bytes]:
    return {str(path.relative_to(root)): path.read_bytes() for path in sorted(root.rglob("*")) if path.is_file()}


class RealSlice:
    """The cohort slice plus a destination root built from the canonical samples and a local index copy."""

    def __init__(self, tmp_path: Path) -> None:
        assert SLICE is not None
        self.samples = json.loads((SLICE / "SAMPLES.json").read_text())
        self.cohort = load_frozen_cohort(
            SLICE / "frozen",
            definition_path=SLICE / "COHORT.json",
            expected_definition_sha256=self.samples["definition_sha256"],
        )
        self.dest = tmp_path / "dest"
        shutil.copytree(SLICE / "canonical", self.dest)
        for species in self.samples["species"]:
            index_dir = self.dest / species / "genome" / "index"
            index_dir.mkdir(parents=True)
            index = index_dir / self.cohort.entry(species).index_name
            try:
                os.link(self.cohort.index_path(species), index)
            except OSError:
                shutil.copyfile(self.cohort.index_path(species), index)
            manifest = json.loads((SLICE / "canonical_reference" / f"{species}.json").read_text())
            for alias in manifest["aliases"]:
                alias["index"] = str(index_dir / Path(alias["index"]).name)
            manifest["index_dir"] = str(index_dir)
            manifest["kallisto_index"] = str(index_dir / Path(manifest["kallisto_index"]).name)
            reference = self.dest / species / "work" / "reference" / "reference_aliases.json"
            reference.parent.mkdir(parents=True)
            reference.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        self.source = ImportSource("october", SLICE / "october", "work")
        self.store = ReceiptStore("october-locked", SLICE / "october-locked", "local-20261019")

    def runs(self, role: str) -> list[tuple[str, str]]:
        return sorted((species, run) for species, runs in self.samples["roles"][role].items() for run in runs)

    def sample(self, root: Path, species: str, run: str) -> Path:
        return root / species / "work" / "quant" / run


@pytest.fixture
def real(tmp_path: Path) -> RealSlice:
    return RealSlice(tmp_path)


def test_real_slice_plan_matches_the_campaign_dry_run(real: RealSlice) -> None:
    before = _tree(real.dest)
    plan = plan_imports(real.cohort, real.dest, [real.source], receipt_stores=[real.store], workers=1)
    decisions = {(d.species, d.accession): d for d in plan.decisions}
    assert all(check.ok for check in plan.references.values())
    for key in real.runs("october_only"):
        decision = decisions[key]
        assert decision.action == ACTION_IMPORT, decision
        assert decision.evidence["kind"] == "receipt" and decision.evidence["index_bytes_verified"] is True
        assert decision.evidence["source_evidence"] == "local_restore_proof"
    for key in real.runs("identical"):
        assert decisions[key].action == ACTION_EXISTS_IDENTICAL
    (differs,) = real.runs("conflict_quantification_differs")
    (same,) = real.runs("conflict_same_quantification")
    for key, relation in ((differs, "different"), (same, "identical")):
        details = decisions[key].details
        assert decisions[key].action == ACTION_EXISTS_DIFFERENT
        assert details["quantification"] == relation
        # Both versions are bound to the cohort index, but only the October one proves the index bytes
        # (receipt); the canonical one carries a restored witness that names the index.
        assert details["bound"] == "both" and details["index_bytes_proven"] == "source_only"
        assert details["destination_evidence"]["kind"] == "restored_witness"
    assert (
        decisions[differs].details["source_abundance_sha256"]
        != decisions[differs].details["destination_abundance_sha256"]
    )
    assert _tree(real.dest) == before


def test_real_slice_imports_october_only_runs_and_rerun_imports_nothing(real: RealSlice, tmp_path: Path) -> None:
    source_before = _tree(SLICE / "october")
    conflicts = real.runs("conflict_quantification_differs") + real.runs("conflict_same_quantification")
    conflict_before = {key: _tree(real.sample(real.dest, *key)) for key in conflicts}

    plan, manifest, summary = execute_imports(
        real.cohort, real.dest, [real.source], receipt_stores=[real.store], manifest_path=tmp_path / "m.jsonl"
    )
    october_only = real.runs("october_only")
    assert summary["imported"] == len(october_only) and summary["failed"] == 0
    assert sorted((d.species, d.accession) for d in plan.decisions if d.action == ACTION_IMPORTED) == october_only
    for species, run in october_only:
        imported = real.sample(real.dest, species, run)
        source = real.sample(SLICE / "october", species, run)
        for name in (QUANT_PROVENANCE_FILENAME, f"{run}_abundance.tsv", f"{run}_run_info.json", f"{run}_abundance.h5"):
            assert (imported / name).read_bytes() == (source / name).read_bytes()
        assert (imported / REFERENCE_BINDING_FILENAME).is_file()
        result = classify_cohort_quantification(imported, real.cohort.expectation(species, run), verify_content=True)
        assert result["status"] in REUSABLE_QUANT_STATUSES and result["evidence"]["kind"] == "binding_witness"
    assert {key: _tree(real.sample(real.dest, *key)) for key in conflicts} == conflict_before
    assert _tree(SLICE / "october") == source_before
    assert len([line for line in manifest.read_text().splitlines() if '"kind": "decision"' in line]) == 6

    tree = _tree(real.dest)
    plan, _, summary = execute_imports(
        real.cohort, real.dest, [real.source], receipt_stores=[real.store], manifest_path=tmp_path / "m2.jsonl"
    )
    assert summary["imported"] == 0
    actions = {(d.species, d.accession): d.action for d in plan.decisions}
    assert all(actions[key] == ACTION_EXISTS_IDENTICAL for key in october_only)
    assert all(actions[key] == ACTION_EXISTS_DIFFERENT for key in conflicts)
    assert _tree(real.dest) == tree


def test_real_slice_adopts_imports_into_a_copy_of_the_canonical_db(real: RealSlice, tmp_path: Path) -> None:
    assert FIXTURE_DIR is not None
    canonical_db = Path(FIXTURE_DIR) / "candb" / "pipeline_progress.db"
    if not canonical_db.is_file():
        pytest.skip("candb/pipeline_progress.db (copy of the canonical progress DB) is not in the fixture dir")
    db_path = tmp_path / "pipeline_progress.db"
    shutil.copyfile(canonical_db, db_path)

    def rows() -> dict[tuple[str, str], tuple[Any, ...]]:
        with sqlite3.connect(db_path) as conn:
            return {(r[0], r[1]): r[2:] for r in conn.execute("SELECT species, srr_id, state, error FROM samples")}

    before = rows()
    october_only = real.runs("october_only")
    assert not any(before.get(key, ("missing",))[0] == "quantified" for key in october_only)
    _, _, summary = execute_imports(
        real.cohort, real.dest, [real.source], receipt_stores=[real.store], adopt_db=db_path, workers=1
    )
    assert summary["adopted"] == len(october_only)
    after = rows()
    for key in october_only:
        assert after[key] == ("quantified", None)
    assert {key: value for key, value in after.items() if key not in october_only} == {
        key: value for key, value in before.items() if key not in october_only
    }
    quantified_before = sum(1 for value in before.values() if value[0] == "quantified")
    assert sum(1 for value in after.values() if value[0] == "quantified") == quantified_before + len(october_only)
    with ProgressDB(db_path) as db:
        assert all(db.get_state(species, run) == "quantified" for species, run in october_only)


def test_real_slice_cohort_refuses_a_changed_table_byte_or_a_wrong_pin(real: RealSlice, tmp_path: Path) -> None:
    assert SLICE is not None
    with pytest.raises(CohortError, match="pinned"):
        load_frozen_cohort(SLICE / "frozen", definition_path=SLICE / "COHORT.json", expected_definition_sha256="0" * 64)
    frozen = tmp_path / "frozen"
    shutil.copytree(SLICE / "frozen", frozen, copy_function=os.link)
    table = frozen / real.cohort.entry("atta_cephalotes").metadata_relpath
    data = bytearray(table.read_bytes())
    table.unlink()  # break the hard link so the fixture copy stays untouched
    data[-3] ^= 0x01
    table.write_bytes(bytes(data))
    with pytest.raises(CohortError, match="differs from its definition"):
        load_frozen_cohort(frozen, definition_path=SLICE / "COHORT.json")

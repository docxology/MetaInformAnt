"""Quantification snapshot export: evidence, idempotency, conflicts, witnesses and Git-tree presence.

Everything runs on real files under ``tmp_path``: a frozen cohort installed by
the production ``install_cohort``, outputs written by ``write_quant_provenance``,
a reference binding written by ``write_reference_binding``, and a real Git
repository for tree-mode presence. The CLI is the production entry point.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from metainformant.rna.engine.cohort import install_cohort, load_frozen_cohort
from metainformant.rna.engine.provenance import (
    QUANT_PROVENANCE_FILENAME,
    REFERENCE_BINDING_FILENAME,
    reference_binding_payload,
    write_reference_binding,
)
from metainformant.rna.engine.quant_snapshot import (
    EXIT_CONFLICT,
    FileDigest,
    compare_destination,
    git_blob_oid,
    lfs_pointer,
    main,
)
from tests.rna.frozen_cohort_support import COHORT_ID, SPECIES, DataRoot, build_frozen_source, sha

SNAPSHOT_FILES = (
    QUANT_PROVENANCE_FILENAME,
    "SRR101_abundance.h5",
    "SRR101_abundance.tsv",
    "SRR101_run_info.json",
)


def _tree(root: Path) -> dict[str, bytes]:
    return {str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()}


class Campaign:
    """A frozen cohort, a data root with one reusable and one non-reusable output, and a destination."""

    def __init__(self, base: Path) -> None:
        build_frozen_source(base / "frozen")
        install_cohort(base / "frozen", base / "cohort", cohort_id=COHORT_ID)
        self.base = base
        self.cohort = base / "cohort"
        self.root = DataRoot(base / "data")
        self.root.write_sample("SRR101")  # current
        other = base / "other.yaml"
        other.write_bytes(b"species_list: [Test_species]\nthreads: 2\n")
        self.root.write_sample("SRR103", config=other)  # configuration differs: not reusable
        self.repo = base / "repo"
        self.dest = self.repo / "quant" / "completed" / COHORT_ID
        self.inventory = self.repo / "quant" / "inventory" / COHORT_ID

    def args(self, *extra: str) -> list[str]:
        return [
            "--data-root",
            str(self.root.base),
            "--config-dir",
            str(self.root.config_dir),
            "--cohort",
            str(self.cohort),
            "--dest",
            str(self.dest),
            "--workers",
            "2",
            *extra,
        ]

    def bind(self, run: str) -> None:
        sample = self.root.sample(run)
        sidecar_bytes = (sample / QUANT_PROVENANCE_FILENAME).read_bytes()
        cohort = load_frozen_cohort(self.cohort)
        payload = reference_binding_payload(
            sidecar_bytes, json.loads(sidecar_bytes), cohort.expectation(SPECIES, run), {"kind": "produced"}
        )
        write_reference_binding(sample, payload)


@pytest.fixture
def campaign(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Campaign:
    monkeypatch.delenv("AMALGKIT_DATA_ROOT", raising=False)
    return Campaign(tmp_path)


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), "-c", "user.name=test", "-c", "user.email=test@example.invalid", *args],
        check=True,
        capture_output=True,
        text=True,
    ).stdout


def test_plan_writes_nothing_and_exports_only_reusable_outputs(
    campaign: Campaign, capsys: pytest.CaptureFixture[str]
) -> None:
    source_before = _tree(campaign.root.base)
    assert main(campaign.args()) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["mode"] == "plan"
    assert report["snapshot_states"] == {"new": 1, "not_reusable": 1}
    assert report["to_publish"]["new"] == 1
    assert [entry["run"] for entry in report["not_exported"]] == ["SRR103"]
    assert report["task_states"] == {"not_reusable": 1, "pending": 2, "present": 1}
    assert not campaign.repo.exists()
    assert _tree(campaign.root.base) == source_before


def test_execute_copies_exact_bytes_writes_inventory_and_is_idempotent(
    campaign: Campaign, capsys: pytest.CaptureFixture[str]
) -> None:
    source_before = _tree(campaign.root.base)
    args = campaign.args("--execute", "--inventory-dir", str(campaign.inventory))
    assert main(args) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["published"] == 1
    published = campaign.dest / SPECIES / "SRR101"
    assert sorted(_tree(published)) == sorted(SNAPSHOT_FILES)
    for name in SNAPSHOT_FILES:
        assert (published / name).read_bytes() == (campaign.root.sample("SRR101") / name).read_bytes()
    assert not (campaign.dest / SPECIES / "SRR103").exists(), "a non-reusable output is never exported"
    assert _tree(campaign.root.base) == source_before, "the data root is only read"

    rows = (campaign.inventory / f"{SPECIES}.tsv").read_text().splitlines()
    assert rows[0] == "run\tstatus\tfile\tsize\tsha256\tgit_oid"
    tsv_row = next(row.split("\t") for row in rows if row.split("\t")[2] == "SRR101_abundance.tsv")
    tsv_bytes = (published / "SRR101_abundance.tsv").read_bytes()
    assert tsv_row[4] == sha(tsv_bytes) and tsv_row[5] == git_blob_oid(tsv_bytes)
    h5_row = next(row.split("\t") for row in rows if row.split("\t")[2] == "SRR101_abundance.h5")
    h5_bytes = (published / "SRR101_abundance.h5").read_bytes()
    assert h5_row[5] == git_blob_oid(lfs_pointer(sha(h5_bytes), len(h5_bytes))), "HDF5 rows carry the LFS pointer id"
    tasks = dict(
        (line.split("\t")[1], line.split("\t")[3])
        for line in (campaign.inventory / "tasks.tsv").read_text().splitlines()[1:]
    )
    assert tasks == {"SRR101": "present", "SRR102": "pending", "SRR103": "not_reusable", "SRR105": "pending"}
    summary = json.loads((campaign.inventory / "SUMMARY.json").read_text())
    assert summary["exported_samples"] == 1 and summary["cohort_id"] == COHORT_ID
    assert "generated_at" not in summary and str(campaign.base) not in json.dumps(summary)

    snapshot_before = _tree(campaign.repo)
    assert main(args) == 0
    again = json.loads(capsys.readouterr().out)
    assert again["published"] == 0 and again["inventory_changed"] == []
    assert again["snapshot_states"]["present"] == 1
    assert _tree(campaign.repo) == snapshot_before, "a second run rewrites nothing"


def test_different_destination_bytes_are_a_conflict_and_never_overwritten(
    campaign: Campaign, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main(campaign.args("--execute")) == 0
    capsys.readouterr()
    table = campaign.dest / SPECIES / "SRR101" / "SRR101_abundance.tsv"
    table.write_text("target_id\tlength\teff_length\test_counts\ttpm\naltered\t1\t1\t1\t1\n")
    altered = table.read_bytes()
    assert main(campaign.args("--execute")) == EXIT_CONFLICT
    report = json.loads(capsys.readouterr().out)
    assert report["conflicts"] == [
        {"species": SPECIES, "run": "SRR101", "reason": "different bytes: SRR101_abundance.tsv"}
    ]
    assert table.read_bytes() == altered


def test_a_new_witness_extends_an_identical_sample_without_rewriting_it(
    campaign: Campaign, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main(campaign.args("--execute")) == 0
    capsys.readouterr()
    published = campaign.dest / SPECIES / "SRR101"
    before = _tree(published)
    campaign.bind("SRR101")
    assert main(campaign.args("--execute")) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["snapshot_states"]["extend"] == 1 and report["to_publish"]["extend"] == 1
    after = _tree(published)
    assert {name: after[name] for name in before} == before
    binding = campaign.root.sample("SRR101") / REFERENCE_BINDING_FILENAME
    assert after[REFERENCE_BINDING_FILENAME] == binding.read_bytes()
    assert main(campaign.args("--execute")) == 0
    assert json.loads(capsys.readouterr().out)["published"] == 0


def test_interrupted_staging_is_removed_on_the_next_execute(
    campaign: Campaign, capsys: pytest.CaptureFixture[str]
) -> None:
    leftover = campaign.dest / SPECIES / ".SRR101.snapshot-staging" / "SRR101_abundance.tsv"
    leftover.parent.mkdir(parents=True)
    leftover.write_text("partial")
    assert main(campaign.args("--execute")) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["removed_stale_staging"] == [f"{SPECIES}/.SRR101.snapshot-staging"]
    assert not leftover.parent.exists()
    assert (campaign.dest / SPECIES / "SRR101" / "SRR101_abundance.tsv").is_file()


@pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")
def test_git_tree_presence_covers_a_sparse_checkout_and_lfs_pointers(
    campaign: Campaign, capsys: pytest.CaptureFixture[str]
) -> None:
    campaign.repo.mkdir()
    _git(campaign.repo, "init", "-q")
    # Commit SRR101 as a repository with Git LFS would store it: the HDF5 file as its pointer.
    source = campaign.root.sample("SRR101")
    committed = campaign.dest / SPECIES / "SRR101"
    committed.mkdir(parents=True)
    for name in SNAPSHOT_FILES:
        data = (source / name).read_bytes()
        (committed / name).write_bytes(lfs_pointer(sha(data), len(data)) if name.endswith(".h5") else data)
    _git(campaign.repo, "add", "-A")
    _git(campaign.repo, "commit", "-q", "-m", "snapshot")
    shutil.rmtree(campaign.dest)  # what a sparse checkout that excludes the snapshot leaves

    tree_args = campaign.args("--git-tree", f"{campaign.repo}:HEAD")
    assert main(tree_args) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["snapshot_states"] == {"not_reusable": 1, "present": 1}
    assert report["to_publish"]["new"] == 0

    campaign.root.write_sample("SRR102")
    assert main([*tree_args, "--execute"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["published"] == 1 and (campaign.dest / SPECIES / "SRR102").is_dir()
    assert not (campaign.dest / SPECIES / "SRR101").exists(), "committed samples are not re-materialized"
    assert main(tree_args) == 0
    assert json.loads(capsys.readouterr().out)["snapshot_states"] == {
        "not_reusable": 1,
        "present": 1,
        "staged": 1,
    }


def test_destination_inside_the_data_root_is_refused(campaign: Campaign) -> None:
    campaign.dest = campaign.root.base / "snapshot"
    with pytest.raises(SystemExit, match="inside the data root"):
        main(campaign.args())


def test_inventory_needs_the_whole_cohort(campaign: Campaign) -> None:
    with pytest.raises(SystemExit, match="whole cohort"):
        main(campaign.args("--species", SPECIES, "--inventory-dir", str(campaign.inventory)))


def _digest(relpath: str, oid: str) -> FileDigest:
    return FileDigest(relpath=relpath, size=1, sha256="0" * 64, git_oid=oid)


@pytest.mark.parametrize(
    ("committed", "on_disk", "tree_mode", "state"),
    [
        ({}, {}, True, "new"),
        ({"a.tsv": "1"}, {}, True, "present"),
        ({}, {"a.tsv": "1"}, True, "staged"),
        ({}, {"a.tsv": "1"}, False, "present"),
        ({"a.tsv": "2"}, {}, True, "conflict"),
        ({"a.tsv": "1"}, {"a.tsv": "2"}, True, "conflict"),
    ],
)
def test_compare_destination_states(
    committed: dict[str, str], on_disk: dict[str, str], tree_mode: bool, state: str
) -> None:
    assert compare_destination([_digest("a.tsv", "1")], committed, on_disk, tree_mode=tree_mode)[0] == state


def test_compare_destination_adds_only_witnesses() -> None:
    expected = [_digest("a.tsv", "1"), _digest(REFERENCE_BINDING_FILENAME, "2")]
    state, _, to_copy = compare_destination(expected, {"a.tsv": "1"}, {}, tree_mode=True)
    assert state == "extend" and [d.relpath for d in to_copy] == [REFERENCE_BINDING_FILENAME]
    state, reason, to_copy = compare_destination(
        [_digest("a.tsv", "1"), _digest("b.json", "3")], {"a.tsv": "1"}, {}, tree_mode=True
    )
    assert state == "conflict" and "b.json" in reason and to_copy == []

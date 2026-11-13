"""Content-addressed, path-independent quantification reuse (quant reuse contract v2).

Every test builds a real data root in ``tmp_path``: a configuration YAML, a
reference manifest, deterministic index bytes, a kallisto-style abundance table
and run_info, and a sidecar written by the production writer. Nothing here
mocks the code under test; "bumped constants" tests replace the pinned Amalgkit
runtime values the same way a version bump would.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from pathlib import Path

import pytest

from metainformant.rna.engine import provenance
from metainformant.rna.engine.durable_quant import (
    DirectoryStore,
    lock_quantification,
    restore_quantification,
)
from metainformant.rna.engine.provenance import (
    EVIDENCE_BINDING_WITNESS,
    EVIDENCE_LOCAL_RESTORE_PROOF,
    EVIDENCE_REFERENCE_MANIFEST,
    EVIDENCE_RESTORED_WITNESS,
    LOCAL_RESTORE_PROOF_FILENAME,
    QUANT_PROVENANCE_FILENAME,
    QUANT_STATUS_CURRENT,
    QUANT_STATUS_INVALID,
    QUANT_STATUS_LEGACY,
    QUANT_STATUS_VERSION_DRIFT,
    REFERENCE_BINDING_FILENAME,
    CohortExpectation,
    build_species_reuse_context,
    classify_cohort_quantification,
    classify_quantification,
    is_current_metadata,
    metadata_runtime_status,
    quantification_contract_id,
    read_reference_binding,
    reference_binding_payload,
    write_metadata_provenance,
    write_quant_provenance,
    write_reference_binding,
)
from metainformant.rna.engine.quant_validation import validate_quantification
from metainformant.rna.engine.reference_digest import IndexDigestCache

COHORT = "cohort-test-20261112"
SPECIES = "test_species"
ACCESSION = "SRR123"
INDEX_BYTES = b"deterministic kallisto index fixture bytes\n" * 64
CONFIG_BYTES = b"species_list: [Test_species]\nthreads: 1\n"
# A hypothetical later Amalgkit build. It must never equal the real pin, or a
# "bump" would be a no-op and every drift assertion below would test nothing.
BUMPED_RUNTIME = {
    "REQUIRED_AMALGKIT_VERSION": "0.16.999",
    "AMALGKIT_RELEASE_TAG": "master@0123456789ab",
    "AMALGKIT_SOURCE_REVISION": "0123456789ab" + "c" * 28,
}
# The runtime recorded by every sidecar on the canonical data root (the pin
# before the 0.16.108 upgrade).
PREVIOUS_PIN_RUNTIME = {
    "amalgkit_version": "0.16.60",
    "amalgkit_release_tag": "v0.16.60",
    "amalgkit_source_revision": "c656a52aacdcee6fd3bf7e8031769ca957204ebc",
}


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _bump_runtime(monkeypatch: pytest.MonkeyPatch) -> None:
    for name, value in BUMPED_RUNTIME.items():
        assert getattr(provenance, name) != value, f"BUMPED_RUNTIME[{name!r}] equals the real pin"
        monkeypatch.setattr(provenance, name, value)


def _record_runtime(sample: Path, runtime: dict[str, str]) -> None:
    """Rewrite a sidecar's runtime fields, as an output of an older Amalgkit build would carry them."""

    sidecar = sample / QUANT_PROVENANCE_FILENAME
    payload = json.loads(sidecar.read_text())
    payload.update(runtime)
    sidecar.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


class Root:
    """A tiny but complete data root: config, index (+ alias), manifest, one sample."""

    def __init__(self, base: Path, accessions: tuple[str, ...] = (ACCESSION,)) -> None:
        self.base = base
        self.data = base / "data"
        self.repo_config = base / "repo" / "config" / "amalgkit_test_species.yaml"
        self.repo_config.parent.mkdir(parents=True)
        self.repo_config.write_bytes(CONFIG_BYTES)
        index_dir = self.data / SPECIES / "genome" / "index"
        index_dir.mkdir(parents=True)
        self.index = index_dir / "Test_species.idx"
        self.index.write_bytes(INDEX_BYTES)
        self.alias = index_dir / "Test_species_subspecies.idx"
        self.alias.write_bytes(INDEX_BYTES)
        self.manifest = self.data / SPECIES / "work" / "reference" / "reference_aliases.json"
        self.manifest.parent.mkdir(parents=True)
        self.manifest.write_text(
            json.dumps(
                {
                    "aliases": [
                        {"index": str(self.index), "metadata_target": "Test species"},
                        {"index": str(self.alias), "metadata_target": "Test species subspecies"},
                    ],
                    "index_dir": str(index_dir),
                    "kallisto_index": str(self.index),
                    "schema": "metainformant.rna.reference_aliases.v1",
                    "species": SPECIES,
                    "status": "complete",
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        for accession in accessions:
            self.write_sample(accession)

    def sample(self, accession: str = ACCESSION) -> Path:
        return self.data / SPECIES / "work" / "quant" / accession

    def write_sample(self, accession: str) -> Path:
        sample = self.sample(accession)
        sample.mkdir(parents=True)
        rows = "".join(f"transcript{i}\t1000\t850.5\t{10 + i}.25\t{1000 + i}.5\n" for i in range(8))
        (sample / f"{accession}_abundance.tsv").write_text("target_id\tlength\teff_length\test_counts\ttpm\n" + rows)
        (sample / f"{accession}_run_info.json").write_text(
            json.dumps(
                {
                    "n_targets": 8,
                    "n_processed": 1000,
                    "n_pseudoaligned": 900,
                    "kallisto_version": "0.52.0",
                    "index_version": 13,
                }
            )
        )
        write_quant_provenance(
            sample,
            species=SPECIES,
            run_accession=accession,
            config_path=self.repo_config,
            command=["amalgkit", "quant", "--out_dir", str(self.data / SPECIES / "work"), "--batch", "1"],
            reference_manifest_path=self.manifest,
        )
        return sample

    def expectation(self, accession: str = ACCESSION, **overrides: str) -> CohortExpectation:
        values = {
            "cohort_id": COHORT,
            "species": SPECIES,
            "run_accession": accession,
            "config_sha256": _sha(CONFIG_BYTES),
            "reference_index_sha256": _sha(INDEX_BYTES),
            "reference_index_name": "Test_species.idx",
        }
        values.update(overrides)
        return CohortExpectation(**values)


@pytest.fixture
def root(tmp_path: Path) -> Root:
    return Root(tmp_path)


def _seal(sample: Path, expectation: CohortExpectation, evidence: dict[str, object] | None = None) -> Path:
    sidecar_bytes = (sample / QUANT_PROVENANCE_FILENAME).read_bytes()
    payload = reference_binding_payload(
        sidecar_bytes,
        json.loads(sidecar_bytes),
        expectation,
        evidence or {"kind": "produced", "index_name": expectation.reference_index_name},
    )
    return write_reference_binding(sample, payload)


# ---------------------------------------------------------------------------
# Configuration is bound by content, not by path
# ---------------------------------------------------------------------------


def test_deleted_config_path_with_identical_bytes_stays_reusable(root: Root) -> None:
    current_config = root.base / "moved-checkout" / "amalgkit_test_species.yaml"
    current_config.parent.mkdir()
    current_config.write_bytes(CONFIG_BYTES)
    root.repo_config.unlink()  # the checkout recorded in every sidecar is gone

    # Today's path-based check is what would have quarantined the canonical outputs.
    assert classify_quantification(root.sample(), ACCESSION)["reason"] == "configuration checksum mismatch"

    context = build_species_reuse_context(SPECIES, current_config)
    legacy = classify_quantification(
        root.sample(), ACCESSION, expected_config_sha256=context.config_sha256, verify_content=True
    )
    assert legacy["status"] == QUANT_STATUS_CURRENT
    assert legacy["evidence"] == {"config": "content", "reference": f"{EVIDENCE_REFERENCE_MANIFEST}:relative"}

    cohort = classify_cohort_quantification(root.sample(), root.expectation(), verify_content=True)
    assert cohort["status"] == QUANT_STATUS_CURRENT
    assert cohort["evidence"] == {
        "kind": EVIDENCE_REFERENCE_MANIFEST,
        "manifest": "relative",
        "index_paths": "recorded",
    }


def test_different_config_bytes_are_invalid(root: Root) -> None:
    other = root.base / "other.yaml"
    other.write_bytes(CONFIG_BYTES + b"# edited\n")
    context = build_species_reuse_context(SPECIES, other)

    legacy = classify_quantification(root.sample(), ACCESSION, expected_config_sha256=context.config_sha256)
    assert legacy["status"] == QUANT_STATUS_INVALID
    assert legacy["reason"] == "configuration differs from current runtime configuration"

    cohort = classify_cohort_quantification(root.sample(), root.expectation(config_sha256=context.config_sha256))
    assert cohort["status"] == QUANT_STATUS_INVALID
    assert cohort["reason"] == "configuration differs from the cohort configuration"


def test_expected_config_sha256_keyword_keeps_default_behavior(root: Root) -> None:
    assert classify_quantification(root.sample(), ACCESSION)["status"] == QUANT_STATUS_CURRENT
    assert classify_quantification(root.sample(), ACCESSION)["evidence"]["config"] == "recorded_path"
    root.repo_config.write_bytes(b"edited in place\n")
    assert classify_quantification(root.sample(), ACCESSION)["status"] == QUANT_STATUS_INVALID
    # Content binding ignores the recorded path entirely.
    assert (
        classify_quantification(root.sample(), ACCESSION, expected_config_sha256=_sha(CONFIG_BYTES))["status"]
        == QUANT_STATUS_CURRENT
    )


# ---------------------------------------------------------------------------
# Reference binding: index bytes, witnesses and proofs
# ---------------------------------------------------------------------------


def test_moved_data_root_rebinds_by_relative_manifest_and_index_bytes(root: Root) -> None:
    moved = root.base / "relocated"
    shutil.copytree(root.data, moved)
    shutil.rmtree(root.data)
    sample = moved / SPECIES / "work" / "quant" / ACCESSION

    cohort = classify_cohort_quantification(sample, root.expectation(), verify_content=True)
    assert cohort["status"] == QUANT_STATUS_CURRENT
    assert cohort["evidence"]["index_paths"] == "rerooted"

    legacy = classify_quantification(sample, ACCESSION, expected_config_sha256=_sha(CONFIG_BYTES))
    assert legacy["status"] == QUANT_STATUS_CURRENT


def test_wrong_or_missing_index_bytes_are_invalid(root: Root) -> None:
    root.alias.write_bytes(INDEX_BYTES[:-1] + b"X")
    result = classify_cohort_quantification(root.sample(), root.expectation())
    assert result["status"] == QUANT_STATUS_INVALID
    assert result["reason"] == "reference index bytes differ from the cohort index"

    root.alias.unlink()
    result = classify_cohort_quantification(root.sample(), root.expectation())
    assert result["status"] == QUANT_STATUS_INVALID
    assert result["reason"] == "reference index unavailable: Test_species_subspecies.idx"


def test_rewritten_manifest_without_other_evidence_is_invalid(root: Root) -> None:
    root.manifest.write_text(root.manifest.read_text().replace('"complete"', '"complete" '))
    result = classify_cohort_quantification(root.sample(), root.expectation())
    assert result["status"] == QUANT_STATUS_INVALID
    assert result["reason"] == "no reference-index binding evidence"


def test_binding_witness_decouples_reuse_from_the_mutable_manifest(root: Root) -> None:
    witness = _seal(root.sample(), root.expectation())
    assert read_reference_binding(root.sample())["reference_index_sha256"] == _sha(INDEX_BYTES)
    root.manifest.unlink()
    root.index.unlink()
    result = classify_cohort_quantification(root.sample(), root.expectation())
    assert result["status"] == QUANT_STATUS_CURRENT
    assert result["evidence"] == {
        "kind": EVIDENCE_BINDING_WITNESS,
        "binding": {"kind": "produced", "index_name": "Test_species.idx"},
    }
    # Downstream tools that know no cohort accept the same witness once the
    # recorded config path is gone.
    root.repo_config.unlink()
    legacy = classify_quantification(root.sample(), ACCESSION)
    assert legacy["status"] == QUANT_STATUS_CURRENT
    assert legacy["evidence"] == {"config": EVIDENCE_BINDING_WITNESS, "reference": EVIDENCE_BINDING_WITNESS}
    assert witness.read_bytes().count(b"/") == 0  # no paths inside the witness


def test_tampered_sidecar_byte_breaks_the_binding_witness(root: Root) -> None:
    _seal(root.sample(), root.expectation())
    sidecar = root.sample() / QUANT_PROVENANCE_FILENAME
    sidecar.write_bytes(sidecar.read_bytes().replace(b"\n}", b"\n }"))
    for result in (
        classify_cohort_quantification(root.sample(), root.expectation()),
        classify_quantification(root.sample(), ACCESSION),
    ):
        assert result["status"] == QUANT_STATUS_INVALID
        assert (
            result["reason"]
            == "invalid reference binding witness: sidecar_sha256 differs from the quantification sidecar"
        )


def test_binding_witness_for_another_cohort_or_index_is_invalid(root: Root) -> None:
    _seal(root.sample(), root.expectation(cohort_id="another-cohort"))
    result = classify_cohort_quantification(root.sample(), root.expectation())
    assert result["status"] == QUANT_STATUS_INVALID
    assert result["reason"] == "reference binding witness cohort_id differs from the cohort"

    other = Root(root.base / "second")
    _seal(other.sample(), other.expectation(reference_index_sha256="a" * 64))
    result = classify_cohort_quantification(other.sample(), other.expectation())
    assert result["reason"] == "reference binding witness reference_index_sha256 differs from the cohort"


def test_write_reference_binding_is_idempotent_and_never_replaces(root: Root) -> None:
    first = _seal(root.sample(), root.expectation())
    before = first.stat()
    assert _seal(root.sample(), root.expectation()) == first
    after = first.stat()
    assert (before.st_ino, before.st_mtime_ns) == (after.st_ino, after.st_mtime_ns)
    original = first.read_bytes()
    with pytest.raises(FileExistsError):
        _seal(root.sample(), root.expectation(), {"kind": "receipt", "receipt_cohort": "x"})
    assert first.read_bytes() == original
    assert not [p for p in root.sample().iterdir() if p.name.startswith(f".{REFERENCE_BINDING_FILENAME}.")]


def test_binding_payload_refuses_absolute_paths_and_mismatched_sidecars(root: Root) -> None:
    sidecar_bytes = (root.sample() / QUANT_PROVENANCE_FILENAME).read_bytes()
    sidecar = json.loads(sidecar_bytes)
    with pytest.raises(ValueError, match="absolute paths"):
        reference_binding_payload(
            sidecar_bytes, sidecar, root.expectation(), {"kind": "receipt", "store": "/Volumes/x"}
        )
    with pytest.raises(ValueError, match="identity"):
        reference_binding_payload(sidecar_bytes, sidecar, root.expectation("SRR999"), {"kind": "produced"})
    with pytest.raises(ValueError, match="configuration"):
        reference_binding_payload(sidecar_bytes, sidecar, root.expectation(config_sha256="b" * 64), {"kind": "x"})
    assert not (root.sample() / REFERENCE_BINDING_FILENAME).exists()


def test_local_restore_proof_is_accepted_only_for_its_contract_and_index(root: Root) -> None:
    root.manifest.unlink()
    contract = quantification_contract_id(json.loads((root.sample() / QUANT_PROVENANCE_FILENAME).read_text()))
    proof = root.sample() / LOCAL_RESTORE_PROOF_FILENAME
    proof.write_text(json.dumps({"contract_id": contract, "reference_index_sha256": _sha(INDEX_BYTES), "x": 1.5}))
    result = classify_cohort_quantification(root.sample(), root.expectation())
    assert result["status"] == QUANT_STATUS_CURRENT
    assert result["evidence"] == {"kind": EVIDENCE_LOCAL_RESTORE_PROOF}

    proof.write_text(json.dumps({"contract_id": "c" * 64, "reference_index_sha256": _sha(INDEX_BYTES)}))
    assert classify_cohort_quantification(root.sample(), root.expectation())["status"] == QUANT_STATUS_INVALID
    proof.write_text(json.dumps({"contract_id": contract, "reference_index_sha256": "d" * 64}))
    assert classify_cohort_quantification(root.sample(), root.expectation())["status"] == QUANT_STATUS_INVALID


def _locked_store(root: Root) -> DirectoryStore:
    store = DirectoryStore(root.base / "store")
    lock_quantification(
        store,
        COHORT,
        SPECIES,
        ACCESSION,
        root.sample(),
        expected_config_sha256=_sha(CONFIG_BYTES),
        expected_reference_index_sha256=_sha(INDEX_BYTES),
    )
    return store


def test_s3_receipt_restore_under_bumped_runtime_is_drift_compatible(
    root: Root, monkeypatch: pytest.MonkeyPatch
) -> None:
    store = _locked_store(root)
    _bump_runtime(monkeypatch)
    destination = root.base / "elsewhere" / SPECIES / "work" / "quant" / ACCESSION
    receipt = restore_quantification(
        store,
        COHORT,
        SPECIES,
        ACCESSION,
        destination,
        expected_config_sha256=_sha(CONFIG_BYTES),
        expected_reference_index_sha256=_sha(INDEX_BYTES),
        verified_config_path=root.repo_config,
    )
    assert receipt["reference_index_sha256"] == _sha(INDEX_BYTES)
    cohort = classify_cohort_quantification(destination, root.expectation(), verify_content=True)
    assert cohort["status"] == QUANT_STATUS_VERSION_DRIFT
    assert cohort["evidence"] == {"kind": EVIDENCE_RESTORED_WITNESS}
    assert cohort["observed_amalgkit_version"] != BUMPED_RUNTIME["REQUIRED_AMALGKIT_VERSION"]
    legacy = classify_quantification(destination, ACCESSION, verify_content=True)
    assert legacy["status"] == QUANT_STATUS_VERSION_DRIFT
    # Re-validation of the restored output also passes under the bumped runtime.
    _, files, rows = validate_quantification(destination, SPECIES, ACCESSION)
    assert rows == 8 and [path.name for path in files] == [
        QUANT_PROVENANCE_FILENAME,
        f"{ACCESSION}_abundance.tsv",
        f"{ACCESSION}_run_info.json",
    ]


def test_validation_requires_a_recorded_runtime(root: Root) -> None:
    sidecar = root.sample() / QUANT_PROVENANCE_FILENAME
    payload = json.loads(sidecar.read_text())
    payload.pop("amalgkit_version")
    sidecar.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="runtime is not recorded"):
        validate_quantification(root.sample(), SPECIES, ACCESSION)


def test_restored_witness_naming_another_index_is_invalid(root: Root) -> None:
    store = _locked_store(root)
    destination = root.base / "restored" / ACCESSION
    restore_quantification(
        store,
        COHORT,
        SPECIES,
        ACCESSION,
        destination,
        expected_reference_index_sha256=_sha(INDEX_BYTES),
        verified_config_path=root.repo_config,
    )
    assert classify_cohort_quantification(destination, root.expectation())["status"] == QUANT_STATUS_CURRENT
    result = classify_cohort_quantification(destination, root.expectation(reference_index_name="Other.idx"))
    assert result["status"] == QUANT_STATUS_INVALID
    assert result["reason"] == "restored reference manifest names a different index"


# ---------------------------------------------------------------------------
# Sidecar integrity and identity
# ---------------------------------------------------------------------------


def test_contract_that_is_not_self_consistent_is_invalid(root: Root) -> None:
    sidecar = root.sample() / QUANT_PROVENANCE_FILENAME
    payload = json.loads(sidecar.read_text())
    payload["quant_contract_id"] = "e" * 64
    sidecar.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    result = classify_cohort_quantification(root.sample(), root.expectation())
    assert result["status"] == QUANT_STATUS_INVALID
    assert result["reason"] == "quantification contract checksum mismatch"


def test_run_identity_mismatches_are_invalid(root: Root) -> None:
    result = classify_cohort_quantification(root.sample(), root.expectation("SRR999"))
    assert result["reason"] == "sample directory name differs from run accession"
    renamed = root.sample().with_name("SRR999")
    root.sample().rename(renamed)
    result = classify_cohort_quantification(renamed, root.expectation("SRR999"))
    assert result["status"] == QUANT_STATUS_INVALID
    assert result["reason"] == "provenance identity differs from the cohort task"


@pytest.mark.parametrize("damage", ["symlink", "escape", "truncated", "missing"])
def test_unsafe_or_partial_quantification_files_are_invalid(root: Root, damage: str) -> None:
    sample = root.sample()
    table = sample / f"{ACCESSION}_abundance.tsv"
    if damage == "symlink":
        outside = root.base / "outside.tsv"
        outside.write_bytes(table.read_bytes())
        table.unlink()
        table.symlink_to(outside)
        expected = "quantification file is not a regular file"
    elif damage == "escape":
        sidecar = sample / QUANT_PROVENANCE_FILENAME
        payload = json.loads(sidecar.read_text())
        payload["quantification_file"] = "../escape.tsv"
        payload["quant_contract_id"] = quantification_contract_id(payload)
        sidecar.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        expected = "quantification file escapes sample directory"
    elif damage == "truncated":
        table.write_bytes(b"target_id\ttpm\n")
        expected = "quantification file is truncated"
    else:
        table.unlink()
        expected = "quantification file is missing"
    result = classify_cohort_quantification(sample, root.expectation())
    assert result["status"] == QUANT_STATUS_INVALID
    assert result["reason"] == expected


def test_verify_content_rehashes_abundance_and_checks_run_counters(root: Root) -> None:
    table = root.sample() / f"{ACCESSION}_abundance.tsv"
    original = table.read_bytes()
    table.write_bytes(original.replace(b"1000.5", b"1000.6"))
    assert classify_cohort_quantification(root.sample(), root.expectation())["status"] == QUANT_STATUS_CURRENT
    result = classify_cohort_quantification(root.sample(), root.expectation(), verify_content=True)
    assert result["reason"] == "quantification checksum mismatch"

    table.write_bytes(original)
    info = root.sample() / f"{ACCESSION}_run_info.json"
    info.write_text(json.dumps({"n_processed": 10, "n_pseudoaligned": 11}))
    result = classify_cohort_quantification(root.sample(), root.expectation(), verify_content=True)
    assert result["reason"] == "run_info n_pseudoaligned exceeds n_processed"


def test_missing_sidecar_is_legacy_unverified(root: Root) -> None:
    (root.sample() / QUANT_PROVENANCE_FILENAME).unlink()
    assert classify_cohort_quantification(root.sample(), root.expectation())["status"] == QUANT_STATUS_LEGACY


def test_recorded_runtime_drift_is_reusable_in_cohort_mode(root: Root) -> None:
    sidecar = root.sample() / QUANT_PROVENANCE_FILENAME
    payload = json.loads(sidecar.read_text())
    payload.update(amalgkit_version="0.16.37", amalgkit_release_tag="v0.16.37", amalgkit_source_revision="old")
    sidecar.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    result = classify_cohort_quantification(root.sample(), root.expectation())
    assert result["status"] == QUANT_STATUS_VERSION_DRIFT
    assert result["observed_amalgkit_version"] == "0.16.37"


# ---------------------------------------------------------------------------
# Per-species precomputation and digest caching
# ---------------------------------------------------------------------------


def test_species_context_hashes_config_and_index_once(root: Root) -> None:
    for accession in ("SRR124", "SRR125", "SRR126"):
        root.write_sample(accession)
    context = build_species_reuse_context(
        SPECIES,
        root.repo_config,
        cohort_id=COHORT,
        reference_index_path=root.index,
    )
    assert context.reference_index_name == "Test_species.idx"
    root.repo_config.unlink()  # never read again: the context carries the content hash
    assert context.classify(root.sample())["status"] == QUANT_STATUS_CURRENT
    entries = len(context.index_digests.values)
    for accession in ("SRR124", "SRR125", "SRR126"):
        assert context.classify(root.sample(accession))["status"] == QUANT_STATUS_CURRENT
    assert len(context.index_digests.values) == entries  # manifest + 2 index files, hashed once

    root.alias.write_bytes(b"rebuilt index")  # a rebuilt index is noticed by the memo
    assert context.classify(root.sample("SRR124"))["status"] == QUANT_STATUS_INVALID

    with pytest.raises(ValueError, match="differs"):
        build_species_reuse_context(
            SPECIES, root.manifest, reference_index_path=root.alias, reference_index_sha256=_sha(INDEX_BYTES)
        )


def test_index_digest_cache_shares_hard_links_and_rehashes_edits(tmp_path: Path) -> None:
    first = tmp_path / "a.idx"
    first.write_bytes(INDEX_BYTES)
    os.link(first, tmp_path / "b.idx")
    cache = IndexDigestCache()
    assert cache.digest(first) == cache.digest(tmp_path / "b.idx") == _sha(INDEX_BYTES)
    assert len(cache.values) == 1
    stat = first.stat()
    first.write_bytes(INDEX_BYTES[:-1] + b"Z")
    os.utime(first, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert cache.digest(tmp_path / "b.idx") == _sha(INDEX_BYTES[:-1] + b"Z")
    assert cache.digest_or_none(tmp_path / "missing.idx") is None


def test_non_cohort_context_uses_content_bound_legacy_classification(root: Root) -> None:
    context = build_species_reuse_context(SPECIES, root.repo_config)
    assert not context.is_cohort
    root.repo_config.unlink()
    result = context.classify(root.sample())
    assert result["status"] == QUANT_STATUS_CURRENT
    assert result["evidence"]["config"] == "content"
    with pytest.raises(ValueError):
        context.expectation(ACCESSION)


# ---------------------------------------------------------------------------
# Sidecar writer: observed runtime, portable fields, durability
# ---------------------------------------------------------------------------


def test_sidecar_records_observed_runtime_portable_fields_and_is_fsynced(
    root: Root, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = json.loads((root.sample() / QUANT_PROVENANCE_FILENAME).read_text())
    assert payload["kallisto_version"] == "0.52.0"
    assert payload["kallisto_index_version"] == 13
    assert payload["config_name"] == "amalgkit_test_species.yaml"
    assert payload["reference_manifest_relpath"] == "../../reference/reference_aliases.json"
    observed = provenance.observed_amalgkit_runtime()
    assert {key: payload[key] for key in observed} == observed
    stripped = {k: v for k, v in payload.items() if k not in observed and not k.startswith(("kallisto", "config_n"))}
    stripped.pop("reference_manifest_relpath")
    assert quantification_contract_id(stripped) == payload["quant_contract_id"]  # new fields stay out of the contract

    synced: list[int] = []
    real_fsync = os.fsync

    def recording_fsync(fd: int) -> None:
        synced.append(fd)
        real_fsync(fd)

    monkeypatch.setattr(os, "fsync", recording_fsync)
    root.write_sample("SRR200")
    assert len(synced) >= 2  # the temp file before os.replace, then the directory entry


def test_bumped_pin_records_the_installed_runtime_as_drift(root: Root, monkeypatch: pytest.MonkeyPatch) -> None:
    _bump_runtime(monkeypatch)
    sample = root.write_sample("SRR300")
    payload = json.loads((sample / QUANT_PROVENANCE_FILENAME).read_text())
    assert payload["amalgkit_version"] != BUMPED_RUNTIME["REQUIRED_AMALGKIT_VERSION"]
    assert classify_cohort_quantification(sample, root.expectation("SRR300"))["status"] == QUANT_STATUS_VERSION_DRIFT


# ---------------------------------------------------------------------------
# Metadata currency is bound to content hashes (MI-05)
# ---------------------------------------------------------------------------


def _metadata_work(tmp_path: Path) -> tuple[Path, Path, Path]:
    work = tmp_path / "work"
    (work / "metadata").mkdir(parents=True)
    (work / "metadata" / "metadata.tsv").write_text("run\tscientific_name\nSRR1\tTest species\n")
    (work / "metadata" / "metadata_selected.tsv").write_text("run\tscientific_name\nSRR1\tTest species\n")
    config = tmp_path / "config.yaml"
    config.write_bytes(CONFIG_BYTES)
    rules = tmp_path / "select_rules.tsv"
    rules.write_text("rule\tvalue\n")
    write_metadata_provenance(work, species=SPECIES, config_path=config, selection_rules_path=rules)
    return work, config, rules


def test_unchanged_metadata_stays_current_after_an_amalgkit_bump(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    work, config, rules = _metadata_work(tmp_path)
    assert metadata_runtime_status(work) == QUANT_STATUS_CURRENT
    before = (work / "metadata" / ".metainformant_metadata_provenance.json").read_bytes()
    _bump_runtime(monkeypatch)
    assert is_current_metadata(work) is True
    assert metadata_runtime_status(work) == QUANT_STATUS_VERSION_DRIFT
    assert is_current_metadata(work, expected_config_sha256=_sha(CONFIG_BYTES)) is True
    assert is_current_metadata(work, expected_config_sha256="f" * 64) is False
    # Rewriting the sidecar over the same tables keeps the runtime that produced them.
    write_metadata_provenance(work, species=SPECIES, config_path=config, selection_rules_path=rules)
    assert (work / "metadata" / ".metainformant_metadata_provenance.json").read_bytes() == before
    (work / "metadata" / "metadata_selected.tsv").write_text("run\tscientific_name\nSRR2\tTest species\n")
    assert is_current_metadata(work) is False


def test_discovery_does_not_requery_ncbi_after_an_amalgkit_bump(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import subprocess

    import yaml

    from metainformant.rna.engine.streaming_orchestrator import (
        StreamingPipelineOrchestrator,
    )

    data_root = tmp_path / "data"
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    index_dir = tmp_path / "index"
    index_dir.mkdir()
    (index_dir / "Test_species.idx").write_bytes(b"kallisto-index")
    config_path = config_dir / "amalgkit_test_species.yaml"
    config_path.write_text(yaml.safe_dump({"species_list": ["Test_species"], "genome": {"index_dir": str(index_dir)}}))
    metadata_path = data_root / SPECIES / "work" / "metadata" / "metadata.tsv"
    calls: list[list[str]] = []

    def amalgkit_stub(command: list[str], **_: object) -> subprocess.CompletedProcess[str]:
        calls.append(list(command))
        metadata_path.parent.mkdir(parents=True, exist_ok=True)
        metadata_path.write_text("run\ttotal_bases\nSRR123\t1000\n")
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr("metainformant.rna.engine.streaming_orchestrator.subprocess.run", amalgkit_stub)
    orchestrator = StreamingPipelineOrchestrator(config_dir, tmp_path / "logs", db_path=tmp_path / "progress.db")
    checks = iter([False, True, True])
    monkeypatch.setattr(orchestrator, "verify_genome_index", lambda *_: next(checks))
    assert len(orchestrator.discover_species_tasks(config_path.name, max_gb=1.0, threads=1)) == 1
    first = [call[1] for call in calls if call[:1] == ["amalgkit"]]
    assert first == ["metadata", "select"]

    _bump_runtime(monkeypatch)
    calls.clear()
    assert len(orchestrator.discover_species_tasks(config_path.name, max_gb=1.0, threads=1)) == 1
    assert [call for call in calls if call[:1] == ["amalgkit"] and call[1] in {"metadata", "select"}] == []


# ---------------------------------------------------------------------------
# The resume path itself: reconcile and the quarantine audit
# ---------------------------------------------------------------------------


def test_resume_keeps_outputs_whose_config_checkout_was_deleted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from metainformant.rna.engine.streaming_orchestrator import (
        StreamingPipelineOrchestrator,
    )

    root = Root(tmp_path, accessions=(ACCESSION, "SRR124"))
    config_dir = tmp_path / "current-checkout"
    config_dir.mkdir()
    current = config_dir / root.repo_config.name
    current.write_bytes(CONFIG_BYTES)
    root.repo_config.unlink()
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(root.data))
    orchestrator = StreamingPipelineOrchestrator(config_dir, tmp_path / "logs", db_path=tmp_path / "progress.db")
    orchestrator.db.init_species(SPECIES, [ACCESSION, "SRR124"])
    quant_dir = root.data / SPECIES / "work" / "quant"

    # Without the content binding the old resume path would promote nothing...
    assert orchestrator.db.reconcile(SPECIES, quant_dir) == 0

    assert orchestrator._register_reuse_context(SPECIES, current, current) is not None
    assert orchestrator.db.reconcile(SPECIES, quant_dir, **orchestrator._reuse_kwargs(SPECIES)) == 2
    orchestrator._reset_non_current_quantified_states(SPECIES)
    assert orchestrator.db.get_state(SPECIES, ACCESSION) == "quantified"
    assert orchestrator.db.get_state(SPECIES, "SRR124") == "quantified"
    assert orchestrator.is_quantified(SPECIES, ACCESSION)

    # ...and an alias runtime YAML keeps the historical path-based check.
    assert orchestrator._register_reuse_context("other_species", current, tmp_path / "alias.yaml") is None


def test_previous_pin_outputs_resume_as_drift_without_requantification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Outputs of the previous pin (0.16.60, as on the canonical root) survive the upgrade to the current pin."""

    from metainformant.rna.engine.streaming_orchestrator import (
        StreamingPipelineOrchestrator,
    )

    assert provenance.pinned_amalgkit_runtime()["amalgkit_version"] != PREVIOUS_PIN_RUNTIME["amalgkit_version"]
    root = Root(tmp_path, accessions=(ACCESSION, "SRR124"))
    for accession in (ACCESSION, "SRR124"):
        _record_runtime(root.sample(accession), PREVIOUS_PIN_RUNTIME)
    config_dir = tmp_path / "current-checkout"
    config_dir.mkdir()
    current = config_dir / root.repo_config.name
    current.write_bytes(CONFIG_BYTES)
    root.repo_config.unlink()  # the recorded checkout is gone too, as for the canonical sidecars

    cohort = classify_cohort_quantification(root.sample(), root.expectation(), verify_content=True)
    assert cohort["status"] == QUANT_STATUS_VERSION_DRIFT
    assert cohort["observed_amalgkit_version"] == "0.16.60"

    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(root.data))
    monkeypatch.delenv("AMALGKIT_REQUANTIFICATION_POLICY", raising=False)
    orchestrator = StreamingPipelineOrchestrator(config_dir, tmp_path / "logs", db_path=tmp_path / "progress.db")
    assert orchestrator.requantification_policy == "preserve"
    orchestrator.db.init_species(SPECIES, [ACCESSION, "SRR124"])
    quant_dir = root.data / SPECIES / "work" / "quant"
    assert orchestrator._register_reuse_context(SPECIES, current, current) is not None
    # The producer's one-pass resume reconcile...
    summary = orchestrator._reconcile_legacy_species(SPECIES, current, quant_dir)
    assert summary["classified"] == 2
    assert sorted(orchestrator.db.get_samples(SPECIES, "quantified")) == [ACCESSION, "SRR124"]
    # ...and the public ProgressDB.reconcile both promote drift outputs.
    for accession in (ACCESSION, "SRR124"):
        orchestrator.db.set_state(SPECIES, accession, "pending")
    assert orchestrator.db.reconcile(SPECIES, quant_dir, **orchestrator._reuse_kwargs(SPECIES)) == 2
    orchestrator._reset_non_current_quantified_states(SPECIES)
    for accession in (ACCESSION, "SRR124"):
        assert orchestrator.db.get_state(SPECIES, accession) == "quantified"
        assert orchestrator.is_quantified(SPECIES, accession)
        recorded = json.loads((root.sample(accession) / QUANT_PROVENANCE_FILENAME).read_text())
        assert recorded["amalgkit_version"] == "0.16.60"  # the resume never rewrites the sidecar


# ---------------------------------------------------------------------------
# Real canonical outputs (copied fixtures; skipped when absent so CI stays hermetic)
# ---------------------------------------------------------------------------

REAL_FIXTURE_ENV = "RNA_REAL_QUANT_FIXTURE_DIR"
REPO_ROOT = Path(__file__).resolve().parents[2]
PROJECT_CONFIG_DIR = REPO_ROOT / "projects" / "hymenoptera_amalgkit" / "config" / "amalgkit"
# Frozen cohort hymenoptera27-20261005 (inventory.json), species apis_mellifera.
APIS_FROZEN_CONFIG_SHA256 = "98b763fd4b16a471b61892df9cd53cc9244ff100e3f20ea44702ef4ca40ef1d1"
APIS_FROZEN_INDEX_SHA256 = "d2e9e3acd83ca0d8e988b901c5cdcfa25e6d47efb351e794e01f708a049c1b17"
APIS_LEGACY_RUNS = ("DRR087067", "SRR12584109")
APIS_WITNESS_RUN = "SRR14374131"


@pytest.fixture
def real_apis_root() -> Path:
    """A copy of ``<data_root>/apis_mellifera/{work/quant,work/reference,genome/index}`` from the canonical root."""

    value = os.environ.get(REAL_FIXTURE_ENV)
    if not value:
        pytest.skip(f"{REAL_FIXTURE_ENV} is not set")
    root = Path(value)
    needed = [root / "apis_mellifera/work/quant" / run for run in (*APIS_LEGACY_RUNS, APIS_WITNESS_RUN)]
    needed += [root / "apis_mellifera/work/reference/reference_aliases.json", root / "apis_mellifera/genome/index"]
    if not all(path.exists() for path in needed):
        pytest.skip(f"{REAL_FIXTURE_ENV} does not hold the copied apis_mellifera fixture")
    if str(root.resolve()).startswith("/Volumes/"):
        pytest.fail("real-output fixtures must be copies, never the data volume itself")
    if not (PROJECT_CONFIG_DIR / "amalgkit_apis_mellifera.yaml").is_file():
        pytest.skip("hymenoptera_amalgkit project configuration is not checked out")
    return root


def _apis_context(config: Path) -> object:
    return build_species_reuse_context(
        "apis_mellifera",
        config,
        cohort_id="hymenoptera27-20261005",
        reference_index_sha256=APIS_FROZEN_INDEX_SHA256,
        reference_index_name="Apis_mellifera.idx",
    )


def test_real_legacy_outputs_are_reusable_under_the_current_repo_config(real_apis_root: Path) -> None:
    config = PROJECT_CONFIG_DIR / "amalgkit_apis_mellifera.yaml"
    context = _apis_context(config)
    assert context.config_sha256 == APIS_FROZEN_CONFIG_SHA256  # checked-in YAML bytes equal the frozen cohort's
    for run in APIS_LEGACY_RUNS:
        sample = real_apis_root / "apis_mellifera/work/quant" / run
        sidecar = json.loads((sample / QUANT_PROVENANCE_FILENAME).read_text())
        assert sidecar["config_path"].startswith("/Volumes/external_drive/Git/")  # the retired, deleted checkout
        assert not (sample / provenance.RESTORED_INPUTS_DIR).exists()
        assert {key: sidecar[key] for key in PREVIOUS_PIN_RUNTIME} == PREVIOUS_PIN_RUNTIME
        result = context.classify(sample, verify_content=True)
        # Produced by the previous pin: reusable as drift under the current one, never re-quantified.
        assert result["status"] == QUANT_STATUS_VERSION_DRIFT, result
        assert result["evidence"] == {
            "kind": EVIDENCE_REFERENCE_MANIFEST,
            "manifest": "relative",
            "index_paths": "rerooted",
        }
        legacy = classify_quantification(sample, run, expected_config_sha256=context.config_sha256)
        assert legacy["status"] == QUANT_STATUS_VERSION_DRIFT


def test_real_restored_witness_output_is_reusable(real_apis_root: Path) -> None:
    context = _apis_context(PROJECT_CONFIG_DIR / "amalgkit_apis_mellifera.yaml")
    sample = real_apis_root / "apis_mellifera/work/quant" / APIS_WITNESS_RUN
    result = context.classify(sample, verify_content=True)
    assert result["observed_amalgkit_version"] == PREVIOUS_PIN_RUNTIME["amalgkit_version"]
    assert result["status"] == QUANT_STATUS_VERSION_DRIFT, result
    assert result["evidence"] == {"kind": EVIDENCE_RESTORED_WITNESS}


def test_real_outputs_survive_a_bump_and_reject_other_config_bytes(
    real_apis_root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _bump_runtime(monkeypatch)
    context = _apis_context(PROJECT_CONFIG_DIR / "amalgkit_apis_mellifera.yaml")
    for run in (*APIS_LEGACY_RUNS, APIS_WITNESS_RUN):
        sample = real_apis_root / "apis_mellifera/work/quant" / run
        assert context.classify(sample)["status"] == QUANT_STATUS_VERSION_DRIFT
    other = _apis_context(PROJECT_CONFIG_DIR / "amalgkit_bombus_terrestris.yaml")
    sample = real_apis_root / "apis_mellifera/work/quant" / APIS_LEGACY_RUNS[0]
    assert other.classify(sample)["status"] == QUANT_STATUS_INVALID

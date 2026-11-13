"""Durable quant receipts exercise real files, hashes and conflicts."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest

from metainformant.rna.engine.durable_quant import lock_quantification, restore_quantification
from metainformant.rna.engine.provenance import (
    QUANT_PROVENANCE_FILENAME,
    sync_sample_outputs,
    write_quant_provenance,
)
from metainformant.rna.engine.quant_storage import DirectoryStore
from metainformant.rna.engine.streaming_orchestrator import (
    build_pipeline_resource_profile,
)


def sample(tmp_path: Path) -> tuple[Path, Path]:
    config = tmp_path / "config.yaml"
    config.write_text("species_list: [Test_species]\n")
    quant = tmp_path / "quant" / "SRR123"
    quant.mkdir(parents=True)
    (quant / "SRR123_abundance.tsv").write_text(
        "target_id\tlength\teff_length\test_counts\ttpm\ntranscript1\t100\t70\t3\t1000000\n"
    )
    (quant / "SRR123_run_info.json").write_text(json.dumps({"n_processed": 100, "n_pseudoaligned": 10}))
    write_quant_provenance(
        quant,
        species="test_species",
        run_accession="SRR123",
        config_path=config,
        command=["amalgkit", "quant", "--batch", "1"],
    )
    return quant, config


def test_explicit_in_flight_limit_is_not_raised() -> None:
    assert build_pipeline_resource_profile(48, 32, max_in_flight=12).max_in_flight == 12


def test_lock_roundtrip_and_idempotence(tmp_path: Path) -> None:
    quant, config = sample(tmp_path)
    store = DirectoryStore(tmp_path / "store")
    receipt = lock_quantification(
        store,
        "cohort",
        "test_species",
        "SRR123",
        quant,
        expected_config_sha256=hashlib.sha256(config.read_bytes()).hexdigest(),
    )
    assert lock_quantification(store, "cohort", "test_species", "SRR123", quant) == receipt
    restored = restore_quantification(store, "cohort", "test_species", "SRR123", tmp_path / "restored")
    assert restored == receipt
    assert (tmp_path / "restored" / "SRR123_abundance.tsv").read_bytes() == (
        quant / "SRR123_abundance.tsv"
    ).read_bytes()


@pytest.mark.parametrize("value", ["nan", "inf", "-1"])
def test_invalid_numeric_values_never_lock(tmp_path: Path, value: str) -> None:
    quant, config = sample(tmp_path)
    text = (quant / "SRR123_abundance.tsv").read_text().replace("\t3\t", f"\t{value}\t")
    (quant / "SRR123_abundance.tsv").write_text(text)
    write_quant_provenance(
        quant,
        species="test_species",
        run_accession="SRR123",
        config_path=config,
        command=["amalgkit", "quant", "--batch", "1"],
    )
    with pytest.raises(ValueError):
        lock_quantification(
            DirectoryStore(tmp_path / "store"),
            "cohort",
            "test_species",
            "SRR123",
            quant,
        )


def test_checksum_corruption_never_locks(tmp_path: Path) -> None:
    quant, _ = sample(tmp_path)
    (quant / "SRR123_abundance.tsv").write_text("corrupt")
    with pytest.raises(ValueError, match="checksum"):
        lock_quantification(
            DirectoryStore(tmp_path / "store"),
            "cohort",
            "test_species",
            "SRR123",
            quant,
        )


def test_existing_receipt_cannot_be_overwritten(tmp_path: Path) -> None:
    quant, config = sample(tmp_path)
    store = DirectoryStore(tmp_path / "store")
    first = lock_quantification(store, "cohort", "test_species", "SRR123", quant)
    (quant / "SRR123_abundance.tsv").write_text((quant / "SRR123_abundance.tsv").read_text().replace("\t3\t", "\t4\t"))
    write_quant_provenance(
        quant,
        species="test_species",
        run_accession="SRR123",
        config_path=config,
        command=["amalgkit", "quant", "--batch", "1"],
    )
    with pytest.raises(FileExistsError):
        lock_quantification(store, "cohort", "test_species", "SRR123", quant)
    assert json.loads(store.get("cohort/receipts/test_species/SRR123.json")) == first


def test_corrupt_stored_blob_refuses_restore(tmp_path: Path) -> None:
    quant, _ = sample(tmp_path)
    store = DirectoryStore(tmp_path / "store")
    receipt = lock_quantification(store, "cohort", "test_species", "SRR123", quant)
    blob = store.root / receipt["files"][0]["key"]
    blob.chmod(0o644)
    blob.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="checksum"):
        restore_quantification(store, "cohort", "test_species", "SRR123", tmp_path / "restored")
    assert not (tmp_path / "restored").exists()


def test_path_traversal_refused(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        DirectoryStore(tmp_path).put("../escape", b"x")


def reference_sample(tmp_path: Path) -> tuple[Path, Path, str]:
    quant, config = sample(tmp_path)
    index = tmp_path / "Test_species.idx"
    index.write_bytes(b"real reference fixture bytes")
    manifest = tmp_path / "reference.json"
    manifest.write_text(
        json.dumps(
            {
                "species": "test_species",
                "status": "complete",
                "kallisto_index": str(index),
            }
        )
    )
    write_quant_provenance(
        quant,
        species="test_species",
        run_accession="SRR123",
        config_path=config,
        command=["amalgkit", "quant"],
        reference_manifest_path=manifest,
    )
    return quant, config, hashlib.sha256(index.read_bytes()).hexdigest()


def test_missing_reference_cannot_receive_bound_receipt(tmp_path: Path) -> None:
    quant, _ = sample(tmp_path)
    store = DirectoryStore(tmp_path / "store")
    with pytest.raises(ValueError, match="lacks a reference"):
        lock_quantification(
            store,
            "cohort",
            "test_species",
            "SRR123",
            quant,
            expected_reference_index_sha256="a" * 64,
        )
    assert not list((tmp_path / "store").rglob("*.json"))


def test_wrong_actual_index_cannot_lock(tmp_path: Path) -> None:
    quant, _, _ = reference_sample(tmp_path)
    with pytest.raises(ValueError, match="differs from frozen"):
        lock_quantification(
            DirectoryStore(tmp_path / "store"),
            "cohort",
            "test_species",
            "SRR123",
            quant,
            expected_reference_index_sha256="a" * 64,
        )


def test_bound_lock_migrates_without_overwriting_recovery_receipt(
    tmp_path: Path,
) -> None:
    from metainformant.rna.engine.durable_quant import bound_receipt_key
    from metainformant.rna.engine.quant_storage import receipt_key

    quant, config, index_hash = reference_sample(tmp_path)
    config_hash = hashlib.sha256(config.read_bytes()).hexdigest()
    store = DirectoryStore(tmp_path / "store")
    lock_quantification(store, "cohort", "test_species", "SRR123", quant)
    original = store.get(receipt_key("cohort", "test_species", "SRR123"))
    with pytest.raises(FileNotFoundError):
        restore_quantification(
            store,
            "cohort",
            "test_species",
            "SRR123",
            tmp_path / "missing",
            expected_reference_index_sha256=index_hash,
        )
    assert not (tmp_path / "missing").exists()
    bound = lock_quantification(
        store,
        "cohort",
        "test_species",
        "SRR123",
        quant,
        expected_config_sha256=config_hash,
        expected_reference_index_sha256=index_hash,
    )
    assert bound["reference_index_sha256"] == index_hash
    assert store.get(receipt_key("cohort", "test_species", "SRR123")) == original
    assert store.get(bound_receipt_key("cohort", "test_species", "SRR123")) != original
    assert (
        restore_quantification(
            store,
            "cohort",
            "test_species",
            "SRR123",
            tmp_path / "fresh",
            expected_config_sha256=config_hash,
            expected_reference_index_sha256=index_hash,
        )
        == bound
    )
    for kwargs in (
        {"expected_reference_index_sha256": "a" * 64},
        {
            "expected_config_sha256": "a" * 64,
            "expected_reference_index_sha256": index_hash,
        },
    ):
        with pytest.raises(ValueError):
            restore_quantification(store, "cohort", "test_species", "SRR123", tmp_path / "bad", **kwargs)
        assert not (tmp_path / "bad").exists()


def test_portable_restore_keeps_original_sidecar_and_validates_local_inputs(
    tmp_path: Path,
) -> None:
    from metainformant.rna.engine.provenance import (
        QUANT_PROVENANCE_FILENAME,
        RESTORED_INPUTS_DIR,
        classify_quantification,
    )

    quant, config, index_hash = reference_sample(tmp_path)
    store = DirectoryStore(tmp_path / "store")
    receipt = lock_quantification(
        store,
        "cohort",
        "test_species",
        "SRR123",
        quant,
        expected_reference_index_sha256=index_hash,
    )
    original = (quant / QUANT_PROVENANCE_FILENAME).read_bytes()
    frozen = tmp_path / "frozen.yaml"
    frozen.write_bytes(config.read_bytes())
    config.write_text("different: config\n")
    assert classify_quantification(quant, "SRR123")["status"] == "invalid"
    destination = tmp_path / "restored_portable"
    kwargs = {
        "expected_config_sha256": receipt["config_sha256"],
        "expected_reference_index_sha256": index_hash,
        "verified_config_path": frozen,
    }
    restore_quantification(store, "cohort", "test_species", "SRR123", destination, **kwargs)
    assert (destination / QUANT_PROVENANCE_FILENAME).read_bytes() == original
    assert classify_quantification(destination, "SRR123", verify_content=True)["status"] == "current"
    restore_quantification(store, "cohort", "test_species", "SRR123", destination, **kwargs)
    (destination / RESTORED_INPUTS_DIR / "config.yaml").write_text("tampered\n")
    assert classify_quantification(destination, "SRR123")["status"] == "invalid"
    with pytest.raises(FileExistsError, match="witness"):
        restore_quantification(store, "cohort", "test_species", "SRR123", destination, **kwargs)


@pytest.mark.parametrize("corruption", ["contract", "reference", "escape"])
def test_portable_restore_witness_rejects_corruption(tmp_path: Path, corruption: str) -> None:
    from metainformant.rna.engine.provenance import (
        RESTORED_INPUTS_DIR,
        classify_quantification,
    )

    quant, config, index_hash = reference_sample(tmp_path)
    store = DirectoryStore(tmp_path / "store")
    lock_quantification(
        store,
        "cohort",
        "test_species",
        "SRR123",
        quant,
        expected_reference_index_sha256=index_hash,
    )
    destination = tmp_path / "restored_portable"
    restore_quantification(
        store,
        "cohort",
        "test_species",
        "SRR123",
        destination,
        expected_reference_index_sha256=index_hash,
        verified_config_path=config,
    )
    witness = destination / RESTORED_INPUTS_DIR
    if corruption == "contract":
        payload = json.loads((witness / "manifest.json").read_text())
        payload["contract_id"] = "a" * 64
        (witness / "manifest.json").write_text(json.dumps(payload))
    elif corruption == "reference":
        (witness / "reference_manifest.json").write_text("{}")
    else:
        outside = tmp_path / "outside.yaml"
        outside.write_bytes(config.read_bytes())
        (witness / "config.yaml").unlink()
        (witness / "config.yaml").symlink_to(outside)
    assert classify_quantification(destination, "SRR123")["status"] == "invalid"


@pytest.mark.parametrize("failure", ["unbound", "wrong_config"])
def test_portable_restore_rejects_unverified_inputs_before_writing(tmp_path: Path, failure: str) -> None:
    quant, config, index_hash = reference_sample(tmp_path)
    store = DirectoryStore(tmp_path / "store")
    kwargs = {} if failure == "unbound" else {"expected_reference_index_sha256": index_hash}
    lock_quantification(store, "cohort", "test_species", "SRR123", quant, **kwargs)
    if failure == "wrong_config":
        config.write_text("wrong: config\n")
    destination = tmp_path / "must_not_exist"
    with pytest.raises(ValueError):
        restore_quantification(
            store,
            "cohort",
            "test_species",
            "SRR123",
            destination,
            verified_config_path=config,
            **kwargs,
        )
    assert not destination.exists()


def test_restored_quant_cannot_checkpoint_under_different_configuration(
    tmp_path: Path,
) -> None:
    from metainformant.rna.engine.provenance import write_downstream_provenance

    quant, config, index_hash = reference_sample(tmp_path)
    store = DirectoryStore(tmp_path / "store")
    lock_quantification(
        store,
        "cohort",
        "test_species",
        "SRR123",
        quant,
        expected_reference_index_sha256=index_hash,
    )
    work = tmp_path / "restored/work"
    restore_quantification(
        store,
        "cohort",
        "test_species",
        "SRR123",
        work / "quant/SRR123",
        expected_reference_index_sha256=index_hash,
        verified_config_path=config,
    )
    metadata = work / "metadata"
    metadata.mkdir()
    for name in ("metadata.tsv", "metadata_selected.tsv"):
        (metadata / name).write_text("run\nSRR123\n")
    stage = work / "merge/Test_species"
    stage.mkdir(parents=True)
    (stage / "Test_species_est_counts.tsv").write_text("feature\tSRR123\nt1\t3\n")
    different = tmp_path / "different.yaml"
    different.write_text("changed: true\n")
    with pytest.raises(ValueError, match="differs from quantified"):
        write_downstream_provenance(
            work,
            species="test_species",
            config_path=different,
            quantified_samples=1,
            steps=["merge"],
            strict=True,
        )
    assert not (work / ".metainformant_downstream_provenance.json").exists()


@pytest.mark.skipif(os.geteuid() == 0, reason="root opens files without read permission")
def test_the_sidecar_is_published_only_after_every_output_was_flushed(tmp_path: Path) -> None:
    quant, config = sample(tmp_path)
    (quant / ".partial").write_text("dot files are not outputs")
    (quant / "subdir").mkdir()
    # The abundance table and the run info; dot files (the sidecar is one, written durably on its own) and
    # directories are skipped.
    assert sync_sample_outputs(quant) == 2

    sidecar = quant / QUANT_PROVENANCE_FILENAME
    sidecar.unlink()
    h5 = quant / "SRR123_abundance.h5"
    h5.write_bytes(b"\x89HDF")
    os.chmod(h5, 0)
    try:
        with pytest.raises(PermissionError):
            write_quant_provenance(
                quant,
                species="test_species",
                run_accession="SRR123",
                config_path=config,
                command=["amalgkit", "quant", "--batch", "1"],
            )
    finally:
        os.chmod(h5, 0o644)
    assert not sidecar.exists()  # never a sidecar over outputs that may not be on the disk

"""Durable quant receipts exercise real files, hashes and conflicts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from metainformant.rna.engine.durable_quant import (
    DirectoryStore,
    lock_quantification,
    restore_quantification,
)
from metainformant.rna.engine.provenance import write_quant_provenance
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
    (quant / "SRR123_run_info.json").write_text(
        json.dumps({"n_processed": 100, "n_pseudoaligned": 10})
    )
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
    assert (
        lock_quantification(store, "cohort", "test_species", "SRR123", quant) == receipt
    )
    restored = restore_quantification(
        store, "cohort", "test_species", "SRR123", tmp_path / "restored"
    )
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
    (quant / "SRR123_abundance.tsv").write_text(
        (quant / "SRR123_abundance.tsv").read_text().replace("\t3\t", "\t4\t")
    )
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
        restore_quantification(
            store, "cohort", "test_species", "SRR123", tmp_path / "restored"
        )
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
    from metainformant.rna.engine.durable_quant import bound_receipt_key, receipt_key
    from metainformant.rna.engine.aws_completion import verify_locked_campaign

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
            restore_quantification(
                store, "cohort", "test_species", "SRR123", tmp_path / "bad", **kwargs
            )
        assert not (tmp_path / "bad").exists()
    inventory = {
        "task_count": 1,
        "species": [
            {
                "species": "test_species",
                "config_sha256": config_hash,
                "index_sha256": index_hash,
                "tasks": [
                    {
                        "accession": "SRR123",
                        "task_id": "test_species/SRR123",
                        "reference_index_sha256": index_hash,
                    }
                ],
            }
        ],
    }
    assert verify_locked_campaign(inventory, store, "cohort", tmp_path / "complete")[
        "all_quant_locked"
    ]


def test_unbound_receipt_never_certifies_frozen_inventory(tmp_path: Path) -> None:
    from metainformant.rna.engine.aws_completion import verify_locked_campaign

    quant, config = sample(tmp_path)
    store = DirectoryStore(tmp_path / "store")
    lock_quantification(store, "cohort", "test_species", "SRR123", quant)
    inventory = {
        "task_count": 1,
        "species": [
            {
                "species": "test_species",
                "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
                "index_sha256": "a" * 64,
                "tasks": [
                    {
                        "accession": "SRR123",
                        "task_id": "test_species/SRR123",
                        "reference_index_sha256": "a" * 64,
                    }
                ],
            }
        ],
    }
    with pytest.raises(FileNotFoundError):
        verify_locked_campaign(inventory, store, "cohort", tmp_path / "complete")
    assert not (tmp_path / "complete" / "quant_completion_certificate.json").exists()

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

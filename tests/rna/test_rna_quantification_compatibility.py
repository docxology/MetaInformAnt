"""Contract-first reuse tests for quantification provenance."""

from __future__ import annotations

import json

from metainformant.rna.engine.provenance import (
    QUANT_STATUS_CURRENT,
    QUANT_STATUS_INVALID,
    QUANT_STATUS_LEGACY,
    QUANT_STATUS_VERSION_DRIFT,
    classify_quantification,
    is_current_quantification,
    is_reusable_quantification,
    write_quant_provenance,
)


def _sample(tmp_path):
    sample_dir = tmp_path / "quant" / "SRR123"
    sample_dir.mkdir(parents=True)
    abundance = sample_dir / "SRR123_abundance.tsv"
    abundance.write_text("target_id\ttpm\nTX1\t1\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text("species_list: [Test_species]\n", encoding="utf-8")
    write_quant_provenance(
        sample_dir,
        species="test_species",
        run_accession="SRR123",
        config_path=config,
        command=["amalgkit", "quant", "--threads", "1"],
        quantification_file=abundance,
    )
    return sample_dir


def test_current_and_version_drift_are_reusable(tmp_path):
    sample_dir = _sample(tmp_path)
    assert classify_quantification(sample_dir, "SRR123")["status"] == QUANT_STATUS_CURRENT
    assert is_current_quantification(sample_dir, "SRR123")
    assert is_reusable_quantification(sample_dir, "SRR123")

    sidecar = sample_dir / ".metainformant_quant_provenance.json"
    payload = json.loads(sidecar.read_text())
    payload["amalgkit_version"] = "0.16.37"
    payload["amalgkit_release_tag"] = "v0.16.37"
    payload["amalgkit_source_revision"] = "legacy-revision"
    sidecar.write_text(json.dumps(payload, indent=2) + "\n")

    assert classify_quantification(sample_dir, "SRR123")["status"] == QUANT_STATUS_VERSION_DRIFT
    assert not is_current_quantification(sample_dir, "SRR123")
    assert is_reusable_quantification(sample_dir, "SRR123")


def test_output_checksum_mismatch_is_invalid(tmp_path):
    sample_dir = _sample(tmp_path)
    (sample_dir / "SRR123_abundance.tsv").write_text("target_id\ttpm\nTX1\t2\n", encoding="utf-8")
    assert classify_quantification(sample_dir, "SRR123")["status"] == QUANT_STATUS_INVALID
    assert not is_reusable_quantification(sample_dir, "SRR123")


def test_missing_provenance_is_legacy_unverified(tmp_path):
    sample_dir = tmp_path / "quant" / "SRR123"
    sample_dir.mkdir(parents=True)
    (sample_dir / "SRR123_abundance.tsv").write_text("target_id\ttpm\nTX1\t1\n", encoding="utf-8")
    assert classify_quantification(sample_dir, "SRR123")["status"] == QUANT_STATUS_LEGACY
    assert not is_reusable_quantification(sample_dir, "SRR123")

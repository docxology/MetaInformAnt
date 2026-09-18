"""Tests for first-class evidence bundles (engine/evidence_bundle).

Deterministic fixtures only: small files under ``tmp_path``, no network,
no external tools. Validator behavior is fail-closed by construction.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metainformant.rna.engine import evidence_bundle as eb


def _write(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


@pytest.fixture()
def campaign_root(tmp_path: Path) -> Path:
    """A small deterministic campaign data root with one file per role."""

    root = tmp_path / "data"
    _write(root / "raw" / "SRR001_download_receipt.json", '{"run": "SRR001"}\n')
    _write(root / "work" / "quarantine_audit.tsv", "species\tsample\tclass\n")
    _write(root / "work" / "speciesA" / "quant_provenance.json", '{"schema": "v1"}\n')
    _write(root / "cross_species" / "species_manifest.tsv", "species\tsamples\nsA\t3\n")
    _write(root / "cross_species" / "results" / "expression_divergence_matrix.tsv", "sA\tsB\n0.1\t0.2\n")
    _write(root / "cross_species" / "results" / "divergence_inference.tsv", "pair\td_stat\nsA_sB\t0.15\n")
    return root


@pytest.fixture()
def built_bundle(tmp_path: Path, campaign_root: Path) -> tuple[Path, dict]:
    """Build a complete, valid bundle from the fixture campaign root."""

    config = _write(tmp_path / "config" / "amalgkit_current.yaml", "min_samples: 3\n")

    builder = eb.EvidenceBundleBuilder(data_root=campaign_root)
    builder.add_input(name="amalgkit_config", path=config)
    builder.add_software(name="amalgkit", version="0.4.1", source_revision="abc1234")
    builder.add_software(name="metainformant", version="1.0.0")
    builder.record_lock(
        name="campaign",
        holder="run_all_species",
        acquired_at_utc="2026-09-17T00:00:00Z",
        released_at_utc="2026-09-17T06:00:00Z",
        heartbeats=[{"at_utc": "2026-09-17T01:00:00Z", "running": 4}],
    )
    builder.record_task_counts(
        terminal={"quantified": 26, "failed": 1},
        unresolved={"pending": 0, "quarantined": 0},
    )
    builder.add_acquisition_receipt(step="download", path=campaign_root / "raw" / "SRR001_download_receipt.json")
    builder.add_recovery_receipt(step="quarantine", path=campaign_root / "work" / "quarantine_audit.tsv")
    builder.add_quantification_receipt(step="quant", path=campaign_root / "work" / "speciesA" / "quant_provenance.json")
    builder.add_descriptive_manifest(path=campaign_root / "cross_species" / "species_manifest.tsv")
    builder.add_descriptive_matrix(
        path=campaign_root / "cross_species" / "results" / "expression_divergence_matrix.tsv"
    )
    builder.add_inference_matrix(path=campaign_root / "cross_species" / "results" / "divergence_inference.tsv")
    builder.record_error("Quantification Failed for speciesA/SRR002", species="speciesA")
    builder.record_command(["amalgkit", "quant", "--config", "amalgkit_current.yaml"], exit_code=0)

    destination = tmp_path / "release" / eb.EVIDENCE_BUNDLE_FILENAME
    builder.build(destination)
    payload = json.loads(destination.read_text(encoding="utf-8"))
    return destination, payload


class TestBundleBuild:
    def test_bundle_is_valid_and_sectioned(self, built_bundle):
        _, payload = built_bundle
        report = eb.validate_bundle(payload)
        assert report.ok, report.failures
        assert report.schema_ok
        assert payload["schema"] == eb.EVIDENCE_BUNDLE_SCHEMA
        assert payload["bundle_id"]

        operational = {a["role"] for a in report.sections["operational"]}
        assert operational == {eb.ROLE_ACQUISITION, eb.ROLE_RECOVERY, eb.ROLE_QUANTIFICATION}
        assert {a["path"] for a in report.sections["descriptive"]} == {
            a["path"] for a in payload["artifacts"] if a["role"] == eb.ROLE_DESCRIPTIVE_ANALYSIS
        }
        assert [a["role"] for a in report.sections["inferential"]] == [eb.ROLE_BIOLOGICAL_INFERENCE]

    def test_task_count_partition_is_recorded(self, built_bundle):
        _, payload = built_bundle
        assert payload["task_counts"]["terminal"] == {"quantified": 26, "failed": 1}
        assert "pending" in payload["task_counts"]["unresolved"]

    def test_error_taxonomy_classifies_durable_failures(self, built_bundle):
        _, payload = built_bundle
        assert payload["errors"][0]["error_class"] == "quantification_failed"

    def test_bundle_id_is_content_deterministic(self, campaign_root):
        matrix = campaign_root / "cross_species" / "results" / "expression_divergence_matrix.tsv"

        first = eb.EvidenceBundleBuilder(data_root=campaign_root)
        first.add_descriptive_matrix(path=matrix)
        second = eb.EvidenceBundleBuilder(data_root=campaign_root)
        second.add_descriptive_matrix(path=matrix)

        assert first.payload()["bundle_id"] == second.payload()["bundle_id"]


class TestAtomicWrite:
    def test_build_leaves_no_staging_files(self, tmp_path, campaign_root):
        builder = eb.EvidenceBundleBuilder(data_root=campaign_root)
        builder.add_descriptive_matrix(
            path=campaign_root / "cross_species" / "results" / "expression_divergence_matrix.tsv"
        )
        destination = tmp_path / "out" / "bundle.json"
        builder.build(destination)
        # The rename is atomic: the destination holds the complete payload and
        # the staging directory holds no temporary leftovers.
        payload = json.loads(destination.read_text(encoding="utf-8"))
        assert payload["schema"] == eb.EVIDENCE_BUNDLE_SCHEMA
        assert list((tmp_path / "out").iterdir()) == [destination]

    def test_rebuild_replaces_destination_atomically(self, tmp_path, campaign_root):
        destination = tmp_path / "out" / "bundle.json"
        first = eb.EvidenceBundleBuilder(data_root=campaign_root)
        first.add_descriptive_matrix(
            path=campaign_root / "cross_species" / "results" / "expression_divergence_matrix.tsv"
        )
        first.build(destination)
        first_id = json.loads(destination.read_text(encoding="utf-8"))["bundle_id"]

        second = eb.EvidenceBundleBuilder(data_root=campaign_root)
        second.add_descriptive_manifest(path=campaign_root / "cross_species" / "species_manifest.tsv")
        second.build(destination)

        replaced = json.loads(destination.read_text(encoding="utf-8"))
        assert replaced["bundle_id"] != first_id
        assert list((tmp_path / "out").iterdir()) == [destination]

    def test_failed_build_writes_nothing(self, tmp_path, campaign_root):
        destination = tmp_path / "out" / "bundle.json"
        builder = eb.EvidenceBundleBuilder(data_root=campaign_root)
        with pytest.raises(eb.MixedRootError):
            builder.add_descriptive_matrix(path=tmp_path / "elsewhere" / "x.tsv")
        assert not destination.exists()
        assert not (tmp_path / "out").exists() or list((tmp_path / "out").iterdir()) == []


class TestNegativeControls:
    def test_stale_bundle_fails_closed(self, built_bundle):
        _, payload = built_bundle
        stale = next(
            a for a in payload["artifacts"] if a["role"] == eb.ROLE_DESCRIPTIVE_ANALYSIS and a["kind"] == "matrix"
        )
        Path(stale["path"]).write_text("sA\tsB\n9.9\t9.9\n", encoding="utf-8")

        report = eb.validate_bundle(payload)
        assert not report.ok
        assert any(f.startswith("stale:") for f in report.failures)
        with pytest.raises(eb.StaleEvidenceError):
            eb.require_valid_bundle(payload)

    def test_missing_artifact_fails_closed(self, built_bundle):
        _, payload = built_bundle
        target = next(a for a in payload["artifacts"] if a["role"] == eb.ROLE_ACQUISITION)
        Path(target["path"]).unlink()

        report = eb.validate_bundle(payload)
        assert not report.ok
        assert any(f.startswith("missing_artifact:") for f in report.failures)

    def test_mixed_root_bundle_fails_closed(self, built_bundle, tmp_path):
        _, payload = built_bundle
        # Tamper: relocate one artifact to a second root and re-sign the body
        # so only the root-mixing invariant is under test.
        other = _write(tmp_path / "other_root" / "divergence.tsv", "sA\tsB\n0.1\t0.2\n")
        for record in payload["artifacts"]:
            if record["role"] == eb.ROLE_BIOLOGICAL_INFERENCE:
                record["path"] = other.as_posix()
        body = {k: v for k, v in payload.items() if k != "bundle_id"}
        payload["bundle_id"] = eb._canonical_digest(body)

        report = eb.validate_bundle(payload)
        assert not report.ok
        assert any(f.startswith("mixed_root:") for f in report.failures)
        with pytest.raises(eb.MixedRootError):
            eb.require_valid_bundle(payload)

    def test_builder_refuses_second_root(self, campaign_root, tmp_path):
        builder = eb.EvidenceBundleBuilder(data_root=campaign_root)
        foreign = _write(tmp_path / "second_root" / "matrix.tsv", "x\n")
        with pytest.raises(eb.MixedRootError):
            builder.add_inference_matrix(path=foreign)

    def test_role_conflated_bundle_fails_closed(self, built_bundle):
        _, payload = built_bundle
        # Tamper: bind one descriptive matrix additionally as inference.
        descriptive = next(
            a for a in payload["artifacts"] if a["role"] == eb.ROLE_DESCRIPTIVE_ANALYSIS and a["kind"] == "matrix"
        )
        payload["artifacts"].append({**descriptive, "role": eb.ROLE_BIOLOGICAL_INFERENCE})
        body = {k: v for k, v in payload.items() if k != "bundle_id"}
        payload["bundle_id"] = eb._canonical_digest(body)

        report = eb.validate_bundle(payload)
        assert not report.ok
        assert any(f.startswith("role_conflation:") for f in report.failures)
        with pytest.raises(eb.RoleConflictError):
            eb.require_valid_bundle(payload)

    def test_builder_refuses_rebinding_a_path_to_another_role(self, campaign_root):
        matrix = campaign_root / "cross_species" / "results" / "expression_divergence_matrix.tsv"
        builder = eb.EvidenceBundleBuilder(data_root=campaign_root)
        builder.add_descriptive_matrix(path=matrix)
        with pytest.raises(eb.RoleConflictError):
            builder.add_inference_matrix(path=matrix)

    def test_state_conflation_rejected_by_builder_and_validator(self, campaign_root):
        builder = eb.EvidenceBundleBuilder(data_root=campaign_root)
        with pytest.raises(eb.EvidenceBundleError):
            builder.record_task_counts(terminal={"quantified": 1}, unresolved={"quantified": 2})
        with pytest.raises(eb.EvidenceBundleError):
            builder.record_task_counts(terminal={"flying": 1})

        payload = {
            "schema": eb.EVIDENCE_BUNDLE_SCHEMA,
            "data_root": campaign_root.as_posix(),
            "task_counts": {"terminal": {"quantified": 1}, "unresolved": {"quantified": 1}},
        }
        report = eb.validate_bundle(payload)
        assert not report.ok
        assert any(f.startswith("state_conflation:") for f in report.failures)

    def test_role_states_are_distinct_constants(self):
        assert len(eb.ALL_ROLES) == 5
        assert not (eb.OPERATIONAL_ROLES & eb.INFERENTIAL_ROLES)
        assert not (eb.DESCRIPTIVE_ROLES & eb.INFERENTIAL_ROLES)
        assert not (eb.OPERATIONAL_ROLES & eb.DESCRIPTIVE_ROLES)
        assert not (eb.TERMINAL_TASK_STATES & eb.UNRESOLVED_TASK_STATES)


class TestValidatorCLI:
    def test_cli_accepts_valid_bundle(self, built_bundle, capsys):
        destination, _ = built_bundle
        assert eb.main([str(destination)]) == 0
        assert "VALID" in capsys.readouterr().out

    def test_cli_rejects_stale_bundle(self, built_bundle, capsys):
        destination, payload = built_bundle
        stale = next(a for a in payload["artifacts"] if a["role"] == eb.ROLE_QUANTIFICATION)
        Path(stale["path"]).write_text('{"schema": "tampered"}\n', encoding="utf-8")
        assert eb.main([str(destination)]) == 1
        assert "INVALID" in capsys.readouterr().out

    def test_cli_rejects_non_bundle_file(self, tmp_path, capsys):
        not_a_bundle = tmp_path / "nope.json"
        not_a_bundle.write_text("{}", encoding="utf-8")
        assert eb.main([str(not_a_bundle)]) == 1

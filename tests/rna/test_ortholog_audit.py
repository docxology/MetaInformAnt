"""Deterministic tests for the ortholog retention audit and artifact manifest (MJ-02).

Covers the per-species retention audit layered on the transcript-level
orthogroup bridge (retention fractions plus explicit below-threshold flags)
and the versioned mapping-artifact manifest schema (schema version, copy
policy, sha256 checksums of consumed inputs, generation timestamp) together
with its fail-closed reader. No network access: every fixture is written to
``tmp_path``.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd
import pytest

from metainformant.rna.analysis.ortholog_mapping import (
    COPY_POLICIES,
    DEFAULT_COPY_POLICY,
    DEFAULT_MIN_RETENTION,
    MAPPING_ARTIFACT_SCHEMA_VERSION,
    MappingArtifactManifest,
    OrthogroupBridgeResult,
    OrthologBridgeError,
    OrthologySourceMetadata,
    audit_species_retention,
    build_orthogroup_bridge,
    read_mapping_artifact_manifest,
    write_mapping_artifact_manifest,
)

ORG_A = "7460_0"
ORG_B = "7461_0"
ORG_C = "7462_0"
SP_A = "Species_A"
SP_B = "Species_B"
SP_C = "Species_C"

TAXON_TO_SPECIES = {ORG_A: SP_A, ORG_B: SP_B}
TAXON_TO_SPECIES_WITH_C = {ORG_A: SP_A, ORG_B: SP_B, ORG_C: SP_C}

EXPRESSION_TIDS = {
    SP_A: {"XM_A1": "XM_A1.1_t1", "XM_A2": "XM_A2.1_t2"},
    SP_B: {"XM_B1": "XM_B1.1_u1"},
    SP_C: {"XM_C1": "XM_C1.1_v1"},
}

# OrthoDB gene id -> protein accession (base); genes :000008/:000009 of A are
# deliberately absent so the partial fixture exercises the unmapped chain.
ORTHODB_PROTEINS = {
    "7460_0:000001": "XP_A1",
    "7460_0:000002": "XP_A2",
    "7461_0:000004": "XP_B1",
}
PROT_TO_RNA = {"XP_A1": "XM_A1", "XP_A2": "XM_A2", "XP_B1": "XM_B1"}

AUDIT_COLUMNS = [
    "species",
    "input_genes",
    "mapped_to_protein",
    "mapped_to_rna",
    "mapped_to_transcript",
    "one_to_one",
    "one_to_many",
    "unmapped",
    "ogs_with_input",
    "ogs_retained",
    "ogs_unmapped",
    "transcript_retention_fraction",
    "orthogroup_retention_fraction",
    "below_threshold",
]

GENERATED_AT = "2026-09-17T12:00:00+00:00"


def _write_orthogroups(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    lines = ["Orthogroup\t" + ORG_A + "\t" + ORG_B]
    lines.extend(f"{og}\t{cell_a}\t{cell_b}" for og, cell_a, cell_b in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _full_fixture(tmp_path: Path) -> Path:
    """Two orthogroups; every input gene maps end-to-end in each species."""
    return _write_orthogroups(
        tmp_path / "orthogroups.tsv",
        [
            ("OG1", "7460_0:000001", "7461_0:000004"),
            ("OG2", "7460_0:000002", ""),
        ],
    )


def _partial_fixture(tmp_path: Path) -> Path:
    """Species A contributes three orthogroups but only OG1 maps; B maps OG1."""
    return _write_orthogroups(
        tmp_path / "orthogroups.tsv",
        [
            ("OG1", "7460_0:000001", "7461_0:000004"),
            ("OG2", "7460_0:000009", ""),
            ("OG3", "7460_0:000008", ""),
        ],
    )


def _bridge(
    og_path: Path,
    *,
    copy_policy: str = DEFAULT_COPY_POLICY,
    strict_duplicates: bool = False,
    duplicate_evidence: list[dict[str, str]] | None = None,
    taxon_to_species: dict[str, str] = TAXON_TO_SPECIES,
    expression_tids: dict[str, dict[str, str]] = EXPRESSION_TIDS,
):
    return build_orthogroup_bridge(
        og_path,
        ORTHODB_PROTEINS,
        PROT_TO_RNA,
        expression_tids,
        taxon_to_species,
        copy_policy=copy_policy,
        strict_duplicates=strict_duplicates,
        duplicate_evidence=duplicate_evidence,
    )


def _source_metadata(tmp_path: Path) -> OrthologySourceMetadata:
    """Versioned source metadata over two small on-disk inputs."""
    genes = tmp_path / "orthodb_genes.gz"
    with gzip.open(genes, "wt") as handle:
        handle.write("7460_0:000001\t7460_0\tXP_A1\n")
    groups = tmp_path / "orthogroups.tsv"
    groups.write_text("Orthogroup\nOG1\n", encoding="utf-8")
    return OrthologySourceMetadata.from_inputs(
        orthodb_release="v12.0",
        gene2refseq_url="https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz",
        retrieved_at="2026-09-17T00:00:00+00:00",
        inputs={"orthodb_genes": genes, "orthogroups": groups},
        taxonomic_scope=[ORG_A, ORG_B],
    )


def _manifest(tmp_path: Path) -> MappingArtifactManifest:
    return MappingArtifactManifest.create(
        source=_source_metadata(tmp_path),
        copy_policy="one-to-one-only",
        generated_at=GENERATED_AT,
    )


# =============================================================================
# Per-species retention audit: fractions and below-threshold flags
# =============================================================================


def test_full_retention_audit_fractions_and_flags(tmp_path: Path) -> None:
    """Fully retained species sit at fraction 1.0 with below_threshold False."""
    audit = audit_species_retention(_bridge(_full_fixture(tmp_path)))
    assert list(audit.columns) == AUDIT_COLUMNS
    by_species = audit.set_index("species")
    assert by_species.loc[SP_A, "ogs_with_input"] == 2
    assert by_species.loc[SP_A, "ogs_retained"] == 2
    assert by_species.loc[SP_A, "ogs_unmapped"] == 0
    assert by_species.loc[SP_A, "transcript_retention_fraction"] == pytest.approx(1.0)
    assert by_species.loc[SP_A, "orthogroup_retention_fraction"] == pytest.approx(1.0)
    assert not by_species.loc[SP_A, "below_threshold"]
    assert by_species.loc[SP_B, "orthogroup_retention_fraction"] == pytest.approx(1.0)
    assert not by_species.loc[SP_B, "below_threshold"]


def test_below_threshold_flag_tracks_orthogroup_fraction(tmp_path: Path) -> None:
    """A species losing two of three orthogroups is explicitly flagged."""
    audit = audit_species_retention(_bridge(_partial_fixture(tmp_path)))
    by_species = audit.set_index("species")

    row_a = by_species.loc[SP_A]
    assert row_a["ogs_with_input"] == 3
    assert row_a["ogs_retained"] == 1
    assert row_a["ogs_unmapped"] == 2
    assert row_a["input_genes"] == 3
    assert row_a["mapped_to_transcript"] == 1
    assert row_a["unmapped"] == 2
    assert row_a["transcript_retention_fraction"] == pytest.approx(1.0 / 3.0)
    assert row_a["orthogroup_retention_fraction"] == pytest.approx(1.0 / 3.0)
    assert bool(row_a["below_threshold"]) is True

    row_b = by_species.loc[SP_B]
    assert row_b["orthogroup_retention_fraction"] == pytest.approx(1.0)
    assert bool(row_b["below_threshold"]) is False


def test_zero_input_species_is_flagged(tmp_path: Path) -> None:
    """A configured species without orthogroup input has fraction 0.0 and is flagged."""
    result = _bridge(
        _full_fixture(tmp_path),
        taxon_to_species=TAXON_TO_SPECIES_WITH_C,
        expression_tids=EXPRESSION_TIDS,
    )
    audit = audit_species_retention(result)
    by_species = audit.set_index("species")
    row_c = by_species.loc[SP_C]
    assert row_c["input_genes"] == 0
    assert row_c["ogs_with_input"] == 0
    assert row_c["ogs_retained"] == 0
    assert row_c["transcript_retention_fraction"] == pytest.approx(0.0)
    assert row_c["orthogroup_retention_fraction"] == pytest.approx(0.0)
    assert bool(row_c["below_threshold"]) is True


@pytest.mark.parametrize("policy", COPY_POLICIES)
def test_retention_audit_invariant_across_copy_policies(tmp_path: Path, policy: str) -> None:
    """Retention accounting is identical under every copy policy from COPY_POLICIES."""
    baseline = audit_species_retention(_bridge(_partial_fixture(tmp_path), copy_policy=DEFAULT_COPY_POLICY))
    result = _bridge(_partial_fixture(tmp_path), copy_policy=policy)
    assert result.copy_policy == policy
    pd.testing.assert_frame_equal(audit_species_retention(result), baseline)


@pytest.mark.parametrize("min_retention,expected", [(0.0, False), (1.0, True)])
def test_min_retention_threshold_boundaries(tmp_path: Path, min_retention: float, expected: bool) -> None:
    """The flag is strict: fraction 1/3 crosses at 0.0 but never at 1.0; 1.0 never flags."""
    audit = audit_species_retention(_bridge(_partial_fixture(tmp_path)), min_retention=min_retention)
    by_species = audit.set_index("species")
    assert bool(by_species.loc[SP_A, "below_threshold"]) is expected
    assert bool(by_species.loc[SP_B, "below_threshold"]) is False
    # The default threshold sits between the boundaries.
    assert DEFAULT_MIN_RETENTION == 0.5


def test_audit_rejects_inconsistent_species_accounting(tmp_path: Path) -> None:
    """An audit and count table covering different species fail closed."""
    result = _bridge(_full_fixture(tmp_path))
    counts = result.orthogroup_counts
    broken = OrthogroupBridgeResult(
        table=result.table,
        retention_audit=result.retention_audit,
        duplicated_evidence=result.duplicated_evidence,
        dropped_orthogroups=result.dropped_orthogroups,
        copy_policy=result.copy_policy,
        orthogroup_counts=counts[counts["species"] != SP_B],
    )
    with pytest.raises(OrthologBridgeError, match="different species"):
        audit_species_retention(broken)


# =============================================================================
# Duplicated-evidence policy derived from COPY_POLICIES
# =============================================================================


@pytest.mark.parametrize("policy", COPY_POLICIES)
def test_strict_duplicates_fail_closed_under_every_policy(tmp_path: Path, policy: str) -> None:
    """strict_duplicates turns seeded evidence into a fail-closed error for all policies."""
    seed = [
        {
            "kind": "gene_multi_transcript",
            "species": SP_A,
            "orthogroup": "",
            "gene_id": "7460_0:000001",
            "transcript_id": "XM_A1.1_t1,XM_A2.1_t2",
            "detail": "seeded evidence",
        }
    ]
    with pytest.raises(OrthologBridgeError, match="duplicated mapping evidence"):
        _bridge(_partial_fixture(tmp_path), copy_policy=policy, strict_duplicates=True, duplicate_evidence=seed)


@pytest.mark.parametrize("policy", COPY_POLICIES)
def test_lenient_duplicates_are_recorded_under_every_policy(tmp_path: Path, policy: str) -> None:
    """Without strict_duplicates the seeded evidence is carried into the audit table."""
    seed = [
        {
            "kind": "gene_multi_transcript",
            "species": SP_A,
            "orthogroup": "",
            "gene_id": "7460_0:000001",
            "transcript_id": "XM_A1.1_t1,XM_A2.1_t2",
            "detail": "seeded evidence",
        }
    ]
    result = _bridge(_partial_fixture(tmp_path), copy_policy=policy, duplicate_evidence=seed)
    assert list(result.duplicated_evidence["kind"]) == ["gene_multi_transcript"]


# =============================================================================
# Versioned mapping-artifact manifest: writer and fail-closed validator
# =============================================================================


def test_manifest_round_trip(tmp_path: Path) -> None:
    """A persisted manifest round-trips with checksums and policy intact."""
    manifest = _manifest(tmp_path)
    assert manifest.schema_version == MAPPING_ARTIFACT_SCHEMA_VERSION == 1
    assert manifest.copy_policy == "one-to-one-only"
    assert manifest.generated_at == GENERATED_AT

    path = tmp_path / "nested" / "mapping_artifact_manifest.json"
    write_mapping_artifact_manifest(manifest, path)
    loaded = read_mapping_artifact_manifest(path)
    assert loaded == manifest

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert set(payload) == {"schema_version", "copy_policy", "generated_at", "source"}
    assert set(payload["source"]["input_sha256"]) == {"orthodb_genes", "orthogroups"}
    assert all(len(digest) == 64 for digest in payload["source"]["input_sha256"].values())


def test_manifest_default_policy_derived_from_constants(tmp_path: Path) -> None:
    """create() defaults to DEFAULT_COPY_POLICY, which must be a member of COPY_POLICIES."""
    manifest = MappingArtifactManifest.create(
        source=_source_metadata(tmp_path),
        generated_at=GENERATED_AT,
    )
    assert DEFAULT_COPY_POLICY in COPY_POLICIES
    assert manifest.copy_policy == DEFAULT_COPY_POLICY


def test_manifest_rejects_malformed_direct_construction(tmp_path: Path) -> None:
    """Direct construction fails closed on version, policy, timestamp, or source drift."""
    source = _source_metadata(tmp_path)
    with pytest.raises(OrthologBridgeError, match="schema version"):
        MappingArtifactManifest(
            schema_version=999,
            copy_policy=DEFAULT_COPY_POLICY,
            generated_at=GENERATED_AT,
            source=source,
        )
    with pytest.raises(OrthologBridgeError, match="schema_version"):
        MappingArtifactManifest(
            schema_version=True,
            copy_policy=DEFAULT_COPY_POLICY,
            generated_at=GENERATED_AT,
            source=source,
        )
    with pytest.raises(OrthologBridgeError, match="copy_policy"):
        MappingArtifactManifest(
            schema_version=MAPPING_ARTIFACT_SCHEMA_VERSION,
            copy_policy="pairwise",
            generated_at=GENERATED_AT,
            source=source,
        )
    with pytest.raises(OrthologBridgeError, match="timezone offset"):
        MappingArtifactManifest(
            schema_version=MAPPING_ARTIFACT_SCHEMA_VERSION,
            copy_policy=DEFAULT_COPY_POLICY,
            generated_at="2026-09-17T12:00:00",
            source=source,
        )
    with pytest.raises(OrthologBridgeError, match="OrthologySourceMetadata"):
        MappingArtifactManifest(
            schema_version=MAPPING_ARTIFACT_SCHEMA_VERSION,
            copy_policy=DEFAULT_COPY_POLICY,
            generated_at=GENERATED_AT,
            source="not-metadata",  # type: ignore[arg-type]
        )


def test_read_manifest_fail_closed(tmp_path: Path) -> None:
    """The reader refuses missing, malformed, or semantically invalid manifests."""
    with pytest.raises(FileNotFoundError):
        read_mapping_artifact_manifest(tmp_path / "absent.json")

    bad_json = tmp_path / "bad.json"
    bad_json.write_text("not json{", encoding="utf-8")
    with pytest.raises(OrthologBridgeError, match="not valid JSON"):
        read_mapping_artifact_manifest(bad_json)

    non_object = tmp_path / "list.json"
    non_object.write_text("[]", encoding="utf-8")
    with pytest.raises(OrthologBridgeError, match="must be a JSON object"):
        read_mapping_artifact_manifest(non_object)

    good_path = tmp_path / "manifest.json"
    write_mapping_artifact_manifest(_manifest(tmp_path), good_path)
    payload = json.loads(good_path.read_text(encoding="utf-8"))
    for name, (mutation, pattern) in {
        "version": ({"schema_version": 999}, "schema version"),
        "policy": ({"copy_policy": "pairwise"}, "copy_policy"),
        "naive": ({"generated_at": "2026-09-17T12:00:00"}, "timezone offset"),
        "garbage": ({"generated_at": "seventeenth of september"}, "ISO 8601"),
    }.items():
        mutated = dict(payload)
        mutated.update(mutation)
        path = tmp_path / f"mutated_{name}.json"
        path.write_text(json.dumps(mutated), encoding="utf-8")
        with pytest.raises(OrthologBridgeError, match=pattern):
            read_mapping_artifact_manifest(path)

    del payload["source"]
    path = tmp_path / "missing_source.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(OrthologBridgeError, match="missing required field"):
        read_mapping_artifact_manifest(path)

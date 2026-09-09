"""Deterministic fixture tests for the ortholog bridge methods (MJ-02).

Covers: copy policies, per-species retention audit, duplicated-evidence
checks, versioned source metadata with fail-closed checksums, and the
orthology presence table. No network access; all fixtures are tiny
OrthoDB-style TSVs written into ``tmp_path``.
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
    OrthogroupBridgeResult,
    OrthologBridgeError,
    OrthologySourceMetadata,
    build_orthogroup_bridge,
    build_transcript_orthogroup_table,
    load_gene2refseq_mapping,
    load_orthodb_proteins,
    orthology_presence_table,
    read_source_manifest,
    write_source_manifest,
)

TAXON_A = "7460"
TAXON_B = "7461"
ORG_A = f"{TAXON_A}_0"
ORG_B = f"{TAXON_B}_0"
SP_A = "Species_A"
SP_B = "Species_B"

TAXON_TO_SPECIES = {ORG_A: SP_A, ORG_B: SP_B}
EXPRESSION_TIDS = {
    SP_A: {"XM_001": "XM_001.1_t1", "XM_002": "XM_002.1_t2"},
    SP_B: {"XM_101": "XM_101.1_u1"},
}

# OrthoDB gene id -> protein accession (base), matching the genes fixture.
ORTHODB_PROTEINS = {
    "7460_0:000001": "XP_001",
    "7460_0:000002": "XP_002",
    "7461_0:000003": "XP_101",
}
PROT_TO_RNA = {"XP_001": "XM_001", "XP_002": "XM_002", "XP_101": "XM_101"}


def _write_genes(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    with gzip.open(path, "wt") as handle:
        for gene_id, org_id, protein in rows:
            handle.write(f"{gene_id}\t{org_id}\t{protein}\n")
    return path


def _write_gene2refseq(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    with gzip.open(path, "wt") as handle:
        handle.write("#tax_id\tGeneID\tsymbol\tRNA\tstatus\tprotein\n")
        for tax_id, rna, protein in rows:
            handle.write(f"{tax_id}\tg\t-\t{rna}\t-\t{protein}\n")
    return path


def _write_orthogroups(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    lines = ["Orthogroup\t" + ORG_A + "\t" + ORG_B]
    lines.extend(f"{og}\t{cell_a}\t{cell_b}" for og, cell_a, cell_b in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


@pytest.fixture()
def og_fixture(tmp_path: Path) -> Path:
    """Two orthogroups: OGX maps one transcript per species, OGY only in A."""
    return _write_orthogroups(
        tmp_path / "orthogroups.tsv",
        [
            ("OGX", "7460_0:000001", "7461_0:000003"),
            ("OGY", "7460_0:000002", ""),
        ],
    )


def test_default_policy_all_joined_and_wrappers_agree(tmp_path: Path, og_fixture: Path) -> None:
    """The default 'all-joined' policy preserves the historical join behaviour."""

    result = build_orthogroup_bridge(og_fixture, ORTHODB_PROTEINS, PROT_TO_RNA, EXPRESSION_TIDS, TAXON_TO_SPECIES)
    assert isinstance(result, OrthogroupBridgeResult)
    assert result.copy_policy == DEFAULT_COPY_POLICY

    table = build_transcript_orthogroup_table(
        og_fixture, ORTHODB_PROTEINS, PROT_TO_RNA, EXPRESSION_TIDS, TAXON_TO_SPECIES
    )
    pd.testing.assert_frame_equal(table, result.table)
    assert list(table.index) == ["OGX", "OGY"]
    assert list(table.columns) == [SP_A, SP_B]
    assert table.loc["OGX", SP_A] == "XM_001.1_t1"
    assert table.loc["OGX", SP_B] == "XM_101.1_u1"
    assert table.loc["OGY", SP_A] == "XM_002.1_t2"
    assert table.loc["OGY", SP_B] == ""


def test_one_to_one_only_drops_and_records_multi_transcript_orthogroups(tmp_path: Path) -> None:
    """Orthogroups with >1 transcript in any included species are dropped and recorded."""

    og_path = _write_orthogroups(
        tmp_path / "orthogroups.tsv",
        [
            ("OG_MULTI", "7460_0:000001,7460_0:000002", "7461_0:000003"),
            ("OG_SINGLE", "7460_0:000001", "7461_0:000003"),
        ],
    )
    result = build_orthogroup_bridge(
        og_path,
        ORTHODB_PROTEINS,
        PROT_TO_RNA,
        EXPRESSION_TIDS,
        TAXON_TO_SPECIES,
        copy_policy="one-to-one-only",
    )
    assert list(result.table.index) == ["OG_SINGLE"]
    assert result.dropped_orthogroups.to_dict("records") == [
        {"orthogroup": "OG_MULTI", "species": SP_A, "transcript_count": 2}
    ]


def test_first_transcript_policy_keeps_single_transcript(tmp_path: Path) -> None:
    """'first-transcript' keeps exactly one deterministic transcript per cell."""

    og_path = _write_orthogroups(
        tmp_path / "orthogroups.tsv",
        [("OGX", "7460_0:000001,7460_0:000002", "")],
    )
    result = build_orthogroup_bridge(
        og_path,
        ORTHODB_PROTEINS,
        PROT_TO_RNA,
        EXPRESSION_TIDS,
        TAXON_TO_SPECIES,
        copy_policy="first-transcript",
    )
    assert result.table.loc["OGX", SP_A] == "XM_001.1_t1"
    # The audit still reports the underlying one-to-many cell.
    audit = result.retention_audit.set_index("species")
    assert audit.loc[SP_A, "one_to_many"] == 1
    assert audit.loc[SP_A, "one_to_one"] == 0


def test_retention_audit_counts_per_species(og_fixture: Path) -> None:
    """The retention audit reports the full per-species mapping chain."""

    audit = build_orthogroup_bridge(
        og_fixture, ORTHODB_PROTEINS, PROT_TO_RNA, EXPRESSION_TIDS, TAXON_TO_SPECIES
    ).retention_audit
    assert list(audit.columns) == [
        "species",
        "input_genes",
        "mapped_to_protein",
        "mapped_to_rna",
        "mapped_to_transcript",
        "one_to_one",
        "one_to_many",
        "unmapped",
    ]
    by_species = audit.set_index("species")
    assert by_species.loc[SP_A].to_dict() == {
        "input_genes": 2,
        "mapped_to_protein": 2,
        "mapped_to_rna": 2,
        "mapped_to_transcript": 2,
        "one_to_one": 2,
        "one_to_many": 0,
        "unmapped": 0,
    }
    assert by_species.loc[SP_B].to_dict() == {
        "input_genes": 1,
        "mapped_to_protein": 1,
        "mapped_to_rna": 1,
        "mapped_to_transcript": 1,
        "one_to_one": 1,
        "one_to_many": 0,
        "unmapped": 0,
    }


def test_invalid_copy_policy_fails_closed(og_fixture: Path) -> None:
    with pytest.raises(OrthologBridgeError, match="copy_policy"):
        build_orthogroup_bridge(
            og_fixture,
            ORTHODB_PROTEINS,
            PROT_TO_RNA,
            EXPRESSION_TIDS,
            TAXON_TO_SPECIES,
            copy_policy="everything",  # type: ignore[arg-type]
        )


def test_transcript_multi_orthogroup_recorded_by_default(tmp_path: Path) -> None:
    """One transcript claimed by two orthogroups is recorded without failing."""

    og_path = _write_orthogroups(
        tmp_path / "orthogroups.tsv",
        [
            ("OG_1", "7460_0:000001", ""),
            ("OG_2", "7460_0:000001", ""),
        ],
    )
    result = build_orthogroup_bridge(og_path, ORTHODB_PROTEINS, PROT_TO_RNA, EXPRESSION_TIDS, TAXON_TO_SPECIES)
    records = result.duplicated_evidence.to_dict("records")
    assert len(records) == 1
    assert records[0]["kind"] == "transcript_multi_orthogroup"
    assert records[0]["transcript_id"] == "XM_001.1_t1"
    assert records[0]["orthogroup"] == "OG_1,OG_2"


def test_strict_duplicates_fails_closed(tmp_path: Path) -> None:
    """--strict-duplicates turns duplicated evidence into a fail-closed error."""

    og_path = _write_orthogroups(
        tmp_path / "orthogroups.tsv",
        [
            ("OG_1", "7460_0:000001", ""),
            ("OG_2", "7460_0:000001", ""),
        ],
    )
    with pytest.raises(OrthologBridgeError, match="duplicated mapping evidence"):
        build_orthogroup_bridge(
            og_path,
            ORTHODB_PROTEINS,
            PROT_TO_RNA,
            EXPRESSION_TIDS,
            TAXON_TO_SPECIES,
            strict_duplicates=True,
        )


def test_load_orthodb_proteins_records_gene_multiple_proteins(tmp_path: Path) -> None:
    """One OrthoDB gene claimed by two protein accessions is recorded, last kept."""

    genes_path = _write_genes(
        tmp_path / "odb12_genes.tab.gz",
        [
            ("7460_0:000001", ORG_A, "XP_001.1"),
            ("7460_0:000001", ORG_A, "XP_009.1"),
        ],
    )
    evidence: list[dict[str, str]] = []
    proteins = load_orthodb_proteins(genes_path, {ORG_A}, duplicate_evidence=evidence)
    assert proteins == {"7460_0:000001": "XP_009"}
    assert len(evidence) == 1
    assert evidence[0]["kind"] == "gene_multiple_proteins"
    assert evidence[0]["gene_id"] == "7460_0:000001"


def test_load_gene2refseq_records_protein_multiple_rnas(tmp_path: Path) -> None:
    """One protein claimed by two RNA accessions is recorded, last kept."""

    gene2refseq_path = _write_gene2refseq(
        tmp_path / "gene2refseq.gz",
        [
            (TAXON_A, "XM_001.1", "XP_001.1"),
            (TAXON_A, "XM_009.1", "XP_001.1"),
        ],
    )
    evidence: list[dict[str, str]] = []
    mapping = load_gene2refseq_mapping(gene2refseq_path, {TAXON_A}, duplicate_evidence=evidence)
    assert mapping == {"XP_001": "XM_009"}
    assert len(evidence) == 1
    assert evidence[0]["kind"] == "protein_multiple_rnas"
    assert evidence[0]["transcript_id"] == "XM_009"


def test_source_metadata_missing_checksum_fails_closed() -> None:
    """Any consumed input without a recorded checksum refuses to run."""

    metadata = OrthologySourceMetadata(
        orthodb_release="odb12",
        gene2refseq_url="https://example.org/gene2refseq.gz",
        retrieved_at="2026-01-01T00:00:00+00:00",
        taxonomic_scope=[TAXON_A],
        input_sha256={"gene2refseq": "0" * 64},
    )
    with pytest.raises(OrthologBridgeError, match="orthogroups"):
        metadata.require_checksums(["gene2refseq", "orthodb_genes", "orthogroups"])


def test_source_metadata_from_inputs_and_verify_detects_modification(tmp_path: Path) -> None:
    """Checksums are recorded from disk and mismatches fail closed."""

    consumed = tmp_path / "orthogroups.tsv"
    _write_orthogroups(consumed, [("OGX", "7460_0:000001", "")])
    metadata = OrthologySourceMetadata.from_inputs(
        orthodb_release="odb12",
        gene2refseq_url="https://example.org/gene2refseq.gz",
        retrieved_at="2026-01-01T00:00:00+00:00",
        inputs={"orthogroups": consumed},
        taxonomic_scope=[TAXON_A, TAXON_B],
    )
    metadata.verify_inputs({"orthogroups": consumed})

    _write_orthogroups(consumed, [("OGY", "7460_0:000002", "")])
    with pytest.raises(OrthologBridgeError, match="does not match its recorded sha256"):
        metadata.verify_inputs({"orthogroups": consumed})

    # A consumed input without a recorded checksum is refused before any
    # digest comparison happens.
    with pytest.raises(OrthologBridgeError, match="lacks a recorded sha256 checksum"):
        metadata.verify_inputs({"gene2refseq": tmp_path / "missing.gz"})


def test_source_metadata_invalid_checksum_format_fails_closed() -> None:
    with pytest.raises(OrthologBridgeError, match="sha256 hex digest"):
        OrthologySourceMetadata(
            orthodb_release="odb12",
            gene2refseq_url="https://example.org/gene2refseq.gz",
            retrieved_at="2026-01-01T00:00:00+00:00",
            taxonomic_scope=[TAXON_A],
            input_sha256={"gene2refseq": "deadbeef"},
        )


def test_source_manifest_round_trip(tmp_path: Path, og_fixture: Path) -> None:
    """A persisted source_manifest.json round-trips and validates checksums."""

    metadata = OrthologySourceMetadata.from_inputs(
        orthodb_release="odb12",
        gene2refseq_url="https://example.org/gene2refseq.gz",
        retrieved_at="2026-01-01T00:00:00+00:00",
        inputs={"orthogroups": og_fixture},
        taxonomic_scope=[TAXON_A, TAXON_B],
    )
    manifest_path = tmp_path / "source_manifest.json"
    write_source_manifest(metadata, manifest_path)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert set(payload["input_sha256"]) == {"orthogroups"}
    assert payload["taxonomic_scope"] == [TAXON_A, TAXON_B]

    restored = read_source_manifest(manifest_path)
    assert restored.to_dict() == metadata.to_dict()
    restored.verify_inputs({"orthogroups": og_fixture})

    broken = tmp_path / "broken.json"
    broken.write_text("{not json", encoding="utf-8")
    with pytest.raises(OrthologBridgeError, match="not valid JSON"):
        read_source_manifest(broken)

    with pytest.raises(FileNotFoundError, match="Source manifest not found"):
        read_source_manifest(tmp_path / "absent.json")


def test_orthology_presence_table(og_fixture: Path) -> None:
    """The presence table is the 0/1 matrix implied by the bridge table."""

    table = build_transcript_orthogroup_table(
        og_fixture, ORTHODB_PROTEINS, PROT_TO_RNA, EXPRESSION_TIDS, TAXON_TO_SPECIES
    )
    presence = orthology_presence_table(table)
    assert presence.loc["OGX", SP_A] == 1
    assert presence.loc["OGX", SP_B] == 1
    assert presence.loc["OGY", SP_A] == 1
    assert presence.loc["OGY", SP_B] == 0
    assert orthology_presence_table(pd.DataFrame()).empty


def test_copy_policy_constants_are_consistent() -> None:
    assert DEFAULT_COPY_POLICY in COPY_POLICIES
    assert set(COPY_POLICIES) == {"one-to-one-only", "first-transcript", "all-joined"}

"""Real XML fixtures and file/hash boundaries for archive source resolution."""

import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import pytest

from metainformant.rna.engine.source_resolution import (
    SCHEMA,
    SourceResolutionError,
    SourceTarget,
    load_source_resolutions,
    parse_ncbi_resolution,
)


def evidence(taxid: int = 7460, spots: int = 10, public: str = "true") -> bytes:
    return (
        "<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE>\n"
        f"<SAMPLE><SAMPLE_NAME><TAXON_ID>{taxid}</TAXON_ID></SAMPLE_NAME></SAMPLE>\n"
        "<EXPERIMENT><DESIGN><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>"
        "</LIBRARY_DESCRIPTOR></DESIGN></EXPERIMENT>\n"
        f'<RUN_SET><RUN accession="SRR123" is_public="{public}" load_done="true" '
        f'total_spots="{spots}" total_bases="1000">\n'
        '<SRAFiles><SRAFile url="https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR123/SRR123" '
        'size="1234" sratoolkit="1"/></SRAFiles>\n'
        "</RUN></RUN_SET></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>"
    ).encode()


def test_source_resolution_uses_authoritative_counts_and_conservative_reservation() -> None:
    # Given a public loaded NCBI RNA-Seq record.
    payload = evidence()
    # When parsed for a frozen target.
    (record,) = parse_ncbi_resolution(payload, [SourceTarget("SRR123", "apis_mellifera", 7460)])
    # Then identity, source hash, and modeled byte reservation are bound to those bytes.
    assert record.evidence_sha256 == hashlib.sha256(payload).hexdigest()
    assert record.total_spots == 10
    assert record.raw_bound_bytes == 1234 + 3 * 1000 + 512 * 10


@pytest.mark.parametrize(
    "payload",
    [
        evidence(taxid=1),
        evidence(spots=0),
        evidence(public="false"),
        evidence().replace(b"RNA-Seq", b"WGS"),
    ],
)
def test_source_resolution_rejects_identity_or_availability_mismatch(
    payload: bytes,
) -> None:
    # Given a wrong or unavailable source.
    with pytest.raises(SourceResolutionError):
        parse_ncbi_resolution(payload, [SourceTarget("SRR123", "apis_mellifera", 7460)])


def test_supplement_cannot_invent_sizes(tmp_path: Path) -> None:
    # Given XML evidence and a supplement with an invented larger spot count.
    payload = evidence()
    (tmp_path / "evidence.xml").write_bytes(payload)
    targets = [SourceTarget("SRR123", "apis_mellifera", 7460)]
    (record,) = parse_ncbi_resolution(payload, targets)
    row = asdict(record)
    row["total_spots"] = 100
    path = tmp_path / "source_resolutions.json"
    path.write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "inventory_sha256": "a" * 64,
                "evidence_file": "evidence.xml",
                "resolutions": [row],
            }
        )
    )
    # When loaded; then the JSON cannot supersede source evidence.
    with pytest.raises(SourceResolutionError, match="differ"):
        load_source_resolutions(path, "a" * 64, targets)


def test_supplement_binds_inventory_and_reparses_source(tmp_path: Path) -> None:
    # Given an exact, hash-bound resolution.
    payload = evidence()
    (tmp_path / "evidence.xml").write_bytes(payload)
    targets = [SourceTarget("SRR123", "apis_mellifera", 7460)]
    records = parse_ncbi_resolution(payload, targets)
    path = tmp_path / "source_resolutions.json"
    path.write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "inventory_sha256": "a" * 64,
                "evidence_file": "evidence.xml",
                "resolutions": [asdict(r) for r in records],
            }
        )
    )
    # When loaded; then it must equal the authoritative parser result.
    assert load_source_resolutions(path, "a" * 64, targets) == records


def test_job_bundle_overlays_resolved_counts_without_rewriting_frozen_metadata(
    tmp_path: Path,
) -> None:
    # Given frozen zero-count metadata and one hash-bound source resolution.
    import csv
    import io
    import tarfile

    from metainformant.rna.engine.aws_completion import _inputs_bundle

    work = tmp_path / "inputs" / "apis_mellifera" / "work"
    (work / "metadata").mkdir(parents=True)
    (work / "index").mkdir()
    metadata = work / "metadata" / "metadata_selected.tsv"
    metadata.write_text("run\ttotal_spots\ttotal_bases\nSRR123\t0\t0\n")
    index = work / "index" / "reference.idx"
    index.write_bytes(b"index-fixture")
    original = metadata.read_bytes()
    species = {
        "species": "apis_mellifera",
        "metadata_sha256": hashlib.sha256(original).hexdigest(),
        "index_name": index.name,
        "index_sha256": hashlib.sha256(index.read_bytes()).hexdigest(),
    }
    task = {
        "task_id": "apis_mellifera/SRR123",
        "accession": "SRR123",
        "species": "apis_mellifera",
        "total_spots": 10,
        "total_bases": 1000,
        "sra_bytes": 1234,
        "source_evidence_sha256": "a" * 64,
    }
    config = tmp_path / "species.yaml"
    config.write_text("species: apis_mellifera\n")
    species.update(config_sha256=hashlib.sha256(config.read_bytes()).hexdigest())
    # When a real job archive is generated.
    bundle, _ = _inputs_bundle(tmp_path, species, [task], tmp_path / "job", config_path=config)
    # Then frozen bytes stay identical and the archive's new metadata is independently hashed.
    assert metadata.read_bytes() == original
    with tarfile.open(bundle) as archive:
        handle = archive.extractfile("data/apis_mellifera/work/metadata/metadata_selected.tsv")
        assert handle is not None
        payload = handle.read()
        (row,) = list(csv.DictReader(io.StringIO(payload.decode()), delimiter="\t"))
        assert row["total_spots"] == "10"
        snapshot_handle = archive.extractfile("snapshot.json")
        assert snapshot_handle is not None
        snapshot = json.load(snapshot_handle)
        record = next(r for r in snapshot["input_files"] if r["path"].endswith("metadata_selected.tsv"))
        assert record["sha256"] == hashlib.sha256(payload).hexdigest()


@pytest.mark.parametrize(
    "doctype",
    [
        b"<!DOCTYPE EXPERIMENT_PACKAGE_SET>",
        b'<!DOCTYPE EXPERIMENT_PACKAGE_SET [<!ENTITY spots "10">]>',
    ],
)
def test_source_xml_rejects_dtd_and_entity_expansion(doctype: bytes) -> None:
    payload = doctype + evidence()
    if b"ENTITY" in doctype:
        payload = payload.replace(b'total_spots="10"', b'total_spots="&spots;"')
    with pytest.raises(SourceResolutionError, match="unsafe or malformed XML"):
        parse_ncbi_resolution(payload, [SourceTarget("SRR123", "apis_mellifera", 7460)])


def test_source_xml_rejects_external_entity_and_malformed_input(tmp_path: Path) -> None:
    external = tmp_path / "private.xml"
    external.write_text("private sentinel")
    payload = (
        f'<!DOCTYPE EXPERIMENT_PACKAGE_SET [<!ENTITY source SYSTEM "{external.as_uri()}">]>'.encode()
        + evidence().replace(b'total_spots="10"', b'total_spots="&source;"')
    )
    for invalid in (payload, b"<EXPERIMENT_PACKAGE_SET>"):
        with pytest.raises(SourceResolutionError, match="unsafe or malformed XML"):
            parse_ncbi_resolution(invalid, [SourceTarget("SRR123", "apis_mellifera", 7460)])
    assert external.read_text() == "private sentinel"

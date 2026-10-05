"""Hash-bound NCBI source supplements for incomplete ENA run records."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Final, Sequence
from urllib.parse import urlparse
from defusedxml import ElementTree as ET
from defusedxml.common import DefusedXmlException

SCHEMA: Final = "metainformant.rna.source_resolution.v1"


class SourceResolutionError(ValueError):
    """An external record cannot safely resolve a frozen target."""

    def __init__(self, accession: str, reason: str) -> None:
        self.accession = accession
        self.reason = reason
        super().__init__(f"{accession}: {reason}")


@dataclass(frozen=True, slots=True)
class SourceTarget:
    """Identity of one already-selected run; resolution cannot expand scope."""

    accession: str
    species: str
    taxid: int


@dataclass(frozen=True, slots=True)
class RunResolution:
    """Public SRA evidence and conservative modeled raw-file reservation."""

    accession: str
    species: str
    taxid: int
    total_spots: int
    total_bases: int
    sra_bytes: int
    source_url: str
    evidence_sha256: str

    def __post_init__(self) -> None:
        for value in (self.taxid, self.total_spots, self.total_bases, self.sra_bytes):
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise SourceResolutionError(
                    self.accession, "counts and sizes must be positive integers"
                )
        parsed = urlparse(self.source_url)
        if (
            (parsed.scheme, parsed.hostname, parsed.path)
            != (
                "https",
                "sra-pub-run-odp.s3.amazonaws.com",
                f"/sra/{self.accession}/{self.accession}",
            )
            or parsed.query
            or parsed.fragment
            or parsed.username
            or parsed.port
        ):
            raise SourceResolutionError(
                self.accession, "unexpected public SRA source URL"
            )
        if len(self.evidence_sha256) != 64 or any(
            c not in "0123456789abcdef" for c in self.evidence_sha256
        ):
            raise SourceResolutionError(self.accession, "invalid evidence hash")

    @property
    def raw_bound_bytes(self) -> int:
        """Reserve SRA plus modeled FASTQ sequence/quality/header overhead.

        This is a conservative scheduling model, not a measured FASTQ size.
        Worker free-space floors and the hard job deadline remain authoritative.
        """
        return self.sra_bytes + 3 * self.total_bases + 512 * self.total_spots


def parse_ncbi_resolution(
    payload: bytes, targets: Sequence[SourceTarget]
) -> tuple[RunResolution, ...]:
    """Require matching taxonomy, RNA-Seq, public loaded runs, and primary SRA files."""
    expected = {t.accession: t for t in targets}
    if not expected or len(expected) != len(targets):
        raise SourceResolutionError("inventory", "targets must be nonempty and unique")
    evidence_hash = hashlib.sha256(payload).hexdigest()
    try:
        root = ET.fromstring(
            payload, forbid_dtd=True, forbid_entities=True, forbid_external=True
        )
    except (ET.ParseError, DefusedXmlException) as exc:
        raise SourceResolutionError("evidence", "unsafe or malformed XML") from exc
    found: dict[str, RunResolution] = {}
    for package in root.findall(".//EXPERIMENT_PACKAGE"):
        taxid = package.findtext("SAMPLE/SAMPLE_NAME/TAXON_ID")
        strategy = package.findtext(
            "EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY"
        )
        for run in package.findall("RUN_SET/RUN"):
            accession = run.get("accession", "")
            if accession not in expected:
                continue
            target = expected[accession]
            if accession in found:
                raise SourceResolutionError(accession, "duplicate run evidence")
            if taxid != str(target.taxid) or strategy != "RNA-Seq":
                raise SourceResolutionError(
                    accession,
                    "taxonomy or library strategy differs from the frozen target",
                )
            if run.get("is_public") != "true" or run.get("load_done") != "true":
                raise SourceResolutionError(accession, "run is not publicly loaded")
            source_url = (
                f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{accession}/{accession}"
            )
            files = [
                f
                for f in run.findall(".//SRAFile")
                if f.get("url") == source_url and f.get("sratoolkit") == "1"
            ]
            if len(files) != 1:
                raise SourceResolutionError(
                    accession, "no unique primary public SRA file"
                )
            try:
                resolution = RunResolution(
                    accession,
                    target.species,
                    target.taxid,
                    int(run.get("total_spots", "0")),
                    int(run.get("total_bases", "0")),
                    int(files[0].get("size", "0")),
                    source_url,
                    evidence_hash,
                )
            except (TypeError, ValueError) as exc:
                raise SourceResolutionError(
                    accession, "invalid quantitative SRA evidence"
                ) from exc
            found[accession] = resolution
    if set(found) != set(expected):
        raise SourceResolutionError(
            "inventory", f"missing run evidence: {sorted(set(expected) - set(found))}"
        )
    return tuple(found[a] for a in sorted(found))


def load_source_resolutions(
    path: Path, inventory_sha256: str, targets: Sequence[SourceTarget]
) -> tuple[RunResolution, ...]:
    """Validate an optional supplement against frozen inventory and XML bytes."""
    if not path.exists():
        return ()
    document = json.loads(path.read_text())
    if document["schema"] != SCHEMA or document["inventory_sha256"] != inventory_sha256:
        raise SourceResolutionError(
            "inventory", "supplement does not bind the frozen inventory"
        )
    evidence_path = path.parent / document["evidence_file"]
    if (
        evidence_path.parent.resolve() != path.parent.resolve()
        or not evidence_path.is_file()
    ):
        raise SourceResolutionError(
            "inventory", "unsafe or missing local evidence file"
        )
    records = tuple(RunResolution(**r) for r in document["resolutions"])
    expected = {t.accession: t for t in targets}
    selected = []
    for record in records:
        target = expected.get(record.accession)
        if target is None or (record.species, record.taxid) != (
            target.species,
            target.taxid,
        ):
            raise SourceResolutionError(
                record.accession, "resolution lies outside the frozen scope"
            )
        selected.append(target)
    # Re-parse the authoritative bytes; JSON cannot invent a size, count, or URL.
    parsed = parse_ncbi_resolution(evidence_path.read_bytes(), selected)
    if records != parsed:
        raise SourceResolutionError(
            "inventory", "resolution records differ from the source evidence"
        )
    return records

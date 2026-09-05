"""Ortholog mapping integration module.

Provides utilities to bridge OrthoDB gene IDs, NCBI protein accessions,
and NCBI RNA accessions to construct transcript-level orthogroup tables.
"""

from __future__ import annotations
import gzip
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Literal, Mapping, Optional, Sequence, Set

import pandas as pd

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def load_gene2refseq_mapping(
    gene2refseq_path: Path,
    target_taxon_ids: Set[str],
    duplicate_evidence: Optional[List[Dict[str, str]]] = None,
) -> Dict[str, str]:
    """Load protein->RNA mappings from NCBI gene2refseq for target taxa.

    Args:
        gene2refseq_path: Path to gene2refseq file (gzipped)
        target_taxon_ids: Set of string taxon IDs to filter for
        duplicate_evidence: Optional list collecting duplicated-evidence records;
            when a protein base maps to several RNA bases, the collision is
            appended here and the last mapping is kept (historical behaviour).

    Returns:
        Dictionary mapping protein_accession_base -> rna_accession_base.
        Base means without version (e.g. XP_001120086 -> XM_001120086).
    """
    prot_to_rna: Dict[str, str] = {}
    logger.info(f"Scanning {gene2refseq_path.name} for {len(target_taxon_ids)} taxa...")

    with gzip.open(gene2refseq_path, "rt") as f:
        f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue

            tax_id = parts[0]
            if tax_id not in target_taxon_ids:
                continue

            rna_acc = parts[3]  # RNA accession with version
            prot_acc = parts[5]  # Protein accession with version

            if rna_acc == "-" or prot_acc == "-":
                continue

            # Strip version numbers for matching
            prot_base = prot_acc.rsplit(".", 1)[0] if "." in prot_acc else prot_acc
            rna_base = rna_acc.rsplit(".", 1)[0] if "." in rna_acc else rna_acc

            existing_rna = prot_to_rna.get(prot_base)
            if existing_rna is not None and existing_rna != rna_base and duplicate_evidence is not None:
                duplicate_evidence.append(
                    {
                        "kind": "protein_multiple_rnas",
                        "species": tax_id,
                        "orthogroup": "",
                        "gene_id": "",
                        "transcript_id": rna_base,
                        "detail": (
                            f"protein {prot_base} maps to RNA {existing_rna} and {rna_base} in taxon "
                            f"{tax_id}; the bridge keeps {rna_base}"
                        ),
                    }
                )
            prot_to_rna[prot_base] = rna_base

    logger.info(f"Loaded {len(prot_to_rna)} protein->RNA mappings")
    return prot_to_rna


def load_orthodb_proteins(
    genes_path: Path,
    target_org_prefixes: Set[str],
    duplicate_evidence: Optional[List[Dict[str, str]]] = None,
) -> Dict[str, str]:
    """Load OrthoDB gene_id -> protein_accession (base, no version).

    Args:
        genes_path: Path to OrthoDB genes file (gzipped)
        target_org_prefixes: Set of organism prefixes (e.g., taxonomy IDs mapped to OrthoDB format)
        duplicate_evidence: Optional list collecting duplicated-evidence records;
            when one OrthoDB gene ID maps to several protein accessions, the
            collision is appended here and the last mapping is kept (historical
            behaviour).

    Returns:
        Dictionary mapping OrthoDB gene_id -> protein_accession_base.
    """
    result: Dict[str, str] = {}
    logger.info(f"Scanning {genes_path.name}...")
    with gzip.open(genes_path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            gene_id = parts[0]
            org_id = parts[1]
            protein_acc = parts[2]
            if org_id in target_org_prefixes:
                # Strip version
                prot_base = protein_acc.rsplit(".", 1)[0] if "." in protein_acc else protein_acc
                existing_protein = result.get(gene_id)
                if existing_protein is not None and existing_protein != prot_base and duplicate_evidence is not None:
                    duplicate_evidence.append(
                        {
                            "kind": "gene_multiple_proteins",
                            "species": org_id,
                            "orthogroup": "",
                            "gene_id": gene_id,
                            "transcript_id": "",
                            "detail": (
                                f"OrthoDB gene {gene_id} maps to proteins {existing_protein} and "
                                f"{prot_base} in {org_id}; the bridge keeps {prot_base}"
                            ),
                        }
                    )
                result[gene_id] = prot_base

    logger.info(f"Loaded {len(result)} OrthoDB gene->protein mappings")
    return result


def _resolve_manifest_matrix_path(row: pd.Series) -> Path:
    """Resolve a manifest matrix path using its optional data-root provenance."""

    path = Path(str(row["finalize_path"])).expanduser()
    if path.is_absolute():
        return path
    data_root = row.get("data_root")
    if data_root is not None and not pd.isna(data_root):
        return Path(str(data_root)).expanduser() / path
    return path


def load_expression_transcript_ids(manifest: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    """Load transcript IDs and build RNA_accession -> full_transcript_id maps.

    Args:
        manifest: DataFrame with columns 'species_title' and 'finalize_path'

    Returns:
        Dictionary mapping species title -> (RNA_accession_base -> full_transcript_id)
    """
    result = {}
    for _, row in manifest.iterrows():
        sp_title = row["species_title"]
        finalize_path = _resolve_manifest_matrix_path(row)
        if not finalize_path.exists():
            continue
        try:
            tid_map = {}
            with open(finalize_path) as f:
                f.readline()
                for line in f:
                    tid = line.split("\t", 1)[0]
                    # Extract RNA accession (XM_, XR_, NM_) without version
                    m = re.search(r"((?:XM|XR|NM)_\d+)", tid)
                    if m:
                        tid_map[m.group(1)] = tid
            result[sp_title] = tid_map
            logger.info(f"{sp_title}: loaded {len(tid_map)} transcript IDs")
        except Exception as e:
            logger.error(f"Failed to load transcript IDs for {sp_title}: {e}")
    return result


# =============================================================================
# Versioned sources, copy policies, and bridge artifacts (MJ-02)
# =============================================================================

CopyPolicy = Literal["one-to-one-only", "first-transcript", "all-joined"]
COPY_POLICIES = ("one-to-one-only", "first-transcript", "all-joined")
DEFAULT_COPY_POLICY: CopyPolicy = "all-joined"

_DUPLICATED_EVIDENCE_COLUMNS = [
    "kind",
    "species",
    "orthogroup",
    "gene_id",
    "transcript_id",
    "detail",
]
_DROPPED_ORTHOGROUP_COLUMNS = ["orthogroup", "species", "transcript_count"]
_RETENTION_AUDIT_COLUMNS = [
    "species",
    "input_genes",
    "mapped_to_protein",
    "mapped_to_rna",
    "mapped_to_transcript",
    "one_to_one",
    "one_to_many",
    "unmapped",
]


class OrthologBridgeError(ValueError):
    """Fail-closed error for bridge inputs, policies, manifests, or duplicates."""


def sha256_file(path: Path) -> str:
    """Compute the lowercase sha256 hex digest of a file in bounded chunks."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class OrthologySourceMetadata:
    """Versioned source record for every input the ortholog bridge consumes.

    The record fails closed: constructing it requires non-empty release and
    provenance strings, a non-empty taxonomic scope, and at least one recorded
    sha256 checksum, and :meth:`require_checksums` refuses any consumed input
    without a recorded checksum before analysis outputs are written.
    """

    orthodb_release: str
    gene2refseq_url: str
    retrieved_at: str
    taxonomic_scope: List[str]
    input_sha256: Dict[str, str]

    def __post_init__(self) -> None:
        for field_name in ("orthodb_release", "gene2refseq_url", "retrieved_at"):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value.strip():
                raise OrthologBridgeError(
                    f"orthology source metadata field '{field_name}' must be a non-empty string"
                )
        if not isinstance(self.taxonomic_scope, (list, tuple)) or not self.taxonomic_scope:
            raise OrthologBridgeError(
                "orthology source metadata requires a non-empty taxonomic_scope list"
            )
        if not isinstance(self.input_sha256, Mapping) or not self.input_sha256:
            raise OrthologBridgeError(
                "orthology source metadata requires at least one recorded input checksum"
            )
        for name, checksum in self.input_sha256.items():
            if not isinstance(name, str) or not name.strip():
                raise OrthologBridgeError("input checksum names must be non-empty strings")
            if (
                not isinstance(checksum, str)
                or len(checksum) != 64
                or any(character not in "0123456789abcdef" for character in checksum)
            ):
                raise OrthologBridgeError(
                    f"recorded checksum for input '{name}' is not a lowercase sha256 hex digest"
                )

    @classmethod
    def from_inputs(
        cls,
        *,
        orthodb_release: str,
        gene2refseq_url: str,
        retrieved_at: str,
        inputs: Mapping[str, Path],
        taxonomic_scope: Sequence[str],
    ) -> OrthologySourceMetadata:
        """Record metadata, computing the sha256 checksum of every input file."""
        checksums: Dict[str, str] = {}
        for name in sorted(inputs):
            path = Path(inputs[name])
            if not path.is_file():
                raise OrthologBridgeError(
                    f"cannot record a checksum for input '{name}': file not found: {path}"
                )
            checksums[name] = sha256_file(path)
        return cls(
            orthodb_release=orthodb_release,
            gene2refseq_url=gene2refseq_url,
            retrieved_at=retrieved_at,
            taxonomic_scope=list(taxonomic_scope),
            input_sha256=checksums,
        )

    def require_checksums(self, consumed_inputs: Iterable[str]) -> Dict[str, str]:
        """Fail closed unless every consumed input has a recorded checksum."""
        missing = sorted({name for name in consumed_inputs if name not in self.input_sha256})
        if missing:
            raise OrthologBridgeError(
                "orthology source metadata lacks a recorded sha256 checksum for consumed "
                "input(s): " + ", ".join(missing)
            )
        return dict(self.input_sha256)

    def verify_inputs(self, paths: Mapping[str, Path]) -> None:
        """Fail closed unless each recorded checksum matches the file on disk."""
        self.require_checksums(paths.keys())
        for name in sorted(paths):
            path = Path(paths[name])
            if not path.is_file():
                raise OrthologBridgeError(f"consumed input '{name}' not found: {path}")
            if sha256_file(path) != self.input_sha256[name]:
                raise OrthologBridgeError(
                    f"consumed input '{name}' ({path}) does not match its recorded sha256 checksum"
                )

    def to_dict(self) -> Dict[str, object]:
        """Return a JSON-serializable representation."""
        return {
            "orthodb_release": self.orthodb_release,
            "gene2refseq_url": self.gene2refseq_url,
            "retrieved_at": self.retrieved_at,
            "taxonomic_scope": list(self.taxonomic_scope),
            "input_sha256": dict(self.input_sha256),
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> OrthologySourceMetadata:
        """Rebuild metadata from a source-manifest payload, failing closed."""
        required = (
            "orthodb_release",
            "gene2refseq_url",
            "retrieved_at",
            "taxonomic_scope",
            "input_sha256",
        )
        missing = [key for key in required if key not in payload]
        if missing:
            raise OrthologBridgeError(
                "source manifest is missing required field(s): " + ", ".join(missing)
            )
        scope = payload["taxonomic_scope"]
        checksums = payload["input_sha256"]
        if not isinstance(scope, (list, tuple)):
            raise OrthologBridgeError("source manifest field 'taxonomic_scope' must be a list")
        if not isinstance(checksums, Mapping):
            raise OrthologBridgeError("source manifest field 'input_sha256' must be a mapping")
        return cls(
            orthodb_release=str(payload["orthodb_release"]),
            gene2refseq_url=str(payload["gene2refseq_url"]),
            retrieved_at=str(payload["retrieved_at"]),
            taxonomic_scope=[str(item) for item in scope],
            input_sha256={str(name): str(value) for name, value in checksums.items()},
        )


def write_source_manifest(metadata: OrthologySourceMetadata, path: Path) -> None:
    """Persist a source_manifest.json next to the bridge outputs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def read_source_manifest(path: Path) -> OrthologySourceMetadata:
    """Load and validate a persisted source_manifest.json."""
    if not path.is_file():
        raise FileNotFoundError(f"Source manifest not found at {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise OrthologBridgeError(f"Source manifest at {path} is not valid JSON: {exc}") from exc
    if not isinstance(payload, Mapping):
        raise OrthologBridgeError(f"Source manifest at {path} must be a JSON object")
    return OrthologySourceMetadata.from_dict(payload)


@dataclass(frozen=True)
class OrthogroupBridgeResult:
    """Artifacts of one bridge build: table plus per-species audits."""

    table: pd.DataFrame
    retention_audit: pd.DataFrame
    duplicated_evidence: pd.DataFrame
    dropped_orthogroups: pd.DataFrame
    copy_policy: CopyPolicy


def build_orthogroup_bridge(
    og_path: Path,
    orthodb_proteins: Dict[str, str],
    prot_to_rna: Dict[str, str],
    expression_tids: Dict[str, Dict[str, str]],
    taxon_to_species: Dict[str, str],
    copy_policy: CopyPolicy = DEFAULT_COPY_POLICY,
    strict_duplicates: bool = False,
    duplicate_evidence: Optional[List[Dict[str, str]]] = None,
) -> OrthogroupBridgeResult:
    """Build the transcript-level orthogroup table with retention and duplicate audits.

    Args:
        og_path: Path to the base OrthoDB orthogroups table (tab-separated).
        orthodb_proteins: Dict mapping OrthoDB gene ID -> protein accession.
        prot_to_rna: Dict mapping protein accession -> RNA accession.
        expression_tids: Dict mapping species -> (RNA accession -> transcript ID).
        taxon_to_species: Dict mapping dataset taxonomy column names to species titles.
        copy_policy: How multi-transcript orthogroup cells are copied:
            ``'all-joined'`` (default) joins every mapped transcript, preserving
            the historical behaviour; ``'first-transcript'`` keeps one
            deterministic transcript per cell; ``'one-to-one-only'`` drops and
            records every orthogroup whose cell for any included species holds
            more than one transcript.
        strict_duplicates: When True, fail closed instead of only recording
            duplicated mapping evidence.
        duplicate_evidence: Optional pre-seeded evidence list (e.g. collected by
            the loaders); build-level evidence is appended to it.

    Returns:
        An :class:`OrthogroupBridgeResult` with the mapped table, a per-species
        retention audit, duplicated-evidence records, and dropped orthogroups.

    Raises:
        OrthologBridgeError: On an unknown copy policy, or on duplicated
            evidence when ``strict_duplicates`` is set.
    """
    if copy_policy not in COPY_POLICIES:
        raise OrthologBridgeError(
            f"copy_policy must be one of {list(COPY_POLICIES)}, got {copy_policy!r}"
        )
    evidence = list(duplicate_evidence) if duplicate_evidence else []

    og_table = pd.read_csv(og_path, sep="\t", index_col=0, dtype=str).fillna("")
    taxon_cols = list(og_table.columns)
    species_names = sorted(set(taxon_to_species.values()))

    logger.info(f"Building transcript-level orthogroup table from {og_table.shape[0]} groups")

    st = {"mapped": 0, "no_protein": 0, "no_rna_map": 0, "no_transcript": 0}
    per_species = {
        species: {
            "input_genes": 0,
            "mapped_to_protein": 0,
            "mapped_to_rna": 0,
            "mapped_to_transcript": 0,
            "one_to_one": 0,
            "one_to_many": 0,
        }
        for species in species_names
    }
    gene_transcripts: Dict[str, Dict[str, Set[str]]] = {}
    transcript_orthogroups: Dict[str, Dict[str, Set[str]]] = {}
    multi_transcript_cells: Dict[str, List[Dict[str, object]]] = {}
    rows: List[Dict[str, object]] = []

    for og_name, row in og_table.iterrows():
        orthogroup = str(og_name)
        new_row: Dict[str, object] = {}
        has_any = False

        for taxon_col in taxon_cols:
            species = taxon_to_species.get(taxon_col)
            if not species or species not in expression_tids:
                continue

            cell = row.get(taxon_col, "")
            if not cell or pd.isna(cell) or not str(cell).strip():
                new_row[species] = ""
                continue

            gene_ids = [gene.strip() for gene in str(cell).split(",") if gene.strip()]
            audit = per_species[species]
            species_gene_tids = gene_transcripts.setdefault(species, {})
            species_tid_ogs = transcript_orthogroups.setdefault(species, {})
            mapped_tids = []

            for gene_id in gene_ids:
                audit["input_genes"] += 1

                # Step 1: OrthoDB gene_id -> protein accession (base)
                prot_base = orthodb_proteins.get(gene_id)
                if not prot_base:
                    st["no_protein"] += 1
                    continue
                audit["mapped_to_protein"] += 1

                # Step 2: protein accession -> RNA accession via gene2refseq
                rna_base = prot_to_rna.get(prot_base)
                if not rna_base:
                    st["no_rna_map"] += 1
                    continue
                audit["mapped_to_rna"] += 1

                # Step 3: RNA accession -> full transcript ID from expression matrix
                tid = expression_tids[species].get(rna_base)
                if tid:
                    mapped_tids.append(tid)
                    st["mapped"] += 1
                    audit["mapped_to_transcript"] += 1
                    species_gene_tids.setdefault(gene_id, set()).add(tid)
                    species_tid_ogs.setdefault(tid, set()).add(orthogroup)
                else:
                    st["no_transcript"] += 1

            if len(mapped_tids) > 1:
                audit["one_to_many"] += 1
                if copy_policy == "one-to-one-only":
                    multi_transcript_cells.setdefault(orthogroup, []).append(
                        {
                            "orthogroup": orthogroup,
                            "species": species,
                            "transcript_count": len(mapped_tids),
                        }
                    )
            elif len(mapped_tids) == 1:
                audit["one_to_one"] += 1

            if mapped_tids:
                if copy_policy == "first-transcript":
                    new_row[species] = mapped_tids[0]
                else:
                    new_row[species] = ",".join(mapped_tids)
                has_any = True
            else:
                new_row[species] = ""

        if has_any:
            new_row["Orthogroup"] = orthogroup
            rows.append(new_row)

    dropped_orthogroups: List[Dict[str, object]] = []
    if copy_policy == "one-to-one-only":
        dropped_set = set(multi_transcript_cells)
        rows = [row for row in rows if row["Orthogroup"] not in dropped_set]
        for orthogroup in sorted(multi_transcript_cells):
            dropped_orthogroups.extend(multi_transcript_cells[orthogroup])

    # Duplicated-evidence checks: (a) one OrthoDB gene mapped to several
    # transcripts of the same species, (b) one transcript claimed by several
    # orthogroups.
    for species in sorted(gene_transcripts):
        for gene_id in sorted(gene_transcripts[species]):
            tids = gene_transcripts[species][gene_id]
            if len(tids) > 1:
                evidence.append(
                    {
                        "kind": "gene_multi_transcript",
                        "species": species,
                        "orthogroup": "",
                        "gene_id": gene_id,
                        "transcript_id": ",".join(sorted(tids)),
                        "detail": f"OrthoDB gene {gene_id} maps to {len(tids)} transcripts in {species}",
                    }
                )
    for species in sorted(transcript_orthogroups):
        for tid in sorted(transcript_orthogroups[species]):
            orthogroups = transcript_orthogroups[species][tid]
            if len(orthogroups) > 1:
                evidence.append(
                    {
                        "kind": "transcript_multi_orthogroup",
                        "species": species,
                        "orthogroup": ",".join(sorted(orthogroups)),
                        "gene_id": "",
                        "transcript_id": tid,
                        "detail": (
                            f"transcript {tid} in {species} is claimed by {len(orthogroups)} orthogroups"
                        ),
                    }
                )

    duplicated_evidence = pd.DataFrame(evidence, columns=_DUPLICATED_EVIDENCE_COLUMNS)
    if strict_duplicates and not duplicated_evidence.empty:
        kinds = sorted(duplicated_evidence["kind"].unique().tolist())
        raise OrthologBridgeError(
            f"duplicated mapping evidence detected ({len(duplicated_evidence)} record(s), kinds: "
            f"{', '.join(kinds)}); refusing to write bridge outputs"
        )

    retention_audit = pd.DataFrame(
        [
            {
                "species": species,
                "input_genes": per_species[species]["input_genes"],
                "mapped_to_protein": per_species[species]["mapped_to_protein"],
                "mapped_to_rna": per_species[species]["mapped_to_rna"],
                "mapped_to_transcript": per_species[species]["mapped_to_transcript"],
                "one_to_one": per_species[species]["one_to_one"],
                "one_to_many": per_species[species]["one_to_many"],
                "unmapped": (
                    per_species[species]["input_genes"]
                    - per_species[species]["mapped_to_transcript"]
                ),
            }
            for species in species_names
        ],
        columns=_RETENTION_AUDIT_COLUMNS,
    )

    logger.info(
        f"Mapping stats: Mapped={st['mapped']}, No Protein={st['no_protein']}, "
        f"No RNA Map={st['no_rna_map']}, No Transcript={st['no_transcript']}"
    )
    logger.info(f"Orthogroups with >=1 mapping: {len(rows)}")

    df = pd.DataFrame(rows)
    if "Orthogroup" in df.columns:
        df.set_index("Orthogroup", inplace=True)

    available = [c for c in species_names if c in df.columns]
    table = df[available]

    return OrthogroupBridgeResult(
        table=table,
        retention_audit=retention_audit,
        duplicated_evidence=duplicated_evidence,
        dropped_orthogroups=pd.DataFrame(
            dropped_orthogroups, columns=_DROPPED_ORTHOGROUP_COLUMNS
        ),
        copy_policy=copy_policy,
    )


def orthology_presence_table(table: pd.DataFrame) -> pd.DataFrame:
    """Return the orthogroup x species 0/1 presence matrix implied by a bridge table."""
    if table.empty:
        return pd.DataFrame(dtype=int)
    return (table.astype(str) != "").astype(int)


def build_transcript_orthogroup_table(
    og_path: Path,
    orthodb_proteins: Dict[str, str],
    prot_to_rna: Dict[str, str],
    expression_tids: Dict[str, Dict[str, str]],
    taxon_to_species: Dict[str, str],
    copy_policy: CopyPolicy = DEFAULT_COPY_POLICY,
) -> pd.DataFrame:
    """Build transcript-level orthogroup table using the full mapping chain.

    Resolves the chain: OrthoDB Gene ID -> NCBI Protein -> NCBI RNA -> Local Transcript ID.

    Args:
        og_path: Path to the base OrthoDB orthogroups table (tab-separated)
        orthodb_proteins: Dict mapping OrthoDB gene ID -> protein accession
        prot_to_rna: Dict mapping protein accession -> RNA accession
        expression_tids: Dict mapping species -> (RNA accession -> transcript ID)
        taxon_to_species: Dict mapping dataset taxonomy column names to species titles
        copy_policy: How multi-transcript orthogroup cells are copied; the default
            ``'all-joined'`` preserves the historical behaviour of joining every
            mapped transcript with commas.

    Returns:
        DataFrame with mapped transcript IDs for each species across orthogroups.
    """
    return build_orthogroup_bridge(
        og_path,
        orthodb_proteins,
        prot_to_rna,
        expression_tids,
        taxon_to_species,
        copy_policy=copy_policy,
    ).table

"""Evidence-derived orthogroup cardinality, retention and duplicate audits."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import pandas as pd

DUPLICATE_PRECEDENCE_VERSION = "1"
DUPLICATE_PRECEDENCE_RULES: dict[str, str] = {
    "fail": "refuse to build on any conflicting mapping evidence",
    "last-recorded": (
        "resolve protein->RNA and gene->protein accession collisions, one gene "
        "mapping to several transcripts, and one transcript claimed by several "
        "orthogroups to the last recorded row of the checksum-pinned inputs; the "
        "bridge keeps every resolvable claim and records each resolution"
    ),
}
ORTHOGROUP_CARDINALITY_CLASSES = ("one_to_one", "one_to_many", "many_to_many")
_DROP_REASON_COLUMNS = [
    "species",
    "genes_seen",
    "genes_retained",
    "dropped_no_protein",
    "dropped_no_rna_map",
    "dropped_no_transcript",
    "dropped_total",
]


def classify_orthogroup_cardinality(
    og_path: Path,
    taxon_to_species: Mapping[str, str],
) -> pd.DataFrame:
    """Label one-to-one / one-to-many / many-to-many orthogroups from input evidence.

    Counts OrthoDB gene membership per included species in the raw
    orthogroups table, mirroring the tokenization the bridge itself applies,
    so every cardinality label is justifiable from the consumed evidence.
    """
    og_table = pd.read_csv(og_path, sep="\t", index_col=0, dtype=str).fillna("")
    species_names = sorted(set(taxon_to_species.values()))
    taxon_cols = [col for col in og_table.columns if col in taxon_to_species]
    counts: dict[str, dict[str, int]] = {}
    for og_name, row in og_table.iterrows():
        og_counts = {species: 0 for species in species_names}
        for taxon_col in taxon_cols:
            species = taxon_to_species[taxon_col]
            cell = row.get(taxon_col, "")
            if not cell or pd.isna(cell) or not str(cell).strip():
                continue
            genes = [gene.strip() for gene in str(cell).split(",") if gene.strip()]
            og_counts[species] += len(genes)
        counts[str(og_name)] = og_counts
    rows: list[dict[str, object]] = []
    for orthogroup in sorted(counts):
        og_counts = counts[orthogroup]
        if not any(og_counts[species] for species in species_names):
            continue
        multi = [species for species in species_names if og_counts[species] > 1]
        if len(multi) >= 2:
            cardinality_class = "many_to_many"
        elif len(multi) == 1:
            cardinality_class = "one_to_many"
        else:
            cardinality_class = "one_to_one"
        row_out: dict[str, object] = {"orthogroup": orthogroup}
        row_out.update({species: og_counts[species] for species in species_names})
        row_out["cardinality_class"] = cardinality_class
        rows.append(row_out)
    return pd.DataFrame(rows, columns=["orthogroup", *species_names, "cardinality_class"])


def audit_gene_drop_reasons(
    og_path: Path,
    orthodb_proteins: Mapping[str, str],
    prot_to_rna: Mapping[str, str],
    expression_tids: Mapping[str, Mapping[str, str]],
    taxon_to_species: Mapping[str, str],
) -> pd.DataFrame:
    """Audit per-species gene retention with typed drop-reason classes.

    Mirrors the mapping chain of ``build_orthogroup_bridge`` so every dropped
    gene carries one of the reason classes ``no_protein``, ``no_rna_map``, or
    ``no_transcript``; drops are never silent.
    """
    og_table = pd.read_csv(og_path, sep="\t", index_col=0, dtype=str).fillna("")
    species_names = sorted(set(taxon_to_species.values()))
    tallies = {
        species: {
            "no_protein": 0,
            "no_rna_map": 0,
            "no_transcript": 0,
            "retained": 0,
            "seen": 0,
        }
        for species in species_names
    }
    for _, row in og_table.iterrows():
        for taxon_col in og_table.columns:
            species = taxon_to_species.get(taxon_col)
            if not species or species not in expression_tids:
                continue
            cell = row.get(taxon_col, "")
            if not cell or pd.isna(cell) or not str(cell).strip():
                continue
            for gene_id in (gene.strip() for gene in str(cell).split(",")):
                if not gene_id:
                    continue
                tally = tallies[species]
                tally["seen"] += 1
                prot_base = orthodb_proteins.get(gene_id)
                if not prot_base:
                    tally["no_protein"] += 1
                    continue
                rna_base = prot_to_rna.get(prot_base)
                if not rna_base:
                    tally["no_rna_map"] += 1
                    continue
                if expression_tids[species].get(rna_base):
                    tally["retained"] += 1
                else:
                    tally["no_transcript"] += 1
    rows = []
    for species in species_names:
        tally = tallies[species]
        dropped_total = tally["no_protein"] + tally["no_rna_map"] + tally["no_transcript"]
        if tally["seen"] != tally["retained"] + dropped_total:
            raise RuntimeError(
                f"drop-reason audit for {species} does not reconcile: {tally['seen']} "
                f"genes seen, {tally['retained']} retained, {dropped_total} dropped"
            )
        rows.append(
            {
                "species": species,
                "genes_seen": tally["seen"],
                "genes_retained": tally["retained"],
                "dropped_no_protein": tally["no_protein"],
                "dropped_no_rna_map": tally["no_rna_map"],
                "dropped_no_transcript": tally["no_transcript"],
                "dropped_total": dropped_total,
            }
        )
    return pd.DataFrame(rows, columns=_DROP_REASON_COLUMNS)


def annotate_duplicate_resolutions(
    evidence: pd.DataFrame,
    og_path: Path,
    orthodb_proteins: Mapping[str, str],
    prot_to_rna: Mapping[str, str],
    expression_tids: Mapping[str, Mapping[str, str]],
    taxon_to_species: Mapping[str, str],
    *,
    precedence: str,
) -> pd.DataFrame:
    """Annotate duplicated-evidence records with the declared precedence resolution.

    Under the versioned ``last-recorded`` rule the last recorded input row is
    authoritative for every conflicted claim; the bridge table keeps every
    resolvable claim and each record names the authoritative orthogroup.
    """
    if evidence.empty or precedence == "fail":
        return evidence
    claims: dict[tuple[str, str], list[str]] = {}
    og_table = pd.read_csv(og_path, sep="\t", index_col=0, dtype=str).fillna("")
    for og_name, row in og_table.iterrows():
        for taxon_col in og_table.columns:
            species = taxon_to_species.get(taxon_col)
            if not species or species not in expression_tids:
                continue
            cell = row.get(taxon_col, "")
            if not cell or pd.isna(cell) or not str(cell).strip():
                continue
            for gene_id in (gene.strip() for gene in str(cell).split(",")):
                if not gene_id:
                    continue
                prot_base = orthodb_proteins.get(gene_id)
                if not prot_base:
                    continue
                rna_base = prot_to_rna.get(prot_base)
                if not rna_base:
                    continue
                tid = expression_tids[species].get(rna_base)
                if tid:
                    claims.setdefault((species, str(tid)), []).append(str(og_name))
    records = []
    for record in evidence.to_dict("records"):
        if record["kind"] == "transcript_multi_orthogroup":
            claimants = claims.get((record["species"], record["transcript_id"]), [])
            authority = claimants[-1] if claimants else ""
            if authority:
                record["detail"] = (
                    f"{record['detail']}; resolved by last-recorded@"
                    f"{DUPLICATE_PRECEDENCE_VERSION} "
                    f"(authoritative orthogroup: {authority})"
                )
                records.append(record)
                continue
        record["detail"] = f"{record['detail']}; resolved by last-recorded@{DUPLICATE_PRECEDENCE_VERSION}"
        records.append(record)
    return pd.DataFrame(records, columns=list(evidence.columns))

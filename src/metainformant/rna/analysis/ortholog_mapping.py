"""Ortholog mapping integration module.

Provides utilities to bridge OrthoDB gene IDs, NCBI protein accessions,
and NCBI RNA accessions to construct transcript-level orthogroup tables.
"""

import gzip
import re
from pathlib import Path
from typing import Dict, Set

import pandas as pd

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def load_gene2refseq_mapping(gene2refseq_path: Path, target_taxon_ids: Set[str]) -> Dict[str, str]:
    """Load protein->RNA mappings from NCBI gene2refseq for target taxa.

    Args:
        gene2refseq_path: Path to gene2refseq file (gzipped)
        target_taxon_ids: Set of string taxon IDs to filter for

    Returns:
        Dictionary mapping protein_accession_base -> rna_accession_base.
        Base means without version (e.g. XP_001120086 -> XM_001120086).
    """
    prot_to_rna = {}
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

            prot_to_rna[prot_base] = rna_base

    logger.info(f"Loaded {len(prot_to_rna)} protein->RNA mappings")
    return prot_to_rna


def load_orthodb_proteins(genes_path: Path, target_org_prefixes: Set[str]) -> Dict[str, str]:
    """Load OrthoDB gene_id -> protein_accession (base, no version).

    Args:
        genes_path: Path to OrthoDB genes file (gzipped)
        target_org_prefixes: Set of organism prefixes (e.g., taxonomy IDs mapped to OrthoDB format)

    Returns:
        Dictionary mapping OrthoDB gene_id -> protein_accession_base.
    """
    result = {}
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
                result[gene_id] = prot_base

    logger.info(f"Loaded {len(result)} OrthoDB gene->protein mappings")
    return result


def load_expression_transcript_ids(manifest: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    """Load transcript IDs and build RNA_accession -> full_transcript_id maps.

    Args:
        manifest: DataFrame with columns 'species_title' and 'curate_path'

    Returns:
        Dictionary mapping species title -> (RNA_accession_base -> full_transcript_id)
    """
    result = {}
    for _, row in manifest.iterrows():
        sp_title = row["species_title"]
        curate_path = Path(row["curate_path"])
        if not curate_path.exists():
            continue
        try:
            tid_map = {}
            with open(curate_path) as f:
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


def build_transcript_orthogroup_table(
    og_path: Path,
    orthodb_proteins: Dict[str, str],
    prot_to_rna: Dict[str, str],
    expression_tids: Dict[str, Dict[str, str]],
    taxon_to_species: Dict[str, str],
) -> pd.DataFrame:
    """Build transcript-level orthogroup table using the full mapping chain.

    Resolves the chain: OrthoDB Gene ID -> NCBI Protein -> NCBI RNA -> Local Transcript ID.

    Args:
        og_path: Path to the base OrthoDB orthogroups table (tab-separated)
        orthodb_proteins: Dict mapping OrthoDB gene ID -> protein accession
        prot_to_rna: Dict mapping protein accession -> RNA accession
        expression_tids: Dict mapping species -> (RNA accession -> transcript ID)
        taxon_to_species: Dict mapping dataset taxonomy column names to species titles

    Returns:
        DataFrame with mapped transcript IDs for each species across orthogroups.
    """
    og_table = pd.read_csv(og_path, sep="\t", index_col=0, dtype=str).fillna("")
    taxon_cols = list(og_table.columns)
    species_names = sorted(set(taxon_to_species.values()))

    logger.info(f"Building transcript-level orthogroup table from {og_table.shape[0]} groups")

    st = {"mapped": 0, "no_protein": 0, "no_rna_map": 0, "no_transcript": 0}
    rows = []

    for og_name, row in og_table.iterrows():
        new_row = {}
        has_any = False

        for taxon_col in taxon_cols:
            sp_title = taxon_to_species.get(taxon_col)
            if not sp_title or sp_title not in expression_tids:
                continue

            cell = row.get(taxon_col, "")
            if not cell or pd.isna(cell) or not str(cell).strip():
                new_row[sp_title] = ""
                continue

            gene_ids = [g.strip() for g in str(cell).split(",") if g.strip()]
            mapped_tids = []

            for gene_id in gene_ids:
                # Step 1: OrthoDB gene_id -> protein accession (base)
                prot_base = orthodb_proteins.get(gene_id)
                if not prot_base:
                    st["no_protein"] += 1
                    continue

                # Step 2: protein accession -> RNA accession via gene2refseq
                rna_base = prot_to_rna.get(prot_base)
                if not rna_base:
                    st["no_rna_map"] += 1
                    continue

                # Step 3: RNA accession -> full transcript ID from expression matrix
                tid = expression_tids[sp_title].get(rna_base)
                if tid:
                    mapped_tids.append(tid)
                    st["mapped"] += 1
                else:
                    st["no_transcript"] += 1

            if mapped_tids:
                new_row[sp_title] = ",".join(mapped_tids)
                has_any = True
            else:
                new_row[sp_title] = ""

        if has_any:
            new_row["Orthogroup"] = og_name
            rows.append(new_row)

    logger.info(
        f"Mapping stats: Mapped={st['mapped']}, No Protein={st['no_protein']}, No RNA Map={st['no_rna_map']}, No Transcript={st['no_transcript']}"
    )
    logger.info(f"Orthogroups with >=1 mapping: {len(rows)}")

    df = pd.DataFrame(rows)
    if "Orthogroup" in df.columns:
        df.set_index("Orthogroup", inplace=True)

    available = [c for c in species_names if c in df.columns]
    return df[available]

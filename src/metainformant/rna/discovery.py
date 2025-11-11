"""Species discovery and configuration generation functions.

This module provides functions to discover species with RNA-seq data
and generate amalgkit YAML configurations.
"""

from __future__ import annotations

from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

from ..core.logging import get_logger

logger = get_logger(__name__)

# Try to import optional dependencies
try:
    from Bio import Entrez

    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    logger.warning("Biopython not available - discovery functions will be limited")

try:
    from ncbi.datasets import GenomeApi
    from ncbi.datasets.openapi import ApiClient

    NCBI_DATASETS_AVAILABLE = True
except ImportError:
    NCBI_DATASETS_AVAILABLE = False
    logger.warning("ncbi-datasets-pylib not available - genome info will be limited")


def search_species_with_rnaseq(
    search_query: str,
    *,
    max_records: int = 10000,
) -> dict[str, dict[str, Any]]:
    """Search NCBI SRA for species with RNA-seq data.

    Args:
        search_query: NCBI Entrez search query
        max_records: Maximum number of records to retrieve

    Returns:
        Dictionary mapping species names to metadata

    Raises:
        ImportError: If Biopython is not available
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("Biopython required for SRA search. Install with: pip install biopython")

    import os

    # Set email from environment
    email = os.environ.get("NCBI_EMAIL") or os.environ.get("ENTREZ_EMAIL")
    if email:
        Entrez.email = email
    else:
        Entrez.email = "your.email@example.com"

    logger.info(f"Searching NCBI SRA: {search_query}")

    try:
        handle = Entrez.esearch(db="sra", term=search_query, retmax=max_records, usehistory="y")
        search_results = Entrez.read(handle)
        handle.close()

        total_records = int(search_results["Count"])
        logger.info(f"Found {total_records} RNA-seq records")

        if total_records == 0:
            return {}

        # Fetch details in batches
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        species_data = defaultdict(
            lambda: {
                "sample_count": 0,
                "run_ids": [],
                "study_ids": set(),
                "taxonomy_id": None,
                "scientific_name": None,
            }
        )

        batch_size = 500
        for start in range(0, min(total_records, max_records), batch_size):
            end = min(start + batch_size, total_records, max_records)
            logger.info(f"Fetching records {start+1} to {end} of {total_records}...")

            try:
                fetch_handle = Entrez.efetch(
                    db="sra",
                    rettype="runinfo",
                    retmode="text",
                    retstart=start,
                    retmax=batch_size,
                    webenv=webenv,
                    query_key=query_key,
                )

                # Parse CSV format
                lines = fetch_handle.read().decode("utf-8").strip().split("\n")
                fetch_handle.close()

                if len(lines) < 2:
                    continue

                # First line is header
                header = lines[0].split(",")
                try:
                    org_idx = header.index("ScientificName")
                    run_idx = header.index("Run")
                    study_idx = header.index("SRAStudy")
                    taxid_idx = header.index("TaxID")
                except ValueError:
                    logger.warning("Could not find required columns in runinfo")
                    continue

                # Parse data lines
                for line in lines[1:]:
                    parts = line.split(",")
                    if len(parts) <= max(org_idx, run_idx, study_idx, taxid_idx):
                        continue

                    scientific_name = parts[org_idx].strip('"')
                    run_id = parts[run_idx].strip('"')
                    study_id = parts[study_idx].strip('"')
                    taxonomy_id = parts[taxid_idx].strip('"')

                    if scientific_name:
                        species_data[scientific_name]["sample_count"] += 1
                        species_data[scientific_name]["scientific_name"] = scientific_name
                        species_data[scientific_name]["taxonomy_id"] = taxonomy_id
                        if run_id:
                            species_data[scientific_name]["run_ids"].append(run_id)
                        if study_id:
                            species_data[scientific_name]["study_ids"].add(study_id)

            except Exception as e:
                logger.error(f"Error fetching batch {start}-{end}: {e}")
                continue

        # Convert sets to lists
        for species_name in species_data:
            species_data[species_name]["study_ids"] = list(species_data[species_name]["study_ids"])

        logger.info(f"Found {len(species_data)} species")
        return dict(species_data)

    except Exception as e:
        logger.error(f"Error searching SRA: {e}")
        return {}


def get_genome_info(taxonomy_id: str, species_name: str) -> dict[str, Any] | None:
    """Get genome assembly information for a species.

    Args:
        taxonomy_id: NCBI taxonomy ID
        species_name: Scientific name

    Returns:
        Genome information dictionary or None
    """
    if not NCBI_DATASETS_AVAILABLE:
        logger.warning(f"Skipping genome lookup for {species_name} (ncbi-datasets-pylib not installed)")
        return None

    logger.info(f"Fetching genome info for {species_name} (TaxID: {taxonomy_id})...")

    try:
        with ApiClient() as api_client:
            genome_api = GenomeApi(api_client)

            # Get assemblies for this taxon
            genome_metadata = genome_api.assembly_descriptors_by_taxon(
                taxon=taxonomy_id,
                filters_reference_only=False,
                filters_assembly_source="RefSeq",
                page_size=10,
            )

            assemblies = genome_metadata.get("assemblies", [])

            if not assemblies:
                # Try GenBank if RefSeq not available
                genome_metadata = genome_api.assembly_descriptors_by_taxon(
                    taxon=taxonomy_id,
                    filters_reference_only=False,
                    filters_assembly_source="GenBank",
                    page_size=10,
                )
                assemblies = genome_metadata.get("assemblies", [])

            if not assemblies:
                logger.warning(f"No genome assemblies found for {species_name}")
                return None

            # Select best assembly
            best_assembly = _select_best_assembly(assemblies)

            if not best_assembly:
                return None

            assembly_info = best_assembly.get("assembly", {})
            assembly_stats = best_assembly.get("assembly_stats", {})
            annotation_info = best_assembly.get("annotation_info", {})

            genome_info = {
                "accession": assembly_info.get("assembly_accession", ""),
                "assembly_name": assembly_info.get("assembly_name", ""),
                "level": assembly_info.get("assembly_level", ""),
                "release_date": assembly_info.get("submission_date", ""),
                "contig_n50": assembly_stats.get("contig_n50", 0),
                "scaffold_n50": assembly_stats.get("scaffold_n50", 0),
                "total_sequence_length": assembly_stats.get("total_sequence_length", 0),
                "number_of_contigs": assembly_stats.get("number_of_contigs", 0),
                "annotation_release": annotation_info.get("release_version", ""),
                "annotation_date": annotation_info.get("release_date", ""),
                "sequencing_tech": assembly_info.get("sequencing_tech", ""),
            }

            return genome_info

    except Exception as e:
        logger.error(f"Error fetching genome info for {species_name}: {e}")
        return None


def _select_best_assembly(assemblies: list[dict[str, Any]]) -> dict[str, Any] | None:
    """Select the best genome assembly from available options.

    Args:
        assemblies: List of assembly metadata dictionaries

    Returns:
        Best assembly dictionary or None
    """
    if not assemblies:
        return None

    # Score each assembly
    scored_assemblies = []
    for assembly in assemblies:
        assembly_info = assembly.get("assembly", {})
        assembly_stats = assembly.get("assembly_stats", {})

        score = 0

        # Prefer RefSeq (GCF) over GenBank (GCA)
        accession = assembly_info.get("assembly_accession", "")
        if accession.startswith("GCF_"):
            score += 1000

        # Prefer chromosome-level assemblies
        level = assembly_info.get("assembly_level", "").lower()
        if level == "chromosome":
            score += 500
        elif level == "scaffold":
            score += 100
        elif level == "contig":
            score += 10

        # Prefer higher N50
        contig_n50 = assembly_stats.get("contig_n50", 0)
        if contig_n50:
            score += min(contig_n50 / 10000, 100)  # Cap at 100 points

        scored_assemblies.append((score, assembly))

    # Return highest scoring assembly
    scored_assemblies.sort(key=lambda x: x[0], reverse=True)
    return scored_assemblies[0][1]


def generate_config_yaml(
    species_name: str,
    species_data: dict[str, Any],
    genome_info: dict[str, Any] | None = None,
    *,
    repo_root: Path | None = None,
) -> str:
    """Generate amalgkit YAML configuration for a species.

    Args:
        species_name: Scientific name
        species_data: RNA-seq metadata
        genome_info: Genome assembly metadata (optional)
        repo_root: Repository root directory for paths (optional)

    Returns:
        YAML configuration string
    """
    species_slug = species_name.replace(" ", "_").lower()
    taxonomy_id = species_data.get("taxonomy_id", "UNKNOWN")
    sample_count = species_data.get("sample_count", 0)

    # Determine work_dir path
    if repo_root:
        work_dir = f"{repo_root}/output/amalgkit/{species_slug}/work"
        log_dir = f"{repo_root}/output/amalgkit/{species_slug}/logs"
    else:
        work_dir = f"output/amalgkit/{species_slug}/work"
        log_dir = f"output/amalgkit/{species_slug}/logs"

    # Build YAML content
    yaml_lines = [
        "# METAINFORMANT Amalgkit Configuration",
        f"# Species: {species_name}",
        f"# NCBI Taxonomy ID: {taxonomy_id}",
    ]

    if genome_info:
        yaml_lines.extend(
            [
                f"# Assembly: {genome_info['accession']} ({genome_info['assembly_name']})",
                f"# Assembly Level: {genome_info['level']}",
                f"# Sequencing: {genome_info.get('sequencing_tech', 'Not specified')}",
            ]
        )
        if genome_info.get("annotation_release"):
            yaml_lines.append(f"# Annotation Release: {genome_info['annotation_release']}")

    yaml_lines.extend(
        [
            f"# RNA-seq Samples Available: {sample_count}",
            f"# Generated: {datetime.now().strftime('%Y-%m-%d')}",
            "",
            "",
            f"work_dir: {work_dir}",
            f"log_dir: {log_dir}",
            "threads: 12",
            "",
            "",
            "auto_install_amalgkit: true",
            "",
            "",
            "# Pre-selection filtering rules",
            "filters:",
            "  require_tissue: false  # Set to true to filter by tissue metadata",
            "",
            "",
            "# Species configuration",
            "species_list:",
            f"  - {species_name.replace(' ', '_')}",
            "",
            "",
            "# Alternative taxonomy specification",
            f"taxon_id: {taxonomy_id}",
            "",
            "",
        ]
    )

    # Add genome configuration if available
    if genome_info:
        accession = genome_info["accession"]
        assembly_name = genome_info["assembly_name"]

        # Construct FTP URL
        accession_parts = accession.split("_")
        if len(accession_parts) >= 2:
            gcf_gca = accession_parts[0]
            numbers = accession_parts[1].split(".")[0]
            num_groups = [numbers[i : i + 3] for i in range(0, min(9, len(numbers)), 3)]
            ftp_path = "/".join([gcf_gca] + num_groups + [f"{accession}_{assembly_name}"])
            ftp_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{ftp_path}/"
        else:
            ftp_url = "# FTP URL construction failed - please verify manually"

        yaml_lines.extend(
            [
                "# Reference genome configuration",
                f"# Based on NCBI assembly: {accession}",
                "genome:",
                f"  accession: {accession}",
                f"  assembly_name: {assembly_name}",
            ]
        )

        if genome_info.get("annotation_release"):
            yaml_lines.append(f"  annotation_release: {genome_info['annotation_release']}")

        yaml_lines.extend(
            [
                f"  dest_dir: output/amalgkit/{species_slug}/genome",
                "",
                "",
                "  # Files to download (based on available NCBI FTP resources)",
                "  include:",
                "    - genome           # Genomic sequences (FASTA)",
                "    - gff3            # Gene annotations (GFF3)",
                "    - rna             # RNA sequences (FASTA)",
                "    - cds             # CDS sequences (FASTA)",
                "    - protein         # Protein sequences (FASTA)",
                "",
                "",
                "  # Direct FTP URL for all genome files",
                f"  ftp_url: {ftp_url}",
                "",
                "",
            ]
        )
    else:
        yaml_lines.extend(
            [
                "# Reference genome configuration",
                "# NOTE: No genome assembly found in NCBI for this species",
                "# You may need to provide a custom genome or use a related species",
                "# genome:",
                "#   accession: GCF_XXXXXXXXX.X",
                "#   assembly_name: AssemblyName",
                "#   dest_dir: output/amalgkit/{species_slug}/genome",
                "",
                "",
            ]
        )

    # Add step configuration
    yaml_lines.extend(
        [
            "# Per-step parameters (merged on top of common params)",
            "steps:",
            "  metadata:",
            f"    out_dir: output/amalgkit/{species_slug}/work",
            f"    search_string: '\"{species_name}\"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]'",
            "    redo: yes",
            "  integrate:",
            f"    fastq_dir: output/amalgkit/{species_slug}/fastq",
            "  config: {}",
            "  select: {}",
            "  getfastq:",
            f"    out_dir: output/amalgkit/{species_slug}/fastq",
            "    threads: 12",
            "    aws: yes",
            "    gcp: yes",
            "    ncbi: yes",
            "    pfd: no",
            "    fastp: yes",
            "  quant:",
            f"    out_dir: output/amalgkit/{species_slug}/quant",
            "    threads: 12",
            "    redo: no",
            "    keep_fastq: no",
            "    build_index: yes",
            "  merge:",
            f"    out: output/amalgkit/{species_slug}/merged/merged_abundance.tsv",
            f"    out_dir: output/amalgkit/{species_slug}/merged",
            "  cstmm:",
            f"    out_dir: output/amalgkit/{species_slug}/cstmm",
            "  curate:",
            f"    out_dir: output/amalgkit/{species_slug}/curate",
            "  csca:",
            f"    out_dir: output/amalgkit/{species_slug}/csca",
            "  sanity: {}",
            "",
        ]
    )

    return "\n".join(yaml_lines)


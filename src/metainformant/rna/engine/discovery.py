"""RNA species discovery and genome configuration generation.

This module provides tools for discovering RNA-seq datasets for specific species
and generating configuration files for amalgkit workflows. It integrates with
NCBI Entrez and NCBI Datasets to find species information and genome assemblies.

Main Functions:
    - generate_config_yaml: Generate YAML config from species discovery data
    - _select_best_assembly: Select optimal genome assembly from candidates
    - search_species_with_rnaseq: Search for species with RNA-seq data (requires Biopython)
    - get_genome_info: Fetch genome metadata from NCBI (requires ncbi-datasets-pylib)

Example:
    >>> from metainformant.rna import discovery
    >>> species_data = {
    ...     "sample_count": 100,
    ...     "taxonomy_id": "9606",
    ...     "run_ids": ["SRR123"],
    ...     "study_ids": ["SRP001"]
    ... }
    >>> yaml_str = discovery.generate_config_yaml("Homo sapiens", species_data)

Environment:
    BioPython and ncbi-datasets-pylib are optional dependencies.
    Functions requiring them will raise ImportError if unavailable.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Optional dependency flags
try:
    from Bio import Entrez

    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    from ncbi_datasets.openapi import openapi_client

    NCBI_DATASETS_AVAILABLE = True
except ImportError:
    NCBI_DATASETS_AVAILABLE = False


def generate_config_yaml(
    species_name: str,
    species_data: Dict[str, Any],
    genome_info: Optional[Dict[str, Any]] = None,
    repo_root: Optional[Path] = None,
) -> str:
    """Generate YAML configuration for amalgkit workflow.

    Creates a YAML configuration string for RNA-seq analysis based on species
    discovery data. Optionally includes genome information and file paths.

    Args:
        species_name: Scientific name of species (e.g., "Homo sapiens")
        species_data: Dictionary containing:
            - sample_count (int): Number of RNA-seq samples
            - taxonomy_id (str): NCBI taxonomy ID
            - run_ids (list[str]): SRA run accessions
            - study_ids (list[str]): SRA study accessions
        genome_info: Optional dict with genome assembly info:
            - accession (str): Assembly accession (e.g., "GCF_000001405.40")
            - assembly_name (str): Assembly name
            - level (str): Assembly level
            - annotation_release (str): Annotation release number
        repo_root: Optional path to repository root for setting work_dir

    Returns:
        YAML content as string formatted for amalgkit configuration.

    Example:
        >>> species_data = {
        ...     "sample_count": 100,
        ...     "taxonomy_id": "9606",
        ...     "run_ids": ["SRR001", "SRR002"],
        ...     "study_ids": ["SRP001"]
        ... }
        >>> yaml = generate_config_yaml("Homo sapiens", species_data)
        >>> assert "Homo sapiens" in yaml
        >>> assert "9606" in yaml
    """
    # Determine work directory
    if repo_root:
        repo_root = Path(repo_root)
        work_dir = str(repo_root / "output" / "amalgkit" / species_name.replace(" ", "_"))
    else:
        work_dir = f"output/amalgkit/{species_name.replace(' ', '_')}"

    # Start building YAML
    yaml_lines = []
    yaml_lines.append("# Amalgkit workflow configuration")
    yaml_lines.append(f"# Generated for {species_name}")
    yaml_lines.append("")
    yaml_lines.append("work_dir: " + work_dir)
    yaml_lines.append("")
    yaml_lines.append("species_list:")
    yaml_lines.append(f"  - name: {species_name}")
    yaml_lines.append(f"    taxonomy_id: {species_data['taxonomy_id']}")
    yaml_lines.append(f"    sample_count: {species_data['sample_count']}")

    # Add run IDs
    if species_data.get("run_ids"):
        yaml_lines.append("    run_ids:")
        for run_id in species_data["run_ids"]:
            yaml_lines.append(f"      - {run_id}")

    # Add study IDs
    if species_data.get("study_ids"):
        yaml_lines.append("    study_ids:")
        for study_id in species_data["study_ids"]:
            yaml_lines.append(f"      - {study_id}")

    # Add genome info if provided
    if genome_info:
        yaml_lines.append("")
        yaml_lines.append("    genome:")
        yaml_lines.append(f"      accession: {genome_info.get('accession', 'unknown')}")
        yaml_lines.append(f"      assembly_name: {genome_info.get('assembly_name', 'unknown')}")
        yaml_lines.append(f"      level: {genome_info.get('level', 'unknown')}")

        # FTP URL construction (NCBI standard pattern)
        accession = genome_info.get("accession", "")
        if accession.startswith("GCF_"):
            # RefSeq FTP URL pattern
            ftp_url = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{accession[:7]}/{accession[8:11]}/{accession[12:]}"
        elif accession.startswith("GCA_"):
            # GenBank FTP URL pattern
            ftp_url = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{accession[:7]}/{accession[8:11]}/{accession[12:]}"
        else:
            ftp_url = ""

        if ftp_url:
            yaml_lines.append(f"      ftp_url: {ftp_url}")

    yaml_lines.append("")

    yaml_content = "\n".join(yaml_lines)
    logger.debug(f"Generated YAML config for {species_name}")

    return yaml_content


def _select_best_assembly(assemblies: list[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Select best genome assembly from list of candidates.

    Uses heuristics to select the optimal assembly:
    1. Prefer RefSeq (GCF_) over GenBank (GCA_)
    2. Prefer chromosome-level over scaffold-level over contig-level
    3. Return first valid option

    Args:
        assemblies: List of assembly dictionaries with structure:
            - assembly:
                - assembly_accession: str (e.g., "GCF_000001405.40")
                - assembly_level: str ("chromosome", "scaffold", "contig")
            - assembly_stats:
                - contig_n50: int

    Returns:
        Best assembly dict from input list or None if list is empty.

    Example:
        >>> assemblies = [
        ...     {
        ...         "assembly": {"assembly_accession": "GCA_001", "assembly_level": "contig"},
        ...         "assembly_stats": {"contig_n50": 1000}
        ...     },
        ...     {
        ...         "assembly": {"assembly_accession": "GCF_001", "assembly_level": "chromosome"},
        ...         "assembly_stats": {"contig_n50": 50000}
        ...     }
        ... ]
        >>> best = _select_best_assembly(assemblies)
        >>> assert best["assembly"]["assembly_accession"].startswith("GCF_")
    """
    if not assemblies:
        return None

    # Level preference ranking
    level_rank = {"chromosome": 0, "scaffold": 1, "contig": 2}

    def assembly_key(assembly: Dict[str, Any]) -> tuple:
        """Create sort key for assembly selection."""
        accession = assembly.get("assembly", {}).get("assembly_accession", "")
        level = assembly.get("assembly", {}).get("assembly_level", "contig")

        # Prefer RefSeq (GCF_) over GenBank (GCA_)
        is_refseq = 0 if accession.startswith("GCF_") else 1

        # Prefer chromosome > scaffold > contig
        level_pref = level_rank.get(level, 3)

        return (is_refseq, level_pref)

    # Sort by preference and return best
    best = min(assemblies, key=assembly_key)
    logger.debug(f"Selected assembly: {best.get('assembly', {}).get('assembly_accession', 'unknown')}")

    return best


def search_species_with_rnaseq(query: str, max_records: int = 10, email: Optional[str] = None) -> Dict[str, Any]:
    """Search for species with RNA-seq data available.

    Searches NCBI Entrez SRA database for RNA-seq studies related to the given query.
    Requires BioPython (Entrez module).

    Args:
        query: Search query (species name, taxonomy ID, etc.)
        max_records: Maximum number of results to return
        email: Email for NCBI Entrez (required by NCBI, uses env var NCBI_EMAIL if not provided)

    Returns:
        Dictionary with:
            - query: Original search query
            - results: List of SRA study records with metadata
            - total: Total number of matching records in database
            - returned: Number of records returned in this query

    Raises:
        ImportError: If BioPython is not installed
        RuntimeError: If NCBI Entrez query fails

    Example:
        >>> if BIOPYTHON_AVAILABLE:
        ...     results = search_species_with_rnaseq("Homo sapiens", max_records=5)
        ...     print(f"Found {results['total']} studies")
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("BioPython required for species search. Install with: uv pip install biopython")

    import os

    # Set email for NCBI (required by NCBI policy)
    entrez_email = email or os.environ.get("NCBI_EMAIL", "metainformant@example.com")
    Entrez.email = entrez_email

    logger.info(f"Searching NCBI SRA for RNA-seq data: {query}")

    # Build search query for RNA-seq data
    # Combine species query with RNA-seq library strategy
    search_query = f'({query}[Organism]) AND ("RNA-Seq"[Strategy] OR "RNA Seq"[Title])'

    try:
        # Search SRA database
        handle = Entrez.esearch(db="sra", term=search_query, retmax=max_records, usehistory="y")
        search_results = Entrez.read(handle)
        handle.close()

        total_count = int(search_results.get("Count", 0))
        id_list = search_results.get("IdList", [])

        if not id_list:
            logger.info(f"No RNA-seq data found for query: {query}")
            return {"query": query, "results": [], "total": 0, "returned": 0}

        # Fetch detailed records for the found IDs
        logger.info(f"Found {total_count} records, fetching details for {len(id_list)}")

        fetch_handle = Entrez.efetch(db="sra", id=",".join(id_list), rettype="xml", retmode="xml")
        fetch_data = fetch_handle.read()
        fetch_handle.close()

        # Parse the XML response
        results = _parse_sra_xml(fetch_data, query)

        return {"query": query, "results": results, "total": total_count, "returned": len(results)}

    except Exception as e:
        logger.error(f"NCBI Entrez query failed: {e}")
        raise RuntimeError(f"Failed to search NCBI SRA: {e}") from e


def _parse_sra_xml(xml_data: bytes, query: str) -> list:
    """Parse SRA XML response into structured records.

    Args:
        xml_data: Raw XML bytes from Entrez efetch
        query: Original search query for reference

    Returns:
        List of dictionaries with study metadata
    """
    import xml.etree.ElementTree as ET

    results = []

    try:
        root = ET.fromstring(xml_data)

        # SRA XML structure has EXPERIMENT_PACKAGE elements
        for package in root.findall(".//EXPERIMENT_PACKAGE"):
            record = {}

            # Extract experiment info
            experiment = package.find(".//EXPERIMENT")
            if experiment is not None:
                record["experiment_accession"] = experiment.get("accession", "")
                record["experiment_alias"] = experiment.get("alias", "")

                # Get title
                title_elem = experiment.find(".//TITLE")
                if title_elem is not None and title_elem.text:
                    record["title"] = title_elem.text

                # Get library info
                library = experiment.find(".//LIBRARY_DESCRIPTOR")
                if library is not None:
                    strategy = library.find("LIBRARY_STRATEGY")
                    if strategy is not None and strategy.text:
                        record["library_strategy"] = strategy.text
                    source = library.find("LIBRARY_SOURCE")
                    if source is not None and source.text:
                        record["library_source"] = source.text

            # Extract study info
            study = package.find(".//STUDY")
            if study is not None:
                record["study_accession"] = study.get("accession", "")
                study_title = study.find(".//STUDY_TITLE")
                if study_title is not None and study_title.text:
                    record["study_title"] = study_title.text

            # Extract sample info
            sample = package.find(".//SAMPLE")
            if sample is not None:
                record["sample_accession"] = sample.get("accession", "")

                # Get organism info
                sample_name = sample.find(".//SAMPLE_NAME")
                if sample_name is not None:
                    taxon_id = sample_name.find("TAXON_ID")
                    if taxon_id is not None and taxon_id.text:
                        record["taxonomy_id"] = taxon_id.text
                    scientific_name = sample_name.find("SCIENTIFIC_NAME")
                    if scientific_name is not None and scientific_name.text:
                        record["organism"] = scientific_name.text

            # Extract run info
            runs = package.findall(".//RUN")
            record["run_accessions"] = [r.get("accession", "") for r in runs if r.get("accession")]

            if record:
                results.append(record)

    except ET.ParseError as e:
        logger.warning(f"Failed to parse SRA XML: {e}")

    return results


def get_genome_info(taxonomy_id: str, species_name: str) -> Optional[Dict[str, Any]]:
    """Fetch genome information for a species from NCBI Datasets.

    Retrieves genome assembly information including accession, assembly name,
    and taxonomic details. Requires ncbi-datasets-pylib.

    Args:
        taxonomy_id: NCBI taxonomy ID as string
        species_name: Scientific species name for reference

    Returns:
        Dictionary with genome info including:
            - accession: Assembly accession (e.g., "GCF_000001405.40")
            - assembly_name: Human-readable assembly name
            - level: Assembly level (chromosome, scaffold, contig)
            - organism: Scientific name
            - taxonomy_id: NCBI taxonomy ID
            - annotation_release: Annotation release info if available
            - ftp_path: FTP download path
        Returns None if no genome found

    Raises:
        ImportError: If ncbi-datasets-pylib is not installed
        RuntimeError: If API query fails

    Example:
        >>> if NCBI_DATASETS_AVAILABLE:
        ...     info = get_genome_info("9606", "Homo sapiens")
        ...     print(f"Assembly: {info['accession']}")
    """
    if not NCBI_DATASETS_AVAILABLE:
        raise ImportError(
            "ncbi-datasets-pylib required for genome info. Install with: uv pip install ncbi-datasets-pylib"
        )

    logger.info(f"Fetching genome info for {species_name} (taxonomy: {taxonomy_id})")

    try:
        from ncbi.datasets import GenomeApi
        from ncbi.datasets.openapi import ApiClient

        # Create API client and genome API instance
        with ApiClient() as api_client:
            genome_api = GenomeApi(api_client)

            # Query by taxonomy ID
            try:
                response = genome_api.genome_tax_id(int(taxonomy_id))
            except Exception:
                # Fallback: try by species name
                logger.info(f"Taxonomy ID query failed, trying species name: {species_name}")
                response = genome_api.genome_tax_name(species_name)

            if not response or not hasattr(response, "assemblies") or not response.assemblies:
                logger.warning(f"No genome assemblies found for {species_name}")
                return None

            # Select best assembly using existing helper
            assemblies = []
            for assembly_data in response.assemblies:
                if hasattr(assembly_data, "assembly"):
                    asm = assembly_data.assembly
                    assemblies.append(
                        {
                            "assembly": {
                                "assembly_accession": getattr(asm, "assembly_accession", ""),
                                "assembly_level": getattr(asm, "assembly_level", ""),
                                "assembly_name": getattr(asm, "assembly_name", ""),
                                "org": {
                                    "sci_name": (
                                        getattr(getattr(asm, "org", None), "sci_name", species_name)
                                        if hasattr(asm, "org")
                                        else species_name
                                    ),
                                    "tax_id": taxonomy_id,
                                },
                            },
                            "assembly_stats": {
                                "contig_n50": (
                                    getattr(getattr(assembly_data, "assembly_stats", None), "contig_n50", 0)
                                    if hasattr(assembly_data, "assembly_stats")
                                    else 0
                                )
                            },
                        }
                    )

            best = _select_best_assembly(assemblies)

            if not best:
                return None

            asm_info = best.get("assembly", {})
            accession = asm_info.get("assembly_accession", "")

            # Construct FTP path from accession
            ftp_path = _construct_ftp_path(accession)

            return {
                "accession": accession,
                "assembly_name": asm_info.get("assembly_name", ""),
                "level": asm_info.get("assembly_level", ""),
                "organism": asm_info.get("org", {}).get("sci_name", species_name),
                "taxonomy_id": taxonomy_id,
                "ftp_path": ftp_path,
            }

    except ImportError as e:
        # Handle case where ncbi.datasets module structure is different
        logger.warning(f"ncbi-datasets API structure issue: {e}")
        return _get_genome_info_fallback(taxonomy_id, species_name)
    except Exception as e:
        logger.error(f"Failed to fetch genome info: {e}")
        raise RuntimeError(f"NCBI Datasets API query failed: {e}") from e


def _construct_ftp_path(accession: str) -> str:
    """Construct NCBI FTP path from assembly accession.

    Args:
        accession: Assembly accession (e.g., "GCF_000001405.40")

    Returns:
        FTP URL for the assembly directory
    """
    if not accession:
        return ""

    # Parse accession: GCF_000001405.40 or GCA_000001405.40
    prefix = accession[:3]  # GCF or GCA
    if prefix not in ("GCF", "GCA"):
        return ""

    # Remove prefix and version
    base = accession[4:].split(".")[0]  # "000001405"

    # Pad to 9 digits if needed
    base = base.zfill(9)

    # Split into 3-character segments
    seg1 = base[0:3]
    seg2 = base[3:6]
    seg3 = base[6:9]

    ftp_base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{seg1}/{seg2}/{seg3}"
    return ftp_base


def _get_genome_info_fallback(taxonomy_id: str, species_name: str) -> Optional[Dict[str, Any]]:
    """Fallback genome info lookup using NCBI Entrez when Datasets API unavailable.

    Args:
        taxonomy_id: NCBI taxonomy ID
        species_name: Scientific species name

    Returns:
        Genome info dict or None
    """
    if not BIOPYTHON_AVAILABLE:
        return None

    import os

    try:
        Entrez.email = os.environ.get("NCBI_EMAIL", "metainformant@example.com")

        # Search NCBI Assembly database
        search_query = f'({species_name}[Organism]) AND "latest refseq"[filter]'
        handle = Entrez.esearch(db="assembly", term=search_query, retmax=5)
        results = Entrez.read(handle)
        handle.close()

        id_list = results.get("IdList", [])
        if not id_list:
            return None

        # Fetch assembly summary
        handle = Entrez.esummary(db="assembly", id=id_list[0])
        summary = Entrez.read(handle)
        handle.close()

        if not summary or "DocumentSummarySet" not in summary:
            return None

        doc_summary = summary["DocumentSummarySet"]["DocumentSummary"][0]

        accession = doc_summary.get("AssemblyAccession", "")
        return {
            "accession": accession,
            "assembly_name": doc_summary.get("AssemblyName", ""),
            "level": doc_summary.get("AssemblyStatus", ""),
            "organism": doc_summary.get("Organism", species_name),
            "taxonomy_id": taxonomy_id,
            "ftp_path": _construct_ftp_path(accession),
        }

    except Exception as e:
        logger.warning(f"Fallback genome lookup failed: {e}")
        return None

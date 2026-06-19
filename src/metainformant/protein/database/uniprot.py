"""UniProt database integration utilities.

This module provides tools for accessing and parsing UniProt protein data,
including sequence retrieval, annotations, and functional information.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import requests

from metainformant.core.utils import logging
from metainformant.protein._network import get_protein_api_timeout

logger = logging.get_logger(__name__)


def fetch_uniprot_record(uniprot_id: str) -> Dict[str, Any]:
    """Fetch complete UniProt record for a protein.

    Args:
        uniprot_id: UniProt accession or entry name

    Returns:
        Dictionary containing UniProt record data

    Raises:
        requests.RequestException: If API request fails

    Example:
        >>> # This would fetch real UniProt data
        >>> # record = fetch_uniprot_record("P12345")
        >>> # "sequence" in record
        >>> # True
    """
    base_url = "https://www.uniprot.org/uniprotkb"
    url = f"{base_url}/{uniprot_id}.json"

    try:
        response = requests.get(url, timeout=get_protein_api_timeout())
        response.raise_for_status()

        data = response.json()

        # Extract relevant information
        record = {
            "accession": data.get("primaryAccession"),
            "entry_name": data.get("uniProtkbId"),
            "protein_name": data.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value"),
            "organism": data.get("organism", {}).get("scientificName"),
            "taxon_id": data.get("organism", {}).get("taxonId"),
            "sequence": data.get("sequence", {}).get("value"),
            "length": data.get("sequence", {}).get("length"),
            "gene_name": data.get("genes", [{}])[0].get("geneName", {}).get("value") if data.get("genes") else None,
            "function": (
                data.get("comments", [{}])[0].get("texts", [{}])[0].get("value") if data.get("comments") else None
            ),
            "subcellular_location": _extract_subcellular_location(data),
            "domains": _extract_domains(data),
            "ptms": _extract_ptms(data),
            "go_terms": _extract_go_terms(data),
            "keywords": _extract_keywords(data),
        }

        logger.info(f"Fetched UniProt record for {uniprot_id}")
        return record

    except requests.RequestException as e:
        logger.error(f"Failed to fetch UniProt record for {uniprot_id}: {e}")
        raise


def _extract_subcellular_location(data: Dict[str, Any]) -> List[str]:
    """Extract subcellular location information."""
    locations = []

    comments = data.get("comments", [])
    for comment in comments:
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            for location in comment.get("subcellularLocations", []):
                location_name = location.get("location", {}).get("value")
                if location_name:
                    locations.append(location_name)

    return locations


def _extract_domains(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Extract domain information."""
    domains = []

    features = data.get("features", [])
    for feature in features:
        if feature.get("type") == "Domain":
            domain = {
                "name": feature.get("description"),
                "start": feature.get("location", {}).get("start", {}).get("value"),
                "end": feature.get("location", {}).get("end", {}).get("value"),
            }
            domains.append(domain)

    return domains


def _extract_ptms(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Extract post-translational modification information."""
    ptms = []

    features = data.get("features", [])
    for feature in features:
        if feature.get("type") in ["Modified residue", "Glycosylation", "Disulfide bond"]:
            ptm = {
                "type": feature.get("type"),
                "description": feature.get("description"),
                "position": feature.get("location", {}).get("position", {}).get("value")
                or feature.get("location", {}).get("start", {}).get("value"),
            }
            ptms.append(ptm)

    return ptms


def _extract_go_terms(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Extract Gene Ontology cross-references from a UniProt JSON record."""
    go_terms = []
    aspect_names = {
        "C": "cellular_component",
        "F": "molecular_function",
        "P": "biological_process",
    }

    for xref in data.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "GO":
            continue

        properties = {prop.get("key"): prop.get("value") for prop in xref.get("properties", [])}
        go_term = properties.get("GoTerm")
        aspect = None
        name = go_term
        if isinstance(go_term, str) and ":" in go_term:
            aspect_code, name = go_term.split(":", 1)
            aspect = aspect_names.get(aspect_code, aspect_code)

        go_terms.append(
            {
                "id": xref.get("id"),
                "name": name,
                "aspect": properties.get("GoAspect") or aspect,
                "evidence": properties.get("GoEvidenceType"),
            }
        )

    return go_terms


def _extract_keywords(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Extract keyword annotations from a UniProt JSON record."""
    keywords = []
    for keyword in data.get("keywords", []):
        if isinstance(keyword, dict):
            keywords.append(
                {
                    "id": keyword.get("id"),
                    "name": keyword.get("name") or keyword.get("value"),
                    "category": keyword.get("category"),
                }
            )
        elif keyword:
            keywords.append({"id": None, "name": str(keyword), "category": None})
    return keywords


def fetch_uniprot_fasta(uniprot_id: str) -> Optional[str]:
    """Fetch protein sequence in FASTA format from UniProt.

    Args:
        uniprot_id: UniProt accession

    Returns:
        FASTA sequence string, or None if not found

    Example:
        >>> # This would fetch real FASTA data
        >>> # fasta = fetch_uniprot_fasta("P12345")
        >>> # fasta.startswith(">")
        >>> # True
    """
    url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}.fasta"

    try:
        response = requests.get(url, timeout=get_protein_api_timeout())
        response.raise_for_status()

        return response.text

    except requests.RequestException as e:
        logger.error(f"Failed to fetch FASTA for {uniprot_id}: {e}")
        return None


def parse_uniprot_fasta_header(header: str) -> Dict[str, str]:
    """Parse UniProt FASTA header to extract metadata.

    Args:
        header: FASTA header line (without >)

    Returns:
        Dictionary with parsed header information

    Example:
        >>> header = "sp|P12345|PROT_HUMAN Protein name OS=Homo sapiens"
        >>> parsed = parse_uniprot_fasta_header(header)
        >>> parsed['accession'] == 'P12345'
        True
    """
    # UniProt FASTA header format: db|accession|entry_name description OS=organism
    parsed = {"database": "", "accession": "", "entry_name": "", "description": "", "organism": "", "gene_name": ""}

    if not header:
        return parsed

    # Split by spaces, but handle OS= and GN= specially
    parts = header.split()

    if len(parts) >= 1:
        # First part: db|accession|entry_name
        id_part = parts[0]
        id_components = id_part.split("|")
        if len(id_components) >= 3:
            parsed["database"] = id_components[0]
            parsed["accession"] = id_components[1]
            parsed["entry_name"] = id_components[2]

    # Extract description and organism
    description_parts = []
    i = 1
    while i < len(parts):
        part = parts[i]
        if part.startswith("OS="):
            # Organism
            organism_parts = []
            while i < len(parts) and not parts[i].startswith("GN=") and not parts[i].startswith("PE="):
                organism_parts.append(parts[i])
                i += 1
            parsed["organism"] = " ".join(organism_parts).replace("OS=", "")
        elif part.startswith("GN="):
            # Gene name
            parsed["gene_name"] = part.replace("GN=", "")
            i += 1
        else:
            description_parts.append(part)
            i += 1

    parsed["description"] = " ".join(description_parts)

    return parsed


def get_uniprot_annotations(uniprot_id: str) -> List[Dict[str, Any]]:
    """Get functional annotations for a UniProt entry.

    Args:
        uniprot_id: UniProt accession

    Returns:
        List of annotation dictionaries

    Example:
        >>> # This would get real annotations
        >>> # annotations = get_uniprot_annotations("P12345")
        >>> # isinstance(annotations, list)
        >>> # True
    """
    try:
        record = fetch_uniprot_record(uniprot_id)
        annotations = []

        for go_term in record.get("go_terms", []):
            annotations.append(
                {
                    "type": "GO",
                    "id": go_term.get("id"),
                    "name": go_term.get("name"),
                    "aspect": go_term.get("aspect"),
                    "evidence": go_term.get("evidence"),
                }
            )

        for keyword in record.get("keywords", []):
            annotations.append(
                {
                    "type": "Keyword",
                    "id": keyword.get("id"),
                    "name": keyword.get("name"),
                    "category": keyword.get("category"),
                }
            )

        return annotations

    except Exception as e:
        logger.error(f"Failed to get annotations for {uniprot_id}: {e}")
        return []


def search_uniprot_proteins(query: str, max_results: int = 100) -> List[Dict[str, Any]]:
    """Search UniProt database for proteins matching query.

    Args:
        query: Search query
        max_results: Maximum number of results

    Returns:
        List of matching protein records

    Example:
        >>> # This would search UniProt
        >>> # results = search_uniprot_proteins("kinase AND human", max_results=10)
        >>> # isinstance(results, list)
        >>> # True
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": query, "format": "json", "size": min(max_results, 500)}  # API limit

    try:
        response = requests.get(base_url, params=params, timeout=get_protein_api_timeout())
        response.raise_for_status()

        data = response.json()

        results = []
        for result in data.get("results", [])[:max_results]:
            protein = {
                "accession": result.get("primaryAccession"),
                "entry_name": result.get("uniProtkbId"),
                "protein_name": result.get("proteinDescription", {})
                .get("recommendedName", {})
                .get("fullName", {})
                .get("value"),
                "organism": result.get("organism", {}).get("scientificName"),
                "sequence_length": result.get("sequence", {}).get("length"),
                "gene_name": (
                    result.get("genes", [{}])[0].get("geneName", {}).get("value") if result.get("genes") else None
                ),
            }
            results.append(protein)

        logger.info(f"Found {len(results)} proteins matching query: {query}")
        return results

    except requests.RequestException as e:
        logger.error(f"UniProt search failed for query '{query}': {e}")
        return []


def get_uniprot_taxonomy_info(taxon_id: int) -> Optional[Dict[str, Any]]:
    """Get taxonomy information from UniProt.

    Args:
        taxon_id: NCBI taxonomy ID

    Returns:
        Taxonomy information dictionary

    Example:
        >>> # This would get taxonomy info
        >>> # taxonomy = get_uniprot_taxonomy_info(9606)
        >>> # taxonomy['scientific_name'] == 'Homo sapiens'
        >>> # True
    """
    url = f"https://www.uniprot.org/taxonomy/{taxon_id}.json"

    try:
        response = requests.get(url, timeout=get_protein_api_timeout())
        response.raise_for_status()

        data = response.json()

        taxonomy = {
            "taxon_id": data.get("taxonId"),
            "scientific_name": data.get("scientificName"),
            "common_name": data.get("commonName"),
            "rank": data.get("rank"),
            "lineage": data.get("lineage", []),
            "parent_id": data.get("parent", {}).get("taxonId"),
        }

        return taxonomy

    except requests.RequestException as e:
        logger.error(f"Failed to get taxonomy info for {taxon_id}: {e}")
        return None


def batch_fetch_uniprot_records(uniprot_ids: List[str]) -> Dict[str, Dict[str, Any]]:
    """Fetch multiple UniProt records in batch.

    Args:
        uniprot_ids: List of UniProt accessions

    Returns:
        Dictionary mapping IDs to records

    Example:
        >>> ids = ["P12345", "P67890"]
        >>> records = batch_fetch_uniprot_records(ids)
        >>> isinstance(records, dict)
        True
    """
    results = {}

    for uniprot_id in uniprot_ids:
        try:
            record = fetch_uniprot_record(uniprot_id)
            results[uniprot_id] = record
        except Exception as e:
            logger.error(f"Failed to fetch record for {uniprot_id}: {e}")
            results[uniprot_id] = None

    return results


def validate_uniprot_accession(accession: str) -> bool:
    """Validate UniProt accession format.

    Args:
        accession: UniProt accession to validate

    Returns:
        True if format is valid

    Example:
        >>> validate_uniprot_accession("P12345")
        True
        >>> validate_uniprot_accession("INVALID")
        False
    """
    import re

    if not isinstance(accession, str):
        return False

    # UniProtKB accession format:
    # - 6 chars: [OPQ][0-9][A-Z0-9]{3}[0-9]
    # - 10 chars: [A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]
    patterns = [
        r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$",
        r"^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]$",
    ]

    return any(re.match(pattern, accession) for pattern in patterns)


def map_ids_uniprot(
    protein_ids: List[str], source_db: str = "auto", target_format: str = "accession"
) -> Dict[str, str]:
    """Map protein identifiers to UniProt identifiers.

    Args:
        protein_ids: List of protein identifiers to map
        source_db: Source database ('auto', 'ensembl', 'refseq', 'pdb', etc.)
        target_format: Target format ('accession', 'entry_name', 'id')

    Returns:
        Dictionary mapping input IDs to UniProt IDs

    Raises:
        ValueError: If parameters are invalid

    Example:
        >>> ids = ["ENSP00000389680", "NP_001005484"]
        >>> mapping = map_ids_uniprot(ids, source_db="ensembl")
        >>> # Returns mapping to UniProt accessions
    """
    if not protein_ids:
        return {}

    if target_format not in ["accession", "entry_name", "id"]:
        raise ValueError(f"Invalid target format: {target_format}")

    # UniProt ID mapping API endpoint
    base_url = "https://www.uniprot.org/uploadlists/"

    # Prepare the request data
    ids_string = " ".join(protein_ids)

    # Determine source database
    if source_db == "auto":
        # Try to auto-detect based on ID patterns
        if protein_ids and protein_ids[0].startswith("ENSP"):
            from_db = "ENSEMBL_PRO_ID"
        elif protein_ids and protein_ids[0].startswith("NP_"):
            from_db = "RefSeq_Protein"
        elif protein_ids and protein_ids[0].startswith("P"):
            from_db = "UniProtKB_AC-ID"
        else:
            from_db = "UniProtKB_AC-ID"  # Default fallback
    else:
        # Map database names to UniProt API codes
        db_mapping = {
            "ensembl": "ENSEMBL_PRO_ID",
            "refseq": "RefSeq_Protein",
            "pdb": "PDB",
            "geneid": "GeneID",
            "uniprot": "UniProtKB_AC-ID",
        }
        from_db = db_mapping.get(source_db.lower())
        if from_db is None:
            raise ValueError(f"Unsupported source database: {source_db}")

    # Map target format to UniProt API codes
    to_db = {"accession": "UniProtKB", "entry_name": "UniProtKB", "id": "UniProtKB"}[target_format]

    try:
        # Make the API request
        response = requests.post(
            base_url,
            data={"from": from_db, "to": to_db, "format": "tab", "query": ids_string},
            timeout=get_protein_api_timeout(default=60.0),
        )
        response.raise_for_status()

        # Parse the tab-separated response
        lines = response.text.strip().split("\n")
        if len(lines) < 2:  # No mappings found
            logger.warning(f"No mappings found for {len(protein_ids)} IDs")
            return {}

        # Skip header line and parse mappings
        mapping = {}
        for line in lines[1:]:
            if "\t" in line:
                parts = line.split("\t")
                if len(parts) >= 2:
                    input_id, uniprot_id = parts[0], parts[1]

                    if target_format == "accession":
                        # Extract just the accession part
                        uniprot_id = uniprot_id.split("_")[0] if "_" in uniprot_id else uniprot_id

                    mapping[input_id] = uniprot_id

        logger.info(f"Successfully mapped {len(mapping)} out of {len(protein_ids)} IDs")
        return mapping

    except requests.RequestException as e:
        logger.error(f"UniProt ID mapping failed: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in ID mapping: {e}")
        raise

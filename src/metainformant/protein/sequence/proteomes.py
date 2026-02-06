"""Protein proteome analysis and taxonomy utilities.

This module provides functions for working with protein proteomes,
taxonomy IDs, and proteome-level analysis.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import io, logging

logger = logging.get_logger(__name__)


def read_taxon_ids(file_path: Union[str, Path]) -> List[int]:
    """Read taxonomy IDs from a file.

    Args:
        file_path: Path to file containing taxonomy IDs

    Returns:
        List of taxonomy IDs as integers

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Taxonomy ID file not found: {file_path}")

    taxon_ids = []

    try:
        with io.open_text_auto(file_path) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    # Extract taxonomy ID (could be just numbers or with prefixes)
                    match = re.search(r"\b(\d+)\b", line)
                    if match:
                        taxon_ids.append(int(match.group(1)))

    except Exception as e:
        raise ValueError(f"Error reading taxonomy IDs from {file_path}: {e}")

    logger.info(f"Read {len(taxon_ids)} taxonomy IDs from {file_path}")
    return taxon_ids


def validate_taxon_ids(taxon_ids: List[str]) -> Tuple[List[str], List[str]]:
    """Validate taxonomy IDs.

    Args:
        taxon_ids: List of taxonomy ID strings

    Returns:
        Tuple of (valid_ids, invalid_ids)
    """
    valid_ids = []
    invalid_ids = []

    for taxon_id in taxon_ids:
        # Basic validation - should be numeric
        if taxon_id.isdigit() and len(taxon_id) >= 1:
            valid_ids.append(taxon_id)
        else:
            invalid_ids.append(taxon_id)

    return valid_ids, invalid_ids


def get_proteome_metadata(taxon_id: str) -> Dict[str, Any]:
    """Get proteome metadata for a taxonomy ID using UniProt Proteomes API.

    Args:
        taxon_id: Taxonomy ID string

    Returns:
        Dictionary with proteome metadata

    Raises:
        requests.RequestException: If API request fails
        ValueError: If no proteome found for the taxonomy ID
    """
    import requests

    # UniProt Proteomes API endpoint
    base_url = "https://www.ebi.ac.uk/proteins/api/proteomes"

    # Search for proteome by taxonomy ID
    params = {"taxid": taxon_id, "size": 1}  # Get first result only

    try:
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()

        data = response.json()

        if not data:
            raise ValueError(f"No proteome found for taxonomy ID {taxon_id}")

        # Take first proteome result
        proteome = data[0]

        # Extract protein count from scores if available
        protein_count = 0
        scores = proteome.get("scores", {})
        if isinstance(scores, dict):
            protein_count = scores.get("proteinCount", 0)

        metadata = {
            "taxon_id": taxon_id,
            "scientific_name": proteome.get("name", "Unknown"),
            "proteome_id": proteome.get("upid"),
            "protein_count": protein_count,
            "busco_complete": (
                proteome.get("scores", {}).get("buscoComplete", 0.0)
                if isinstance(proteome.get("scores"), dict)
                else 0.0
            ),
            "annotation_level": (
                proteome.get("annotationScore", {}).get("category")
                if isinstance(proteome.get("annotationScore"), dict)
                else "unknown"
            ),
            "is_reference": proteome.get("isReferenceProteome", False),
            "genome_assembly": proteome.get("genomeAssembly", {}),
            "genome_annotation": proteome.get("genomeAnnotation", {}),
            "strain": proteome.get("strain"),
            "description": proteome.get("description"),
            "superregnum": proteome.get("superregnum"),
            "source": proteome.get("source"),
        }

        logger.info(f"Retrieved proteome metadata for taxon {taxon_id}: {metadata['proteome_id']}")
        return metadata

    except requests.RequestException as e:
        logger.error(f"Failed to fetch proteome metadata for taxon {taxon_id}: {e}")
        raise
    except (KeyError, IndexError) as e:
        logger.error(f"Unexpected API response format for taxon {taxon_id}: {e}")
        raise ValueError(f"Invalid API response for taxonomy ID {taxon_id}")


def download_proteome_fasta(taxon_id: str, output_path: Union[str, Path], include_isoforms: bool = False) -> bool:
    """Download proteome FASTA file for a taxonomy ID using UniProt API.

    Args:
        taxon_id: Taxonomy ID string
        output_path: Output file path
        include_isoforms: Whether to include protein isoforms

    Returns:
        True if download successful, False otherwise

    Raises:
        requests.RequestException: If API request fails
    """
    import requests

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Downloading proteome FASTA for taxon {taxon_id} to {output_path}")

    try:
        # First get proteome metadata to find the proteome ID
        metadata = get_proteome_metadata(taxon_id)
        proteome_id = metadata.get("proteome_id")

        if not proteome_id:
            logger.error(f"No proteome ID found for taxon {taxon_id}")
            return False

        # UniProt proteome download URL
        base_url = "https://rest.uniprot.org/uniprotkb/stream"
        params = {"query": f"proteome:{proteome_id}", "format": "fasta"}

        # Note: includeIsoform parameter doesn't seem to be supported in this endpoint
        # All canonical sequences are returned

        response = requests.get(base_url, params=params, timeout=60)
        response.raise_for_status()

        # Write FASTA content to file
        with io.open_text_auto(output_path, "w") as f:
            f.write(response.text)

        logger.info(f"Downloaded proteome FASTA for taxon {taxon_id} ({len(response.text)} characters)")
        return True

    except requests.RequestException as e:
        logger.error(f"Failed to download proteome FASTA for taxon {taxon_id}: {e}")
        raise
    except Exception as e:
        logger.error(f"Error processing proteome download for taxon {taxon_id}: {e}")
        return False


def proteome_statistics(fasta_path: Union[str, Path]) -> Dict[str, Any]:
    """Calculate statistics for a proteome FASTA file.

    Args:
        fasta_path: Path to proteome FASTA file

    Returns:
        Dictionary with proteome statistics
    """
    try:
        from metainformant.protein.sequence.sequences import molecular_weight, read_fasta

        sequences = read_fasta(fasta_path)

        stats = {
            "protein_count": len(sequences),
            "total_residues": 0,
            "average_length": 0.0,
            "molecular_weights": [],
        }

        for seq in sequences.values():
            length = len(seq)
            stats["total_residues"] += length
            stats["molecular_weights"].append(molecular_weight(seq))

        if sequences:
            stats["average_length"] = stats["total_residues"] / len(sequences)
            stats["average_molecular_weight"] = sum(stats["molecular_weights"]) / len(stats["molecular_weights"])

        logger.info(f"Calculated statistics for proteome with {stats['protein_count']} proteins")

        return stats

    except Exception as e:
        logger.error(f"Error calculating proteome statistics: {e}")
        return {}


def compare_proteomes(proteome1_path: Union[str, Path], proteome2_path: Union[str, Path]) -> Dict[str, Any]:
    """Compare two proteomes.

    Args:
        proteome1_path: Path to first proteome FASTA
        proteome2_path: Path to second proteome FASTA

    Returns:
        Dictionary with comparison results
    """
    stats1 = proteome_statistics(proteome1_path)
    stats2 = proteome_statistics(proteome2_path)

    comparison = {
        "proteome1": stats1,
        "proteome2": stats2,
        "protein_count_ratio": stats2.get("protein_count", 0) / max(stats1.get("protein_count", 1), 1),
        "size_difference": stats2.get("total_residues", 0) - stats1.get("total_residues", 0),
    }

    return comparison

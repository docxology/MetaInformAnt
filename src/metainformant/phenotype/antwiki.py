"""AntWiki data loading and processing.

This module provides functionality to load, validate, and process AntWiki
phenotype data from JSON files, including morphological traits, behavioral
characteristics, and taxonomic information.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Any, Optional, Iterator, Set, Tuple
from collections import defaultdict
import json

from metainformant.core import logging, errors, validation, io

logger = logging.get_logger(__name__)


class AntWikiRecord:
    """Represents a single AntWiki species record."""

    def __init__(self, data: Dict[str, Any]):
        """Initialize from AntWiki JSON data.

        Args:
            data: Raw AntWiki record data

        Raises:
            ValueError: If required fields are missing
        """
        self.raw_data = data
        self.species_name = data.get('species_name', '')
        self.genus = data.get('genus', '')
        self.subfamily = data.get('subfamily', '')
        self.tribe = data.get('tribe', '')

        # Validate required fields
        if not self.species_name:
            raise ValueError("AntWiki record missing species_name")
        if not self.genus:
            raise ValueError(f"AntWiki record for {self.species_name} missing genus")

        # Extract phenotype data
        self.phenotypes = self._extract_phenotypes(data)
        self.morphology = data.get('morphology', {})
        self.behavior = data.get('behavior', {})
        self.ecology = data.get('ecology', {})
        self.distribution = data.get('distribution', {})

        # Metadata
        self.last_updated = data.get('last_updated')
        self.data_source = data.get('data_source', 'antwiki')
        self.confidence_score = data.get('confidence_score', 1.0)

    def _extract_phenotypes(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Extract phenotype information from raw data."""
        phenotypes = {}

        # Size-related phenotypes
        if 'morphology' in data:
            morph = data['morphology']
            phenotypes.update({
                'body_length': morph.get('body_length_mm'),
                'head_width': morph.get('head_width_mm'),
                'worker_size': morph.get('worker_size_category'),
                'queen_size': morph.get('queen_size_category'),
            })

        # Color phenotypes
        if 'color' in data:
            color_data = data['color']
            phenotypes.update({
                'head_color': color_data.get('head'),
                'thorax_color': color_data.get('thorax'),
                'gaster_color': color_data.get('gaster'),
                'leg_color': color_data.get('legs'),
            })

        # Behavioral phenotypes
        if 'behavior' in data:
            behavior = data['behavior']
            phenotypes.update({
                'foraging_strategy': behavior.get('foraging_strategy'),
                'nest_type': behavior.get('nest_type'),
                'colony_size': behavior.get('colony_size_category'),
                'activity_pattern': behavior.get('activity_pattern'),
            })

        # Remove None values
        return {k: v for k, v in phenotypes.items() if v is not None}

    @property
    def full_species_name(self) -> str:
        """Get full scientific name."""
        return f"{self.genus} {self.species_name}"

    @property
    def taxonomic_path(self) -> List[str]:
        """Get taxonomic classification path."""
        path = []
        if self.tribe:
            path.append(self.tribe)
        if self.subfamily:
            path.append(self.subfamily)
        path.append(self.genus)
        path.append(self.species_name)
        return path

    def get_phenotype_value(self, phenotype_name: str) -> Any:
        """Get a specific phenotype value."""
        return self.phenotypes.get(phenotype_name)

    def has_phenotype(self, phenotype_name: str) -> bool:
        """Check if record has a specific phenotype."""
        return phenotype_name in self.phenotypes

    def get_morphological_traits(self) -> Dict[str, Any]:
        """Get all morphological traits."""
        return self.morphology

    def get_behavioral_traits(self) -> Dict[str, Any]:
        """Get all behavioral traits."""
        return self.behavior

    def to_dict(self) -> Dict[str, Any]:
        """Convert record to dictionary."""
        return {
            'species_name': self.species_name,
            'genus': self.genus,
            'subfamily': self.subfamily,
            'tribe': self.tribe,
            'phenotypes': self.phenotypes,
            'morphology': self.morphology,
            'behavior': self.behavior,
            'ecology': self.ecology,
            'distribution': self.distribution,
            'last_updated': self.last_updated,
            'data_source': self.data_source,
            'confidence_score': self.confidence_score,
        }


def load_antwiki_json(path: str | Path, validate: bool = True) -> List[AntWikiRecord]:
    """Load AntWiki data from JSON file.

    Args:
        path: Path to JSON file containing AntWiki data
        validate: Whether to validate records during loading

    Returns:
        List of AntWikiRecord objects

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If JSON is invalid or records fail validation
    """
    path = validation.validate_path_exists(Path(path))

    logger.info(f"Loading AntWiki data from {path}")

    try:
        # Load JSON data
        data = io.load_json(path)

        # Handle different JSON structures
        if isinstance(data, dict):
            # Single record or wrapped data
            if 'records' in data:
                records_data = data['records']
            elif 'species' in data:
                records_data = data['species']
            else:
                # Assume it's a single record
                records_data = [data]
        elif isinstance(data, list):
            # List of records
            records_data = data
        else:
            raise ValueError("Invalid AntWiki JSON structure")

        # Convert to AntWikiRecord objects
        records = []
        for i, record_data in enumerate(records_data):
            try:
                record = AntWikiRecord(record_data)
                if validate:
                    _validate_antwiki_record(record)
                records.append(record)
            except Exception as e:
                if validate:
                    logger.warning(f"Skipping invalid record {i}: {e}")
                    continue
                else:
                    raise ValueError(f"Invalid record {i}: {e}") from e

        logger.info(f"Successfully loaded {len(records)} AntWiki records")
        return records

    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON in {path}: {e}")
        raise ValueError(f"Invalid JSON format: {e}") from e
    except Exception as e:
        logger.error(f"Error loading AntWiki data from {path}: {e}")
        raise errors.FileIOError(f"Failed to load AntWiki data: {e}") from e


def _validate_antwiki_record(record: AntWikiRecord) -> None:
    """Validate an AntWikiRecord object.

    Args:
        record: Record to validate

    Raises:
        ValueError: If validation fails
    """
    # Check required fields
    if not record.species_name:
        raise ValueError("Missing species_name")

    if not record.genus:
        raise ValueError("Missing genus")

    # Check taxonomic consistency
    if record.genus and not record.genus[0].isupper():
        raise ValueError(f"Genus should start with capital letter: {record.genus}")

    # Check phenotype data structure
    if not isinstance(record.phenotypes, dict):
        raise ValueError("Phenotypes should be a dictionary")

    # Check confidence score range
    if not (0.0 <= record.confidence_score <= 1.0):
        raise ValueError(f"Confidence score out of range: {record.confidence_score}")


def save_antwiki_json(records: List[AntWikiRecord], path: str | Path) -> None:
    """Save AntWiki records to JSON file.

    Args:
        path: Output path for JSON file
        records: List of AntWikiRecord objects to save
    """
    path = Path(path)
    io.ensure_directory(path.parent)

    # Convert records to dictionaries
    data = {
        'metadata': {
            'total_records': len(records),
            'export_timestamp': str(path.stat().st_mtime) if path.exists() else None,
            'data_source': 'antwiki',
        },
        'records': [record.to_dict() for record in records]
    }

    io.dump_json(data, path, indent=2)
    logger.info(f"Saved {len(records)} AntWiki records to {path}")


def filter_antwiki_records(records: List[AntWikiRecord],
                          genus: Optional[str] = None,
                          subfamily: Optional[str] = None,
                          min_confidence: float = 0.0,
                          required_phenotypes: Optional[List[str]] = None) -> List[AntWikiRecord]:
    """Filter AntWiki records based on criteria.

    Args:
        records: List of records to filter
        genus: Filter by genus name
        subfamily: Filter by subfamily name
        min_confidence: Minimum confidence score
        required_phenotypes: List of phenotype names that must be present

    Returns:
        Filtered list of records
    """
    filtered = []

    for record in records:
        # Apply filters
        if genus and record.genus != genus:
            continue

        if subfamily and record.subfamily != subfamily:
            continue

        if record.confidence_score < min_confidence:
            continue

        if required_phenotypes:
            has_all_phenotypes = all(record.has_phenotype(p) for p in required_phenotypes)
            if not has_all_phenotypes:
                continue

        filtered.append(record)

    logger.info(f"Filtered {len(records)} records to {len(filtered)} based on criteria")
    return filtered


def get_phenotype_distribution(records: List[AntWikiRecord],
                              phenotype_name: str) -> Dict[str, Any]:
    """Get distribution statistics for a phenotype across records.

    Args:
        records: List of AntWiki records
        phenotype_name: Name of phenotype to analyze

    Returns:
        Dictionary with distribution statistics
    """
    values = []
    value_counts = defaultdict(int)

    for record in records:
        value = record.get_phenotype_value(phenotype_name)
        if value is not None:
            values.append(value)
            value_counts[str(value)] += 1

    if not values:
        return {"phenotype": phenotype_name, "total_records": len(records), "values_found": 0}

    # Calculate statistics
    stats = {
        "phenotype": phenotype_name,
        "total_records": len(records),
        "values_found": len(values),
        "coverage": len(values) / len(records),
        "unique_values": len(set(str(v) for v in values)),
        "value_counts": dict(value_counts),
    }

    # Add numeric statistics if applicable
    numeric_values = [v for v in values if isinstance(v, (int, float))]
    if numeric_values:
        stats.update({
            "mean": sum(numeric_values) / len(numeric_values),
            "min": min(numeric_values),
            "max": max(numeric_values),
            "median": sorted(numeric_values)[len(numeric_values) // 2],
        })

    return stats


def find_similar_species(records: List[AntWikiRecord],
                        target_record: AntWikiRecord,
                        phenotype_weights: Optional[Dict[str, float]] = None,
                        top_k: int = 10) -> List[Tuple[AntWikiRecord, float]]:
    """Find species similar to a target record based on phenotypes.

    Args:
        records: List of all AntWiki records
        target_record: Target record to find similar species for
        phenotype_weights: Optional weights for different phenotypes
        top_k: Number of most similar records to return

    Returns:
        List of (record, similarity_score) tuples, sorted by similarity
    """
    if phenotype_weights is None:
        phenotype_weights = {
            'body_length': 0.3,
            'head_width': 0.3,
            'worker_size': 0.2,
            'foraging_strategy': 0.1,
            'nest_type': 0.1,
        }

    similarities = []

    for record in records:
        if record.species_name == target_record.species_name:
            continue  # Skip self

        similarity = 0.0
        total_weight = 0.0

        for phenotype, weight in phenotype_weights.items():
            target_value = target_record.get_phenotype_value(phenotype)
            record_value = record.get_phenotype_value(phenotype)

            if target_value is not None and record_value is not None:
                if target_value == record_value:
                    similarity += weight
                total_weight += weight

        if total_weight > 0:
            normalized_similarity = similarity / total_weight
            similarities.append((record, normalized_similarity))

    # Sort by similarity (descending) and return top_k
    similarities.sort(key=lambda x: x[1], reverse=True)
    return similarities[:top_k]


def create_phenotype_matrix(records: List[AntWikiRecord],
                           phenotype_names: Optional[List[str]] = None) -> Tuple[List[str], List[str], List[List[Any]]]:
    """Create a phenotype matrix for statistical analysis.

    Args:
        records: List of AntWiki records
        phenotype_names: Specific phenotypes to include (None for all)

    Returns:
        Tuple of (species_names, phenotype_names, matrix)
    """
    if not records:
        return [], [], []

    # Get all phenotype names if not specified
    if phenotype_names is None:
        all_phenotypes = set()
        for record in records:
            all_phenotypes.update(record.phenotypes.keys())
        phenotype_names = sorted(all_phenotypes)

    species_names = [record.full_species_name for record in records]
    matrix = []

    for record in records:
        row = []
        for phenotype in phenotype_names:
            value = record.get_phenotype_value(phenotype)
            row.append(value)
        matrix.append(row)

    return species_names, phenotype_names, matrix


def generate_antwiki_report(records: List[AntWikiRecord],
                           output_path: Optional[str | Path] = None) -> str:
    """Generate a summary report of AntWiki data.

    Args:
        output_path: Optional path to save the report
        records: List of AntWiki records to summarize

    Returns:
        Formatted report string
    """
    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append("ANTWIKI DATA SUMMARY REPORT")
    report_lines.append("=" * 60)
    report_lines.append("")

    # Basic statistics
    report_lines.append(f"Total Records: {len(records)}")
    report_lines.append(f"Unique Genera: {len(set(r.genus for r in records))}")
    report_lines.append(f"Unique Subfamilies: {len(set(r.subfamily for r in records if r.subfamily))}")
    report_lines.append("")

    # Genus distribution
    genus_counts = defaultdict(int)
    for record in records:
        genus_counts[record.genus] += 1

    report_lines.append("Top Genera:")
    for genus, count in sorted(genus_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
        report_lines.append(f"  {genus}: {count} species")
    report_lines.append("")

    # Phenotype coverage
    all_phenotypes = set()
    for record in records:
        all_phenotypes.update(record.phenotypes.keys())

    report_lines.append(f"Total Phenotypes: {len(all_phenotypes)}")
    report_lines.append("Phenotype Coverage:")

    for phenotype in sorted(all_phenotypes):
        count = sum(1 for r in records if r.has_phenotype(phenotype))
        coverage = count / len(records) * 100
        report_lines.append(f"  {phenotype}: {count}/{len(records)} ({coverage:.1f}%)")

    report_lines.append("")

    # Quality metrics
    confidence_scores = [r.confidence_score for r in records if r.confidence_score is not None]
    if confidence_scores:
        avg_confidence = sum(confidence_scores) / len(confidence_scores)
        report_lines.append(f"Average Confidence Score: {avg_confidence:.2f}")
        report_lines.append(f"High Confidence Records (>0.8): {sum(1 for s in confidence_scores if s > 0.8)}")

    report = "\n".join(report_lines)

    if output_path:
        output_path = Path(output_path)
        with open(output_path, 'w') as f:
            f.write(report)
        logger.info(f"AntWiki report saved to {output_path}")

    return report

"""Sample selection step for RNA-seq workflow.

This step selects appropriate SRA samples based on quality criteria,
experimental design, and biological relevance for downstream analysis.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict, List

from metainformant.core import io, logging
from metainformant.rna.steps import StepResult

logger = logging.get_logger(__name__)


def run_step(step_params: Dict[str, Any]) -> StepResult:
    """Execute the sample selection step.

    Args:
        step_params: Parameters for the selection step
            - work_dir: Working directory
            - integrated_metadata: Path to integrated metadata file
            - species_list: List of species to process
            - selection_criteria: Criteria for sample selection
            - min_reads: Minimum number of reads required
            - max_samples_per_species: Maximum samples per species
            - preferred_platforms: Preferred sequencing platforms

    Returns:
        StepResult with selection results
    """
    start_time = time.time()

    try:
        work_dir = Path(step_params.get('work_dir', '.'))
        integrated_metadata_file = step_params.get('integrated_metadata',
                                                 work_dir / "integration" / "integrated_metadata.json")
        species_list = step_params.get('species_list', [])
        selection_criteria = step_params.get('selection_criteria', {})
        min_reads = step_params.get('min_reads', 1000000)
        max_samples_per_species = step_params.get('max_samples_per_species', 50)
        preferred_platforms = step_params.get('preferred_platforms', ['ILLUMINA'])

        # Load integrated metadata
        if isinstance(integrated_metadata_file, str):
            integrated_metadata_file = Path(integrated_metadata_file)

        if not integrated_metadata_file.exists():
            raise FileNotFoundError(f"Integrated metadata file not found: {integrated_metadata_file}")

        integrated_metadata = io.load_json(integrated_metadata_file)

        # Create selection directory
        selection_dir = work_dir / "selection"
        selection_dir.mkdir(parents=True, exist_ok=True)

        selected_samples = {}

        for species in species_list:
            if species not in integrated_metadata:
                logger.warning(f"Species {species} not found in integrated metadata")
                continue

            logger.info(f"Selecting samples for {species}")

            species_samples = integrated_metadata[species]
            selected = select_species_samples(
                species_samples,
                selection_criteria,
                min_reads,
                max_samples_per_species,
                preferred_platforms
            )

            selected_samples[species] = selected

            logger.info(f"Selected {len(selected)} samples for {species}")

        # Save selection results
        selection_file = selection_dir / "selected_samples.json"
        io.dump_json(selected_samples, selection_file)

        # Create sample list for downstream steps
        sample_list = create_selected_sample_list(selected_samples)
        sample_list_file = selection_dir / "selected_sample_list.txt"

        with open(sample_list_file, 'w') as f:
            for sample in sample_list:
                f.write(f"{sample}\n")

        # Generate selection summary
        summary = create_selection_summary(selected_samples, integrated_metadata)
        summary_file = selection_dir / "selection_summary.json"
        io.dump_json(summary, summary_file)

        execution_time = time.time() - start_time

        return StepResult(
            success=True,
            outputs={
                'selection_dir': str(selection_dir),
                'selected_samples': str(selection_file),
                'sample_list': str(sample_list_file),
                'summary': str(summary_file),
                'total_selected': sum(len(samples) for samples in selected_samples.values())
            },
            metadata={
                'selection_criteria': selection_criteria,
                'min_reads': min_reads,
                'max_samples_per_species': max_samples_per_species,
                'preferred_platforms': preferred_platforms,
                'execution_time': execution_time
            },
            execution_time=execution_time
        )

    except Exception as e:
        execution_time = time.time() - start_time
        logger.error(f"Selection step failed: {e}")

        return StepResult(
            success=False,
            outputs={},
            metadata={},
            error_message=str(e),
            execution_time=execution_time
        )


def select_species_samples(species_samples: Dict[str, Any],
                          selection_criteria: Dict[str, Any],
                          min_reads: int, max_samples: int,
                          preferred_platforms: List[str]) -> Dict[str, Any]:
    """Select samples for a specific species based on criteria.

    Args:
        species_samples: Sample metadata for the species
        selection_criteria: Selection criteria dictionary
        min_reads: Minimum read count
        max_samples: Maximum samples to select
        preferred_platforms: Preferred sequencing platforms

    Returns:
        Dictionary of selected samples
    """
    # Score and rank samples
    scored_samples = []

    for sample_id, sample_info in species_samples.items():
        score = calculate_sample_score(
            sample_info,
            selection_criteria,
            min_reads,
            preferred_platforms
        )

        scored_samples.append((sample_id, sample_info, score))

    # Sort by score (descending)
    scored_samples.sort(key=lambda x: x[2], reverse=True)

    # Select top samples
    selected = {}
    for i, (sample_id, sample_info, score) in enumerate(scored_samples):
        if i >= max_samples:
            break

        selected[sample_id] = {
            'original_info': sample_info,
            'selection_score': score,
            'rank': i + 1
        }

    return selected


def calculate_sample_score(sample_info: Dict[str, Any],
                          criteria: Dict[str, Any],
                          min_reads: int,
                          preferred_platforms: List[str]) -> float:
    """Calculate selection score for a sample.

    Args:
        sample_info: Sample metadata
        criteria: Selection criteria
        min_reads: Minimum read count
        preferred_platforms: Preferred platforms

    Returns:
        Selection score (higher is better)
    """
    score = 0.0

    original_metadata = sample_info.get('original_metadata', {})

    # Read count score
    spots = original_metadata.get('spots', 0)
    if spots >= min_reads:
        score += min(spots / 1000000, 10)  # Cap at 10 points for 1M+ reads
    else:
        return 0.0  # Reject samples with too few reads

    # Platform preference score
    platform = original_metadata.get('platform', '').upper()
    if platform in [p.upper() for p in preferred_platforms]:
        score += 5.0
    elif platform:
        score += 2.0  # Some preference for known platforms

    # Library layout preference (paired-end preferred)
    layout = original_metadata.get('library_layout', '').upper()
    if layout == 'PAIRED':
        score += 3.0
    elif layout == 'SINGLE':
        score += 1.0

    # Data completeness score
    completeness = sample_info.get('validation', {}).get('completeness_score', 0.0)
    score += completeness * 2.0

    # Recency preference (newer data preferred, but not strongly)
    # Could add date-based scoring here

    return score


def create_selected_sample_list(selected_samples: Dict[str, Dict[str, Any]]) -> List[str]:
    """Create a list of selected sample accessions.

    Args:
        selected_samples: Selected samples dictionary

    Returns:
        List of sample accessions
    """
    sample_list = []

    for species_samples in selected_samples.values():
        for sample_info in species_samples.values():
            accession = sample_info.get('original_info', {}).get('accession', '')
            if accession:
                sample_list.append(accession)

    return sample_list


def create_selection_summary(selected_samples: Dict[str, Dict[str, Any]],
                           original_metadata: Dict[str, Any]) -> Dict[str, Any]:
    """Create a summary of selection results.

    Args:
        selected_samples: Selected samples
        original_metadata: Original metadata for comparison

    Returns:
        Selection summary
    """
    summary = {
        'total_species': len(selected_samples),
        'total_selected_samples': sum(len(samples) for samples in selected_samples.values()),
        'selection_efficiency': 0.0,
        'average_score': 0.0,
        'platform_distribution': {},
        'species_breakdown': {}
    }

    total_original = sum(len(samples) for samples in original_metadata.values())
    total_selected = summary['total_selected_samples']

    if total_original > 0:
        summary['selection_efficiency'] = total_selected / total_original

    total_score = 0.0

    for species, samples in selected_samples.items():
        species_summary = {
            'selected_count': len(samples),
            'original_count': len(original_metadata.get(species, {})),
            'average_score': 0.0,
            'platform_counts': {}
        }

        species_scores = []

        for sample_info in samples.values():
            score = sample_info.get('selection_score', 0.0)
            species_scores.append(score)
            total_score += score

            # Platform distribution
            platform = sample_info.get('original_info', {}).get('original_metadata', {}).get('platform', 'Unknown')
            summary['platform_distribution'][platform] = summary['platform_distribution'].get(platform, 0) + 1
            species_summary['platform_counts'][platform] = species_summary['platform_counts'].get(platform, 0) + 1

        if species_scores:
            species_summary['average_score'] = sum(species_scores) / len(species_scores)

        summary['species_breakdown'][species] = species_summary

    if total_selected > 0:
        summary['average_score'] = total_score / total_selected

    return summary


def validate_selection(selected_samples: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Validate selected samples for downstream compatibility.

    Args:
        selected_samples: Selected samples dictionary

    Returns:
        Validation results
    """
    validation = {
        'is_valid': True,
        'issues': [],
        'warnings': []
    }

    for species, samples in selected_samples.items():
        if not samples:
            validation['issues'].append(f"No samples selected for {species}")
            validation['is_valid'] = False
            continue

        for sample_id, sample_info in samples.items():
            accession = sample_info.get('original_info', {}).get('accession')

            if not accession:
                validation['issues'].append(f"Sample {sample_id} missing accession")
                validation['is_valid'] = False

            # Check for data availability
            has_local = sample_info.get('original_info', {}).get('local_fastq_files')
            has_sra = sample_info.get('original_info', {}).get('sra_available')

            if not has_local and not has_sra:
                validation['warnings'].append(f"Sample {sample_id} has no data source")

    return validation



run_select = run_step






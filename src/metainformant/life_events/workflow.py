"""Workflow functions for life events analysis.

This module provides high-level workflow functions for analyzing life event sequences,
training models, and generating predictions.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from metainformant.core import io, logging

logger = logging.get_logger(__name__)

# Optional imports for ML functionality
try:
    from . import embeddings
    from . import models
    from . import config
    EMBEDDINGS_AVAILABLE = True
except ImportError:
    EMBEDDINGS_AVAILABLE = False
    logger.warning("Life events ML components not available")


def analyze_life_course(
    sequences: List[Any],
    outcomes: Optional[List[str]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    config_obj: Optional[Any] = None,
    **kwargs: Any
) -> Dict[str, Any]:
    """Analyze life event sequences with optional outcome prediction.

    This function provides a complete workflow for life course analysis including:
    - Sequence statistics and summary
    - Optional embedding learning and model training
    - Optional outcome prediction and evaluation

    Args:
        sequences: List of EventSequence objects
        outcomes: Optional list of outcome labels for supervised learning
        output_dir: Directory to save results (defaults to output/life_events/)
        config_obj: Optional LifeEventsWorkflowConfig object
        **kwargs: Additional configuration options

    Returns:
        Dictionary containing analysis results
    """
    if not sequences:
        raise ValueError("No sequences provided for analysis")

    # Set up output directory
    if output_dir is None:
        output_dir = Path("output") / "life_events"
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    results = {
        "n_sequences": len(sequences),
        "output_dir": str(output_dir),
        "config_used": {},
        "sequence_stats": {},
        "model_results": {},
        "predictions": {},
        "visualizations": []
    }

    # Basic sequence statistics
    results["sequence_stats"] = _compute_sequence_stats(sequences)

    # If outcomes provided, do supervised analysis
    if outcomes is not None:
        if len(outcomes) != len(sequences):
            raise ValueError("Number of outcomes must match number of sequences")

        if not EMBEDDINGS_AVAILABLE:
            logger.warning("ML components not available, skipping supervised analysis")
            results["model_results"] = {"error": "ML components not available"}
        else:
            # Learn embeddings
            logger.info("Learning event embeddings...")
            embedding_results = embeddings.learn_event_embeddings(
                sequences, output_dir=output_dir, **kwargs
            )

            # Train prediction model
            logger.info("Training prediction model...")
            model_results = models.train_event_predictor(
                sequences, outcomes, embedding_results=embedding_results,
                output_dir=output_dir, **kwargs
            )

            # Generate predictions
            logger.info("Generating predictions...")
            predictions = models.predict_outcomes(
                sequences, model_results["model"], embedding_results=embedding_results
            )

            results["model_results"] = model_results
            results["predictions"] = predictions

            # Save model if requested
            if kwargs.get("save_model", True):
                model_path = output_dir / "trained_model.json"
                models.save_model(model_results["model"], model_path)
                results["model_path"] = str(model_path)

    # Generate visualizations if matplotlib available
    try:
        from . import visualization
        vis_results = _generate_analysis_visualizations(
            sequences, outcomes, output_dir, results
        )
        results["visualizations"] = vis_results
    except ImportError:
        logger.warning("Visualization not available")
        results["visualizations"] = []

    # Save results
    results_file = output_dir / "analysis_results.json"
    io.dump_json(results, results_file)

    logger.info(f"Analysis complete. Results saved to {results_file}")
    return results


def _compute_sequence_stats(sequences: List[Any]) -> Dict[str, Any]:
    """Compute basic statistics for event sequences."""
    if not sequences:
        return {}

    stats = {
        "total_sequences": len(sequences),
        "total_events": sum(len(seq.events) for seq in sequences),
        "avg_events_per_sequence": sum(len(seq.events) for seq in sequences) / len(sequences),
        "sequence_lengths": [len(seq.events) for seq in sequences],
        "unique_event_types": set(),
        "unique_domains": set(),
    }

    # Collect unique event types and domains
    for seq in sequences:
        for event in seq.events:
            stats["unique_event_types"].add(event.event_type)
            if hasattr(event, 'domain'):
                stats["unique_domains"].add(event.domain)

    stats["unique_event_types"] = list(stats["unique_event_types"])
    stats["unique_domains"] = list(stats["unique_domains"])
    stats["n_unique_event_types"] = len(stats["unique_event_types"])
    stats["n_unique_domains"] = len(stats["unique_domains"])

    return stats


def _generate_analysis_visualizations(
    sequences: List[Any],
    outcomes: Optional[List[str]],
    output_dir: Path,
    results: Dict[str, Any]
) -> List[Dict[str, str]]:
    """Generate analysis visualizations."""
    visualizations = []

    try:
        from . import visualization

        # Sequence length distribution
        fig = visualization.plot_sequence_length_distribution(sequences)
        if fig:
            length_plot = output_dir / "sequence_lengths.png"
            fig.savefig(length_plot, dpi=300, bbox_inches='tight')
            visualizations.append({
                "type": "sequence_length_distribution",
                "file": str(length_plot)
            })

        # Domain distribution if domains available
        if any(hasattr(event, 'domain') for seq in sequences for event in seq.events):
            fig = visualization.plot_domain_distribution(sequences)
            if fig:
                domain_plot = output_dir / "domain_distribution.png"
                fig.savefig(domain_plot, dpi=300, bbox_inches='tight')
                visualizations.append({
                    "type": "domain_distribution",
                    "file": str(domain_plot)
                })

        # Outcome distribution if outcomes available
        if outcomes:
            fig = visualization.plot_outcome_distribution(outcomes)
            if fig:
                outcome_plot = output_dir / "outcome_distribution.png"
                fig.savefig(outcome_plot, dpi=300, bbox_inches='tight')
                visualizations.append({
                    "type": "outcome_distribution",
                    "file": str(outcome_plot)
                })

    except Exception as e:
        logger.warning(f"Error generating visualizations: {e}")

    return visualizations


def compare_populations(
    sequences1: List[Any],
    sequences2: List[Any],
    group_names: tuple[str, str] = ("Group1", "Group2"),
    output_dir: Optional[Union[str, Path]] = None,
    **kwargs: Any
) -> Dict[str, Any]:
    """Compare two populations of life event sequences.

    Args:
        sequences1: Sequences from first population
        sequences2: Sequences from second population
        group_names: Names for the two groups
        output_dir: Directory to save results
        **kwargs: Additional analysis options

    Returns:
        Dictionary containing comparison results
    """
    if not sequences1 or not sequences2:
        raise ValueError("Both sequence groups must be non-empty")

    # Set up output directory
    if output_dir is None:
        output_dir = Path("output") / "life_events" / "comparisons"
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    results = {
        "group1_name": group_names[0],
        "group2_name": group_names[1],
        "group1_stats": _compute_sequence_stats(sequences1),
        "group2_stats": _compute_sequence_stats(sequences2),
        "comparison": {},
        "visualizations": []
    }

    # Statistical comparisons
    results["comparison"] = _compare_sequence_groups(sequences1, sequences2, group_names)

    # Generate comparison visualizations
    try:
        from . import visualization
        fig = visualization.plot_population_comparison(
            sequences1, sequences2, group_names=group_names
        )
        if fig:
            comp_plot = output_dir / "population_comparison.png"
            fig.savefig(comp_plot, dpi=300, bbox_inches='tight')
            results["visualizations"].append({
                "type": "population_comparison",
                "file": str(comp_plot)
            })
    except Exception as e:
        logger.warning(f"Error generating comparison visualizations: {e}")

    # Save results
    results_file = output_dir / "comparison_results.json"
    io.dump_json(results, results_file)

    logger.info(f"Population comparison complete. Results saved to {results_file}")
    return results


def _compare_sequence_groups(
    sequences1: List[Any],
    sequences2: List[Any],
    group_names: tuple[str, str]
) -> Dict[str, Any]:
    """Compare two groups of sequences statistically."""
    # Basic statistical tests
    lengths1 = [len(seq.events) for seq in sequences1]
    lengths2 = [len(seq.events) for seq in sequences2]

    comparison = {
        "length_comparison": {
            "group1_mean": sum(lengths1) / len(lengths1) if lengths1 else 0,
            "group2_mean": sum(lengths2) / len(lengths2) if lengths2 else 0,
            "group1_median": sorted(lengths1)[len(lengths1)//2] if lengths1 else 0,
            "group2_median": sorted(lengths2)[len(lengths2)//2] if lengths2 else 0,
        }
    }

    # Event type overlap
    types1 = set()
    types2 = set()

    for seq in sequences1:
        for event in seq.events:
            types1.add(event.event_type)

    for seq in sequences2:
        for event in seq.events:
            types2.add(event.event_type)

    comparison["event_types"] = {
        "group1_unique": len(types1),
        "group2_unique": len(types2),
        "overlap": len(types1 & types2),
        "jaccard_similarity": len(types1 & types2) / len(types1 | types2) if (types1 | types2) else 0
    }

    return comparison

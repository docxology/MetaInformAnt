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

def compare_populations(group1: List[EventSequence], group2: List[EventSequence],
                      output_dir: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """Compare two populations of life event sequences.

    Args:
        group1: First group of EventSequence objects
        group2: Second group of EventSequence objects
        output_dir: Directory to save comparison results (optional)

    Returns:
        Dictionary containing comparison statistics and plots
    """
    from pathlib import Path

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Analyze group 1
    group1_stats = analyze_life_course(group1, outcomes=None)
    group1_stats["n_sequences"] = len(group1)

    # Analyze group 2
    group2_stats = analyze_life_course(group2, outcomes=None)
    group2_stats["n_sequences"] = len(group2)

    # Compute comparison metrics
    comparison = {
        "sequence_length": {
            "group1_mean": np.mean([len(seq.events) for seq in group1]),
            "group2_mean": np.mean([len(seq.events) for seq in group2]),
            "difference": abs(np.mean([len(seq.events) for seq in group1]) -
                              np.mean([len(seq.events) for seq in group2]))
        },
        "event_diversity": {
            "group1_unique_events": len(set(event.event_type
                                            for seq in group1
                                            for event in seq.events)),
            "group2_unique_events": len(set(event.event_type
                                            for seq in group2
                                            for event in seq.events))
        },
        "temporal_range": {
            "group1_days": (max(event.timestamp for seq in group1 for event in seq.events) -
                           min(event.timestamp for seq in group1 for event in seq.events)).days,
            "group2_days": (max(event.timestamp for seq in group2 for event in seq.events) -
                           min(event.timestamp for seq in group2 for event in seq.events)).days
        }
    }

    # Generate comparison plots if output directory provided
    if output_dir:
        try:
            from .visualization import plot_population_comparison
            plot_path = output_dir / "population_comparison.png"
            plot_population_comparison(group1, group2,
                                     group1_label="Group 1", group2_label="Group 2",
                                     output_path=plot_path)
        except ImportError:
            logger.warning("Population comparison plot not available")

    return {
        "group1": group1_stats,
        "group2": group2_stats,
        "comparison": comparison,
        "output_dir": str(output_dir) if output_dir else None
    }

def intervention_analysis(sequences: List[EventSequence],
                        intervention_time: float,
                        output_dir: Optional[Union[str, Path]] = None,
                        outcomes: Optional[List[Any]] = None) -> Dict[str, Any]:
    """Analyze the effects of an intervention on life event sequences.

    Args:
        sequences: List of EventSequence objects
        intervention_time: Timestamp of the intervention
        output_dir: Directory to save analysis results (optional)
        outcomes: Optional list of outcomes to compare before/after intervention

    Returns:
        Dictionary containing pre/post intervention analysis
    """
    from pathlib import Path
    import datetime

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Convert timestamp to datetime if needed
    if isinstance(intervention_time, (int, float)):
        intervention_datetime = datetime.datetime.fromtimestamp(intervention_time)
    else:
        intervention_datetime = intervention_time

    # Split sequences into pre and post intervention
    pre_sequences = []
    post_sequences = []

    for seq in sequences:
        pre_events = [event for event in seq.events if event.timestamp < intervention_datetime]
        post_events = [event for event in seq.events if event.timestamp >= intervention_datetime]

        if pre_events:
            pre_sequences.append(EventSequence(person_id=f"{seq.person_id}_pre", events=pre_events))
        if post_events:
            post_sequences.append(EventSequence(person_id=f"{seq.person_id}_post", events=post_events))

    # Analyze pre-intervention
    pre_analysis = analyze_life_course(pre_sequences, outcomes[:len(pre_sequences)] if outcomes else None)

    # Analyze post-intervention
    post_analysis = analyze_life_course(post_sequences, outcomes[len(pre_sequences):] if outcomes and len(outcomes) > len(pre_sequences) else None)

    # Calculate differences
    differences = {}
    for key in pre_analysis.keys():
        if key in post_analysis and isinstance(pre_analysis[key], (int, float)):
            if isinstance(post_analysis[key], (int, float)):
                differences[f"{key}_change"] = post_analysis[key] - pre_analysis[key]

    # Generate intervention effect plots if output directory provided
    if output_dir:
        try:
            from .visualization import plot_intervention_effects
            plot_path = output_dir / "intervention_effects.png"
            plot_intervention_effects(pre_sequences, post_sequences,
                                    intervention_time=intervention_datetime,
                                    output_path=plot_path)
        except ImportError:
            logger.warning("Intervention effects plot not available")

    return {
        "pre_intervention": pre_analysis,
        "post_intervention": post_analysis,
        "differences": differences,
        "intervention_time": intervention_time,
        "n_pre_sequences": len(pre_sequences),
        "n_post_sequences": len(post_sequences),
        "output_dir": str(output_dir) if output_dir else None
    }

def event_importance(sequences: List[EventSequence], method: str = 'frequency',
                    normalize: bool = True) -> Dict[str, float]:
    """Calculate event importance scores using various methods.

    Args:
        sequences: List of EventSequence objects
        method: Method for calculating importance ('frequency', 'temporal', 'transition')
        normalize: Whether to normalize scores to [0,1] range

    Returns:
        Dictionary mapping event types to importance scores

    Examples:
        >>> sequences = [EventSequence(...), ...]
        >>> importance = event_importance(sequences, method='frequency')
        >>> print(f"Most important event: {max(importance, key=importance.get)}")
    """
    if not sequences:
        return {}

    event_counts = {}
    event_positions = {}
    event_transitions = {}

    for seq in sequences:
        for i, event in enumerate(seq.events):
            event_type = event.event_type

            # Count frequency
            event_counts[event_type] = event_counts.get(event_type, 0) + 1

            # Track positions
            if event_type not in event_positions:
                event_positions[event_type] = []
            event_positions[event_type].append(i / max(1, len(seq.events) - 1))  # Normalized position

            # Track transitions
            if i > 0:
                prev_event = seq.events[i-1].event_type
                transition_key = f"{prev_event}->{event_type}"
                event_transitions[transition_key] = event_transitions.get(transition_key, 0) + 1

    if method == 'frequency':
        # Simple frequency-based importance
        importance_scores = event_counts.copy()

    elif method == 'temporal':
        # Importance based on temporal positioning
        importance_scores = {}
        for event_type, positions in event_positions.items():
            # Events that appear early get higher importance
            avg_position = sum(positions) / len(positions)
            importance_scores[event_type] = 1.0 - avg_position  # Earlier = more important

    elif method == 'transition':
        # Importance based on being transition hubs
        importance_scores = {}
        for event_type in event_counts.keys():
            # Count incoming and outgoing transitions
            incoming = sum(count for key, count in event_transitions.items()
                          if key.endswith(f"->{event_type}"))
            outgoing = sum(count for key, count in event_transitions.items()
                          if key.startswith(f"{event_type}->"))
            importance_scores[event_type] = incoming + outgoing

    if normalize and importance_scores:
        max_score = max(importance_scores.values())
        min_score = min(importance_scores.values())

        if max_score > min_score:
            for event_type in importance_scores:
                importance_scores[event_type] = (importance_scores[event_type] - min_score) / (max_score - min_score)
        else:
            # All scores are the same, set to 0.5
            for event_type in importance_scores:
                importance_scores[event_type] = 0.5

    return importance_scores

"""Workflow functions for life events analysis.

This module provides high-level workflow functions for analyzing life event sequences,
training models, and generating predictions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from metainformant.core import io
from metainformant.core.utils import logging
from metainformant.life_events.core.config import LifeEventsWorkflowConfig, load_life_events_config
from metainformant.life_events.core.events import EventSequence
from metainformant.life_events.core.utils import convert_sequences_to_tokens, get_event_statistics
from metainformant.life_events.models.embeddings import learn_event_embeddings
from metainformant.life_events.models.predictor import EventSequencePredictor

logger = logging.get_logger(__name__)


def analyze_life_course(
    sequences: List[Any],
    outcomes: Optional[List[str]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    config_obj: Optional[Any] = None,
    config_path: Optional[Union[str, Path]] = None,
    **kwargs: Any,
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

    config = config_obj
    if config is None and config_path is not None:
        config = load_life_events_config(config_path)

    embedding_config: dict[str, Any] = {}
    model_config: dict[str, Any] = {}
    workflow_config: dict[str, Any] = {}
    if isinstance(config, LifeEventsWorkflowConfig):
        embedding_config = dict(config.embedding)
        model_config = dict(config.model)
        workflow_config = dict(config.workflow)
        if output_dir is None:
            output_dir = config.work_dir
    elif isinstance(config, dict):
        embedding_config = dict(config.get("embedding", {}))
        model_config = dict(config.get("model", {}))
        workflow_config = dict(config.get("workflow", {}))
        for key in ("embedding_dim", "window_size", "epochs", "learning_rate"):
            if key in config:
                embedding_config.setdefault(key, config[key])
        for key in ("model_type", "task_type", "random_state"):
            if key in config:
                model_config.setdefault(key, config[key])

    if output_dir is None:
        output_dir = Path("output") / "life_events"
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    embedding_dim = int(kwargs.get("embedding_dim", embedding_config.get("embedding_dim", 100)))
    window_size = int(kwargs.get("window_size", embedding_config.get("window_size", 5)))
    epochs = int(kwargs.get("epochs", embedding_config.get("epochs", 5)))
    random_state = int(kwargs.get("random_state", model_config.get("random_state", 42)))
    model_type = str(kwargs.get("model_type", model_config.get("model_type", "embedding")))
    task_type = str(kwargs.get("task_type", model_config.get("task_type", "classification")))

    token_sequences = convert_sequences_to_tokens(sequences)
    embeddings_result = learn_event_embeddings(
        token_sequences,
        embedding_dim=embedding_dim,
        window_size=window_size,
        epochs=epochs,
        random_state=random_state,
    )

    embeddings_path = output_dir / "event_embeddings.json"
    io.dump_json({event: vector.tolist() for event, vector in embeddings_result.items()}, embeddings_path)

    stats = get_event_statistics(sequences)

    results = {
        "n_sequences": len(sequences),
        "n_events": stats.get("total_events", 0),
        "output_dir": str(output_dir),
        "config_used": {
            "embedding": embedding_config,
            "model": model_config,
            "workflow": workflow_config,
        },
        "sequence_stats": stats,
        "embeddings": str(embeddings_path),
        "embedding_dim": embedding_dim,
        "model_type": model_type,
        "task_type": task_type,
        "predictions": [],
        "visualizations": [],
    }

    if outcomes is not None:
        if len(outcomes) != len(sequences):
            raise ValueError("Number of outcomes must match number of sequences")

        predictor = EventSequencePredictor(
            model_type=model_type,
            task_type=task_type,
            embedding_dim=embedding_dim,
            random_state=random_state,
        )
        predictor.fit(token_sequences, outcomes, event_embeddings=embeddings_result)
        predictions = predictor.predict(token_sequences)
        results["predictions"] = predictions.tolist() if hasattr(predictions, "tolist") else list(predictions)

        if kwargs.get("save_model", workflow_config.get("save_model", True)):
            model_path = output_dir / "trained_model.json"
            predictor.save_model(model_path)
            results["model"] = str(model_path)
            results["model_path"] = str(model_path)

    results["visualizations"] = _generate_analysis_visualizations(sequences, outcomes, output_dir, results)

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
            if hasattr(event, "domain"):
                stats["unique_domains"].add(event.domain)

    stats["unique_event_types"] = list(stats["unique_event_types"])
    stats["unique_domains"] = list(stats["unique_domains"])
    stats["n_unique_event_types"] = len(stats["unique_event_types"])
    stats["n_unique_domains"] = len(stats["unique_domains"])

    return stats


def _generate_analysis_visualizations(
    sequences: List[Any], outcomes: Optional[List[str]], output_dir: Path, results: Dict[str, Any]
) -> List[Dict[str, str]]:
    """Generate analysis visualizations."""
    visualizations = []

    try:
        from . import visualization

        # Sequence length distribution
        fig = visualization.plot_sequence_length_distribution(sequences)
        if fig:
            length_plot = output_dir / "sequence_lengths.png"
            fig.savefig(length_plot, dpi=300, bbox_inches="tight")
            visualizations.append({"type": "sequence_length_distribution", "file": str(length_plot)})

        # Domain distribution if domains available
        if any(hasattr(event, "domain") for seq in sequences for event in seq.events):
            fig = visualization.plot_domain_distribution(sequences)
            if fig:
                domain_plot = output_dir / "domain_distribution.png"
                fig.savefig(domain_plot, dpi=300, bbox_inches="tight")
                visualizations.append({"type": "domain_distribution", "file": str(domain_plot)})

        # Outcome distribution if outcomes available
        if outcomes:
            fig = visualization.plot_outcome_distribution(outcomes)
            if fig:
                outcome_plot = output_dir / "outcome_distribution.png"
                fig.savefig(outcome_plot, dpi=300, bbox_inches="tight")
                visualizations.append({"type": "outcome_distribution", "file": str(outcome_plot)})

    except Exception as e:
        logger.warning(f"Error generating visualizations: {e}")

    return visualizations


def compare_populations(
    sequences1: List[Any],
    sequences2: List[Any],
    group_names: tuple[str, str] = ("Group1", "Group2"),
    output_dir: Optional[Union[str, Path]] = None,
    **kwargs: Any,
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

    group1_stats = _compute_sequence_stats(sequences1)
    group2_stats = _compute_sequence_stats(sequences2)

    results = {
        "group1_name": group_names[0],
        "group2_name": group_names[1],
        "group1_stats": group1_stats,
        "group2_stats": group2_stats,
        # Backward-compatible aliases used by older tests and scripts.
        "group1": {"n_sequences": len(sequences1), **group1_stats},
        "group2": {"n_sequences": len(sequences2), **group2_stats},
        "comparison": {},
        "visualizations": [],
        "output_dir": str(output_dir),
    }

    # Statistical comparisons
    results["comparison"] = _compare_sequence_groups(sequences1, sequences2, group_names)

    # Generate comparison visualizations
    try:
        from . import visualization

        fig = visualization.plot_population_comparison(sequences1, sequences2, group_names=group_names)
        if fig:
            comp_plot = output_dir / "population_comparison.png"
            fig.savefig(comp_plot, dpi=300, bbox_inches="tight")
            results["visualizations"].append({"type": "population_comparison", "file": str(comp_plot)})
    except Exception as e:
        logger.warning(f"Error generating comparison visualizations: {e}")

    # Save results
    results_file = output_dir / "comparison_results.json"
    io.dump_json(results, results_file)

    logger.info(f"Population comparison complete. Results saved to {results_file}")
    return results


def _compare_sequence_groups(
    sequences1: List[Any], sequences2: List[Any], group_names: tuple[str, str]
) -> Dict[str, Any]:
    """Compare two groups of sequences statistically."""
    # Basic statistical tests
    lengths1 = [len(seq.events) for seq in sequences1]
    lengths2 = [len(seq.events) for seq in sequences2]

    comparison = {
        "length_comparison": {
            "group1_mean": sum(lengths1) / len(lengths1) if lengths1 else 0,
            "group2_mean": sum(lengths2) / len(lengths2) if lengths2 else 0,
            "group1_median": sorted(lengths1)[len(lengths1) // 2] if lengths1 else 0,
            "group2_median": sorted(lengths2)[len(lengths2) // 2] if lengths2 else 0,
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
        "jaccard_similarity": len(types1 & types2) / len(types1 | types2) if (types1 | types2) else 0,
    }

    common_event_types = sorted(types1 & types2)
    unique_to_group1 = sorted(types1 - types2)
    unique_to_group2 = sorted(types2 - types1)
    comparison.update(
        {
            "common_event_types": common_event_types,
            "unique_to_group1": unique_to_group1,
            "unique_to_group2": unique_to_group2,
            "sequence_length": comparison["length_comparison"],
            "event_diversity": {
                "group1_unique_events": len(types1),
                "group2_unique_events": len(types2),
            },
        }
    )

    return comparison


def intervention_analysis(
    sequences: List[EventSequence],
    intervention_time: float,
    output_dir: Optional[Union[str, Path]] = None,
    outcomes: Optional[List[Any]] = None,
    pre_intervention_outcomes: Optional[Any] = None,
    post_intervention_outcomes: Optional[Any] = None,
) -> Dict[str, Any]:
    """Analyze the effects of an intervention on life event sequences.

    Args:
        sequences: List of EventSequence objects
        intervention_time: Timestamp of the intervention
        output_dir: Directory to save analysis results (optional)
        outcomes: Optional list of outcomes to compare before/after intervention

    Returns:
        Dictionary containing pre/post intervention analysis
    """
    import datetime
    from pathlib import Path

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

    pre_outcomes = pre_intervention_outcomes
    post_outcomes = post_intervention_outcomes
    if outcomes is not None and pre_outcomes is None and post_outcomes is None:
        pre_outcomes = outcomes[: len(pre_sequences)]
        post_outcomes = outcomes[len(pre_sequences) :] if len(outcomes) > len(pre_sequences) else None

    pre_output = output_dir / "pre" if output_dir else None
    post_output = output_dir / "post" if output_dir else None

    # Analyze pre-intervention
    pre_analysis = analyze_life_course(pre_sequences, pre_outcomes, output_dir=pre_output) if pre_sequences else {}

    # Analyze post-intervention
    post_analysis = analyze_life_course(post_sequences, post_outcomes, output_dir=post_output) if post_sequences else {}

    # Calculate differences
    differences = {}
    for key in pre_analysis.keys():
        if key in post_analysis and isinstance(pre_analysis[key], (int, float)):
            if isinstance(post_analysis[key], (int, float)):
                differences[f"{key}_change"] = post_analysis[key] - pre_analysis[key]

    outcome_change = {}
    if pre_outcomes is not None and post_outcomes is not None:
        import numpy as np

        pre_arr = np.asarray(pre_outcomes, dtype=float)
        post_arr = np.asarray(post_outcomes, dtype=float)
        if len(pre_arr) and len(post_arr):
            outcome_change = {
                "mean": float(post_arr.mean() - pre_arr.mean()),
                "pre_mean": float(pre_arr.mean()),
                "post_mean": float(post_arr.mean()),
            }

    # Generate intervention effect plots if output directory provided
    if output_dir:
        try:
            from .visualization import plot_intervention_effects

            plot_path = output_dir / "intervention_effects.png"
            plot_intervention_effects(pre_sequences, post_sequences, output_path=plot_path)
        except ImportError:
            logger.warning("Intervention effects plot not available")

    return {
        "pre_intervention": pre_analysis,
        "post_intervention": post_analysis,
        "differences": differences,
        "outcome_change": outcome_change,
        "intervention_time": intervention_time,
        "n_pre_sequences": len(pre_sequences),
        "n_post_sequences": len(post_sequences),
        "output_dir": str(output_dir) if output_dir else None,
    }


def event_importance(
    sequences: List[EventSequence], method: str = "frequency", normalize: bool = True
) -> Dict[str, float]:
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
                prev_event = seq.events[i - 1].event_type
                transition_key = f"{prev_event}->{event_type}"
                event_transitions[transition_key] = event_transitions.get(transition_key, 0) + 1

    if method == "frequency":
        # Simple frequency-based importance
        importance_scores = event_counts.copy()

    elif method == "temporal":
        # Importance based on temporal positioning
        importance_scores = {}
        for event_type, positions in event_positions.items():
            # Events that appear early get higher importance
            avg_position = sum(positions) / len(positions)
            importance_scores[event_type] = 1.0 - avg_position  # Earlier = more important

    elif method == "transition":
        # Importance based on being transition hubs
        importance_scores = {}
        for event_type in event_counts.keys():
            # Count incoming and outgoing transitions
            incoming = sum(count for key, count in event_transitions.items() if key.endswith(f"->{event_type}"))
            outgoing = sum(count for key, count in event_transitions.items() if key.startswith(f"{event_type}->"))
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

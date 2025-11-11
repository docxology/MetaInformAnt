"""End-to-end analysis pipelines for life course analysis.

This module provides complete workflows for analyzing event sequences,
from raw data to predictions and comparisons.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
from numpy.typing import NDArray

from ..core import config, io, paths
from ..core.errors import error_context
from ..core.logging import log_with_metadata, setup_logger
from .config import LifeEventsWorkflowConfig, load_life_events_config
from .events import EventDatabase, EventSequence
from .embeddings import learn_event_embeddings, sequence_embeddings
from .models import EventSequencePredictor


def analyze_life_course(
    sequences: List[EventSequence],
    outcomes: Optional[NDArray] = None,
    config_path: Optional[str | Path] = None,
    config_obj: Optional[LifeEventsWorkflowConfig | Dict[str, Any]] = None,
    output_dir: Optional[str | Path] = None
) -> Dict[str, Any]:
    """Complete workflow from raw event sequences to predictions.
    
    Performs:
    1. Event sequence loading and validation
    2. Embedding learning
    3. Model training (if outcomes provided)
    4. Predictions and analysis
    
    Args:
        sequences: List of event sequences (must not be empty)
        outcomes: Optional outcome labels/values for training (must match sequence length if provided)
        config_path: Optional path to configuration file
        output_dir: Output directory (defaults to output/life_events/)
        
    Returns:
        Dictionary with analysis results including embeddings path, model type, and predictions
        
    Raises:
        ValueError: If sequences is empty or outcomes length doesn't match sequences length
        
    Examples:
        >>> from metainformant.life_events import EventSequence, Event, analyze_life_course
        >>> from datetime import datetime
        >>> seq = EventSequence("p1", [Event("degree", datetime(2010, 1, 1), "education")])
        >>> results = analyze_life_course([seq], output_dir="output/test")
        >>> "embeddings" in results
        True
    """
    logger = setup_logger(__name__)
    
    # Log workflow start with metadata
    log_with_metadata(
        logger,
        "Starting life course analysis",
        {
            "num_sequences": len(sequences),
            "has_outcomes": outcomes is not None,
            "output_dir": str(output_dir) if output_dir else "default",
        },
    )
    
    # Validate inputs
    if not sequences:
        with error_context("Life course analysis validation failed"):
            raise ValueError("sequences list cannot be empty")
    
    if outcomes is not None and len(outcomes) != len(sequences):
        with error_context("Life course analysis validation failed"):
            raise ValueError(
                f"outcomes length ({len(outcomes)}) must match sequences length ({len(sequences)})"
            )
    
    if output_dir is None:
        output_dir = Path("output/life_events")
    else:
        output_dir = Path(output_dir)
    paths.ensure_directory(output_dir)
    
    # Load config if provided
    cfg: Dict[str, Any] = {}
    if config_obj is not None:
        if isinstance(config_obj, LifeEventsWorkflowConfig):
            # Convert config object to dict for compatibility
            cfg = {
                "embedding_dim": config_obj.embedding.get("embedding_dim", 100),
                "window_size": config_obj.embedding.get("window_size", 5),
                "epochs": config_obj.embedding.get("epochs", 10),
                "model_type": config_obj.model.get("model_type", "embedding"),
                "task_type": config_obj.model.get("task_type", "classification"),
                "random_state": config_obj.model.get("random_state", 42),
            }
            # Merge in other config sections
            cfg.update(config_obj.embedding)
            cfg.update(config_obj.model)
            cfg.update(config_obj.workflow)
            if output_dir is None and config_obj.work_dir:
                output_dir = config_obj.work_dir
        elif isinstance(config_obj, dict):
            cfg = config_obj
    elif config_path:
        try:
            # Try loading as LifeEventsWorkflowConfig
            config_obj = load_life_events_config(config_path)
            cfg = {
                "embedding_dim": config_obj.embedding.get("embedding_dim", 100),
                "window_size": config_obj.embedding.get("window_size", 5),
                "epochs": config_obj.embedding.get("epochs", 10),
                "model_type": config_obj.model.get("model_type", "embedding"),
                "task_type": config_obj.model.get("task_type", "classification"),
                "random_state": config_obj.model.get("random_state", 42),
            }
            cfg.update(config_obj.embedding)
            cfg.update(config_obj.model)
            cfg.update(config_obj.workflow)
            if output_dir is None and config_obj.work_dir:
                output_dir = config_obj.work_dir
        except Exception:
            # Fallback to old config loading
            cfg = config.load_config(config_path)
    
    # Convert sequences to token format
    sequences_tokens = []
    for seq in sequences:
        tokens = [f"{e.domain}:{e.event_type}" for e in seq.events]
        sequences_tokens.append(tokens)
    
    logger.info(f"Processing {len(sequences)} event sequences")
    
    # Learn embeddings
    embedding_dim = cfg.get("embedding_dim", 100)
    window_size = cfg.get("window_size", 5)
    epochs = cfg.get("epochs", 10)
    
    logger.info("Learning event embeddings...")
    embeddings = learn_event_embeddings(
        sequences_tokens,
        embedding_dim=embedding_dim,
        window_size=window_size,
        epochs=epochs
    )
    
    # Save embeddings
    embeddings_file = output_dir / "embeddings.json"
    embeddings_dict = {k: v.tolist() for k, v in embeddings.items()}
    io.dump_json(embeddings_dict, embeddings_file)
    logger.info(f"Saved embeddings to {embeddings_file}")
    
    results = {
        "n_sequences": len(sequences),
        "n_events": sum(len(seq.events) for seq in sequences),
        "embeddings": str(embeddings_file),
        "embedding_dim": embedding_dim,
    }
    
    # Train model if outcomes provided
    if outcomes is not None:
        logger.info("Training prediction model...")
        
        task_type = cfg.get("task_type", "classification")
        if len(np.unique(outcomes)) > 10 or np.issubdtype(outcomes.dtype, np.floating):
            task_type = "regression"
        
        predictor = EventSequencePredictor(
            model_type=cfg.get("model_type", "embedding"),
            task_type=task_type,
            embedding_dim=embedding_dim,
            random_state=cfg.get("random_state", 42)
        )
        
        predictor.fit(sequences_tokens, outcomes, event_embeddings=embeddings)
        
        # Save model
        model_file = output_dir / "model.json"
        predictor.save_model(model_file)
        logger.info(f"Saved model to {model_file}")
        results["model"] = str(model_file)
        
        # Make predictions
        predictions = predictor.predict(sequences_tokens)
        
        results["model_type"] = predictor.model_type
        results["task_type"] = predictor.task_type
        results["predictions"] = predictions.tolist()
        
        if task_type == "classification":
            # Calculate accuracy
            accuracy = np.mean(predictions == outcomes)
            results["accuracy"] = float(accuracy)
            logger.info(f"Classification accuracy: {accuracy:.3f}")
        else:
            # Calculate R²
            try:
                from sklearn.metrics import r2_score
                r2 = r2_score(outcomes, predictions)
                results["r2_score"] = float(r2)
                logger.info(f"Regression R²: {r2:.3f}")
            except ImportError:
                # sklearn not available, use simple correlation
                correlation = np.corrcoef(outcomes, predictions)[0, 1]
                results["correlation"] = float(correlation)
                logger.info(f"Regression correlation: {correlation:.3f}")
    
    # Save results
    results_file = output_dir / "analysis_results.json"
    io.dump_json(results, results_file)
    
    # Log workflow completion with metadata
    log_with_metadata(
        logger,
        "Life course analysis completed",
        {
            "num_sequences": len(sequences),
            "has_model": "model" in results,
            "has_embeddings": "embeddings" in results,
            "results_file": str(results_file),
        },
    )
    logger.info(f"Analysis complete. Results saved to {results_file}")
    
    return results


def compare_populations(
    sequences_group1: List[EventSequence],
    sequences_group2: List[EventSequence],
    output_dir: Optional[str | Path] = None
) -> Dict[str, Any]:
    """Compare event patterns across two population groups.
    
    Args:
        sequences_group1: Event sequences for group 1 (must not be empty)
        sequences_group2: Event sequences for group 2 (must not be empty)
        output_dir: Output directory
        
    Returns:
        Dictionary with comparison statistics including event type overlaps and domain distributions
        
    Raises:
        ValueError: If either group is empty
        
    Examples:
        >>> from metainformant.life_events import EventSequence, Event, compare_populations
        >>> from datetime import datetime
        >>> group1 = [EventSequence("p1", [Event("degree", datetime(2010, 1, 1), "education")])]
        >>> group2 = [EventSequence("p2", [Event("job_change", datetime(2015, 1, 1), "occupation")])]
        >>> comparison = compare_populations(group1, group2)
        >>> "comparison" in comparison
        True
    """
    logger = setup_logger(__name__)
    
    # Validate inputs
    if not sequences_group1:
        raise ValueError("sequences_group1 cannot be empty")
    if not sequences_group2:
        raise ValueError("sequences_group2 cannot be empty")
    
    if output_dir is None:
        output_dir = Path("output/life_events/comparison")
    else:
        output_dir = Path(output_dir)
    paths.ensure_directory(output_dir)
    
    # Convert to databases
    db1 = EventDatabase(sequences=sequences_group1)
    db2 = EventDatabase(sequences=sequences_group2)
    
    # Get statistics
    stats1 = db1.get_statistics()
    stats2 = db2.get_statistics()
    
    # Compare event types
    event_types1 = set(stats1["event_types"].keys())
    event_types2 = set(stats2["event_types"].keys())
    
    common_types = event_types1 & event_types2
    unique_to_group1 = event_types1 - event_types2
    unique_to_group2 = event_types2 - event_types1
    
    # Compare domains
    domains1 = set(stats1["domains"].keys())
    domains2 = set(stats2["domains"].keys())
    
    comparison = {
        "group1": {
            "n_sequences": stats1["n_sequences"],
            "n_events": stats1["n_events"],
            "avg_events_per_sequence": stats1["avg_events_per_sequence"],
            "domains": stats1["domains"],
            "event_types": stats1["event_types"],
        },
        "group2": {
            "n_sequences": stats2["n_sequences"],
            "n_events": stats2["n_events"],
            "avg_events_per_sequence": stats2["avg_events_per_sequence"],
            "domains": stats2["domains"],
            "event_types": stats2["event_types"],
        },
        "comparison": {
            "common_event_types": sorted(list(common_types)),
            "unique_to_group1": sorted(list(unique_to_group1)),
            "unique_to_group2": sorted(list(unique_to_group2)),
            "common_domains": sorted(list(domains1 & domains2)),
        },
    }
    
    # Save comparison
    comparison_file = output_dir / "comparison.json"
    io.dump_json(comparison, comparison_file)
    logger.info(f"Comparison saved to {comparison_file}")
    
    return comparison


def intervention_analysis(
    sequences: List[EventSequence],
    intervention_time: float,
    pre_intervention_outcomes: Optional[NDArray] = None,
    post_intervention_outcomes: Optional[NDArray] = None,
    output_dir: Optional[str | Path] = None
) -> Dict[str, Any]:
    """Analyze effects of interventions on life courses.
    
    Compares pre- and post-intervention event patterns and outcomes.
    
    Args:
        sequences: Event sequences (must not be empty)
        intervention_time: Time point of intervention (timestamp or numeric)
        pre_intervention_outcomes: Outcomes before intervention (must match sequence length if provided)
        post_intervention_outcomes: Outcomes after intervention (must match sequence length if provided)
        output_dir: Output directory
        
    Returns:
        Dictionary with intervention analysis results including pre/post statistics and outcome changes
        
    Raises:
        ValueError: If sequences is empty or outcome arrays don't match sequence length
        
    Examples:
        >>> from metainformant.life_events import EventSequence, Event, intervention_analysis
        >>> from datetime import datetime
        >>> seq = EventSequence("p1", [Event("degree", datetime(2010, 1, 1), "education")])
        >>> results = intervention_analysis([seq], intervention_time=datetime(2015, 1, 1).timestamp())
        >>> "pre_intervention" in results
        True
    """
    logger = setup_logger(__name__)
    
    # Validate inputs
    if not sequences:
        raise ValueError("sequences list cannot be empty")
    
    if pre_intervention_outcomes is not None and len(pre_intervention_outcomes) != len(sequences):
        raise ValueError(
            f"pre_intervention_outcomes length ({len(pre_intervention_outcomes)}) must match "
            f"sequences length ({len(sequences)})"
        )
    
    if post_intervention_outcomes is not None and len(post_intervention_outcomes) != len(sequences):
        raise ValueError(
            f"post_intervention_outcomes length ({len(post_intervention_outcomes)}) must match "
            f"sequences length ({len(sequences)})"
        )
    
    if output_dir is None:
        output_dir = Path("output/life_events/intervention")
    else:
        output_dir = Path(output_dir)
    paths.ensure_directory(output_dir)
    
    # Split sequences into pre and post intervention
    pre_sequences = []
    post_sequences = []
    
    for seq in sequences:
        def get_timestamp(event):
            if isinstance(event.timestamp, float):
                return event.timestamp
            return event.timestamp.timestamp()
        
        pre_events = [e for e in seq.events if get_timestamp(e) < intervention_time]
        post_events = [e for e in seq.events if get_timestamp(e) >= intervention_time]
        
        pre_sequences.append(EventSequence(
            person_id=seq.person_id,
            events=pre_events,
            metadata=seq.metadata
        ))
        post_sequences.append(EventSequence(
            person_id=seq.person_id,
            events=post_events,
            metadata=seq.metadata
        ))
    
    # Compare pre and post
    pre_db = EventDatabase(sequences=pre_sequences)
    post_db = EventDatabase(sequences=post_sequences)
    
    pre_stats = pre_db.get_statistics()
    post_stats = post_db.get_statistics()
    
    results = {
        "intervention_time": intervention_time,
        "pre_intervention": pre_stats,
        "post_intervention": post_stats,
    }
    
    # Compare outcomes if provided
    if pre_intervention_outcomes is not None and post_intervention_outcomes is not None:
        outcome_diff = post_intervention_outcomes - pre_intervention_outcomes
        results["outcome_change"] = {
            "mean": float(np.mean(outcome_diff)),
            "std": float(np.std(outcome_diff)),
            "median": float(np.median(outcome_diff)),
        }
        
        # Statistical test (simple t-test)
        try:
            from scipy import stats
            t_stat, p_value = stats.ttest_rel(pre_intervention_outcomes, post_intervention_outcomes)
            results["outcome_change"]["t_statistic"] = float(t_stat)
            results["outcome_change"]["p_value"] = float(p_value)
        except ImportError:
            # scipy not available, skip statistical test
            results["outcome_change"]["t_statistic"] = None
            results["outcome_change"]["p_value"] = None
    
    # Save results
    results_file = output_dir / "intervention_analysis.json"
    io.dump_json(results, results_file)
    logger.info(f"Intervention analysis saved to {results_file}")
    
    return results


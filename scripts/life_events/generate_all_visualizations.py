#!/usr/bin/env python3
"""Generate all comprehensive visualizations for life events analysis.

This script generates the complete suite of visualizations for event sequences.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, paths
from metainformant.core.logging import setup_logger
from metainformant.life_events import (
    EventSequencePredictor,
    convert_sequences_to_tokens,
    event_importance,
    learn_event_embeddings,
    load_sequences_from_json,
    plot_domain_timeline,
    plot_embedding_clusters,
    plot_event_cooccurrence,
    plot_event_frequency_heatmap,
    plot_outcome_distribution,
    plot_prediction_accuracy,
    plot_sequence_similarity,
    plot_temporal_patterns,
    plot_transition_network,
)

logger = setup_logger(__name__, level="INFO")


def main():
    """Generate all visualizations."""
    sequences_file = Path("output/life_events/full_analysis/synthetic_sequences.json")
    outcomes_file = Path("output/life_events/full_analysis/synthetic_outcomes.json")
    model_file = Path("output/life_events/full_analysis/model.json")
    embeddings_file = Path("output/life_events/full_analysis/embeddings.json")
    output_dir = Path("output/life_events/full_analysis/visualizations")
    
    paths.ensure_directory(output_dir)
    
    logger.info("Loading data...")
    sequences = load_sequences_from_json(sequences_file)
    logger.info(f"✅ Loaded {len(sequences)} sequences")
    
    outcomes = None
    if outcomes_file.exists():
        outcomes_data = io.load_json(outcomes_file)
        outcomes = np.array(outcomes_data.get("outcomes", []))
        logger.info(f"✅ Loaded {len(outcomes)} outcomes")
    
    embeddings = None
    if embeddings_file.exists():
        embeddings_data = io.load_json(embeddings_file)
        embeddings = {k: np.array(v) for k, v in embeddings_data.items()}
        logger.info(f"✅ Loaded embeddings for {len(embeddings)} events")
    
    sequences_tokens = convert_sequences_to_tokens(sequences)
    
    logger.info("Generating comprehensive visualizations...")
    
    # 1. Event co-occurrence
    try:
        logger.info("  Generating event co-occurrence heatmap...")
        plot_event_cooccurrence(
            sequences,
            output_path=output_dir / "event_cooccurrence.png",
            top_n=20
        )
        logger.info("  ✅ Event co-occurrence generated")
    except Exception as e:
        logger.warning(f"  ⚠️  Event co-occurrence failed: {e}")
    
    # 2. Outcome distribution
    if outcomes is not None:
        try:
            logger.info("  Generating outcome distribution...")
            plot_outcome_distribution(
                outcomes,
                output_path=output_dir / "outcome_distribution.png",
                plot_type="histogram"
            )
            logger.info("  ✅ Outcome distribution generated")
        except Exception as e:
            logger.warning(f"  ⚠️  Outcome distribution failed: {e}")
    
    # 3. Sequence similarity
    try:
        logger.info("  Generating sequence similarity matrix...")
        plot_sequence_similarity(
            sequences,
            embeddings=embeddings,
            output_path=output_dir / "sequence_similarity.png"
        )
        logger.info("  ✅ Sequence similarity generated")
    except Exception as e:
        logger.warning(f"  ⚠️  Sequence similarity failed: {e}")
    
    # 4. Transition network
    try:
        logger.info("  Generating transition network...")
        plot_transition_network(
            sequences,
            output_path=output_dir / "transition_network.png",
            top_n=15
        )
        logger.info("  ✅ Transition network generated")
    except ImportError:
        logger.warning("  ⚠️  Transition network skipped (networkx not available)")
    except Exception as e:
        logger.warning(f"  ⚠️  Transition network failed: {e}")
    
    # 5. Domain timeline
    try:
        logger.info("  Generating domain timeline...")
        plot_domain_timeline(
            sequences,
            output_path=output_dir / "domain_timeline.png",
            max_sequences=10
        )
        logger.info("  ✅ Domain timeline generated")
    except Exception as e:
        logger.warning(f"  ⚠️  Domain timeline failed: {e}")
    
    # 6. Prediction accuracy
    if outcomes is not None and model_file.exists():
        try:
            logger.info("  Generating prediction accuracy plots...")
            predictor = EventSequencePredictor.load_model(model_file)
            predictions = predictor.predict(sequences_tokens)
            
            plot_prediction_accuracy(
                outcomes,
                predictions,
                task_type="classification" if len(np.unique(outcomes)) <= 10 else "regression",
                output_path=output_dir / "prediction_accuracy.png"
            )
            logger.info("  ✅ Prediction accuracy generated")
        except Exception as e:
            logger.warning(f"  ⚠️  Prediction accuracy failed: {e}")
    
    # 7. Temporal patterns
    if outcomes is not None:
        try:
            logger.info("  Generating temporal patterns...")
            # Create importance scores based on outcomes
            importance_scores = {i: float(abs(outcomes[i])) for i in range(len(outcomes))}
            plot_temporal_patterns(
                sequences,
                importance_scores=importance_scores,
                output_path=output_dir / "temporal_patterns.png"
            )
            logger.info("  ✅ Temporal patterns generated")
        except Exception as e:
            logger.warning(f"  ⚠️  Temporal patterns failed: {e}")
    
    # 8. Event frequency heatmap
    try:
        logger.info("  Generating event frequency heatmap...")
        plot_event_frequency_heatmap(
            sequences,
            output_path=output_dir / "event_frequency_heatmap.png",
            time_bins=10
        )
        logger.info("  ✅ Event frequency heatmap generated")
    except Exception as e:
        logger.warning(f"  ⚠️  Event frequency heatmap failed: {e}")
    
    # 9. Embedding clusters
    if embeddings is not None:
        try:
            logger.info("  Generating embedding clusters...")
            plot_embedding_clusters(
                embeddings,
                clusters=None,
                method="umap",
                output_path=output_dir / "embedding_clusters.png"
            )
            logger.info("  ✅ Embedding clusters generated")
        except Exception as e:
            logger.warning(f"  ⚠️  Embedding clusters failed: {e}")
    
    logger.info("=" * 60)
    logger.info("✅ Comprehensive visualization suite complete!")
    logger.info(f"   Output directory: {output_dir}")
    logger.info("=" * 60)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())


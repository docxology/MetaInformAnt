#!/usr/bin/env python3
"""Life events analysis workflow orchestrator.

This script provides comprehensive orchestration for life course analysis workflows,
including synthetic data generation, embedding learning, model training, predictions,
and comprehensive visualization.

Usage:
    python3 scripts/life_events/run_life_events_analysis.py --synthetic --n-sequences 100
    python3 scripts/life_events/run_life_events_analysis.py --input data/life_events/sequences.json --config config/life_events_template.yaml
    python3 scripts/life_events/run_life_events_analysis.py --help
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Optional

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, paths
from metainformant.core.utils.logging import setup_logger

logger = setup_logger(__name__, level="INFO")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Life events analysis workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate synthetic data and run complete analysis
  %(prog)s --synthetic --n-sequences 100 --generate-outcomes

  # Analyze existing sequences with config file
  %(prog)s --input data/life_events/sequences.json --config config/life_events_template.yaml

  # Run with custom embedding parameters
  %(prog)s --synthetic --embedding-dim 200 --window-size 10 --epochs 20

  # Generate visualizations only
  %(prog)s --input data/life_events/sequences.json --visualize-only
        """,
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Input event sequences file (JSON format)",
    )
    parser.add_argument(
        "--synthetic",
        action="store_true",
        help="Generate synthetic event sequences instead of loading from file",
    )
    parser.add_argument(
        "--n-sequences",
        type=int,
        default=50,
        help="Number of sequences to generate (default: 50)",
    )
    parser.add_argument(
        "--min-events",
        type=int,
        default=5,
        help="Minimum events per sequence (default: 5)",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=30,
        help="Maximum events per sequence (default: 30)",
    )
    parser.add_argument(
        "--generate-outcomes",
        action="store_true",
        help="Generate outcome labels for synthetic data",
    )
    parser.add_argument(
        "--outcome-relationship",
        type=str,
        default="complex",
        choices=["random", "health_focused", "education_focused", "complex"],
        help="How outcomes relate to events (default: complex)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        help="Configuration file path (YAML, TOML, or JSON)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/life_events"),
        help="Output directory (default: output/life_events)",
    )
    parser.add_argument(
        "--embedding-dim",
        type=int,
        help="Embedding dimension (overrides config)",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        help="Context window size (overrides config)",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        help="Training epochs (overrides config)",
    )
    parser.add_argument(
        "--model-type",
        type=str,
        choices=["embedding", "simple", "lstm"],
        help="Model type (overrides config)",
    )
    parser.add_argument(
        "--visualize-only",
        action="store_true",
        help="Only generate visualizations (requires existing model)",
    )
    parser.add_argument(
        "--all-visualizations",
        action="store_true",
        help="Generate all visualization types",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without executing",
    )
    return parser.parse_args()


def generate_synthetic_data(
    n_sequences: int,
    min_events: int,
    max_events: int,
    generate_outcomes: bool,
    outcome_relationship: str,
    random_state: int = 42,
) -> tuple[list, Optional[Any]]:
    """Generate synthetic life event sequences.

    Args:
        n_sequences: Number of sequences to generate
        min_events: Minimum events per sequence
        max_events: Maximum events per sequence
        generate_outcomes: Whether to generate outcomes
        outcome_relationship: Outcome relationship pattern
        random_state: Random seed

    Returns:
        Tuple of (sequences, outcomes)
    """
    logger.info(f"Generating {n_sequences} synthetic event sequences...")

    from metainformant.life_events import generate_synthetic_life_events

    sequences, outcomes = generate_synthetic_life_events(
        n_sequences=n_sequences,
        min_events_per_sequence=min_events,
        max_events_per_sequence=max_events,
        generate_outcomes=generate_outcomes,
        outcome_relationship=outcome_relationship,
        random_state=random_state,
    )

    logger.info(f"✅ Generated {len(sequences)} sequences")
    if outcomes is not None:
        logger.info(f"✅ Generated outcomes for {len(outcomes)} sequences")

    return sequences, outcomes


def load_sequences(input_path: Path) -> list:
    """Load event sequences from file.

    Args:
        input_path: Path to input file

    Returns:
        List of EventSequence objects
    """
    logger.info(f"Loading event sequences from {input_path}")

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    from metainformant.life_events import load_sequences_from_json

    sequences = load_sequences_from_json(input_path)
    logger.info(f"✅ Loaded {len(sequences)} sequences")

    return sequences


def run_analysis(
    sequences: list,
    outcomes: Optional[Any],
    config_path: Optional[Path],
    output_dir: Path,
    embedding_dim: Optional[int],
    window_size: Optional[int],
    epochs: Optional[int],
    model_type: Optional[str],
) -> dict[str, Any]:
    """Run complete life course analysis workflow.

    Args:
        sequences: List of event sequences
        outcomes: Optional outcome labels/values
        config_path: Optional configuration file path
        output_dir: Output directory
        embedding_dim: Optional embedding dimension override
        window_size: Optional window size override
        epochs: Optional epochs override
        model_type: Optional model type override

    Returns:
        Analysis results dictionary
    """
    logger.info("Running complete life course analysis workflow...")

    import numpy as np

    from metainformant.life_events import analyze_life_course, load_life_events_config

    # Prepare configuration
    config_obj = None
    if config_path and config_path.exists():
        config_obj = load_life_events_config(config_path)
        logger.info(f"✅ Loaded configuration from {config_path}")

    # Override config with command-line arguments
    if config_obj and isinstance(config_obj, dict):
        config_dict = config_obj
    elif config_obj:
        config_dict = {
            "embedding_dim": config_obj.embedding.get("embedding_dim", 100),
            "window_size": config_obj.embedding.get("window_size", 5),
            "epochs": config_obj.embedding.get("epochs", 10),
            "model_type": config_obj.model.get("model_type", "embedding"),
            "task_type": config_obj.model.get("task_type", "classification"),
            "random_state": config_obj.model.get("random_state", 42),
        }
        config_dict.update(config_obj.embedding)
        config_dict.update(config_obj.model)
    else:
        config_dict = {}

    # Apply command-line overrides
    if embedding_dim is not None:
        config_dict["embedding_dim"] = embedding_dim
    if window_size is not None:
        config_dict["window_size"] = window_size
    if epochs is not None:
        config_dict["epochs"] = epochs
    if model_type is not None:
        config_dict["model_type"] = model_type

    # Convert outcomes to numpy array if needed
    outcomes_array = None
    if outcomes is not None:
        if not isinstance(outcomes, np.ndarray):
            outcomes_array = np.array(outcomes)
        else:
            outcomes_array = outcomes

    # Run workflow
    results = analyze_life_course(
        sequences,
        outcomes=outcomes_array,
        config_obj=config_dict if config_dict else None,
        output_dir=output_dir,
    )

    logger.info("✅ Analysis workflow complete")
    logger.info(f"   Model type: {results.get('model_type', 'N/A')}")
    logger.info(f"   Embedding dimension: {results.get('embedding_dim', 'N/A')}")
    if "accuracy" in results:
        logger.info(f"   Accuracy: {results['accuracy']:.3f}")
    if "r2_score" in results:
        logger.info(f"   R² score: {results['r2_score']:.3f}")

    return results


def generate_visualizations(
    sequences: list,
    results: dict[str, Any],
    output_dir: Path,
    all_viz: bool = False,
) -> None:
    """Generate all visualizations.

    Args:
        sequences: List of event sequences
        results: Analysis results dictionary
        output_dir: Output directory
        all_viz: Whether to generate all visualization types
    """
    logger.info("Generating visualizations...")

    viz_dir = output_dir / "visualizations"
    paths.ensure_directory(viz_dir)

    try:
        from metainformant.life_events import (
            EventSequencePredictor,
            convert_sequences_to_tokens,
            event_importance,
            learn_event_embeddings,
            plot_event_embeddings,
            plot_event_timeline,
            plot_prediction_importance,
        )

        # 1. Timeline visualization
        if sequences:
            logger.info("  Generating timeline visualization...")
            timeline_path = viz_dir / "event_timeline_first_sequence.png"
            plot_event_timeline(sequences[0], output_path=timeline_path)
            logger.info(f"  ✅ Timeline saved: {timeline_path}")

            if all_viz and len(sequences) > 1:
                # Generate timeline for a few more sequences
                for i in range(1, min(4, len(sequences))):
                    timeline_path = viz_dir / f"event_timeline_sequence_{i+1}.png"
                    plot_event_timeline(sequences[i], output_path=timeline_path)

        # 2. Embedding visualization
        if "embeddings" in results:
            logger.info("  Generating embedding visualization...")
            embeddings_path = Path(results["embeddings"])
            if embeddings_path.exists():
                embeddings_data = io.load_json(embeddings_path)
                import numpy as np

                embeddings = {k: np.array(v) for k, v in embeddings_data.items()}

                # Try different methods
                for method in ["pca", "umap"]:
                    try:
                        embedding_path = viz_dir / f"event_embeddings_{method}.png"
                        plot_event_embeddings(
                            embeddings,
                            method=method,
                            n_components=2,
                            output_path=embedding_path,
                        )
                        logger.info(f"  ✅ Embeddings ({method}) saved: {embedding_path}")
                        break
                    except Exception as e:
                        logger.warning(f"  ⚠️  {method} failed: {e}")
                        continue

        # 3. Importance visualization
        if "model" in results:
            logger.info("  Generating importance visualization...")
            try:
                model_path = Path(results["model"])
                if model_path.exists():
                    predictor = EventSequencePredictor.load_model(model_path)
                    sequences_tokens = convert_sequences_to_tokens(sequences)

                    # Get embeddings
                    if hasattr(predictor, "event_embeddings"):
                        embeddings = predictor.event_embeddings
                    else:
                        embeddings = learn_event_embeddings(
                            sequences_tokens,
                            embedding_dim=results.get("embedding_dim", 100),
                            random_state=42,
                        )

                    # Compute importance
                    importance = event_importance(
                        predictor,
                        sequences_tokens,
                        embeddings,
                        method="permutation",
                    )

                    importance_path = viz_dir / "prediction_importance.png"
                    plot_prediction_importance(
                        importance,
                        top_n=20,
                        output_path=importance_path,
                    )
                    logger.info(f"  ✅ Importance plot saved: {importance_path}")
            except Exception as e:
                logger.warning(f"  ⚠️  Importance visualization failed: {e}")

        logger.info("✅ All visualizations generated")

    except ImportError as e:
        logger.warning(f"Visualization dependencies not available: {e}")
    except Exception as e:
        logger.error(f"Error generating visualizations: {e}", exc_info=True)


def main():
    """Main workflow function."""
    args = parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    if args.dry_run:
        logger.info("DRY RUN MODE - No files will be created")
        logger.info(f"Would process: {'synthetic' if args.synthetic else args.input}")
        logger.info(f"Would output to: {args.output}")
        return 0

    try:
        # Setup output directory
        output_dir = Path(args.output)
        paths.ensure_directory(output_dir)
        logger.info(f"Output directory: {output_dir}")

        # Load or generate sequences
        sequences = None
        outcomes = None

        if args.synthetic:
            sequences, outcomes = generate_synthetic_data(
                n_sequences=args.n_sequences,
                min_events=args.min_events,
                max_events=args.max_events,
                generate_outcomes=args.generate_outcomes,
                outcome_relationship=args.outcome_relationship,
            )

            # Save synthetic data
            sequences_file = output_dir / "synthetic_sequences.json"
            sequences_data = [seq.to_dict() for seq in sequences]
            io.dump_json(sequences_data, sequences_file)
            logger.info(f"✅ Synthetic sequences saved: {sequences_file}")

            if outcomes is not None:
                outcomes_file = output_dir / "synthetic_outcomes.json"
                io.dump_json(
                    {
                        "outcomes": outcomes.tolist() if hasattr(outcomes, "tolist") else list(outcomes),
                        "n_sequences": len(outcomes),
                    },
                    outcomes_file,
                )
                logger.info(f"✅ Synthetic outcomes saved: {outcomes_file}")

        elif args.input:
            sequences = load_sequences(args.input)
            # Try to load outcomes if they exist alongside input
            outcomes_path = args.input.parent / f"{args.input.stem}_outcomes.json"
            if outcomes_path.exists():
                outcomes_data = io.load_json(outcomes_path)
                import numpy as np

                outcomes = np.array(outcomes_data.get("outcomes", []))
                logger.info(f"✅ Loaded outcomes from {outcomes_path}")
        else:
            logger.error("Must specify either --synthetic or --input")
            return 1

        if not sequences:
            logger.error("No sequences to process")
            return 1

        # Run analysis workflow
        if not args.visualize_only:
            results = run_analysis(
                sequences,
                outcomes,
                args.config,
                output_dir,
                args.embedding_dim,
                args.window_size,
                args.epochs,
                args.model_type,
            )

            # Save workflow summary
            summary_path = output_dir / "workflow_summary.json"
            summary = {
                "n_sequences": len(sequences),
                "n_events": sum(len(seq.events) for seq in sequences),
                "results": results,
                "config_overrides": {
                    "embedding_dim": args.embedding_dim,
                    "window_size": args.window_size,
                    "epochs": args.epochs,
                    "model_type": args.model_type,
                },
            }
            io.dump_json(summary, summary_path, indent=2)
            logger.info(f"✅ Workflow summary saved: {summary_path}")
        else:
            # Load existing results
            results_path = output_dir / "workflow_summary.json"
            if results_path.exists():
                summary = io.load_json(results_path)
                results = summary.get("results", {})
                logger.info("✅ Loaded existing workflow results")
            else:
                logger.error("No existing results found for visualization-only mode")
                return 1

        # Generate visualizations
        generate_visualizations(sequences, results, output_dir, args.all_visualizations)

        logger.info("=" * 60)
        logger.info("✅ Life events analysis workflow complete!")
        logger.info(f"   Output directory: {output_dir}")
        logger.info(f"   Sequences processed: {len(sequences)}")
        logger.info("=" * 60)

        return 0

    except Exception as e:
        logger.error(f"Workflow failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())

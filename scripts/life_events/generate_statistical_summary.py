#!/usr/bin/env python3
"""Generate comprehensive statistical summary of life events analysis.

This script produces detailed statistical summaries including:
- Event frequency statistics
- Sequence statistics
- Temporal patterns
- Domain distributions
- Outcome statistics
- Model performance metrics
"""

import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, paths
from metainformant.core.utils.logging import setup_logger
from metainformant.life_events import EventSequence, load_sequences_from_json

logger = setup_logger(__name__, level="INFO")


def calculate_sequence_statistics(sequences: List[EventSequence]) -> Dict[str, Any]:
    """Calculate comprehensive sequence statistics."""
    sequence_lengths = [len(seq.events) for seq in sequences]

    all_events = []
    events_by_domain = defaultdict(list)
    events_by_type = defaultdict(int)

    for seq in sequences:
        for event in seq.events:
            all_events.append(event)
            events_by_domain[event.domain].append(event)
            events_by_type[f"{event.domain}:{event.event_type}"] += 1

    # Calculate temporal spans
    temporal_spans = []
    for seq in sequences:
        if len(seq.events) > 1:
            dates = [event.timestamp for event in seq.events]
            span_days = (max(dates) - min(dates)).days
            temporal_spans.append(span_days)

    return {
        "n_sequences": len(sequences),
        "total_events": len(all_events),
        "unique_event_types": len(events_by_type),
        "sequence_lengths": {
            "mean": float(np.mean(sequence_lengths)),
            "median": float(np.median(sequence_lengths)),
            "std": float(np.std(sequence_lengths)),
            "min": int(np.min(sequence_lengths)),
            "max": int(np.max(sequence_lengths)),
        },
        "temporal_spans_days": {
            "mean": float(np.mean(temporal_spans)) if temporal_spans else None,
            "median": float(np.median(temporal_spans)) if temporal_spans else None,
            "std": float(np.std(temporal_spans)) if temporal_spans else None,
            "min": int(np.min(temporal_spans)) if temporal_spans else None,
            "max": int(np.max(temporal_spans)) if temporal_spans else None,
        },
        "events_by_domain": {domain: len(events) for domain, events in events_by_domain.items()},
        "top_event_types": dict(Counter(events_by_type).most_common(20)),
    }


def calculate_temporal_statistics(sequences: List[EventSequence]) -> Dict[str, Any]:
    """Calculate temporal pattern statistics."""
    all_dates = []
    events_by_month = defaultdict(int)
    events_by_year = defaultdict(int)
    events_by_weekday = defaultdict(int)

    for seq in sequences:
        for event in seq.events:
            all_dates.append(event.timestamp)
            events_by_month[event.timestamp.month] += 1
            events_by_year[event.timestamp.year] += 1
            events_by_weekday[event.timestamp.weekday()] += 1

    # Calculate inter-event intervals
    intervals = []
    for seq in sequences:
        if len(seq.events) > 1:
            dates = sorted([event.timestamp for event in seq.events])
            for i in range(1, len(dates)):
                interval_days = (dates[i] - dates[i - 1]).days
                intervals.append(interval_days)

    return {
        "date_range": {
            "start": min(all_dates).isoformat() if all_dates else None,
            "end": max(all_dates).isoformat() if all_dates else None,
        },
        "events_by_month": dict(sorted(events_by_month.items())),
        "events_by_year": dict(sorted(events_by_year.items())),
        "events_by_weekday": {
            ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"][d]: count
            for d, count in sorted(events_by_weekday.items())
        },
        "inter_event_intervals_days": {
            "mean": float(np.mean(intervals)) if intervals else None,
            "median": float(np.median(intervals)) if intervals else None,
            "std": float(np.std(intervals)) if intervals else None,
            "min": int(np.min(intervals)) if intervals else None,
            "max": int(np.max(intervals)) if intervals else None,
        },
    }


def calculate_cooccurrence_statistics(sequences: List[EventSequence]) -> Dict[str, Any]:
    """Calculate event co-occurrence statistics."""
    cooccurrence = defaultdict(int)
    domain_cooccurrence = defaultdict(int)

    for seq in sequences:
        event_types = [f"{e.domain}:{e.event_type}" for e in seq.events]
        domains = [e.domain for e in seq.events]

        # Event type co-occurrence
        for i, e1 in enumerate(event_types):
            for e2 in event_types[i + 1 :]:
                pair = tuple(sorted([e1, e2]))
                cooccurrence[pair] += 1

        # Domain co-occurrence
        unique_domains = set(domains)
        for d1 in unique_domains:
            for d2 in unique_domains:
                if d1 != d2:
                    pair = tuple(sorted([d1, d2]))
                    domain_cooccurrence[pair] += 1

    return {
        "top_event_cooccurrences": {
            f"{k[0]}|{k[1]}": v for k, v in sorted(cooccurrence.items(), key=lambda x: x[1], reverse=True)[:30]
        },
        "domain_cooccurrences": {f"{k[0]}|{k[1]}": v for k, v in domain_cooccurrence.items()},
    }


def main():
    """Generate comprehensive statistical summary."""
    sequences_file = Path("output/life_events/full_analysis/synthetic_sequences.json")
    outcomes_file = Path("output/life_events/full_analysis/synthetic_outcomes.json")
    analysis_file = Path("output/life_events/full_analysis/analysis_results.json")
    output_file = Path("output/life_events/full_analysis/statistical_summary.json")

    logger.info("Loading data...")
    sequences = load_sequences_from_json(sequences_file)
    logger.info(f"✅ Loaded {len(sequences)} sequences")

    outcomes = None
    if outcomes_file.exists():
        outcomes_data = io.load_json(outcomes_file)
        outcomes = outcomes_data.get("outcomes", [])
        logger.info(f"✅ Loaded {len(outcomes)} outcomes")

    analysis_results = None
    if analysis_file.exists():
        analysis_results = io.load_json(analysis_file)
        logger.info("✅ Loaded analysis results")

    logger.info("Calculating statistics...")

    # Sequence statistics
    sequence_stats = calculate_sequence_statistics(sequences)
    logger.info("  ✅ Sequence statistics calculated")

    # Temporal statistics
    temporal_stats = calculate_temporal_statistics(sequences)
    logger.info("  ✅ Temporal statistics calculated")

    # Co-occurrence statistics
    cooccurrence_stats = calculate_cooccurrence_statistics(sequences)
    logger.info("  ✅ Co-occurrence statistics calculated")

    # Outcome statistics
    outcome_stats = None
    if outcomes:
        outcomes_array = np.array(outcomes)
        outcome_stats = {
            "n_outcomes": len(outcomes),
            "mean": float(np.mean(outcomes_array)),
            "median": float(np.median(outcomes_array)),
            "std": float(np.std(outcomes_array)),
            "min": float(np.min(outcomes_array)),
            "max": float(np.max(outcomes_array)),
            "unique_values": int(len(np.unique(outcomes_array))),
        }
        logger.info("  ✅ Outcome statistics calculated")

    # Model performance
    model_stats = None
    if analysis_results:
        predictions = analysis_results.get("predictions", [])
        if predictions:
            pred_array = np.array(predictions)
            if outcomes:
                true_array = np.array(outcomes)
                correlation = float(np.corrcoef(true_array, pred_array)[0, 1])
                mse = float(np.mean((true_array - pred_array) ** 2))
                mae = float(np.mean(np.abs(true_array - pred_array)))

                model_stats = {
                    "correlation": correlation,
                    "mse": mse,
                    "mae": mae,
                    "rmse": float(np.sqrt(mse)),
                    "predictions_mean": float(np.mean(pred_array)),
                    "predictions_std": float(np.std(pred_array)),
                }
            else:
                model_stats = {
                    "predictions_mean": float(np.mean(pred_array)),
                    "predictions_std": float(np.std(pred_array)),
                }
            logger.info("  ✅ Model statistics calculated")

    # Compile summary
    summary = {
        "sequence_statistics": sequence_stats,
        "temporal_statistics": temporal_stats,
        "cooccurrence_statistics": cooccurrence_stats,
        "outcome_statistics": outcome_stats,
        "model_statistics": model_stats,
    }

    # Save summary
    io.dump_json(summary, output_file)
    logger.info("=" * 60)
    logger.info("✅ Statistical summary complete!")
    logger.info(f"   Output file: {output_file}")
    logger.info("=" * 60)

    # Print key statistics
    logger.info("\n📊 Key Statistics:")
    logger.info(f"   Sequences: {sequence_stats['n_sequences']}")
    logger.info(f"   Total events: {sequence_stats['total_events']}")
    logger.info(f"   Unique event types: {sequence_stats['unique_event_types']}")
    logger.info(f"   Mean sequence length: {sequence_stats['sequence_lengths']['mean']:.1f}")
    logger.info(
        f"   Mean temporal span: {sequence_stats['temporal_spans_days']['mean']:.1f} days"
        if sequence_stats["temporal_spans_days"]["mean"]
        else "   Mean temporal span: N/A"
    )

    if outcome_stats:
        logger.info(f"   Outcome mean: {outcome_stats['mean']:.3f}")
        logger.info(f"   Outcome std: {outcome_stats['std']:.3f}")

    if model_stats and "correlation" in model_stats:
        logger.info(f"   Prediction correlation: {model_stats['correlation']:.3f}")
        logger.info(f"   RMSE: {model_stats['rmse']:.3f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

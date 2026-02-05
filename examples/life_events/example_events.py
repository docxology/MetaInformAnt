#!/usr/bin/env python3
"""Life events analysis example.

This example demonstrates life course event analysis using METAINFORMANT's life_events toolkit.

Usage:
    python examples/life_events/example_events.py

Output:
    output/examples/life_events/event_analysis.json
"""

from __future__ import annotations

from pathlib import Path

from metainformant.core import io


def main():
    """Demonstrate life events analysis."""
    # Setup output directory
    output_dir = Path("output/examples/life_events")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Life Events Example ===")

    # Simulated life event sequences
    event_sequences = {
        "individual_1": ["birth", "school", "work", "marriage", "retirement"],
        "individual_2": ["birth", "school", "education", "work", "promotion"],
        "individual_3": ["birth", "school", "work", "divorce", "remarriage"],
    }

    # Basic sequence analysis
    analysis_results = {}

    for individual, events in event_sequences.items():
        analysis_results[individual] = {
            "events": events,
            "sequence_length": len(events),
            "unique_events": len(set(events)),
            "transition_points": len(events) - 1,
            "life_stage": "early" if "school" in events[:2] else "standard",
        }

    print(f"✓ Analyzed {len(event_sequences)} life event sequences")

    # Save results
    results_file = output_dir / "event_analysis.json"
    io.dump_json(
        {"life_events_analysis": {"sequences_analyzed": len(event_sequences), "results": analysis_results}},
        results_file,
        indent=2,
    )

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Life Events Example Complete ===")


if __name__ == "__main__":
    main()

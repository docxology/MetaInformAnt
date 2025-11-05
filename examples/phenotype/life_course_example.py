#!/usr/bin/env python3
"""Example script demonstrating life course phenotype extraction.

This script shows how to:
- Extract phenotypes from event sequences
- Aggregate temporal phenotypes
- Map events to trait categories
- Handle errors gracefully

Usage:
    python3 examples/phenotype/life_course_example.py
"""

import sys
from datetime import datetime
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

try:
    from metainformant.life_events import Event, EventSequence
    from metainformant.phenotype import (
        aggregate_temporal_phenotypes,
        extract_phenotypes_from_events,
        map_events_to_traits,
    )
    from metainformant.core.errors import ValidationError
    from metainformant.core.io import write_json
    from metainformant.core.paths import expand_and_resolve
    
    LIFE_EVENTS_AVAILABLE = True
except ImportError as e:
    print(f"Error: life_events module not available: {e}")
    print("Please ensure metainformant is installed with life_events support.")
    sys.exit(1)


def create_example_sequences():
    """Create example event sequences for demonstration."""
    sequences = [
        EventSequence(
            person_id="person_001",
            events=[
                Event("diabetes", datetime(2020, 1, 15), "health", {"severity": "moderate"}),
                Event("bachelors", datetime(2010, 6, 1), "education", {"degree": "BS"}),
                Event("job_change", datetime(2015, 3, 20), "occupation", {"company": "TechCorp"}),
                Event("promotion", datetime(2018, 9, 10), "occupation", {"title": "Senior Engineer"}),
            ]
        ),
        EventSequence(
            person_id="person_002",
            events=[
                Event("diagnosis", datetime(2019, 5, 10), "health", {"condition": "hypertension"}),
                Event("masters", datetime(2012, 5, 15), "education", {"degree": "MS"}),
                Event("job_change", datetime(2014, 1, 5), "occupation"),
                Event("move", datetime(2016, 8, 1), "address", {"city": "Boston"}),
            ]
        ),
        EventSequence(
            person_id="person_003",
            events=[
                Event("treatment", datetime(2021, 2, 20), "health"),
                Event("degree", datetime(2008, 5, 1), "education", {"degree": "PhD"}),
            ]
        ),
    ]
    return sequences


def main():
    """Main example function."""
    print("Life Course Phenotype Extraction Example")
    print("=" * 60)
    print()
    
    if not LIFE_EVENTS_AVAILABLE:
        print("Error: life_events module not available")
        return 1
    
    # Create example sequences
    print("Creating example event sequences...")
    sequences = create_example_sequences()
    print(f"✓ Created {len(sequences)} sequences")
    print()
    
    # Extract phenotypes from each sequence
    print("Extracting phenotypes from sequences:")
    print("-" * 60)
    all_phenotypes = []
    
    for sequence in sequences:
        try:
            phenotypes = extract_phenotypes_from_events(sequence)
            all_phenotypes.append(phenotypes)
            
            print(f"\nPerson: {phenotypes['person_id']}")
            print(f"  Total events: {phenotypes['total_events']}")
            print(f"  Domains: {', '.join(phenotypes['domains'])}")
            print(f"  Domain counts: {phenotypes['domain_counts']}")
            
            if 'health_events' in phenotypes:
                print(f"  Health events: {phenotypes['health_events']}")
                print(f"    Conditions: {phenotypes.get('health_conditions', [])}")
            
            if 'education_events' in phenotypes:
                print(f"  Education events: {phenotypes['education_events']}")
                print(f"    Achievements: {phenotypes.get('education_achievements', [])}")
            
            if 'occupation_events' in phenotypes:
                print(f"  Occupation events: {phenotypes['occupation_events']}")
                print(f"    Changes: {phenotypes.get('occupation_changes', [])}")
            
            if 'event_span_years' in phenotypes:
                print(f"  Event span: {phenotypes['event_span_years']:.2f} years")
                
        except ValidationError as e:
            print(f"✗ Validation error for {sequence.person_id}: {e}")
        except Exception as e:
            print(f"✗ Error processing {sequence.person_id}: {e}")
    
    print()
    print()
    
    # Map events to traits
    print("Mapping events to trait categories:")
    print("-" * 60)
    
    for sequence in sequences:
        try:
            traits = map_events_to_traits(sequence)
            
            print(f"\nPerson: {sequence.person_id}")
            for trait_name, trait_data in traits.items():
                if trait_data['count'] > 0:
                    print(f"  {trait_name}: {trait_data['count']} events")
                    print(f"    Events: {', '.join(trait_data['events'])}")
                    
        except Exception as e:
            print(f"✗ Error mapping traits for {sequence.person_id}: {e}")
    
    print()
    print()
    
    # Aggregate temporal phenotypes
    print("Aggregating temporal phenotypes:")
    print("-" * 60)
    
    try:
        aggregated = aggregate_temporal_phenotypes(sequences, time_window_years=5.0)
        
        print(f"\nAggregate Statistics:")
        print(f"  Total events: {aggregated['aggregates']['total_events']}")
        print(f"  Total people: {aggregated['aggregates']['total_people']}")
        print(f"  Time span: {aggregated['aggregates']['time_span_years']:.2f} years")
        print(f"  Number of time windows: {len(aggregated['time_windows'])}")
        
        print(f"\nTime Windows:")
        for i, window in enumerate(aggregated['time_windows'][:5], 1):  # Show first 5
            start = datetime.fromtimestamp(window['start_time']).strftime('%Y-%m-%d')
            end = datetime.fromtimestamp(window['end_time']).strftime('%Y-%m-%d')
            print(f"  Window {i}: {start} to {end}")
            print(f"    Events: {window['n_events']}, People: {window['n_people']}")
            print(f"    Domains: {window['domain_counts']}")
        
        if len(aggregated['time_windows']) > 5:
            print(f"  ... and {len(aggregated['time_windows']) - 5} more windows")
            
    except ValidationError as e:
        print(f"✗ Validation error: {e}")
        return 1
    except Exception as e:
        print(f"✗ Error aggregating: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    print()
    print("✓ Example completed successfully")
    return 0


if __name__ == "__main__":
    sys.exit(main())


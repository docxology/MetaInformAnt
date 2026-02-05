#!/usr/bin/env python3
"""Example script demonstrating AntWiki JSON data loading with error handling.

This script shows how to:
- Load AntWiki phenotype data from JSON files
- Handle errors gracefully
- Use core utilities for path handling
- Validate data structure
- Process phenotype entries

Usage:
    python3 examples/phenotype/load_antwiki_example.py [--input data/phenotype/antwiki.json]
"""

import argparse
import sys
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.errors import IOError as CoreIOError
from metainformant.core.errors import ValidationError
from metainformant.core.io import write_json
from metainformant.core.paths import expand_and_resolve
from metainformant.phenotype.antwiki import load_antwiki_json


def main():
    """Main example function."""
    parser = argparse.ArgumentParser(description="Load and display AntWiki phenotype data")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/phenotype/antwiki_species.json"),
        help="Path to AntWiki JSON file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output path for processed data",
    )
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="Disable data validation",
    )
    args = parser.parse_args()

    # Expand and resolve input path
    input_path = expand_and_resolve(args.input)

    print(f"Loading AntWiki data from: {input_path}")
    print("-" * 60)

    try:
        # Load data with validation (unless disabled)
        data = load_antwiki_json(input_path, validate=not args.no_validate)

        print(f"✓ Successfully loaded {len(data)} entries")
        print()

        # Display summary
        print("Summary:")
        print(f"  Total entries: {len(data)}")

        # Count entries with different fields
        with_species = sum(1 for e in data if "species" in e)
        with_taxon = sum(1 for e in data if "taxon" in e)
        with_measurements = sum(1 for e in data if "measurements" in e)
        with_traits = sum(1 for e in data if "traits" in e)

        print(f"  Entries with 'species': {with_species}")
        print(f"  Entries with 'taxon': {with_taxon}")
        print(f"  Entries with 'measurements': {with_measurements}")
        print(f"  Entries with 'traits': {with_traits}")
        print()

        # Display first few entries
        print("Sample entries:")
        for i, entry in enumerate(data[:3], 1):
            species = entry.get("species") or entry.get("taxon", "unknown")
            measurements = entry.get("measurements", {})
            traits = entry.get("traits", [])

            print(f"  {i}. {species}")
            print(f"     Measurements: {len(measurements)}")
            print(f"     Traits: {len(traits)}")
            if traits:
                print(f"     Example traits: {', '.join(traits[:3])}")
            print()

        # Save processed data if output specified
        if args.output:
            output_path = expand_and_resolve(args.output)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # Create summary
            summary = {
                "total_entries": len(data),
                "entries": [
                    {
                        "species": e.get("species") or e.get("taxon", "unknown"),
                        "n_measurements": len(e.get("measurements", {})),
                        "n_traits": len(e.get("traits", [])),
                    }
                    for e in data
                ],
            }

            write_json(summary, output_path)
            print(f"✓ Saved summary to: {output_path}")

        return 0

    except FileNotFoundError as e:
        print(f"✗ Error: File not found: {e}")
        print(f"  Please ensure the file exists at: {input_path}")
        return 1

    except ValidationError as e:
        print(f"✗ Validation error: {e}")
        print(f"  The data structure is invalid. Try --no-validate to skip validation.")
        return 1

    except CoreIOError as e:
        print(f"✗ I/O error: {e}")
        print(f"  Failed to read or parse the JSON file.")
        return 1

    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

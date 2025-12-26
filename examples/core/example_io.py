#!/usr/bin/env python3
"""File I/O operations example.

This example demonstrates METAINFORMANT's file input/output utilities for JSON, CSV, and JSONL formats with gzip compression support.

Usage:
    python examples/core/example_io.py

Output:
    output/examples/core/io_example.{json,csv,jsonl,json.gz,jsonl.gz}
"""

from __future__ import annotations

import csv
from pathlib import Path
from metainformant.core import io

def main():
    """Demonstrate file I/O operations."""
    # Setup output directory
    output_dir = Path("output/examples/core")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT I/O Operations Example ===")

    # Sample data for demonstration
    sample_data = {
        "analysis": {
            "type": "dna_sequence_analysis",
            "timestamp": "2024-12-26T10:00:00Z",
            "version": "1.0"
        },
        "sequences": [
            {"id": "seq1", "sequence": "ATCGATCG", "length": 8, "gc_content": 0.5},
            {"id": "seq2", "sequence": "GCTAGCTA", "length": 8, "gc_content": 0.5},
            {"id": "seq3", "sequence": "TTTTAAAA", "length": 8, "gc_content": 0.0}
        ],
        "statistics": {
            "total_sequences": 3,
            "average_length": 8.0,
            "average_gc": 0.333
        }
    }

    # 1. JSON file operations
    print("\n1. JSON file operations...")

    json_file = output_dir / "io_example.json"
    io.dump_json(sample_data, json_file, indent=2)
    print(f"✓ Saved JSON data to: {json_file}")

    # Load it back
    loaded_json = io.load_json(json_file)
    print(f"✓ Loaded JSON data: {len(loaded_json)} top-level keys")

    # 2. Compressed JSON file operations
    print("\n2. Compressed JSON file operations...")

    json_gz_file = output_dir / "io_example.json.gz"
    io.dump_json(sample_data, json_gz_file, indent=2)
    print(f"✓ Saved compressed JSON data to: {json_gz_file}")

    # Load compressed file
    loaded_json_gz = io.load_json(json_gz_file)
    print(f"✓ Loaded compressed JSON data: {len(loaded_json_gz)} top-level keys")

    # 3. CSV file operations
    print("\n3. CSV file operations...")

    # Convert sequence data to CSV format
    csv_data = []
    for seq in sample_data["sequences"]:
        csv_data.append({
            "id": seq["id"],
            "sequence": seq["sequence"],
            "length": seq["length"],
            "gc_content": seq["gc_content"]
        })

    csv_file = output_dir / "io_example.csv"
    io.write_csv(csv_data, csv_file)
    print(f"✓ Saved CSV data to: {csv_file}")

    # Load CSV back
    loaded_csv = io.read_csv(csv_file)
    print(f"✓ Loaded CSV data: {len(loaded_csv)} rows")

    # 4. JSONL (JSON Lines) file operations
    print("\n4. JSONL file operations...")

    # Prepare JSONL data (one JSON object per line)
    jsonl_data = [
        {"step": 1, "operation": "load_data", "status": "success", "duration_ms": 150},
        {"step": 2, "operation": "validate_data", "status": "success", "duration_ms": 50},
        {"step": 3, "operation": "analyze_sequences", "status": "success", "duration_ms": 500},
        {"step": 4, "operation": "generate_report", "status": "success", "duration_ms": 200}
    ]

    jsonl_file = output_dir / "io_example.jsonl"
    io.write_jsonl(jsonl_data, jsonl_file)
    print(f"✓ Saved JSONL data to: {jsonl_file}")

    # Load JSONL back (iterator)
    loaded_jsonl = list(io.read_jsonl(jsonl_file))
    print(f"✓ Loaded JSONL data: {len(loaded_jsonl)} records")

    # 5. Compressed JSONL operations
    print("\n5. Compressed JSONL operations...")

    jsonl_gz_file = output_dir / "io_example.jsonl.gz"
    io.write_jsonl(jsonl_data, jsonl_gz_file)
    print(f"✓ Saved compressed JSONL data to: {jsonl_gz_file}")

    # Load compressed JSONL
    loaded_jsonl_gz = list(io.read_jsonl(jsonl_gz_file))
    print(f"✓ Loaded compressed JSONL data: {len(loaded_jsonl_gz)} records")

    # 6. Demonstrate atomic writes (prevent corruption)
    print("\n6. Atomic write operations...")

    # All dump_json and write_jsonl operations use atomic writes by default
    # This prevents file corruption if the process is interrupted
    print("✓ All write operations use atomic writes (prevent corruption)")

    # 7. File size comparison
    print("\n7. File size comparison...")

    import os
    sizes = {
        "JSON": os.path.getsize(json_file),
        "JSON.gz": os.path.getsize(json_gz_file),
        "CSV": os.path.getsize(csv_file),
        "JSONL": os.path.getsize(jsonl_file),
        "JSONL.gz": os.path.getsize(jsonl_gz_file)
    }

    print("File sizes (bytes):")
    for format_name, size in sizes.items():
        print(f"  {format_name:8}: {size:4}")

    compression_ratio = sizes["JSON.gz"] / sizes["JSON"]
    print(".1f")

    # 8. Create summary report
    print("\n8. Creating summary report...")

    summary = {
        "io_operations_demo": {
            "timestamp": "2024-12-26T10:00:00Z",
            "files_created": {
                "json": str(json_file.relative_to(output_dir)),
                "json_gz": str(json_gz_file.relative_to(output_dir)),
                "csv": str(csv_file.relative_to(output_dir)),
                "jsonl": str(jsonl_file.relative_to(output_dir)),
                "jsonl_gz": str(jsonl_gz_file.relative_to(output_dir))
            },
            "file_sizes_bytes": sizes,
            "compression_ratio": compression_ratio,
            "features_demonstrated": [
                "JSON reading/writing",
                "Gzip compression support",
                "CSV data handling",
                "JSONL (JSON Lines) format",
                "Atomic write operations",
                "Error handling and validation"
            ],
            "data_summary": {
                "sequences_processed": len(sample_data["sequences"]),
                "workflow_steps": len(jsonl_data),
                "total_data_points": sum(len(seq["sequence"]) for seq in sample_data["sequences"])
            }
        }
    }

    summary_file = output_dir / "io_operations_summary.json"
    io.dump_json(summary, summary_file, indent=2)
    print(f"✓ Summary saved to: {summary_file}")

    print("\n=== I/O Operations Example Complete ===")
    print(f"All example files created in: {output_dir}")
    print("\nKey takeaways:")
    print("- Use dump_json/load_json for structured data")
    print("- Use write_jsonl/read_jsonl for streaming large datasets")
    print("- Use .gz extension for automatic compression")
    print("- All writes are atomic to prevent corruption")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Path management and security example.

This example demonstrates METAINFORMANT's path handling utilities for safe file operations, path validation, and directory management.

Usage:
    python examples/core/example_paths.py

Output:
    output/examples/core/paths_example.json
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from metainformant.core import paths, io

def main():
    """Demonstrate path management functionality."""
    # Setup output directory
    output_dir = Path("output/examples/core")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Path Management Example ===")

    # 1. Path expansion and resolution
    print("\n1. Path expansion and resolution...")

    test_paths = [
        "~/Documents/project",  # User home expansion
        "./relative/path",      # Relative path
        "/absolute/path",       # Absolute path
        "output/examples/../core"  # Path with ..
    ]

    expanded_paths = {}
    for test_path in test_paths:
        expanded = paths.expand_and_resolve(test_path)
        expanded_paths[test_path] = str(expanded)
        print(f"  {test_path} -> {expanded}")

    # 2. Path containment checking
    print("\n2. Path containment validation...")

    containment_tests = [
        ("output/examples/core", "output/examples"),
        ("output/examples/core", "output/different"),
        ("/tmp/test", "/tmp"),
        ("/var/log", "/tmp")
    ]

    containment_results = {}
    for path, parent in containment_tests:
        is_contained = paths.is_within(path, parent)
        containment_results[f"{path} in {parent}"] = is_contained
        status = "✓ INSIDE" if is_contained else "✗ OUTSIDE"
        print(f"  {path} {status} {parent}")

    # 3. Safe path validation
    print("\n3. Safe path validation...")

    safe_tests = [
        "data/sequences.fasta",      # Safe relative path
        "output/results.json",       # Safe output path
        "../../../etc/passwd",       # Path traversal attack
        "/etc/passwd",               # Absolute system path
        "data;rm -rf /",             # Command injection
        "data|cat /etc/passwd",      # Pipe injection
        "data && evil_command",      # Command chaining
        "data$HOME",                 # Environment variable
    ]

    safety_results = {}
    for test_path in safe_tests:
        is_safe = paths.is_safe_path(test_path)
        safety_results[test_path] = is_safe
        status = "✓ SAFE" if is_safe else "✗ UNSAFE"
        print(f"  {status}: {test_path}")

    # 4. Directory creation and file preparation
    print("\n4. Directory creation and file preparation...")

    # Create nested directory structure
    test_dirs = [
        "output/examples/core/analysis",
        "output/examples/core/visualizations",
        "output/examples/core/temp"
    ]

    for dir_path in test_dirs:
        created_dir = io.ensure_directory(dir_path)
        print(f"  ✓ Created directory: {created_dir}")

    # Prepare file paths (ensure parent directories exist)
    test_files = [
        "output/examples/core/analysis/sequences.json",
        "output/examples/core/visualizations/plot.png",
        "output/examples/core/temp/cache.tmp"
    ]

    for file_path in test_files:
        file_path_obj = Path(file_path)
        paths.prepare_file_path(file_path_obj)
        # Create empty file to demonstrate
        file_path_obj.touch()
        print(f"  ✓ Prepared file path: {file_path}")

    # 5. File discovery and extension handling
    print("\n5. File discovery and extension handling...")

    # Create some test files
    test_file_dir = output_dir / "test_files"
    test_file_dir.mkdir(exist_ok=True)

    test_files_to_create = [
        "sequences.fasta",
        "variants.vcf",
        "expression.tsv",
        "results.json",
        "report.html",
        "data.csv"
    ]

    for filename in test_files_to_create:
        (test_file_dir / filename).touch()

    # Find files by extension
    extensions_to_find = ["fasta", "json", "csv", "tsv"]
    found_files = {}

    for ext in extensions_to_find:
        files = paths.find_files_by_extension(test_file_dir, ext)
        found_files[ext] = [str(f.relative_to(test_file_dir)) for f in files]
        print(f"  Found {len(files)} .{ext} files: {found_files[ext]}")

    # 6. Filename sanitization
    print("\n6. Filename sanitization...")

    unsafe_names = [
        "file with spaces.txt",
        "file<with>brackets.txt",
        "file:with:colons.txt",
        "file|with|pipes.txt",
        "file;with;semicolons.txt",
        "file?with?questions.txt",
        "file*with*asterisks.txt",
        'file"with"quotes.txt',
        "file/with/slashes.txt",
        "file\\with\\backslashes.txt"
    ]

    sanitized_names = {}
    for unsafe_name in unsafe_names:
        safe_name = paths.sanitize_filename(unsafe_name)
        sanitized_names[unsafe_name] = safe_name
        print(f"  {unsafe_name} -> {safe_name}")

    # 7. Temporary file management
    print("\n7. Temporary file management...")

    # Create temporary files with proper cleanup
    temp_files = []
    temp_dir = Path("output/examples/core/temp")

    for i in range(3):
        temp_file = paths.create_temp_file(
            suffix=f"_example_{i}.tmp",
            prefix="metainformant_",
            directory=temp_dir
        )
        temp_files.append(str(temp_file.relative_to(output_dir)))
        # Write something to demonstrate
        temp_file.write_text(f"Temporary file {i} content")

    print(f"  ✓ Created {len(temp_files)} temporary files in {temp_dir}")

    # 8. Directory size calculation
    print("\n8. Directory size calculation...")

    # Calculate sizes
    sizes = {}
    dirs_to_check = [
        "output/examples/core",
        "output/examples/core/analysis",
        "output/examples/core/temp"
    ]

    for dir_path in dirs_to_check:
        try:
            size = paths.get_directory_size(dir_path)
            sizes[dir_path] = size
            print(f"  {dir_path}: {size} bytes")
        except Exception as e:
            sizes[dir_path] = f"Error: {e}"

    # 9. Create comprehensive summary
    print("\n9. Creating comprehensive summary...")

    summary = {
        "path_management_demo": {
            "timestamp": "2024-12-26T10:00:00Z",
            "features_demonstrated": [
                "Path expansion and resolution",
                "Path containment validation",
                "Safe path validation (security)",
                "Directory creation and file preparation",
                "File discovery by extension",
                "Filename sanitization",
                "Temporary file management",
                "Directory size calculation"
            ],
            "path_operations": {
                "expanded_paths": expanded_paths,
                "containment_checks": containment_results,
                "safety_validation": safety_results,
                "sanitized_filenames": sanitized_names
            },
            "file_system_operations": {
                "directories_created": len(test_dirs),
                "files_prepared": len(test_files),
                "temp_files_created": len(temp_files),
                "files_found_by_extension": found_files,
                "directory_sizes": sizes
            },
            "security_features": [
                "Path traversal attack prevention",
                "Command injection detection",
                "Safe filename sanitization",
                "Containment validation"
            ],
            "output_structure": {
                "main_output_dir": str(output_dir),
                "test_files_dir": str(test_file_dir.relative_to(output_dir)),
                "analysis_dir": "analysis",
                "visualizations_dir": "visualizations",
                "temp_dir": "temp"
            }
        }
    }

    summary_file = output_dir / "paths_example.json"
    io.dump_json(summary, summary_file, indent=2)

    print(f"✓ Comprehensive summary saved to: {summary_file}")

    print("\n=== Path Management Example Complete ===")
    print(f"All example operations completed in: {output_dir}")
    print("\nKey takeaways:")
    print("- Use expand_and_resolve() for robust path handling")
    print("- Use is_within() to validate path containment")
    print("- Use is_safe_path() to prevent security issues")
    print("- Use ensure_directory() and prepare_file_path() for setup")
    print("- Use sanitize_filename() for safe file naming")

if __name__ == "__main__":
    main()

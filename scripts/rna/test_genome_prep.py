#!/usr/bin/env python3
"""Test script for genome_prep module functionality.

This script tests the genome preparation functions to ensure they work correctly.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add src directory to path for imports
script_dir = Path(__file__).parent
repo_root = script_dir.parent.parent
src_dir = repo_root / "src"
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

from metainformant.rna.genome_prep import (
    find_rna_fasta_in_genome_dir,
    get_expected_index_path,
)

def test_find_rna_fasta():
    """Test finding RNA FASTA in genome directory."""
    print("Testing find_rna_fasta_in_genome_dir...")
    
    # Test with a species that has a downloaded genome
    test_cases = [
        ("output/amalgkit/cfloridanus/genome", "GCF_003227725.1"),
        ("output/amalgkit/atta_cephalotes/genome", "GCF_000143395.1"),
        ("output/amalgkit/ooceraea_biroi/genome", "GCF_000956425.1"),
    ]
    
    repo_root = Path(__file__).parent.parent.parent
    for genome_dir_str, accession in test_cases:
        genome_dir = repo_root / genome_dir_str
        if genome_dir.exists():
            print(f"  Checking {genome_dir_str}...")
            rna_fasta = find_rna_fasta_in_genome_dir(genome_dir, accession)
            if rna_fasta:
                print(f"    ✓ Found: {rna_fasta}")
            else:
                print(f"    ✗ Not found")
        else:
            print(f"  Skipping {genome_dir_str} (not found)")

def test_expected_index_path():
    """Test expected index path generation."""
    print("\nTesting get_expected_index_path...")
    
    test_cases = [
        ("output/amalgkit/camponotus_floridanus/work", "Camponotus_floridanus"),
        ("output/amalgkit/pogonomyrmex_barbatus/work", "Pogonomyrmex_barbatus"),
    ]
    
    repo_root = Path(__file__).parent.parent.parent
    for work_dir_str, species_name in test_cases:
        work_dir = repo_root / work_dir_str
        index_path = get_expected_index_path(work_dir, species_name)
        expected_name = f"{species_name}_transcripts.idx"
        print(f"  {species_name}: {index_path.name}")
        assert index_path.name == expected_name, f"Expected {expected_name}, got {index_path.name}"
        print(f"    ✓ Path generation correct")

def main():
    """Run all tests."""
    print("=" * 80)
    print("Genome Prep Module Tests")
    print("=" * 80)
    
    try:
        test_find_rna_fasta()
        test_expected_index_path()
        print("\n" + "=" * 80)
        print("All tests completed successfully!")
        print("=" * 80)
        return 0
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())


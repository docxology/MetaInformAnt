#!/usr/bin/env python3
"""Test runner for all METAINFORMANT examples.

This script discovers and runs all example files in the examples/ directory headlessly,
validates their outputs, and generates a comprehensive test report.

Usage:
    python scripts/test_examples.py [--verbose] [--domain DOMAIN] [--continue-on-error]

Arguments:
    --verbose: Enable verbose output
    --domain: Run examples for specific domain only (e.g., dna, rna, gwas)
    --continue-on-error: Continue testing even if individual examples fail
"""

# Backward compatibility wrapper - imports from the new modular package
from test_examples.main import main

if __name__ == "__main__":
    main()

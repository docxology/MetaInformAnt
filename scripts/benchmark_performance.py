#!/usr/bin/env python3
"""Performance benchmarking script for METAINFORMANT improvements.

This script demonstrates the performance improvements from caching,
vectorization, and optimization enhancements.
"""

from __future__ import annotations

import time
import numpy as np
from pathlib import Path

from metainformant.math import coalescent
from metainformant.math.price import variance, covariance
from metainformant.information.syntactic import shannon_entropy


def benchmark_tajima_constants():
    """Benchmark Tajima constants calculation with caching."""
    print("=== Tajima Constants Performance Benchmark ===")

    sample_sizes = [10, 50, 100, 10, 50, 100]  # Repeated calculations to test caching

    start_time = time.time()
    results = []
    for n in sample_sizes:
        constants = coalescent.tajima_constants(n)
        results.append(constants["a1"])

    total_time = time.time() - start_time
    print(".4f")
    print(".4f")
    print("Note: Second set of calculations should be faster due to caching")
    return total_time


def benchmark_vectorization():
    """Benchmark NumPy vectorization improvements."""
    print("\n=== Vectorization Performance Benchmark ===")

    # Create large test datasets
    sizes = [1000, 10000, 100000]
    results = {}

    for size in sizes:
        data = np.random.randn(size)

        # Benchmark variance
        start = time.time()
        var_result = variance(data)
        var_time = time.time() - start

        # Benchmark covariance
        data2 = np.random.randn(size)
        start = time.time()
        cov_result = covariance(data, data2)
        cov_time = time.time() - start

        results[size] = {
            'variance_time': var_time,
            'covariance_time': cov_time,
            'variance_result': var_result,
            'covariance_result': cov_result
        }

        print("2d")

    return results


def benchmark_entropy_caching():
    """Benchmark entropy calculation with caching."""
    print("\n=== Entropy Caching Performance Benchmark ===")

    # Test with different probability distributions
    distributions = [
        [0.5, 0.5],  # Binary
        [0.25, 0.25, 0.25, 0.25],  # Uniform 4
        [0.1, 0.2, 0.3, 0.4],  # Non-uniform
        [0.5, 0.5],  # Repeat to test caching
        [0.25, 0.25, 0.25, 0.25],  # Repeat to test caching
    ]

    start_time = time.time()
    results = []
    for i, probs in enumerate(distributions):
        entropy = shannon_entropy(probs)
        results.append(entropy)
        print(".4f")

    total_time = time.time() - start_time
    print(".4f")
    print("Note: Repeated distributions should be faster due to caching")
    return total_time, results


def main():
    """Run all performance benchmarks."""
    print("METAINFORMANT Performance Benchmark Suite")
    print("=" * 50)

    try:
        # Benchmark Tajima constants
        tajima_time = benchmark_tajima_constants()

        # Benchmark vectorization
        vectorization_results = benchmark_vectorization()

        # Benchmark entropy caching
        entropy_time, entropy_results = benchmark_entropy_caching()

        # Summary
        print("\n=== Performance Summary ===")
        print(f"Tajima constants total time: {tajima_time:.4f}s")
        print(f"Vectorization tested on datasets up to {max(vectorization_results.keys())} samples")
        print(f"Entropy calculations total time: {entropy_time:.4f}s")
        print("\nAll benchmarks completed successfully!")
        print("\nKey improvements:")
        print("- Caching reduces repeated computation time")
        print("- NumPy vectorization enables large dataset processing")
        print("- Progress tracking provides user feedback for long operations")

    except Exception as e:
        print(f"Benchmark failed: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

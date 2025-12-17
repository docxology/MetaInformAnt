#!/bin/bash
# METAINFORMANT Unified Test Runner
# Consolidates functionality from run_tests.sh, uv_test.sh, and uv_test_optimized.sh
# Supports all test modes and execution options

set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"

# Default values (can be overridden by environment variables)
COVERAGE=${COVERAGE:-false}
PARALLEL=${PARALLEL:-false}
VERBOSE=${VERBOSE:-false}
OUTPUT_DIR=${OUTPUT_DIR:-"output"}
TEST_PATTERN=${TEST_PATTERN:-""}
DISABLE_WARNINGS=${DISABLE_WARNINGS:-false}

# Test mode (will be set by command line args)
TEST_MODE=""

usage() {
    cat << EOF
METAINFORMANT Unified Test Runner - REAL IMPLEMENTATIONS ONLY

This test suite uses REAL external APIs, databases, and file systems.
NO mocking, faking, or stubbing is allowed per project policy.

Usage: $0 [OPTIONS] [--mode MODE]

Options:
    -h, --help              Show this help message
    --mode MODE             Test mode (see Modes section below)
    --coverage              Enable coverage reporting
    --parallel              Enable parallel test execution
    --verbose               Verbose output (-vv)
    --pattern PATTERN       Run tests matching pattern (-k PATTERN)
    --clean                 Clean previous test artifacts
    --disable-warnings      Disable pytest warnings
    --output-dir DIR        Set output directory (default: output)

Modes:
    ultra-fast      Core functionality only (~5s)
    fast            Essential tests (~15s) [DEFAULT]
    coverage        Full coverage analysis (~1-2min)
    coverage-fast   Fast coverage on core modules (~20s)
    parallel        Parallel execution with coverage
    network         Include network-dependent tests (REAL API calls)
    external        Include external CLI tool tests
    integration     Integration tests (~30s)
    smoke           Basic smoke tests (~3s)
    network-only    Network analysis tests only
    ml-only         Machine learning tests only
    benchmark       Performance benchmarks
    all             Full test suite (~2-5min)

Examples:
    $0                           # Run fast tests
    $0 --mode coverage           # Full coverage analysis
    $0 --mode network --parallel # Parallel network tests
    $0 --pattern "test_core_*"   # Run core tests only
    $0 --clean                   # Clean and run default tests

Environment Variables:
    COVERAGE=true        Enable coverage
    PARALLEL=true        Enable parallel execution
    VERBOSE=true         Verbose output
    OUTPUT_DIR=path      Set output directory
    TEST_PATTERN=pattern Set test pattern

EOF
}

# Clean previous test artifacts
clean_artifacts() {
    print_status "INFO" "Cleaning previous test artifacts..."
    rm -rf "${OUTPUT_DIR}/coverage_html" "${OUTPUT_DIR}/coverage_html_fast"
    rm -rf "${OUTPUT_DIR}/coverage.xml" "${OUTPUT_DIR}/.coverage"
    rm -rf ".pytest_cache" "__pycache__" "*/__pycache__" "*.pyc"
    find . -name "*.pyc" -delete 2>/dev/null || true
    find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
    print_status "OK" "Cleaned test artifacts"
}

# Setup test environment
setup_test_env() {
    # Setup basic environment (from _common.sh)
    setup_environment

    # Create output directories
    mkdir -p "$OUTPUT_DIR"
    export COVERAGE_FILE="$OUTPUT_DIR/.coverage"
}

# Build pytest arguments based on mode and options
build_pytest_args() {
    local args=()

    # Basic configuration
    if [[ "$VERBOSE" == "true" ]]; then
        args+=("-vv" "--tb=short")
    else
        args+=("-v" "--tb=short")
    fi

    # Coverage options
    if [[ "$COVERAGE" == "true" ]] || [[ "$TEST_MODE" == "coverage" ]] || [[ "$TEST_MODE" == "coverage-fast" ]]; then
        args+=("--cov=src/metainformant")
        if [[ "$TEST_MODE" == "coverage-fast" ]]; then
            args+=("--cov=src/metainformant/core" "--cov=src/metainformant/simulation")
            args+=("--cov-report=html:$OUTPUT_DIR/coverage_html_fast")
        else
            args+=("--cov-report=html:$OUTPUT_DIR/coverage_html")
            args+=("--cov-report=xml:$OUTPUT_DIR/coverage.xml")
        fi
        args+=("--cov-report=term-missing")
        args+=("--cov-fail-under=85")
    fi

    # Parallel execution
    if [[ "$PARALLEL" == "true" ]] || [[ "$TEST_MODE" == "parallel" ]]; then
        # Check if pytest-xdist is available
        if uv run python -c "import pytest_xdist" >/dev/null 2>&1; then
            args+=("-n" "auto")
        else
            print_status "WARN" "pytest-xdist not available, running sequentially"
        fi
    fi

    # Warnings
    if [[ "$DISABLE_WARNINGS" == "true" ]]; then
        args+=("--disable-warnings")
    fi

    # Pattern matching
    if [[ -n "$TEST_PATTERN" ]]; then
        args+=("-k" "$TEST_PATTERN")
    fi

    # Mode-specific options
    case "$TEST_MODE" in
        "ultra-fast")
            args+=("--timeout=5" "-x" "-q" "--disable-warnings")
            ;;
        "fast")
            args+=("--timeout=10" "-x" "-q")
            ;;
        "coverage-fast")
            args+=("--timeout=15" "--cov-fail-under=0")
            ;;
        "integration")
            args+=("--timeout=30" "--disable-warnings")
            ;;
        "smoke"|"network-only"|"ml-only")
            args+=("--timeout=5")
            ;;
        "network"|"external")
            args+=("--timeout=20")
            ;;
        "benchmark")
            args+=("--benchmark-json=$OUTPUT_DIR/benchmarks/results.json")
            ;;
    esac

    # Additional pytest options
    args+=("--strict-markers")
    args+=("--durations=10")  # Show 10 slowest tests

    echo "${args[@]}"
}

# Get test files for specific mode
get_test_files() {
    case "$TEST_MODE" in
        "ultra-fast")
            echo "tests/test_core_functionality.py"
            ;;
        "fast")
            echo "tests/test_core_functionality.py tests/test_core_comprehensive.py tests/test_dna_comprehensive.py"
            ;;
        "coverage-fast")
            echo "tests/test_core_functionality.py tests/test_core_*.py"
            ;;
        "integration")
            echo "tests/test_integration_comprehensive.py"
            ;;
        "smoke")
            echo "tests/test_core_functionality.py::TestCoreIO::test_json_operations tests/test_core_functionality.py::TestDNAAnalysis::test_gc_content_calculation tests/test_core_functionality.py::TestSimulation::test_sequence_generation"
            ;;
        "network-only")
            echo "tests/test_networks_*.py"
            ;;
        "ml-only")
            echo "tests/test_ml_*.py"
            ;;
        "network")
            echo "tests/ -m network"
            ;;
        "external")
            echo "tests/ -m external_tool"
            ;;
        "benchmark")
            echo "tests/ -m benchmark"
            ;;
        *)
            echo "tests/"
            ;;
    esac
}

# Sync dependencies for test mode
sync_test_deps() {
    local test_type="fast"
    case "$TEST_MODE" in
        "network"|"network-only")
            test_type="network"
            ;;
        "external")
            test_type="external"
            ;;
        "all"|"coverage"|"parallel")
            test_type="all"
            ;;
    esac

    sync_dependencies "$test_type"
}

# Run the tests
run_test_suite() {
    local pytest_cmd
    pytest_cmd=$(get_pytest_cmd)

    local pytest_args
    pytest_args=$(build_pytest_args)

    local test_files
    test_files=$(get_test_files)

    # Show what we're running
    echo "Mode: $TEST_MODE"
    echo "Command: $pytest_cmd"
    echo "Args: $pytest_args"
    echo "Files: $test_files"
    echo ""

    # Set PYTHONPATH
    export PYTHONPATH="${PWD}/src:${PYTHONPATH:-}"

    # Run the tests
    if [[ "$TEST_MODE" == "benchmark" ]]; then
        # Special handling for benchmarks
        if $pytest_cmd $test_files --benchmark-json="$OUTPUT_DIR/benchmarks/results.json"; then
            print_status "OK" "Benchmark tests completed successfully"
            return 0
        else
            print_status "ERROR" "Some benchmark tests failed"
            return 1
        fi
    else
        # Regular test execution
        if $pytest_cmd $pytest_args $test_files; then
            print_status "OK" "Tests completed successfully"
            return 0
        else
            print_status "ERROR" "Some tests failed"
            return 1
        fi
    fi
}

# Generate coverage report
generate_coverage_report() {
    if [[ "$COVERAGE" == "true" ]] || [[ "$TEST_MODE" == "coverage" ]] || [[ "$TEST_MODE" == "coverage-fast" ]] || [[ "$TEST_MODE" == "parallel" ]]; then
        if [[ -f "$OUTPUT_DIR/.coverage" ]]; then
            print_status "INFO" "Generating coverage reports..."

            # Generate HTML report
            coverage html --directory="$OUTPUT_DIR/coverage_html" --data-file="$OUTPUT_DIR/.coverage" 2>/dev/null || true

            # Generate console summary
            echo -e "${BLUE}Coverage Summary:${NC}"
            coverage report --data-file="$OUTPUT_DIR/.coverage" --show-missing 2>/dev/null || true

            print_status "OK" "Coverage reports generated in $OUTPUT_DIR/coverage_html"
        fi
    fi
}

# Main function
main() {
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                usage
                exit 0
                ;;
            --mode)
                TEST_MODE="${2:-fast}"
                shift 2
                ;;
            --coverage)
                COVERAGE=true
                shift
                ;;
            --parallel)
                PARALLEL=true
                shift
                ;;
            --verbose)
                VERBOSE=true
                shift
                ;;
            --pattern)
                TEST_PATTERN="$2"
                shift 2
                ;;
            --clean)
                clean_artifacts
                exit 0
                ;;
            --disable-warnings)
                DISABLE_WARNINGS=true
                shift
                ;;
            --output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            *)
                print_status "ERROR" "Unknown option: $1"
                usage
                exit 1
                ;;
        esac
    done

    # Set default mode
    TEST_MODE="${TEST_MODE:-fast}"

    # Check prerequisites
    check_uv

    # Setup environment
    setup_test_env

    # Sync dependencies
    sync_test_deps

    # Run tests
    print_status "INFO" "Running tests in mode: $TEST_MODE"
    if run_test_suite; then
        generate_coverage_report
        echo ""
        print_status "OK" "All tests in mode '$TEST_MODE' completed successfully!"
        if [[ "$COVERAGE" == "true" ]] || [[ "$TEST_MODE" == "coverage" ]] || [[ "$TEST_MODE" == "coverage-fast" ]]; then
            echo -e "${BLUE}📊 View coverage report: $OUTPUT_DIR/coverage_html/index.html${NC}"
        fi
        exit 0
    else
        echo ""
        print_status "ERROR" "Some tests in mode '$TEST_MODE' failed. Check output above for details."
        exit 1
    fi
}

# Run main function
main "$@"

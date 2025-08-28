#!/usr/bin/env bash
# METAINFORMANT Test Runner Script
# Comprehensive testing with coverage, parallel execution, and reporting

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
COVERAGE=${COVERAGE:-true}
PARALLEL=${PARALLEL:-true}
VERBOSE=${VERBOSE:-false}
FAST_ONLY=${FAST_ONLY:-false}
NETWORK_TESTS=${NETWORK_TESTS:-false}
OUTPUT_DIR=${OUTPUT_DIR:-"output"}
TEST_PATTERN=${TEST_PATTERN:-""}

# Functions
print_header() {
    echo -e "${BLUE}============================================${NC}"
    echo -e "${BLUE}  METAINFORMANT Test Suite Runner${NC}"
    echo -e "${BLUE}  REAL IMPLEMENTATIONS ONLY - NO MOCKS${NC}"
    echo -e "${BLUE}============================================${NC}"
}

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

METAINFORMANT Test Runner - REAL IMPLEMENTATIONS ONLY
This test suite uses REAL external APIs, databases, and file systems.
NO mocking, faking, or stubbing is allowed per project policy.

Options:
    --help                Show this help message
    --fast                Run only fast tests (skip slow/network tests)
    --no-coverage         Skip coverage reporting
    --no-parallel         Run tests sequentially
    --network             Include network-dependent tests (REAL API calls)
    --verbose             Verbose output
    --pattern PATTERN     Run tests matching pattern
    --clean               Clean previous test artifacts
    --report-only         Generate coverage report from existing data

Examples:
    $0                    # Run all tests with real implementations
    $0 --fast             # Skip network-dependent real API tests
    $0 --pattern "test_core_*"  # Run core module tests (no external deps)
    $0 --network          # Include real API tests (requires connectivity)
    $0 --clean            # Clean and run all real implementation tests

Network Requirements:
    - UniProt API: https://rest.uniprot.org (real HTTP requests)
    - PDB Downloads: https://files.rcsb.org (actual file downloads)
    - AlphaFold: https://alphafold.ebi.ac.uk (real model fetching)
    - NCBI Entrez: requires NCBI_EMAIL environment variable

Environment Variables:
    COVERAGE=false        Disable coverage
    PARALLEL=false        Disable parallel execution
    NETWORK_TESTS=true    Enable real network API tests
    OUTPUT_DIR=path       Set output directory
    NCBI_EMAIL=email      Required for NCBI Entrez real API calls
EOF
}

clean_artifacts() {
    echo -e "${YELLOW}Cleaning previous test artifacts...${NC}"
    rm -rf "${OUTPUT_DIR}/coverage_html" "${OUTPUT_DIR}/coverage.xml" "${OUTPUT_DIR}/.coverage"
    rm -rf ".pytest_cache" "__pycache__" "*/__pycache__" "*.pyc"
    find . -name "*.pyc" -delete 2>/dev/null || true
    find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
    echo -e "${GREEN}✓ Cleaned test artifacts${NC}"
}

check_dependencies() {
    echo -e "${YELLOW}Checking test dependencies...${NC}"
    
    # Check if pytest is available
    if ! command -v pytest &> /dev/null; then
        echo -e "${RED}✗ pytest not found. Install with: pip install pytest${NC}"
        exit 1
    fi
    
    # Check if coverage is needed and available
    if [[ "$COVERAGE" == "true" ]] && ! command -v pytest-cov &> /dev/null; then
        echo -e "${YELLOW}Warning: pytest-cov not found. Coverage will be disabled.${NC}"
        COVERAGE=false
    fi
    
    echo -e "${GREEN}✓ Dependencies checked${NC}"
}

setup_output_dir() {
    mkdir -p "$OUTPUT_DIR"
    export COVERAGE_FILE="$OUTPUT_DIR/.coverage"
}

build_pytest_args() {
    local args=()
    
    # Basic configuration
    args+=("-v")
    
    # Coverage options
    if [[ "$COVERAGE" == "true" ]]; then
        args+=("--cov=src/metainformant")
        args+=("--cov-report=html:$OUTPUT_DIR/coverage_html")
        args+=("--cov-report=xml:$OUTPUT_DIR/coverage.xml")
        args+=("--cov-report=term-missing")
        args+=("--cov-fail-under=85")
    fi
    
    # Parallel execution
    if [[ "$PARALLEL" == "true" ]] && command -v pytest-xdist &> /dev/null; then
        args+=("-n" "auto")
    fi
    
    # Verbose output
    if [[ "$VERBOSE" == "true" ]]; then
        args+=("-vv" "--tb=short")
    fi
    
    # Test selection
    if [[ "$FAST_ONLY" == "true" ]]; then
        args+=("-m" "not slow and not network and not external_tool")
    elif [[ "$NETWORK_TESTS" == "false" ]]; then
        args+=("-m" "not network")
    fi
    
    # Pattern matching
    if [[ -n "$TEST_PATTERN" ]]; then
        args+=("-k" "$TEST_PATTERN")
    fi
    
    # Additional pytest options
    args+=("--strict-markers")
    args+=("--tb=short")
    args+=("--durations=10")  # Show 10 slowest tests
    
    echo "${args[@]}"
}

run_tests() {
    local pytest_args
    pytest_args=($(build_pytest_args))
    
    echo -e "${YELLOW}Running tests with options: ${pytest_args[*]}${NC}"
    echo ""
    
    # Set environment variables
    export PYTHONPATH="${PWD}/src:${PYTHONPATH:-}"
    
    # Run the tests
    if pytest "${pytest_args[@]}" tests/; then
        echo -e "${GREEN}✓ Tests passed successfully${NC}"
        return 0
    else
        echo -e "${RED}✗ Some tests failed${NC}"
        return 1
    fi
}

generate_coverage_report() {
    if [[ "$COVERAGE" == "true" ]] && [[ -f "$OUTPUT_DIR/.coverage" ]]; then
        echo -e "${YELLOW}Generating coverage reports...${NC}"
        
        # Generate HTML report
        coverage html --directory="$OUTPUT_DIR/coverage_html" --data-file="$OUTPUT_DIR/.coverage"
        
        # Generate console summary
        echo -e "${BLUE}Coverage Summary:${NC}"
        coverage report --data-file="$OUTPUT_DIR/.coverage" --show-missing
        
        echo -e "${GREEN}✓ Coverage reports generated in $OUTPUT_DIR/coverage_html${NC}"
    fi
}

run_specific_test_suites() {
    echo -e "${BLUE}Running domain-specific test suites...${NC}"
    
    local suites=(
        "core:tests/test_core_*.py"
        "dna:tests/test_dna_*.py"
        "rna:tests/test_rna_*.py"
        "protein:tests/test_protein_*.py"
        "math:tests/test_math_*.py"
        "simulation:tests/test_simulation.py"
    )
    
    local results=()
    
    for suite in "${suites[@]}"; do
        local name="${suite%%:*}"
        local pattern="${suite##*:}"
        
        echo -e "${YELLOW}Running $name tests...${NC}"
        
        if pytest -v --tb=short $pattern 2>/dev/null; then
            results+=("$name:PASS")
            echo -e "${GREEN}✓ $name tests passed${NC}"
        else
            results+=("$name:FAIL")
            echo -e "${RED}✗ $name tests failed${NC}"
        fi
    done
    
    echo -e "${BLUE}Test Suite Summary:${NC}"
    for result in "${results[@]}"; do
        local name="${result%%:*}"
        local status="${result##*:}"
        if [[ "$status" == "PASS" ]]; then
            echo -e "  ${GREEN}✓ $name${NC}"
        else
            echo -e "  ${RED}✗ $name${NC}"
        fi
    done
}

run_performance_tests() {
    echo -e "${YELLOW}Running performance benchmarks...${NC}"
    
    # Run a subset of tests with timing information
    pytest --durations=0 --tb=no -q tests/test_simulation.py tests/test_math_*.py
}

generate_test_report() {
    local report_file="$OUTPUT_DIR/test_report.md"
    
    cat > "$report_file" << EOF
# METAINFORMANT Test Report

Generated: $(date)

## Test Environment
- Python: $(python --version)
- Pytest: $(pytest --version | head -n1)
- Coverage: $(coverage --version 2>/dev/null || echo "Not available")

## Test Execution Summary

EOF
    
    if [[ -f "$OUTPUT_DIR/coverage.xml" ]]; then
        echo "## Coverage Summary" >> "$report_file"
        echo "" >> "$report_file"
        # Extract coverage percentage from XML report
        if command -v xmllint &> /dev/null; then
            local coverage_pct
            coverage_pct=$(xmllint --xpath "//coverage/@line-rate" "$OUTPUT_DIR/coverage.xml" 2>/dev/null | sed 's/line-rate="\([^"]*\)"/\1/' || echo "0")
            echo "- Line Coverage: $(echo "$coverage_pct * 100" | bc -l | cut -d. -f1)%" >> "$report_file"
        fi
        echo "- Detailed report: [HTML Coverage Report](coverage_html/index.html)" >> "$report_file"
        echo "" >> "$report_file"
    fi
    
    echo -e "${GREEN}✓ Test report generated: $report_file${NC}"
}

main() {
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --help)
                print_usage
                exit 0
                ;;
            --fast)
                FAST_ONLY=true
                shift
                ;;
            --no-coverage)
                COVERAGE=false
                shift
                ;;
            --no-parallel)
                PARALLEL=false
                shift
                ;;
            --network)
                NETWORK_TESTS=true
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
            --report-only)
                setup_output_dir
                generate_coverage_report
                generate_test_report
                exit 0
                ;;
            --performance)
                setup_output_dir
                run_performance_tests
                exit 0
                ;;
            --suites)
                run_specific_test_suites
                exit 0
                ;;
            *)
                echo -e "${RED}Unknown option: $1${NC}"
                print_usage
                exit 1
                ;;
        esac
    done
    
    print_header
    setup_output_dir
    check_dependencies
    
    # Run the main test suite
    if run_tests; then
        generate_coverage_report
        generate_test_report
        
        echo ""
        echo -e "${GREEN}🎉 All tests completed successfully!${NC}"
        echo -e "${BLUE}📊 View detailed coverage report: $OUTPUT_DIR/coverage_html/index.html${NC}"
        echo -e "${BLUE}📋 Test report: $OUTPUT_DIR/test_report.md${NC}"
        exit 0
    else
        echo ""
        echo -e "${RED}❌ Some tests failed. Check the output above for details.${NC}"
        exit 1
    fi
}

# Run main function with all arguments
main "$@"

#!/bin/bash
# METAINFORMANT Unified Verification Script
# Merges verify_uv_setup.sh and verify_test_deps.sh with mode support
# Verifies UV setup, test dependencies, and environment configuration

set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"

# Default values
VERIFY_MODE="all"
FIX_MODE=0
TEST_TYPE="all"

usage() {
    cat << EOF
METAINFORMANT Unified Verification Script

Verifies UV setup, test dependencies, and environment configuration.

Usage: $0 [OPTIONS] [--mode MODE]

Options:
    -h, --help              Show this help message
    --mode MODE             Verification mode (see Modes section below)
    --fix                   Attempt to fix missing dependencies automatically
    --test-type TYPE        Test type for dependency verification (fast, network, external, all)

Modes:
    setup                   Verify UV installation and environment setup
    deps                    Verify test dependencies and packages
    all                     Full verification (setup + deps) [DEFAULT]

Examples:
    $0                           # Full verification
    $0 --mode setup             # Check UV setup only
    $0 --mode deps --fix        # Check and fix test dependencies
    $0 --mode deps --test-type fast  # Check fast test dependencies only

EOF
}

# Status tracking
ISSUES_FOUND=0
FIXES_APPLIED=0

# Override print_status to track issues
print_status() {
    local status="$1"
    local message="$2"
    case "$status" in
        "OK")
            echo -e "${GREEN}✅${NC} $message"
            ;;
        "WARN")
            echo -e "${YELLOW}⚠️${NC}  $message"
            ;;
        "ERROR")
            echo -e "${RED}❌${NC} $message"
            ISSUES_FOUND=1
            ;;
        "INFO")
            echo -e "${BLUE}ℹ️${NC}  $message"
            ;;
        "FIXED")
            echo -e "${GREEN}🔧${NC} $message"
            FIXES_APPLIED=1
            ;;
        *)
            echo "$message"
            ;;
    esac
}

# Check UV availability (from setup verification)
check_uv() {
    print_status "INFO" "Checking UV availability..."

    if ! command -v uv >/dev/null 2>&1; then
        print_status "ERROR" "uv is not installed or not in PATH"
        echo "       Install uv from: https://github.com/astral-sh/uv"
        echo "       curl -LsSf https://astral.sh/uv/install.sh | sh"
        return 1
    fi

    local uv_version
    uv_version=$(uv --version 2>/dev/null || echo "unknown")
    print_status "OK" "uv available: $uv_version"
    return 0
}

# Check Python availability
check_python() {
    print_status "INFO" "Checking Python environment..."

    if ! uv run python --version >/dev/null 2>&1; then
        print_status "ERROR" "Python not available in uv environment"
        if [[ $FIX_MODE -eq 1 ]]; then
            print_status "INFO" "Attempting to install Python 3.11..."
            uv python install 3.11
            if uv run python --version >/dev/null 2>&1; then
                print_status "FIXED" "Python 3.11 installed"
            else
                print_status "ERROR" "Failed to install Python"
                return 1
            fi
        else
            return 1
        fi
    fi

    local python_version
    python_version=$(uv run python --version 2>&1 || echo "unknown")
    print_status "OK" "Python available: $python_version"
    return 0
}

# Detect filesystem type and configuration
check_filesystem() {
    print_status "INFO" "Detecting filesystem type..."

    if detect_filesystem; then
        print_status "OK" "Filesystem: standard (supports symlinks)"
    else
        print_status "WARN" "Filesystem: FAT (limited symlink support)"
    fi

    # Check UV cache configuration
    configure_uv_cache
}

# Check virtual environment
check_venv() {
    print_status "INFO" "Checking virtual environment..."

    local venv_dir
    venv_dir=$(get_venv_dir)

    if [[ -d "$venv_dir" ]]; then
        if [[ -f "$venv_dir/bin/python3" ]]; then
            local python_version
            python_version=$("$venv_dir/bin/python3" --version 2>&1 || echo "unknown")
            print_status "OK" "Virtual environment found at: $venv_dir ($python_version)"
        else
            print_status "WARN" "Directory exists but Python not found: $venv_dir"
        fi
    else
        print_status "WARN" "Virtual environment not found (expected: $venv_dir)"
        if [[ $FIX_MODE -eq 1 ]]; then
            print_status "INFO" "Creating virtual environment..."
            uv venv "$venv_dir"
            if [[ -d "$venv_dir" ]]; then
                print_status "FIXED" "Virtual environment created"
            else
                print_status "ERROR" "Failed to create virtual environment"
            fi
        fi
    fi
}

# Check pytest and core testing dependencies
check_pytest_deps() {
    local test_type="$1"
    print_status "INFO" "Checking pytest and core testing dependencies..."

    # Core pytest
    if ! uv run python -c "import pytest; print('pytest available')" >/dev/null 2>&1; then
        print_status "ERROR" "pytest not available"
        if [[ $FIX_MODE -eq 1 ]]; then
            print_status "INFO" "Syncing test-fast dependencies..."
            uv sync --group test-fast
            if uv run python -c "import pytest" >/dev/null 2>&1; then
                print_status "FIXED" "pytest installed"
            else
                print_status "ERROR" "Failed to install pytest"
                return 1
            fi
        else
            return 1
        fi
    fi

    print_status "OK" "pytest available"

    # Check additional plugins based on test type
    local plugins_to_check=("pytest_cov")
    case "$test_type" in
        "network")
            plugins_to_check+=("requests")
            ;;
        "external")
            plugins_to_check+=()
            ;;
        "all")
            plugins_to_check+=("pytest_xdist" "pytest_benchmark" "pytest_asyncio" "pytest_clarity")
            ;;
    esac

    for plugin in "${plugins_to_check[@]}"; do
        if ! uv run python -c "import $plugin" >/dev/null 2>&1; then
            print_status "ERROR" "$plugin not available"
            if [[ $FIX_MODE -eq 1 ]]; then
                case "$test_type" in
                    "fast")
                        uv sync --group test-fast
                        ;;
                    "network")
                        uv sync --group test-network
                        ;;
                    "external")
                        uv sync --group test-external
                        ;;
                    "all")
                        uv sync --group test-all
                        ;;
                esac
                if uv run python -c "import $plugin" >/dev/null 2>&1; then
                    print_status "FIXED" "$plugin installed"
                else
                    print_status "ERROR" "Failed to install $plugin"
                fi
            fi
        else
            print_status "OK" "$plugin available"
        fi
    done

    return 0
}

# Check optional scientific dependencies
check_scientific_deps() {
    local test_type="$1"
    print_status "INFO" "Checking optional scientific dependencies..."

    local deps_to_check=("numpy" "pandas" "scipy" "matplotlib")

    case "$test_type" in
        "all")
            deps_to_check+=("scikit-learn" "networkx" "biopython" "seaborn")
            ;;
    esac

    for dep in "${deps_to_check[@]}"; do
        if ! uv run python -c "import $dep" >/dev/null 2>&1; then
            print_status "WARN" "$dep not available (optional)"
            if [[ $FIX_MODE -eq 1 ]]; then
                print_status "INFO" "Installing scientific dependencies..."
                uv sync
                if uv run python -c "import $dep" >/dev/null 2>&1; then
                    print_status "FIXED" "$dep installed"
                fi
            fi
        else
            print_status "OK" "$dep available"
        fi
    done

    return 0
}

# Check external CLI tools
check_external_tools() {
    local test_type="$1"
    print_status "INFO" "Checking external CLI tools..."

    local tools_to_check=("muscle")

    case "$test_type" in
        "external"|"all")
            tools_to_check+=("amalgkit" "seqkit" "sra-tools")
            ;;
    esac

    for tool in "${tools_to_check[@]}"; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            print_status "WARN" "$tool not available on PATH"
            case "$tool" in
                "amalgkit")
                    echo "       Install from: https://github.com/jnarun/Amalgkit"
                    ;;
                "muscle")
                    echo "       Install from: https://github.com/rcedgar/muscle"
                    ;;
                "seqkit")
                    echo "       Install from: https://github.com/shenwei356/seqkit"
                    ;;
                "sra-tools")
                    echo "       Install from: https://github.com/ncbi/sra-tools"
                    ;;
            esac
        else
            print_status "OK" "$tool available"
        fi
    done

    return 0
}

# Check test data and directories
check_test_data() {
    print_status "INFO" "Checking test data and directories..."

    # Check test data directory
    if [[ ! -d "tests/data" ]]; then
        print_status "ERROR" "tests/data directory not found"
        if [[ $FIX_MODE -eq 1 ]]; then
            mkdir -p tests/data
            print_status "FIXED" "tests/data directory created"
        fi
    else
        print_status "OK" "tests/data directory exists"
    fi

    # Check output directory
    if [[ ! -d "output" ]]; then
        mkdir -p output
        print_status "FIXED" "output directory created"
    else
        print_status "OK" "output directory exists"
    fi

    # Check for basic test files
    local basic_files=("pyproject.toml" "tests/conftest.py")
    for file in "${basic_files[@]}"; do
        if [[ ! -f "$file" ]]; then
            print_status "ERROR" "$file not found"
        else
            print_status "OK" "$file exists"
        fi
    done

    return 0
}

# Check network connectivity for network tests
check_network_connectivity() {
    local test_type="$1"
    if [[ "$test_type" != "network" && "$test_type" != "all" ]]; then
        return 0
    fi

    print_status "INFO" "Checking network connectivity..."

    # Test basic connectivity
    if ! uv run python -c "
import requests
try:
    response = requests.get('https://httpbin.org/status/200', timeout=10)
    print('Network available' if response.status_code == 200 else 'Network issue')
except:
    print('Network unavailable')
" >/dev/null 2>&1; then
        print_status "WARN" "Network connectivity check failed"
        echo "       Network tests may be skipped"
    else
        print_status "OK" "Network connectivity available"
    fi

    return 0
}

# Setup verification mode
verify_setup() {
    print_status "INFO" "Running setup verification..."

    check_uv
    check_python
    check_filesystem
    check_venv

    # Test basic UV commands
    if [[ -d "$(get_venv_dir)" ]]; then
        if uv pip list --python "$(get_venv_dir)/bin/python3" >/dev/null 2>&1; then
            print_status "OK" "uv pip list command works"
        else
            print_status "WARN" "uv pip list command failed"
        fi
    fi
}

# Dependencies verification mode
verify_deps() {
    print_status "INFO" "Running dependencies verification..."

    check_pytest_deps "$TEST_TYPE"
    check_scientific_deps "$TEST_TYPE"
    check_external_tools "$TEST_TYPE"
    check_test_data
    check_network_connectivity "$TEST_TYPE"
}

# Main verification function
main() {
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                usage
                exit 0
                ;;
            --mode)
                VERIFY_MODE="${2:-all}"
                shift 2
                ;;
            --fix)
                FIX_MODE=1
                shift
                ;;
            --test-type)
                TEST_TYPE="${2:-all}"
                shift 2
                ;;
            *)
                print_status "ERROR" "Unknown option: $1"
                usage
                exit 1
                ;;
        esac
    done

    echo "🔍 METAINFORMANT Verification ($VERIFY_MODE mode)"
    echo "=============================================="
    echo ""

    # Setup basic environment
    setup_environment

    # Run verification based on mode
    case "$VERIFY_MODE" in
        "setup")
            verify_setup
            ;;
        "deps")
            verify_deps
            ;;
        "all"|*)
            verify_setup
            echo ""
            verify_deps
            ;;
    esac

    echo ""
    echo "📊 Verification Summary"
    echo "======================"

    if [[ $ISSUES_FOUND -eq 0 ]]; then
        print_status "OK" "All required components are available!"
        if [[ $FIXES_APPLIED -eq 1 ]]; then
            echo ""
            print_status "INFO" "Some issues were automatically fixed"
        fi
        exit 0
    else
        echo ""
        print_status "ERROR" "Some components are missing or misconfigured"
        if [[ $FIX_MODE -eq 1 ]]; then
            echo ""
            print_status "INFO" "Run again with --fix to attempt automatic fixes"
        else
            echo ""
            echo "To fix issues automatically:"
            echo "  bash scripts/package/verify.sh --fix"
            echo ""
            echo "Or manually:"
            echo "  bash scripts/package/setup.sh --with-all"
        fi
        exit 1
    fi
}

# Run main function
main "$@"


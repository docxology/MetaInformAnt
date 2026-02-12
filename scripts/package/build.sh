#!/bin/bash
# METAINFORMANT Package Builder
# Builds distributable Python packages using uv

set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"

# Default values
CLEAN_BUILD=false
WHEEL_ONLY=false
SDIST_ONLY=false
CHECK_PACKAGE=false
OUTPUT_DIR="dist"
BUILD_ALL=true

usage() {
    cat << EOF
METAINFORMANT Package Builder - Build distributable Python packages

Usage: $0 [OPTIONS]

Options:
    -h, --help          Show this help message
    --clean             Clean previous build artifacts before building
    --wheel-only        Build only wheel package (.whl)
    --sdist-only        Build only source distribution (.tar.gz)
    --check             Run package validation checks (requires twine)
    --output-dir DIR    Specify output directory (default: dist/)
    --verbose           Enable verbose output

Examples:
    $0                          # Build both wheel and sdist
    $0 --clean                  # Clean and build all
    $0 --wheel-only             # Build only wheel
    $0 --check                  # Build and validate packages

Output:
    Creates distribution packages in dist/ directory
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            exit 0
            ;;
        --clean)
            CLEAN_BUILD=true
            shift
            ;;
        --wheel-only)
            WHEEL_ONLY=true
            BUILD_ALL=false
            shift
            ;;
        --sdist-only)
            SDIST_ONLY=true
            BUILD_ALL=false
            shift
            ;;
        --check)
            CHECK_PACKAGE=true
            shift
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --verbose)
            set -x
            shift
            ;;
        *)
            print_status "ERROR" "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

print_status "INFO" "Starting METAINFORMANT package build"

# Setup environment
setup_environment

# Clean previous builds if requested
if [[ "$CLEAN_BUILD" == "true" ]]; then
    print_status "INFO" "Cleaning previous build artifacts"
    rm -rf dist/ build/ *.egg-info/
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Build packages
if [[ "$BUILD_ALL" == "true" ]]; then
    print_status "INFO" "Building wheel and source distribution"
    uv build --out-dir "$OUTPUT_DIR"
elif [[ "$WHEEL_ONLY" == "true" ]]; then
    print_status "INFO" "Building wheel package only"
    uv build --wheel --out-dir "$OUTPUT_DIR"
elif [[ "$SDIST_ONLY" == "true" ]]; then
    print_status "INFO" "Building source distribution only"
    uv build --sdist --out-dir "$OUTPUT_DIR"
fi

# List built packages
if [[ -d "$OUTPUT_DIR" ]]; then
    print_status "INFO" "Build artifacts created:"
    ls -la "$OUTPUT_DIR"/*.whl "$OUTPUT_DIR"/*.tar.gz 2>/dev/null || true
fi

# Validate packages if requested
if [[ "$CHECK_PACKAGE" == "true" ]]; then
    print_status "INFO" "Validating built packages"

    # Check if twine is available
    if command -v twine &> /dev/null; then
        print_status "INFO" "Running twine checks"
        twine check "$OUTPUT_DIR"/*.whl "$OUTPUT_DIR"/*.tar.gz
        print_status "OK" "Package validation completed"
    else
        print_status "WARNING" "twine not available, skipping package validation"
        print_status "INFO" "Install twine: uv add --dev twine"
    fi
fi

print_status "SUCCESS" "Package build completed successfully"

# Print next steps
cat << EOF

Build completed! Next steps:
1. Test installation: uv pip install $OUTPUT_DIR/metainformant-*.whl
2. Test functionality: python -c "import metainformant; print('Import successful')"
3. Upload to PyPI: twine upload $OUTPUT_DIR/*

To verify:
- Create venv: uv venv .test_venv && source .test_venv/bin/activate
- Install: uv pip install $OUTPUT_DIR/metainformant-*.whl
- Test: python -c 'import metainformant; print(metainformant.__version__)'
- Clean: rm -rf .test_venv
EOF

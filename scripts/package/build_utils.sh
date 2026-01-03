#!/bin/bash
# METAINFORMANT Build Utilities
# Shared utility functions for build scripts

set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"

# Build utility functions

# Clean build artifacts
clean_build_artifacts() {
    print_status "INFO" "Cleaning build artifacts..."

    # Remove standard build directories
    rm -rf dist/ build/ *.egg-info/

    # Remove documentation build artifacts
    rm -rf docs/_build/

    # Remove test artifacts
    rm -rf .pytest_cache/ .coverage output/coverage/

    # Remove other build artifacts
    rm -rf .uv-cache/archive-v0/

    print_status "SUCCESS" "Build artifacts cleaned"
}

# Validate package structure
validate_package() {
    local dist_dir="${1:-dist}"

    if [[ ! -d "$dist_dir" ]]; then
        print_status "ERROR" "Distribution directory '$dist_dir' not found"
        return 1
    fi

    print_status "INFO" "Validating package structure..."

    # Check for required files
    local has_wheel=false
    local has_sdist=false

    if compgen -G "$dist_dir"/*.whl >/dev/null 2>&1; then
        has_wheel=true
        print_status "INFO" "Found wheel packages: $(ls "$dist_dir"/*.whl)"
    fi

    if compgen -G "$dist_dir"/*.tar.gz >/dev/null 2>&1; then
        has_sdist=true
        print_status "INFO" "Found source distributions: $(ls "$dist_dir"/*.tar.gz)"
    fi

    if [[ "$has_wheel" == "false" ]] && [[ "$has_sdist" == "false" ]]; then
        print_status "ERROR" "No wheel or source distribution found in $dist_dir"
        return 1
    fi

    # Validate with twine if available
    if command -v twine &> /dev/null; then
        print_status "INFO" "Running twine check..."
        if twine check "$dist_dir"/*.whl "$dist_dir"/*.tar.gz; then
            print_status "SUCCESS" "Package validation passed"
        else
            print_status "ERROR" "Package validation failed"
            return 1
        fi
    else
        print_status "WARNING" "twine not available, skipping validation"
    fi

    return 0
}

# Extract version from pyproject.toml
get_version() {
    local pyproject_file="${1:-pyproject.toml}"

    if [[ ! -f "$pyproject_file" ]]; then
        print_status "ERROR" "pyproject.toml not found"
        return 1
    fi

    # Extract version using grep and sed
    local version
    version=$(grep '^version = ' "$pyproject_file" | sed 's/version = "\(.*\)"/\1/')

    if [[ -z "$version" ]]; then
        print_status "ERROR" "Could not extract version from pyproject.toml"
        return 1
    fi

    echo "$version"
}

# Check build requirements
check_build_requirements() {
    print_status "INFO" "Checking build requirements..."

    local missing_deps=()

    # Check for uv
    if ! command -v uv &> /dev/null; then
        missing_deps+=("uv")
    fi

    # Check for build dependencies in pyproject.toml
    if [[ ! -f "pyproject.toml" ]]; then
        print_status "ERROR" "pyproject.toml not found"
        return 1
    fi

    # Check Python version
    local python_version
    python_version=$(python --version 2>&1 | awk '{print $2}' | cut -d. -f1-2)

    if [[ -z "$python_version" ]]; then
        print_status "ERROR" "Could not determine Python version"
        return 1
    fi

    # Check minimum Python version from pyproject.toml
    local min_python
    min_python=$(grep 'requires-python' pyproject.toml | sed 's/.*">=\([0-9.]*\)".*/\1/')

    if [[ -n "$min_python" ]]; then
        if ! python -c "import sys; sys.exit(0 if tuple(map(int, '$python_version'.split('.'))) >= tuple(map(int, '$min_python'.split('.'))) else 1)"; then
            print_status "ERROR" "Python $python_version is below minimum required version $min_python"
            return 1
        fi
    fi

    if [[ ${#missing_deps[@]} -gt 0 ]]; then
        print_status "ERROR" "Missing build dependencies: ${missing_deps[*]}"
        return 1
    fi

    print_status "SUCCESS" "Build requirements satisfied"
    return 0
}

# Test package installation
test_package_installation() {
    local package_path="$1"

    if [[ ! -f "$package_path" ]]; then
        print_status "ERROR" "Package file not found: $package_path"
        return 1
    fi

    print_status "INFO" "Testing package installation..."

    # Create temporary virtual environment
    local temp_venv
    temp_venv=$(mktemp -d)

    # Setup cleanup trap
    trap "rm -rf '$temp_venv'" EXIT

    # Create and activate venv using uv
    if ! uv venv "$temp_venv" >/dev/null 2>&1; then
        print_status "ERROR" "Failed to create virtual environment"
        return 1
    fi
    source "$temp_venv/bin/activate"

    # Install package using uv (UV-exclusive policy)
    if uv pip install "$package_path"; then
        print_status "SUCCESS" "Package installed successfully"
    else
        print_status "ERROR" "Package installation failed"
        deactivate
        return 1
    fi

    # Test basic import
    if python -c "import metainformant; print(f'Import successful: {metainformant.__version__}')"; then
        print_status "SUCCESS" "Package import test passed"
    else
        print_status "ERROR" "Package import test failed"
        deactivate
        return 1
    fi

    # Test CLI if available
    if command -v metainformant &> /dev/null; then
        if metainformant --help >/dev/null 2>&1; then
            print_status "SUCCESS" "CLI test passed"
        else
            print_status "WARNING" "CLI test failed (but package import worked)"
        fi
    fi

    # Deactivate and cleanup
    deactivate

    print_status "SUCCESS" "Package installation test completed"
    return 0
}

# Get package metadata
get_package_metadata() {
    local package_file="$1"

    if [[ ! -f "$package_file" ]]; then
        print_status "ERROR" "Package file not found: $package_file"
        return 1
    fi

    # Extract metadata using Python
    python -c "
import sys
sys.path.insert(0, '.')
from metainformant.core.io import dump_json
import json

# Try to extract metadata from wheel or sdist
import zipfile
import tarfile
import tempfile
import os

metadata = {}

try:
    if '$package_file'.endswith('.whl'):
        with zipfile.ZipFile('$package_file', 'r') as zf:
            # Look for metadata files
            for name in zf.namelist():
                if name.endswith('METADATA'):
                    with zf.open(name) as f:
                        content = f.read().decode('utf-8')
                        # Parse basic metadata
                        for line in content.split('\n'):
                            if line.startswith('Name:'):
                                metadata['name'] = line.split(':', 1)[1].strip()
                            elif line.startswith('Version:'):
                                metadata['version'] = line.split(':', 1)[1].strip()
                            elif line.startswith('Summary:'):
                                metadata['summary'] = line.split(':', 1)[1].strip()
                        break
    elif '$package_file'.endswith('.tar.gz'):
        with tarfile.open('$package_file', 'r:gz') as tf:
            # Look for PKG-INFO in tar
            for member in tf.getmembers():
                if member.name.endswith('PKG-INFO'):
                    with tf.extractfile(member) as f:
                        content = f.read().decode('utf-8')
                        for line in content.split('\n'):
                            if line.startswith('Name:'):
                                metadata['name'] = line.split(':', 1)[1].strip()
                            elif line.startswith('Version:'):
                                metadata['version'] = line.split(':', 1)[1].strip()
                            elif line.startswith('Summary:'):
                                metadata['summary'] = line.split(':', 1)[1].strip()
                        break

    print(json.dumps(metadata, indent=2))

except Exception as e:
    print(f'Error extracting metadata: {e}')
    sys.exit(1)
"
}

# Validate build environment
validate_build_environment() {
    print_status "INFO" "Validating build environment..."

    # Check if we're in the right directory
    if [[ ! -f "pyproject.toml" ]]; then
        print_status "ERROR" "pyproject.toml not found. Are you in the project root?"
        return 1
    fi

    if [[ ! -d "src/metainformant" ]]; then
        print_status "ERROR" "Source directory not found. Expected src/metainformant/"
        return 1
    fi

    # Check build requirements
    if ! check_build_requirements; then
        return 1
    fi

    print_status "SUCCESS" "Build environment validation passed"
    return 0
}

# Print build summary
print_build_summary() {
    local dist_dir="${1:-dist}"

    echo
    echo "========================================"
    echo "        BUILD SUMMARY"
    echo "========================================"

    if [[ -d "$dist_dir" ]]; then
        echo "Built packages in $dist_dir/:"
        ls -la "$dist_dir"/*.whl "$dist_dir"/*.tar.gz 2>/dev/null || echo "No packages found"
    else
        echo "No distribution directory found"
    fi

    echo
    echo "Next steps:"
    echo "1. Test installation: pip install $dist_dir/metainformant-*.whl"
    echo "2. Test functionality: python -c \"import metainformant; print('OK')\""
    echo "3. Upload to PyPI: twine upload $dist_dir/*"
    echo "4. Build docs: bash scripts/package/uv_docs.sh build"
    echo
}

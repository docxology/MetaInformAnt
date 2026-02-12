#!/bin/bash
# METAINFORMANT Release Preparation Script
# Prepare releases with version bumping, changelog generation, and validation

set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
source "$SCRIPT_DIR/build_utils.sh"

# Default values
RELEASE_TYPE="patch"
DRY_RUN=false
SKIP_TESTS=false
SKIP_BUILD=false
SKIP_DOCS=false
CREATE_TAG=false
PUSH_CHANGES=false

usage() {
    cat << EOF
METAINFORMANT Release Preparation Script

Usage: $0 [OPTIONS] [VERSION]

Arguments:
    VERSION     Specific version to release (e.g., 1.2.3), or leave empty to auto-bump

Options:
    -t, --type TYPE      Release type: major, minor, patch (default: patch)
    --dry-run           Show what would be done without making changes
    --skip-tests        Skip running tests
    --skip-build        Skip building packages
    --skip-docs         Skip building documentation
    --create-tag        Create git tag after successful preparation
    --push              Push changes and tag to remote (implies --create-tag)
    -h, --help          Show this help message

Examples:
    $0                          # Bump patch version and prepare release
    $0 --type minor            # Bump minor version
    $0 1.2.3                   # Release specific version
    $0 --dry-run               # Show what would be done
    $0 --push                  # Prepare and push release

Release Process:
1. Validate current state (clean git, tests pass)
2. Bump version in pyproject.toml
3. Update changelog
4. Run tests and build validation
5. Build packages and documentation
6. Create git tag (optional)
7. Push changes (optional)
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--type)
            RELEASE_TYPE="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --skip-tests)
            SKIP_TESTS=true
            shift
            ;;
        --skip-build)
            SKIP_BUILD=true
            shift
            ;;
        --skip-docs)
            SKIP_DOCS=true
            shift
            ;;
        --create-tag)
            CREATE_TAG=true
            shift
            ;;
        --push)
            PUSH_CHANGES=true
            CREATE_TAG=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            print_status "ERROR" "Unknown option: $1"
            usage
            exit 1
            ;;
        *)
            # Treat remaining argument as version
            if [[ -z "${SPECIFIC_VERSION:-}" ]]; then
                SPECIFIC_VERSION="$1"
            else
                print_status "ERROR" "Multiple version arguments provided"
                exit 1
            fi
            shift
            ;;
    esac
done

# Setup dry run mode
if [[ "$DRY_RUN" == "true" ]]; then
    print_status "INFO" "DRY RUN MODE - No actual changes will be made"
    export DRY_RUN=true
fi

print_status "INFO" "Starting METAINFORMANT release preparation"

# Validate environment
print_status "INFO" "Validating release environment..."
if ! validate_build_environment; then
    exit 1
fi

# Check git status
if [[ -n "$(git status --porcelain)" ]]; then
    print_status "ERROR" "Working directory is not clean. Please commit or stash changes."
    git status --short
    exit 1
fi

# Get current version
CURRENT_VERSION=$(get_version)
print_status "INFO" "Current version: $CURRENT_VERSION"

# Determine new version
if [[ -n "${SPECIFIC_VERSION:-}" ]]; then
    NEW_VERSION="$SPECIFIC_VERSION"
    print_status "INFO" "Using specific version: $NEW_VERSION"
else
    # Parse current version and bump
    IFS='.' read -ra VERSION_PARTS <<< "$CURRENT_VERSION"
    case "$RELEASE_TYPE" in
        "major")
            NEW_VERSION="$((VERSION_PARTS[0] + 1)).0.0"
            ;;
        "minor")
            NEW_VERSION="${VERSION_PARTS[0]}.$((VERSION_PARTS[1] + 1)).0"
            ;;
        "patch")
            NEW_VERSION="${VERSION_PARTS[0]}.${VERSION_PARTS[1]}.$((VERSION_PARTS[2] + 1))"
            ;;
        *)
            print_status "ERROR" "Invalid release type: $RELEASE_TYPE"
            exit 1
            ;;
    esac
    print_status "INFO" "Bumping $RELEASE_TYPE version: $CURRENT_VERSION → $NEW_VERSION"
fi

# Update version in pyproject.toml
if [[ "$DRY_RUN" != "true" ]]; then
    print_status "INFO" "Updating version in pyproject.toml..."
    sed -i.bak "s/version = \"$CURRENT_VERSION\"/version = \"$NEW_VERSION\"/" pyproject.toml
    rm pyproject.toml.bak
else
    print_status "INFO" "[DRY RUN] Would update version in pyproject.toml: $CURRENT_VERSION → $NEW_VERSION"
fi

# Update changelog
if [[ -f "CHANGELOG.md" ]]; then
    print_status "INFO" "Updating changelog..."
    if [[ "$DRY_RUN" != "true" ]]; then
        # Add new version entry to changelog
        TODAY=$(date +%Y-%m-%d)
        sed -i.bak "1a ## [$NEW_VERSION] - $TODAY\n\n### Added\n- Release preparation improvements\n\n### Changed\n- Version bumped to $NEW_VERSION\n\n" CHANGELOG.md
        rm CHANGELOG.md.bak
    else
        print_status "INFO" "[DRY RUN] Would update CHANGELOG.md with version $NEW_VERSION"
    fi
else
    print_status "WARNING" "CHANGELOG.md not found, skipping changelog update"
fi

# Run tests
if [[ "$SKIP_TESTS" != "true" ]]; then
    print_status "INFO" "Running tests..."
    if ! bash scripts/package/test.sh --mode fast; then
        print_status "ERROR" "Tests failed, aborting release"
        exit 1
    fi
    print_status "SUCCESS" "Tests passed"
else
    print_status "INFO" "Skipping tests"
fi

# Build packages
if [[ "$SKIP_BUILD" != "true" ]]; then
    print_status "INFO" "Building packages..."
    if ! bash scripts/package/build.sh --check; then
        print_status "ERROR" "Package build failed, aborting release"
        exit 1
    fi
    print_status "SUCCESS" "Packages built successfully"
else
    print_status "INFO" "Skipping package build"
fi

# Build documentation
if [[ "$SKIP_DOCS" != "true" ]]; then
    print_status "INFO" "Building documentation..."
    if ! bash scripts/package/uv_docs.sh build; then
        print_status "WARNING" "Documentation build failed, but continuing with release"
    else
        print_status "SUCCESS" "Documentation built successfully"
    fi
else
    print_status "INFO" "Skipping documentation build"
fi

# Validate release artifacts
print_status "INFO" "Validating release artifacts..."
if [[ -d "dist/" ]]; then
    if ! validate_package "dist/"; then
        print_status "ERROR" "Package validation failed"
        exit 1
    fi
else
    print_status "WARNING" "No dist/ directory found"
fi

# Test package installation
if [[ -d "dist/" ]] && compgen -G "dist/*.whl" >/dev/null 2>&1; then
    print_status "INFO" "Testing package installation..."
    WHEEL_FILE=$(ls dist/*.whl | head -1)
    if ! test_package_installation "$WHEEL_FILE"; then
        print_status "ERROR" "Package installation test failed"
        exit 1
    fi
    print_status "SUCCESS" "Package installation test passed"
fi

# Commit changes
if [[ "$DRY_RUN" != "true" ]]; then
    print_status "INFO" "Committing release changes..."
    git add pyproject.toml CHANGELOG.md
    git commit -m "Release $NEW_VERSION

- Bump version to $NEW_VERSION
- Update changelog
- Prepare for release"
else
    print_status "INFO" "[DRY RUN] Would commit release changes"
fi

# Create git tag
if [[ "$CREATE_TAG" == "true" ]]; then
    TAG_NAME="v$NEW_VERSION"
    if [[ "$DRY_RUN" != "true" ]]; then
        print_status "INFO" "Creating git tag: $TAG_NAME"
        git tag -a "$TAG_NAME" -m "Release $NEW_VERSION"
    else
        print_status "INFO" "[DRY RUN] Would create git tag: $TAG_NAME"
    fi
fi

# Push changes
if [[ "$PUSH_CHANGES" == "true" ]]; then
    if [[ "$DRY_RUN" != "true" ]]; then
        print_status "INFO" "Pushing changes and tags..."
        git push origin main
        git push origin "$TAG_NAME"
    else
        print_status "INFO" "[DRY RUN] Would push changes and tag"
    fi
fi

# Print release summary
echo
echo "========================================"
echo "      RELEASE PREPARATION COMPLETE"
echo "========================================"
echo
echo "Version: $NEW_VERSION"
echo "Previous: $CURRENT_VERSION"
echo "Release Type: $RELEASE_TYPE"
echo
if [[ -d "dist/" ]]; then
    echo "Built packages:"
    ls -la dist/ 2>/dev/null || echo "No packages found"
    echo
fi

echo "Next steps:"
echo "1. Review changes: git log --oneline -5"
echo "2. Test manually: uv pip install dist/metainformant-*.whl"
echo "3. Create GitHub release with tag v$NEW_VERSION"
if [[ "$PUSH_CHANGES" != "true" ]]; then
    echo "4. Push changes: git push origin main && git push origin v$NEW_VERSION"
fi
echo "5. Monitor CI/CD pipelines"
echo
print_status "SUCCESS" "Release preparation completed for version $NEW_VERSION"

#!/bin/bash
# METAINFORMANT Documentation Builder
# Enhanced Sphinx documentation build with multiple formats and validation

set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"

# Default values
BUILD_FORMAT="html"
SERVE_PORT=8000
CHECK_LINKS=false
VALIDATE=false
WATCH=false
OUTPUT_DIR="docs/_build"
SOURCE_DIR="docs"

usage() {
    cat << EOF
METAINFORMANT Documentation Builder

Usage: $0 [OPTIONS] [COMMAND]

Commands:
    build       Build documentation (default)
    serve       Serve built documentation
    clean       Clean build artifacts
    watch       Auto-rebuild on file changes
    check       Validate documentation and links

Options:
    -f, --format FORMAT    Output format: html, pdf, epub, latex (default: html)
    -p, --port PORT        Port for serve command (default: 8000)
    -o, --output DIR       Output directory (default: docs/_build)
    --links                Check external links during build
    --validate             Run additional validation checks
    --help                 Show this help message

Examples:
    $0                      # Build HTML documentation
    $0 build --format pdf  # Build PDF documentation
    $0 serve               # Serve documentation at http://localhost:8000
    $0 watch               # Auto-rebuild on changes
    $0 check --links       # Build with link checking
    $0 clean               # Clean build artifacts

Output formats:
    html    - Web documentation (default)
    pdf     - PDF document (requires LaTeX)
    epub    - EPUB e-book format
    latex   - LaTeX source files
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--format)
            BUILD_FORMAT="$2"
            shift 2
            ;;
        -p|--port)
            SERVE_PORT="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --links)
            CHECK_LINKS=true
            shift
            ;;
        --validate)
            VALIDATE=true
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        -*)
            print_status "ERROR" "Unknown option: $1"
            usage
            exit 1
            ;;
        *)
            COMMAND="$1"
            shift
            break
            ;;
    esac
done

COMMAND="${COMMAND:-build}"

# Setup environment
setup_environment

print_status "INFO" "METAINFORMANT Documentation Builder"

case "$COMMAND" in
    "build")
        print_status "INFO" "Building $BUILD_FORMAT documentation..."

        # Build command with options
        BUILD_CMD="uv run sphinx-build -b $BUILD_FORMAT"

        # Add link checking if requested
        if [[ "$CHECK_LINKS" == "true" ]]; then
            BUILD_CMD="$BUILD_CMD -W --keep-going -E"
            print_status "INFO" "Link checking enabled"
        fi

        # Add validation if requested
        if [[ "$VALIDATE" == "true" ]]; then
            print_status "INFO" "Running additional validation"
            # Check for required tools
            if ! command -v uv &> /dev/null; then
                print_status "ERROR" "uv not found. Install uv first."
                exit 1
            fi
        fi

        # Execute build
        $BUILD_CMD "$SOURCE_DIR" "$OUTPUT_DIR/$BUILD_FORMAT"

        # Report results
        if [[ "$BUILD_FORMAT" == "html" ]]; then
            DOC_INDEX="$OUTPUT_DIR/html/index.html"
            if [[ -f "$DOC_INDEX" ]]; then
                print_status "SUCCESS" "Documentation built successfully"
                print_status "INFO" "Open in browser: file://$(pwd)/$DOC_INDEX"
            else
                print_status "ERROR" "Documentation build failed"
                exit 1
            fi
        else
            print_status "SUCCESS" "$BUILD_FORMAT documentation built in $OUTPUT_DIR/$BUILD_FORMAT/"
        fi
        ;;

    "serve")
        DOC_DIR="$OUTPUT_DIR/html"
        if [[ ! -d "$DOC_DIR" ]]; then
            print_status "ERROR" "Documentation not built yet. Run: $0 build"
            exit 1
        fi

        print_status "INFO" "Serving documentation at http://localhost:$SERVE_PORT"
        print_status "INFO" "Press Ctrl+C to stop"

        # Use Python's built-in server
        cd "$DOC_DIR"
        uv run python -m http.server "$SERVE_PORT"
        ;;

    "watch")
        print_status "INFO" "Watching for changes... (Ctrl+C to stop)"

        # Check if fswatch is available
        if command -v fswatch &> /dev/null; then
            fswatch -o "$SOURCE_DIR" | while read -r num; do
                print_status "INFO" "Changes detected, rebuilding..."
                $0 build --format "$BUILD_FORMAT"
                print_status "INFO" "Rebuild complete"
            done
        else
            print_status "WARNING" "fswatch not found. Install fswatch for auto-rebuilding."
            print_status "INFO" "Manual rebuild: $0 build"
            exit 1
        fi
        ;;

    "check")
        print_status "INFO" "Running documentation validation..."

        # Check if documentation builds without errors
        if $0 build --format html >/dev/null 2>&1; then
            print_status "SUCCESS" "Documentation builds without errors"
        else
            print_status "ERROR" "Documentation build failed"
            exit 1
        fi

        # Check external links if requested
        if [[ "$CHECK_LINKS" == "true" ]]; then
            print_status "INFO" "Checking external links..."
            uv run sphinx-build -b linkcheck "$SOURCE_DIR" "$OUTPUT_DIR/linkcheck"
            print_status "SUCCESS" "Link check completed"
        fi

        # Additional validation
        if [[ "$VALIDATE" == "true" ]]; then
            print_status "INFO" "Running additional validation checks..."

            # Check for broken internal references
            if uv run sphinx-build -b html -W "$SOURCE_DIR" "$OUTPUT_DIR/validate" >/dev/null 2>&1; then
                print_status "SUCCESS" "No broken internal references"
            else
                print_status "WARNING" "Found broken internal references"
            fi
        fi
        ;;

    "clean")
        print_status "INFO" "Cleaning documentation build artifacts..."
        rm -rf "$OUTPUT_DIR"
        rm -rf docs/_build/
        print_status "SUCCESS" "Documentation build directory cleaned"
        ;;

    *)
        print_status "ERROR" "Unknown command: $COMMAND"
        usage
        exit 1
        ;;
esac

#!/bin/bash
# METAINFORMANT Orchestrator Script
# Provides interactive menu access and direct script execution

set -euo pipefail

# Get script directory and repo root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"

# Source common utilities if available
if [ -f "$REPO_ROOT/scripts/package/_common.sh" ]; then
    source "$REPO_ROOT/scripts/package/_common.sh"
fi

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print usage
print_usage() {
    cat <<EOF
METAINFORMANT Orchestrator

Usage:
    $0 [OPTIONS] [CATEGORY] [SCRIPT] [ARGS...]

Options:
    --help, -h          Show this help message
    --menu, -m         Launch interactive menu
    --list, -l         List available scripts
    --category, -c     Show scripts in a category

Examples:
    $0                    # Launch interactive menu (default)
    $0 --menu             # Launch interactive menu explicitly
    $0 --help             # Show this help message
    $0 --list             # List all available scripts
    $0 rna                # Show RNA scripts menu
    $0 rna run_workflow.py --config config.yaml  # Execute script directly

Categories:
    core, rna, gwas, dna, protein, networks, multiomics, math, ml,
    singlecell, quality, simulation, visualization, epigenome, ecology,
    ontology, phenotype, life_events, package

EOF
}

# Function to check if uv is available (required)
check_python_env() {
    if ! command -v uv >/dev/null 2>&1; then
        echo "❌ Error: uv package manager is required"
        echo "   METAINFORMANT uses uv for all Python package management and execution"
        echo "   Install uv from: https://github.com/astral.sh/uv"
        echo "   curl -LsSf https://astral.sh/uv/install.sh | sh"
        return 1
    fi

    # Verify uv can run Python
    if ! uv run python --version >/dev/null 2>&1; then
        echo "❌ Error: uv cannot run Python"
        echo "   Try: uv python install 3.11"
        return 1
    fi

    return 0
}

# Function to launch interactive menu
launch_menu() {
    check_python_env || return 1

    # Launch the dedicated Python menu script
    local menu_script="$REPO_ROOT/scripts/menu/run_menu.py"

    if [ ! -f "$menu_script" ]; then
        echo "❌ Error: Menu launcher script not found: $menu_script"
        return 1
    fi

    # Execute menu script with uv
    uv run python "$menu_script"
}

# Function to list available scripts
list_scripts() {
    check_python_env || return 1

    uv run python -c "
from pathlib import Path
import sys
sys.path.insert(0, str(Path('$REPO_ROOT/src')))
from metainformant.menu import discover_scripts

scripts = discover_scripts(Path('$REPO_ROOT'))
for category, script_list in sorted(scripts.items()):
    print(f'\n{category.upper()}:')
    for script in script_list:
        print(f'  {script.name:<40} {script.description}')
"
}

# Function to show category menu
show_category() {
    local category="$1"
    check_python_env || return 1

    uv run python -c "
from pathlib import Path
import sys
sys.path.insert(0, str(Path('$REPO_ROOT/src')))
from metainformant.menu import discover_scripts, generate_menu_from_scripts
from metainformant.menu.display import show_menu

scripts = discover_scripts(Path('$REPO_ROOT'))
menus = generate_menu_from_scripts(scripts)

menu_id = f'menu_{category}'
if menu_id in menus:
    menu = menus[menu_id]
    show_menu(menu.items, menu.title)
else:
    print(f'❌ Category \"{category}\" not found')
    sys.exit(1)
"
}

# Function to execute script directly
execute_direct() {
    local category="$1"
    local script_name="$2"
    shift 2
    local script_args=("$@")
    
    # Find script
    local script_path="$REPO_ROOT/scripts/$category/$script_name"
    
    if [ ! -f "$script_path" ]; then
        # Try with .py extension
        script_path="$REPO_ROOT/scripts/$category/${script_name}.py"
    fi
    
    if [ ! -f "$script_path" ]; then
        # Try with .sh extension
        script_path="$REPO_ROOT/scripts/$category/${script_name}.sh"
    fi
    
    if [ ! -f "$script_path" ]; then
        echo "Error: Script not found: $category/$script_name"
        echo "Available scripts in $category:"
        ls -1 "$REPO_ROOT/scripts/$category/" 2>/dev/null | grep -E '\.(py|sh)$' || echo "  (none)"
        return 1
    fi
    
    # Execute script
    if [[ "$script_path" == *.py ]]; then
        uv run python "$script_path" "${script_args[@]}"
    elif [[ "$script_path" == *.sh ]]; then
        bash "$script_path" "${script_args[@]}"
    else
        echo "❌ Error: Unsupported script type: $script_path"
        return 1
    fi
}

# Main argument parsing
main() {
    # Handle help first (explicit help overrides default behavior)
    if [ $# -gt 0 ] && ([ "$1" = "--help" ] || [ "$1" = "-h" ]); then
        print_usage
        exit 0
    fi
    
    # Handle menu launch
    if [ $# -gt 0 ] && ([ "$1" = "--menu" ] || [ "$1" = "-m" ]); then
        launch_menu
        exit $?
    fi

    # Handle list
    if [ $# -gt 0 ] && ([ "$1" = "--list" ] || [ "$1" = "-l" ]); then
        list_scripts
        exit $?
    fi

    # Handle category display
    if [ $# -gt 0 ] && ([ "$1" = "--category" ] || [ "$1" = "-c" ]); then
        if [ $# -lt 2 ]; then
            echo "Error: Category name required"
            print_usage
            exit 1
        fi
        show_category "$2"
        exit $?
    fi
    
    # Direct execution: category script [args...]
    if [ $# -ge 2 ]; then
        execute_direct "$@"
        exit $?
    fi
    
    # Default: launch menu
    launch_menu
    exit $?
}

# Run main function
main "$@"


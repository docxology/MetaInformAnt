#!/usr/bin/env python3
"""Interactive menu launcher for METAINFORMANT.

This script provides an interactive text-based menu system for accessing
scripts and workflows throughout the METAINFORMANT repository.

Usage:
    python3 scripts/menu/run_menu.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

# Add src to path for imports
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from metainformant.core.utils.logging import get_logger
from metainformant.menu import discover_scripts, generate_menu_from_scripts
from metainformant.menu.display import clear_screen, format_breadcrumb, get_choice, show_menu
from metainformant.menu.executor import execute_script, prompt_for_args
from metainformant.menu.navigation import MenuSystem

logger = get_logger(__name__)

# Optional dependency groups and their key packages
OPTIONAL_DEPENDENCY_GROUPS = {
    "database": ["psycopg2"],
    "networks": ["community", "cdlib"],  # python-louvain installs as 'community'
    "scientific": ["scipy", "sklearn"],
    "ml": ["optuna", "joblib"],
    "singlecell": ["scanpy", "leidenalg"],
    "visualization": ["plotly", "bokeh", "altair"],
    "external-tools": ["pysam", "dendropy"],
    "scraping": ["cloudscraper"],
}


def check_uv_available() -> bool:
    """Check if uv is available in the system PATH.

    Returns:
        True if uv is available, False otherwise
    """
    try:
        result = subprocess.run(["uv", "--version"], capture_output=True, text=True, timeout=5)
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


def check_optional_dependencies() -> dict[str, bool]:
    """Check which optional dependency groups are available.

    Returns:
        Dictionary mapping group name to availability status
    """
    availability = {}

    for group_name, packages in OPTIONAL_DEPENDENCY_GROUPS.items():
        group_available = True
        for package in packages:
            try:
                # Use uv run to check if package is importable
                result = subprocess.run(
                    ["uv", "run", "python", "-c", f"import {package}"],
                    capture_output=True,
                    text=True,
                    timeout=10,
                    cwd=REPO_ROOT,
                )
                if result.returncode != 0:
                    group_available = False
                    break
            except (subprocess.TimeoutExpired, FileNotFoundError):
                group_available = False
                break

        availability[group_name] = group_available

    return availability


def install_optional_dependencies(missing_groups: list[str]) -> bool:
    """Install missing optional dependency groups using uv.

    Args:
        missing_groups: List of dependency group names to install

    Returns:
        True if installation was successful, False otherwise
    """
    if not missing_groups:
        return True

    try:
        # Build uv sync command with extras
        cmd = ["uv", "sync"]
        for group in missing_groups:
            cmd.extend(["--extra", group])

        print(f"📦 Installing dependencies: {' '.join(missing_groups)}")
        print("   This may take a few minutes...")

        result = subprocess.run(cmd, cwd=REPO_ROOT, timeout=300)  # 5 minute timeout

        return result.returncode == 0
    except subprocess.TimeoutExpired:
        print("❌ Dependency installation timed out")
        return False
    except Exception as e:
        print(f"❌ Error during dependency installation: {e}")
        return False


def ensure_dependencies() -> bool:
    """Ensure all dependencies are available, offering to install missing ones.

    Returns:
        True if all dependencies are available or successfully installed,
        False if user declined or installation failed
    """
    print("🔍 Checking dependencies...")

    # Check uv availability
    if not check_uv_available():
        print("❌ uv package manager is required but not found")
        print("   Install uv from: https://github.com/astral.sh/uv")
        print("   curl -LsSf https://astral.sh/uv/install.sh | sh")
        return False

    # Check optional dependencies
    dep_status = check_optional_dependencies()

    available_groups = [group for group, available in dep_status.items() if available]
    missing_groups = [group for group, available in dep_status.items() if not available]

    if available_groups:
        print(f"✅ Available dependency groups: {', '.join(available_groups)}")

    if missing_groups:
        print("⚠️  Missing optional dependency groups:")
        for group in missing_groups:
            packages = OPTIONAL_DEPENDENCY_GROUPS[group]
            print(f"   • {group} ({', '.join(packages)})")

        print()
        while True:
            try:
                response = input("💡 Install missing optional dependencies? [y/N]: ").strip().lower()
                if response in ("", "n", "no"):
                    print("ℹ️  Continuing without optional dependencies")
                    print("   Some features may be disabled (warnings will appear)")
                    return True
                elif response in ("y", "yes"):
                    if install_optional_dependencies(missing_groups):
                        print("✅ All dependencies installed successfully!")
                        print()
                        return True
                    else:
                        print("❌ Failed to install dependencies")
                        response = input("💡 Continue anyway? [y/N]: ").strip().lower()
                        return response in ("y", "yes")
                else:
                    print("Please enter 'y' or 'n'")
            except (KeyboardInterrupt, EOFError):
                print("\nℹ️  Continuing without optional dependencies")
                return True
    else:
        print("✅ All optional dependencies are available")

    return True


def run_interactive_menu(menu_system: MenuSystem) -> None:
    """Run the main interactive menu loop.

    Args:
        menu_system: Initialized menu system
    """
    try:
        while True:
            # Clear screen for clean display
            clear_screen()

            # Get current menu
            current_menu = menu_system.get_current_menu()

            # Display breadcrumb
            breadcrumb = format_breadcrumb(menu_system.history.get_path())
            if breadcrumb != "Home":
                print(f"📍 {breadcrumb}")
                print()

            # Display menu
            show_menu(current_menu.items, current_menu.title)

            # Show navigation hints
            if menu_system.current_menu_id == "root":
                print("\n💡 Select a category to browse available scripts")
            else:
                print(
                    f"\n💡 Select a script to run or '0' to go back to {menu_system.history.get_path()[-2] if len(menu_system.history.get_path()) > 1 else 'main menu'}"
                )

            # Get user choice
            try:
                choice = get_choice("Select option: ").strip()
            except (KeyboardInterrupt, EOFError):
                print("\n\n👋 Goodbye!")
                return

            # Handle back/exit
            if choice == "0":
                if menu_system.current_menu_id == "root":
                    print("👋 Goodbye!")
                    return
                menu_system.go_back()
                continue

            # Handle menu selection
            try:
                choice_num = int(choice)
                if 1 <= choice_num <= len(current_menu.items):
                    item = current_menu.items[choice_num - 1]

                    if not item.enabled:
                        print(f"❌ {item.label} is currently disabled.")
                        input("Press Enter to continue...")
                        continue

                    # Handle action
                    action_result = handle_menu_action(item, menu_system)
                    if action_result is None:
                        # Wait for user input before continuing
                        input("Press Enter to continue...")
                else:
                    print("❌ Invalid choice. Please enter a number from the menu.")
                    input("Press Enter to continue...")
            except ValueError:
                print("❌ Invalid choice. Please enter a number.")
                input("Press Enter to continue...")
            except Exception as e:
                print(f"❌ Error processing choice: {e}")
                input("Press Enter to continue...")

    except KeyboardInterrupt:
        print("\n\n👋 Goodbye!")
        return


def handle_menu_action(item, menu_system: MenuSystem) -> str | None:
    """Handle menu item action.

    Args:
        item: MenuItem to handle
        menu_system: Menu system instance

    Returns:
        None if user should be prompted to continue, otherwise a result string
    """
    if isinstance(item.action, str):
        if item.action.startswith("submenu:"):
            # Navigate to submenu
            submenu_id = item.action.split(":", 1)[1]
            menu_system.navigate_to(submenu_id)
            return None  # Continue menu loop

        elif item.action.startswith("script:"):
            # Execute script
            script_path_str = item.action.split(":", 1)[1]
            script_path = Path(script_path_str)

            if script_path.exists():
                clear_screen()
                print(f"🚀 Executing: {script_path.name}")
                print(f"📄 Description: {item.description}")
                print(f"📁 Location: {script_path.relative_to(REPO_ROOT)}")
                print("─" * 60)

                try:
                    exit_code = execute_script(script_path)
                    print("─" * 60)
                    if exit_code == 0:
                        print("✅ Script completed successfully")
                        print("💡 Press Enter to return to menu...")
                    else:
                        print(f"⚠️  Script exited with code: {exit_code}")
                        print("💡 Press Enter to return to menu...")
                except Exception as e:
                    print(f"❌ Error executing script: {e}")
                    print("💡 Press Enter to return to menu...")

                return None  # Prompt user to continue
            else:
                print(f"❌ Script not found: {script_path}")
                print("💡 Press Enter to return to menu...")
                return None

        else:
            print(f"❓ Unknown action: {item.action}")
            return None

    elif callable(item.action):
        # Execute callable action
        try:
            result = item.action()
            if isinstance(result, str):
                print(result)
            return None
        except Exception as e:
            print(f"❌ Error executing action: {e}")
            return None

    else:
        print(f"❓ No action defined for {item.label}")
        return None


def main() -> int:
    """Main entry point."""
    print("🔬 METAINFORMANT Interactive Script Menu")
    print("🧬 Access bioinformatics workflows and tools")
    print("=" * 60)

    try:
        # Check dependencies first
        if not ensure_dependencies():
            return 1

        print("🔍 Discovering available scripts...")
        scripts = discover_scripts(REPO_ROOT)
        menus = generate_menu_from_scripts(scripts)

        if not menus:
            print("❌ No scripts found in scripts/ directory.")
            print("   This may indicate an installation issue.")
            print("   Try running: bash scripts/package/setup.sh")
            return 1

        if "root" not in menus:
            print("❌ Menu structure could not be generated.")
            return 1

        # Create menu system
        menu_system = MenuSystem(menus=menus, current_menu_id="root")

        # Display welcome message
        total_scripts = sum(len(s) for s in scripts.values())
        print(f"✅ Found {len(scripts)} script categories with {total_scripts} total scripts")
        print()
        print("🎮 Navigation:")
        print("   • Use number keys (1-9) to select options")
        print("   • Press '0' to go back or exit")
        print("   • Press Ctrl+C to quit at any time")
        print()
        print("📂 Categories include: RNA-seq, GWAS, DNA analysis, ML, visualization, and more")
        print()
        input("🎯 Press Enter to explore available tools...")

        # Run interactive menu
        run_interactive_menu(menu_system)

        return 0

    except KeyboardInterrupt:
        print("\n\n👋 Goodbye!")
        return 0
    except Exception as e:
        print(f"❌ Error launching menu: {e}")
        logger.exception("Menu launch failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())

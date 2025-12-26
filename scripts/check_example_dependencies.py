#!/usr/bin/env python3
"""Check dependencies for METAINFORMANT examples.

This script validates that all required dependencies are available before running examples,
helping prevent runtime failures due to missing packages.

Usage:
    python scripts/check_example_dependencies.py [--domain DOMAIN] [--example PATH] [--fix]

Arguments:
    --domain: Check dependencies for specific domain only
    --example: Check dependencies for specific example file
    --fix: Attempt to install missing dependencies (requires pip)
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Set


class DependencyChecker:
    """Check dependencies for METAINFORMANT examples."""

    def __init__(self, dependencies_file: Path | None = None):
        self.dependencies_file = dependencies_file or Path("examples/dependencies.json")
        self.dependencies = self._load_dependencies()

    def _load_dependencies(self) -> Dict[str, any]:
        """Load dependency manifest."""
        if not self.dependencies_file.exists():
            raise FileNotFoundError(f"Dependencies file not found: {self.dependencies_file}")

        with open(self.dependencies_file, 'r') as f:
            return json.load(f)

    def check_all_dependencies(self, domain_filter: str | None = None,
                              example_filter: str | None = None) -> Dict[str, any]:
        """Check all dependencies."""
        results = {
            "total_examples": 0,
            "checked_examples": 0,
            "missing_required": [],
            "missing_optional": [],
            "external_tools": [],
            "example_status": {}
        }

        # Check global dependencies
        global_deps = self._check_dependencies(self.dependencies["global_dependencies"])
        results["global_deps"] = global_deps

        # Check domain dependencies
        if domain_filter:
            domains_to_check = [domain_filter] if domain_filter in self.dependencies["domain_dependencies"] else []
        else:
            domains_to_check = list(self.dependencies["domain_dependencies"].keys())

        for domain in domains_to_check:
            if domain not in self.dependencies["domain_dependencies"]:
                continue

            domain_deps = self._check_dependencies(self.dependencies["domain_dependencies"][domain])
            results[f"{domain}_deps"] = domain_deps

            # Check examples in this domain
            examples_dir = Path("examples") / domain
            if examples_dir.exists():
                for example_file in examples_dir.glob("example_*.py"):
                    example_path = f"{domain}/{example_file.name}"
                    results["total_examples"] += 1

                    if example_filter and example_path != example_filter:
                        continue

                    results["checked_examples"] += 1
                    example_status = self._check_example_dependencies(example_path)
                    results["example_status"][example_path] = example_status

                    # Aggregate issues
                    if example_status["missing_required"]:
                        results["missing_required"].extend(example_status["missing_required"])
                    if example_status["missing_optional"]:
                        results["missing_optional"].extend(example_status["missing_optional"])

        # Check external tools
        for tool_name, tool_info in self.dependencies.get("external_tools", {}).items():
            tool_status = self._check_external_tool(tool_name, tool_info)
            results["external_tools"].append(tool_status)

        # Remove duplicates
        results["missing_required"] = list(set(results["missing_required"]))
        results["missing_optional"] = list(set(results["missing_optional"]))

        return results

    def _check_dependencies(self, deps_config: Dict[str, List[str]]) -> Dict[str, any]:
        """Check a set of dependencies."""
        required = deps_config.get("required", [])
        optional = deps_config.get("optional", [])

        missing_required = []
        missing_optional = []

        for dep in required:
            if not self._is_package_available(dep):
                missing_required.append(dep)

        for dep in optional:
            if not self._is_package_available(dep):
                missing_optional.append(dep)

        return {
            "required": required,
            "optional": optional,
            "missing_required": missing_required,
            "missing_optional": missing_optional,
            "all_available": len(missing_required) == 0
        }

    def _check_example_dependencies(self, example_path: str) -> Dict[str, any]:
        """Check dependencies for a specific example."""
        if example_path in self.dependencies.get("example_dependencies", {}):
            deps_config = self.dependencies["example_dependencies"][example_path]
            return self._check_dependencies(deps_config)
        else:
            # Fallback to domain-level dependencies
            domain = example_path.split("/")[0]
            if domain in self.dependencies.get("domain_dependencies", {}):
                return self._check_dependencies(self.dependencies["domain_dependencies"][domain])
            else:
                return {
                    "required": [],
                    "optional": [],
                    "missing_required": [],
                    "missing_optional": [],
                    "all_available": True
                }

    def _is_package_available(self, package_name: str) -> bool:
        """Check if a Python package is available."""
        try:
            # Handle submodules (e.g., metainformant.core)
            if "." in package_name:
                module_parts = package_name.split(".")
                module = importlib.import_module(module_parts[0])
                for part in module_parts[1:]:
                    module = getattr(module, part)
                return True
            else:
                importlib.import_module(package_name)
                return True
        except ImportError:
            return False

    def _check_external_tool(self, tool_name: str, tool_info: Dict[str, any]) -> Dict[str, any]:
        """Check if an external tool is available."""
        check_command = tool_info.get("check_command", f"which {tool_name}")

        try:
            result = subprocess.run(
                check_command,
                shell=True,
                capture_output=True,
                text=True,
                timeout=10
            )
            available = result.returncode == 0
        except (subprocess.TimeoutExpired, subprocess.SubprocessError):
            available = False

        return {
            "tool": tool_name,
            "available": available,
            "description": tool_info.get("description", ""),
            "examples": tool_info.get("examples", []),
            "install_hint": tool_info.get("install_hint", "")
        }

    def print_report(self, results: Dict[str, any]) -> None:
        """Print a human-readable dependency report."""
        print("METAINFORMANT Examples Dependency Check")
        print("=" * 50)

        # Global dependencies
        global_deps = results.get("global_deps", {})
        if global_deps.get("missing_required"):
            print("❌ Missing global required dependencies:")
            for dep in global_deps["missing_required"]:
                print(f"  - {dep}")

        # Example status
        if results["checked_examples"] > 0:
            print(f"\nChecked {results['checked_examples']} examples:")

            all_good = True
            for example, status in results["example_status"].items():
                if not status["all_available"]:
                    all_good = False
                    missing = status["missing_required"] + status["missing_optional"]
                    print(f"❌ {example}: missing {', '.join(missing)}")

            if all_good:
                print("✅ All example dependencies available")

        # External tools
        if results.get("external_tools"):
            print("\nExternal tools:")
            for tool in results["external_tools"]:
                status = "✅" if tool["available"] else "❌"
                print(f"{status} {tool['tool']}: {tool['description']}")
                if not tool["available"] and tool.get("install_hint"):
                    print(f"   Hint: {tool['install_hint']}")

        # Summary
        missing_required = len(results.get("missing_required", []))
        missing_optional = len(results.get("missing_optional", []))

        print("\nSummary:")
        print(f"Examples checked: {results['checked_examples']}")
        if missing_required > 0:
            print(f"❌ Missing required dependencies: {missing_required}")
        if missing_optional > 0:
            print(f"⚠️  Missing optional dependencies: {missing_optional}")
        if missing_required == 0 and missing_optional == 0:
            print("✅ All dependencies satisfied")

    def attempt_fix(self, results: Dict[str, any]) -> bool:
        """Attempt to fix missing dependencies."""
        missing_deps = set(results.get("missing_required", []) + results.get("missing_optional", []))

        if not missing_deps:
            print("✅ No missing dependencies to fix")
            return True

        print(f"Attempting to install {len(missing_deps)} missing dependencies...")

        # Filter to pip-installable packages (exclude metainformant modules)
        pip_packages = [dep for dep in missing_deps if not dep.startswith("metainformant")]

        if not pip_packages:
            print("⚠️  No pip-installable packages to fix")
            return True

        try:
            # Use uv pip install for consistency
            cmd = [sys.executable, "-m", "pip", "install"] + pip_packages
            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode == 0:
                print("✅ Dependencies installed successfully")
                return True
            else:
                print(f"❌ Failed to install dependencies: {result.stderr}")
                return False

        except Exception as e:
            print(f"❌ Error installing dependencies: {e}")
            return False


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Check METAINFORMANT example dependencies")
    parser.add_argument("--domain", help="Check specific domain only")
    parser.add_argument("--example", help="Check specific example file")
    parser.add_argument("--fix", action="store_true", help="Attempt to install missing dependencies")
    parser.add_argument("--quiet", action="store_true", help="Suppress detailed output")

    args = parser.parse_args()

    # Create checker
    checker = DependencyChecker()

    try:
        # Check dependencies
        results = checker.check_all_dependencies(
            domain_filter=args.domain,
            example_filter=args.example
        )

        # Print report
        if not args.quiet:
            checker.print_report(results)

        # Attempt fix if requested
        if args.fix:
            if checker.attempt_fix(results):
                # Re-check after fix
                results = checker.check_all_dependencies(
                    domain_filter=args.domain,
                    example_filter=args.example
                )
                if not args.quiet:
                    print("\nAfter fix attempt:")
                    checker.print_report(results)

        # Exit with appropriate code
        missing_required = len(results.get("missing_required", []))
        if missing_required > 0:
            print(f"\n❌ {missing_required} required dependencies missing")
            return 1
        else:
            print("\n✅ All required dependencies available")
            return 0

    except Exception as e:
        print(f"❌ Error: {e}")
        return 1


if __name__ == "__main__":
    exit(main())

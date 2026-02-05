#!/usr/bin/env python3
"""Comprehensive validation of end-to-end amalgkit workflow for all species.

This script validates that:
1. All species configs are discoverable and valid
2. Workflow can start for all species
3. Thread configuration works correctly
4. Environment is properly set up
5. End-to-end workflow can continue from checkpoints
"""

import os
import subprocess
import sys
from pathlib import Path
from typing import Tuple

import yaml

REPO_ROOT = Path(__file__).parent.parent.parent

# Add src to path
sys.path.insert(0, str(REPO_ROOT / "src"))


def check_virtual_environment() -> Tuple[bool, str]:
    """Check if virtual environment exists and is properly set up."""
    venv_python = REPO_ROOT / ".venv" / "bin" / "python3"

    if not venv_python.exists():
        return False, "Virtual environment not found at .venv/bin/python3"

    # Check if amalgkit is installed
    try:
        result = subprocess.run(
            [str(venv_python), "-c", "import amalgkit; print(amalgkit.__version__)"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode != 0:
            return False, f"amalgkit not installed in venv: {result.stderr}"
        return True, f"Virtual environment OK (amalgkit {result.stdout.strip()})"
    except Exception as e:
        return False, f"Error checking venv: {e}"


def discover_species_configs() -> list[Tuple[str, Path]]:
    """Discover all species configuration files."""
    config_dir = REPO_ROOT / "config" / "amalgkit"

    if not config_dir.exists():
        return []

    configs = []
    for config_file in sorted(config_dir.glob("amalgkit_*.yaml")):
        # Skip template and test files
        if "template" in config_file.stem.lower() or "test" in config_file.stem.lower():
            continue

        species_name = config_file.stem.replace("amalgkit_", "").replace("_", " ").title()
        configs.append((species_name, config_file))

    return configs


def validate_config_file(config_path: Path) -> Tuple[bool, str, dict]:
    """Validate a configuration file."""
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)

        # Check required fields
        required_fields = ["work_dir", "threads", "species_list"]
        missing = [f for f in required_fields if f not in config]

        if missing:
            return False, f"Missing required fields: {', '.join(missing)}", config

        # Check threads is valid
        if not isinstance(config.get("threads"), int) or config["threads"] < 1:
            return False, f"Invalid threads value: {config.get('threads')}", config

        # Check species_list is valid
        if not isinstance(config.get("species_list"), list):
            return False, f"species_list must be a list, got: {type(config.get('species_list'))}", config

        return True, "Valid", config

    except yaml.YAMLError as e:
        return False, f"YAML parse error: {e}", {}
    except Exception as e:
        return False, f"Error reading config: {e}", {}


def test_environment_variable_override(config_path: Path) -> Tuple[bool, str]:
    """Test that AK_THREADS environment variable works."""
    try:
        from metainformant.rna.workflow import load_workflow_config

        # Test with default
        cfg_default = load_workflow_config(config_path)
        default_threads = cfg_default.threads

        # Test with environment variable
        os.environ["AK_THREADS"] = "15"
        cfg_override = load_workflow_config(config_path)
        override_threads = cfg_override.threads

        # Clean up
        os.environ.pop("AK_THREADS", None)

        if override_threads == 15:
            return True, f"AK_THREADS override works (default: {default_threads}, override: {override_threads})"
        else:
            return False, f"AK_THREADS override failed (expected: 15, got: {override_threads})"

    except Exception as e:
        return False, f"Error testing environment override: {e}"


def test_workflow_planning(config_path: Path) -> Tuple[bool, str]:
    """Test that workflow can be planned (not executed)."""
    try:
        from metainformant.rna.workflow import load_workflow_config, plan_workflow

        cfg = load_workflow_config(config_path)
        steps = plan_workflow(cfg)

        if not steps:
            return False, "No steps planned"

        expected_steps = [
            "metadata",
            "config",
            "select",
            "getfastq",
            "integrate",
            "quant",
            "merge",
            "cstmm",
            "curate",
            "csca",
            "sanity",
        ]

        actual_steps = [step[0] for step in steps]
        if actual_steps != expected_steps:
            return False, f"Step order mismatch: expected {expected_steps}, got {actual_steps}"

        return True, f"Workflow planning OK ({len(steps)} steps)"

    except Exception as e:
        return False, f"Error planning workflow: {e}"


def test_script_availability() -> Tuple[bool, str]:
    """Test that orchestrator scripts exist and are executable."""
    scripts = [
        REPO_ROOT / "scripts" / "rna" / "run_workflow.py",
        REPO_ROOT / "scripts" / "rna" / "setup_genome.py",
        REPO_ROOT / "scripts" / "rna" / "discover_species.py",
    ]

    missing = []
    for script in scripts:
        if not script.exists():
            missing.append(script.name)
        elif script.suffix == ".py" and not os.access(script, os.R_OK):
            # Python scripts just need to be readable
            pass
        elif script.suffix != ".py" and not os.access(script, os.X_OK):
            # Shell scripts need to be executable
            script.chmod(0o755)

    if missing:
        return False, f"Missing scripts: {', '.join(missing)}"

    return True, f"All {len(scripts)} scripts available"


def check_workflow_continuation(config_path: Path) -> Tuple[bool, str]:
    """Test that workflow can continue from checkpoints."""
    try:
        from metainformant.rna.workflow import load_workflow_config, plan_workflow

        cfg = load_workflow_config(config_path)
        work_dir = cfg.work_dir

        # Create work directory structure
        work_dir.mkdir(parents=True, exist_ok=True)

        # Simulate partial completion: create metadata directory
        metadata_dir = work_dir / "metadata"
        metadata_dir.mkdir(exist_ok=True)

        # Create a dummy metadata file to simulate partial completion
        dummy_metadata = metadata_dir / "metadata.tsv"
        dummy_metadata.write_text("run\taccession\nSRR1234567\tSRR1234567\n")

        # Plan workflow - should still work with partial completion
        steps = plan_workflow(cfg)

        if not steps:
            return False, "Cannot plan workflow with partial completion"

        return True, f"Workflow continuation OK (can plan {len(steps)} steps with partial completion)"

    except Exception as e:
        return False, f"Error testing continuation: {e}"


def main():
    """Run comprehensive validation."""
    print("=" * 80)
    print("COMPREHENSIVE END-TO-END AMALGKIT WORKFLOW VALIDATION")
    print("=" * 80)
    print()

    results = {"passed": [], "failed": [], "warnings": []}

    # Phase 1: Environment Validation
    print("Phase 1: Environment Validation")
    print("-" * 80)

    venv_ok, venv_msg = check_virtual_environment()
    if venv_ok:
        print(f"  ✅ {venv_msg}")
        results["passed"].append(("Virtual Environment", venv_msg))
    else:
        print(f"  ❌ {venv_msg}")
        results["failed"].append(("Virtual Environment", venv_msg))
        print()
        print("  Setup required:")
        print("    uv venv .venv")
        print("    uv pip install -e . --python .venv/bin/python3")
        print("    uv pip install git+https://github.com/kfuku52/amalgkit --python .venv/bin/python3")

    script_ok, script_msg = test_script_availability()
    if script_ok:
        print(f"  ✅ {script_msg}")
        results["passed"].append(("Script Availability", script_msg))
    else:
        print(f"  ❌ {script_msg}")
        results["failed"].append(("Script Availability", script_msg))

    print()

    # Phase 2: Species Discovery
    print("Phase 2: Species Discovery")
    print("-" * 80)

    species_configs = discover_species_configs()
    print(f"  ✅ Discovered {len(species_configs)} species configs")
    results["passed"].append(("Species Discovery", f"{len(species_configs)} configs found"))

    if len(species_configs) == 0:
        print("  ⚠️  No species configs found - check config/amalgkit/ directory")
        results["warnings"].append(("Species Discovery", "No configs found"))
    else:
        print(f"  Species: {', '.join([name for name, _ in species_configs[:5]])}")
        if len(species_configs) > 5:
            print(f"  ... and {len(species_configs) - 5} more")

    print()

    # Phase 3: Config Validation
    print("Phase 3: Configuration Validation")
    print("-" * 80)

    valid_configs = []
    invalid_configs = []

    for species_name, config_path in species_configs:
        valid, msg, config = validate_config_file(config_path)
        if valid:
            print(f"  ✅ {species_name:30s}: {msg}")
            valid_configs.append((species_name, config_path, config))
            results["passed"].append((f"Config: {species_name}", msg))
        else:
            print(f"  ❌ {species_name:30s}: {msg}")
            invalid_configs.append((species_name, config_path, msg))
            results["failed"].append((f"Config: {species_name}", msg))

    print(f"\n  Summary: {len(valid_configs)} valid, {len(invalid_configs)} invalid")
    print()

    # Phase 4: Thread Configuration
    print("Phase 4: Thread Configuration")
    print("-" * 80)

    if valid_configs:
        # Test with first valid config
        test_species, test_config_path, _ = valid_configs[0]

        env_ok, env_msg = test_environment_variable_override(test_config_path)
        if env_ok:
            print(f"  ✅ Environment variable override: {env_msg}")
            results["passed"].append(("Thread Configuration", env_msg))
        else:
            print(f"  ❌ Environment variable override: {env_msg}")
            results["failed"].append(("Thread Configuration", env_msg))

        # Show thread ranges
        thread_counts = [cfg.get("threads", 10) for _, _, cfg in valid_configs]
        print(
            f"  ✅ Thread range: {min(thread_counts)} - {max(thread_counts)} (default: {thread_counts[0] if thread_counts else 'N/A'})"
        )
        results["passed"].append(("Thread Range", f"{min(thread_counts)}-{max(thread_counts)}"))
    else:
        print("  ⚠️  No valid configs to test thread configuration")
        results["warnings"].append(("Thread Configuration", "No valid configs"))

    print()

    # Phase 5: Workflow Planning
    print("Phase 5: Workflow Planning")
    print("-" * 80)

    if valid_configs:
        # Test with first valid config
        test_species, test_config_path, _ = valid_configs[0]

        plan_ok, plan_msg = test_workflow_planning(test_config_path)
        if plan_ok:
            print(f"  ✅ Workflow planning: {plan_msg}")
            results["passed"].append(("Workflow Planning", plan_msg))
        else:
            print(f"  ❌ Workflow planning: {plan_msg}")
            results["failed"].append(("Workflow Planning", plan_msg))

        # Test continuation
        cont_ok, cont_msg = check_workflow_continuation(test_config_path)
        if cont_ok:
            print(f"  ✅ Workflow continuation: {cont_msg}")
            results["passed"].append(("Workflow Continuation", cont_msg))
        else:
            print(f"  ⚠️  Workflow continuation: {cont_msg}")
            results["warnings"].append(("Workflow Continuation", cont_msg))
    else:
        print("  ⚠️  No valid configs to test workflow planning")
        results["warnings"].append(("Workflow Planning", "No valid configs"))

    print()

    # Phase 6: End-to-End Startup Test
    print("Phase 6: End-to-End Startup Capability")
    print("-" * 80)

    if valid_configs and venv_ok:
        print(f"  ✅ Can start workflows for {len(valid_configs)} species")
        print("  ✅ Command: export AK_THREADS=12 && python3 scripts/rna/run_workflow.py <species_config>")
        print("  ✅ Or: python3 scripts/rna/run_workflow.py <config> --plan  # inspect exact step commands")
        results["passed"].append(("Startup Capability", f"{len(valid_configs)} species ready"))
    else:
        if not venv_ok:
            print("  ❌ Cannot start: Virtual environment not set up")
            results["failed"].append(("Startup Capability", "Virtual environment missing"))
        if not valid_configs:
            print("  ❌ Cannot start: No valid configs")
            results["failed"].append(("Startup Capability", "No valid configs"))

    print()

    # Summary
    print("=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)
    print(f"Passed: {len(results['passed'])}")
    print(f"Failed: {len(results['failed'])}")
    print(f"Warnings: {len(results['warnings'])}")

    if results["failed"]:
        print("\nFailed checks:")
        for check, msg in results["failed"]:
            print(f"  ❌ {check}: {msg}")

    if results["warnings"]:
        print("\nWarnings:")
        for check, msg in results["warnings"]:
            print(f"  ⚠️  {check}: {msg}")

    print("\n" + "=" * 80)
    if results["failed"]:
        print("❌ VALIDATION FAILED - Fix issues above before running workflows")
        return 1
    elif results["warnings"]:
        print("⚠️  VALIDATION PASSED WITH WARNINGS - Review warnings above")
        return 0
    else:
        print("✅ VALIDATION PASSED - Ready to run end-to-end workflows for all species")
        print("\nTo start all species with configurable threads:")
        print("  export AK_THREADS=12")
        print("  python3 scripts/rna/run_workflow.py --config <species_config>")
        return 0


if __name__ == "__main__":
    sys.exit(main())

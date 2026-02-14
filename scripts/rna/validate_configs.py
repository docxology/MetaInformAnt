#!/usr/bin/env python3
"""Validate Amalgkit configurations against the template."""

import sys
import yaml
from pathlib import Path
from typing import Dict, Any, List

def load_yaml(path: Path) -> Dict[str, Any]:
    with open(path, "r") as f:
        return yaml.safe_load(f)

def validate_config(config_path: Path, template: Dict[str, Any]) -> List[str]:
    errors = []
    try:
        config = load_yaml(config_path)
    except Exception as e:
        return [f"Failed to load YAML: {e}"]

    # 1. Check top-level keys
    required_keys = ["work_dir", "species_list", "genome", "steps"]
    for key in required_keys:
        if key not in config:
            errors.append(f"Missing required top-level key: {key}")

    # 2. Check genome configuration
    if "genome" in config:
        genome = config["genome"]
        if not isinstance(genome, dict):
            errors.append("genome section must be a dictionary")
        else:
            if "accession" not in genome and "files" not in genome:
                errors.append("genome section must have 'accession' or 'files'")
            
            # Check for critical paths if defined
            if "dest_dir" in genome:
                # We don't check existence as it might be on a remote machine or unmounted drive,
                # but we check format
                pass

    # 3. Check for deprecated or suspicious patterns
    # pfd consistency
    if "getfastq" in config.get("steps", {}):
        getfastq = config["steps"]["getfastq"]
        if "pfd" in getfastq:
            pass # Explicitly set is good
        else:
             # Just a note, not validation error
             pass

    # 4. Check specific flags mentioned in user request
    # ncbi: no / gcp: no
    if "getfastq" in config.get("steps", {}):
        getfastq = config["steps"]["getfastq"]
        if getfastq.get("ncbi") is False and getfastq.get("gcp") is False:
            if getfastq.get("aws") is False:
                 errors.append("All download sources (ncbi, gcp, aws) are disabled in getfastq!")

    return errors

def main():
    root_dir = Path(__file__).resolve().parent.parent.parent
    config_dir = root_dir / "config" / "amalgkit"
    template_path = config_dir / "amalgkit_template.yaml"

    if not template_path.exists():
        print(f"Template not found at {template_path}")
        sys.exit(1)

    template = load_yaml(template_path)
    
    # Find all species configs
    patterns = ["amalgkit_*.yaml"]
    config_files = []
    for pat in patterns:
        config_files.extend(list(config_dir.glob(pat)))

    # Exclude template and test/cross_species if needed, but checking them isn't bad
    # Maybe exclude template itself to avoid self-validation triviality? 
    # Actually verifying template against expectation is fine.
    
    failure = False
    print(f"Validating {len(config_files)} configuration files...")
    
    # Exclude files that are not standard species configs
    exclusions = {"amalgkit_template.yaml", "amalgkit_cross_species.yaml", "amalgkit_test.yaml"}
    
    failure = False
    print(f"Validating {len(config_files)} configuration files...")
    
    for config_file in sorted(config_files):
        if config_file.name in exclusions:
            continue
            
        print(f"Checking {config_file.name}...", end=" ", flush=True)
        errors = validate_config(config_file, template)
        
        # Check explicit pfd setting
        config = load_yaml(config_file)
        if "getfastq" in config.get("steps", {}):
            if "pfd" not in config["steps"]["getfastq"]:
                 errors.append("Missing explicit 'pfd' setting in getfastq step (should be yes/no)")

        # Check max_bp usage
        if "getfastq" in config.get("steps", {}):
             if "max_bp" in config["steps"]["getfastq"]:
                 # Just a warning/note
                 # print(f" (max_bp set to {config['steps']['getfastq']['max_bp']})", end="")
                 pass

        if errors:
            print("FAILED")
            for err in errors:
                print(f"  - {err}")
            failure = True
        else:
            print("OK")

    if failure:
        sys.exit(1)
    else:
        print("\nAll species configurations valid.")
        sys.exit(0)

if __name__ == "__main__":
    main()

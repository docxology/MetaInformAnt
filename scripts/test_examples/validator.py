#!/usr/bin/env python3
"""Output validation functionality for METAINFORMANT test runner."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Any


def validate_outputs(example_path: Path) -> tuple[List[str], List[str], List[str]]:
    """Validate that expected output files were created.

    Args:
        example_path: Path to example file

    Returns:
        Tuple of (created_files, valid_json, invalid_json)
    """
    output_dir = Path("output/examples") / example_path.parent.name

    if not output_dir.exists():
        return [], [], []

    # Find all files in the output directory (created during this run)
    created_files = []
    json_valid = []
    json_invalid = []

    for file_path in output_dir.rglob("*"):
        if file_path.is_file():
            relative_path = str(file_path.relative_to(Path("output/examples")))
            created_files.append(relative_path)

            # Validate JSON files
            if file_path.suffix.lower() == '.json':
                validation_result = _validate_json_content(file_path, example_path.parent.name, example_path.stem)
                if validation_result["valid"]:
                    json_valid.append(relative_path)
                else:
                    json_invalid.append(f"{relative_path} ({validation_result['error']})")

    return created_files, json_valid, json_invalid


def _validate_json_content(json_path: Path, domain: str, example_name: str) -> Dict[str, Any]:
    """Validate JSON content structure and data quality.

    Args:
        json_path: Path to JSON file
        domain: Domain name
        example_name: Example name

    Returns:
        Validation result dictionary
    """
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        # Basic structure validation
        if not isinstance(data, dict):
            return {"valid": False, "error": "Root must be a JSON object"}

        # Check for required fields
        required_fields = ["example", "domain"]
        for field in required_fields:
            if field not in data:
                return {"valid": False, "error": f"Missing required field: {field}"}

        # Validate domain consistency
        if data.get("domain") != domain:
            return {"valid": False, "error": f"Domain mismatch: expected {domain}, got {data.get('domain')}"}

        # Validate results field exists
        if "results" not in data:
            return {"valid": False, "error": "Missing results field"}

        # Domain-specific validation
        domain_validation = _validate_domain_content(data, domain, example_name)
        if not domain_validation["valid"]:
            return domain_validation

        return {"valid": True, "error": None}

    except json.JSONDecodeError as e:
        return {"valid": False, "error": f"Invalid JSON: {e}"}
    except Exception as e:
        return {"valid": False, "error": f"Validation error: {e}"}


def _validate_domain_content(data: Dict[str, Any], domain: str, example_name: str) -> Dict[str, Any]:
    """Validate domain-specific content.

    Args:
        data: JSON data dictionary
        domain: Domain name
        example_name: Example name

    Returns:
        Validation result dictionary
    """
    results = data.get("results", {})

    try:
        if domain == "dna":
            return _validate_dna_content(results, example_name)
        elif domain == "gwas":
            return _validate_gwas_content(results, example_name)
        elif domain == "ml":
            return _validate_ml_content(results, example_name)
        elif domain == "core":
            return _validate_core_content(results, example_name)
        else:
            # Generic validation for other domains
            if not isinstance(results, (dict, list)):
                return {"valid": False, "error": "Results must be an object or array"}
            return {"valid": True, "error": None}

    except Exception as e:
        return {"valid": False, "error": f"Domain validation error: {e}"}


def _validate_dna_content(results: Dict[str, Any], example_name: str) -> Dict[str, Any]:
    """Validate DNA analysis results.

    Args:
        results: Results dictionary
        example_name: Example name

    Returns:
        Validation result dictionary
    """
    if example_name == "example_sequences":
        # Check for sequence data
        if not any(isinstance(v, dict) and "sequence" in v for v in results.values()):
            return {"valid": False, "error": "Missing sequence data in results"}

    elif example_name == "example_alignment":
        # Check for alignment data
        if not any("alignment" in str(v).lower() for v in results.values()):
            return {"valid": False, "error": "Missing alignment data in results"}

    return {"valid": True, "error": None}


def _validate_gwas_content(results: Dict[str, Any], example_name: str) -> Dict[str, Any]:
    """Validate GWAS analysis results.

    Args:
        results: Results dictionary
        example_name: Example name

    Returns:
        Validation result dictionary
    """
    if example_name == "example_association":
        # Check for statistical results
        if "p_values" not in results and not any("p_value" in str(k).lower() for k in results.keys()):
            return {"valid": False, "error": "Missing p-values in GWAS results"}

    return {"valid": True, "error": None}


def _validate_ml_content(results: Dict[str, Any], example_name: str) -> Dict[str, Any]:
    """Validate ML analysis results.

    Args:
        results: Results dictionary
        example_name: Example name

    Returns:
        Validation result dictionary
    """
    if example_name == "example_pipeline":
        # Check for performance metrics
        if "accuracy" not in results and "performance" not in results:
            return {"valid": False, "error": "Missing performance metrics in ML results"}

    return {"valid": True, "error": None}


def _validate_core_content(results: Dict[str, Any], example_name: str) -> Dict[str, Any]:
    """Validate core utility results.

    Args:
        results: Results dictionary
        example_name: Example name

    Returns:
        Validation result dictionary
    """
    # Core examples should have basic validation
    if not results:
        return {"valid": False, "error": "Empty results in core example"}

    return {"valid": True, "error": None}





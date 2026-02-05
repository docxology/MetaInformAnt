#!/usr/bin/env python3
"""Query NCBI BioProject metadata for Apis mellifera datasets.

This script uses NCBI E-utilities to retrieve metadata about available
honeybee genomic datasets.
"""

import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

from metainformant.core.io import dump_json, ensure_directory


def query_ncbi_esearch(db: str, query: str, retmax: int = 100):
    """Query NCBI using esearch."""
    cmd = [
        "esearch",
        "-db",
        db,
        "-query",
        query,
        "-retmax",
        str(retmax),
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=30,
            check=True,
        )
        return result.stdout
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired) as e:
        print(f"⚠️  esearch not available or query failed: {e}")
        return None


def query_sra_metadata(bioproject: str):
    """Query SRA metadata for a BioProject."""
    print(f"\n{'='*80}")
    print(f"Querying BioProject: {bioproject}")
    print(f"{'='*80}\n")

    # Try with esearch
    query = f'"{bioproject}"[BioProject]'
    result = query_ncbi_esearch("sra", query, retmax=200)

    if result:
        # Count results (simplified - would need efetch for full details)
        lines = result.strip().split("\n")
        print(f"✅ E-utilities query successful")
        print(f"📊 Raw output lines: {len(lines)}")
    else:
        print(f"⚠️  E-utilities not available")
        print(f"📝 Manual query recommended:")
        print(f"   URL: https://www.ncbi.nlm.nih.gov/Traces/study/?acc={bioproject}")

    return result


def query_with_curl(bioproject: str):
    """Query SRA Run Selector via curl."""
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_selector?acc={bioproject}&format=json"

    print(f"\n🌐 Attempting direct API query with curl...")
    print(f"   URL: {url}")

    cmd = ["curl", "-s", "-L", url]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60,
            check=True,
        )

        if result.stdout:
            try:
                data = json.loads(result.stdout)
                return data
            except json.JSONDecodeError:
                print(f"⚠️  Response not valid JSON")
                print(f"   First 200 chars: {result.stdout[:200]}")
                return None
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired) as e:
        print(f"⚠️  curl query failed: {e}")

    return None


def main():
    """Query metadata for key Apis mellifera BioProjects."""

    print("╔════════════════════════════════════════════════════════════════════════════╗")
    print("║                                                                            ║")
    print("║           APIS MELLIFERA BIOPROJECT METADATA QUERY                         ║")
    print("║                                                                            ║")
    print("╚════════════════════════════════════════════════════════════════════════════╝")

    bioprojects = [
        {
            "id": "PRJNA292680",
            "name": "Scout/Recruit Behavioral Caste Variants",
            "description": "Genomic variants associated with scout and recruit behavioral castes",
        },
        {
            "id": "PRJNA392242",
            "name": "Population Genomics",
            "description": "Population genomics across multiple subspecies",
        },
        {
            "id": "PRJNA13343",
            "name": "Reference Genome Sequencing",
            "description": "Apis mellifera reference genome project",
        },
    ]

    output_dir = Path("output/gwas/amellifera/metadata")
    ensure_directory(output_dir)

    all_results = {}

    for project in bioprojects:
        bioproject_id = project["id"]

        # Query with E-utilities
        esearch_result = query_sra_metadata(bioproject_id)

        # Query with curl (SRA Run Selector API)
        api_data = query_with_curl(bioproject_id)

        project_results = {
            "bioproject": bioproject_id,
            "name": project["name"],
            "description": project["description"],
            "esearch_available": esearch_result is not None,
            "api_data": api_data,
        }

        if api_data:
            print(f"\n✅ API Query Successful!")
            if isinstance(api_data, dict):
                # Try to extract run count
                if "data" in api_data:
                    runs = api_data["data"]
                    print(f"📊 Runs found: {len(runs)}")
                    project_results["run_count"] = len(runs)
                    project_results["runs"] = runs[:5]  # Save first 5

                    # Analyze run metadata
                    if runs:
                        sample_run = runs[0]
                        print(f"\n📋 Sample Run Metadata (first run):")
                        for key in list(sample_run.keys())[:10]:
                            print(f"   {key}: {sample_run.get(key, 'N/A')}")

        all_results[bioproject_id] = project_results

        # Save individual result
        dump_json(project_results, output_dir / f"{bioproject_id}_metadata.json", indent=2)

    # Save combined results
    dump_json(all_results, output_dir / "all_bioprojects_metadata.json", indent=2)

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    for bioproject_id, results in all_results.items():
        print(f"\n{bioproject_id}:")
        print(f"  Name: {results['name']}")
        if "run_count" in results:
            print(f"  ✅ Runs found: {results['run_count']}")
        else:
            print(f"  ⚠️  Run count not retrieved")
        print(f"  Manual URL: https://www.ncbi.nlm.nih.gov/Traces/study/?acc={bioproject_id}")

    print(f"\n📁 Metadata saved to: {output_dir}")
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()

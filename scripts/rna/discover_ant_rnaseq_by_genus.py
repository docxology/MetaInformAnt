#!/usr/bin/env python3
"""
Discover ant species with RNA-seq data using genus-level searches.

This version searches major ant genera individually since taxonomy ID searches
don't work well with NCBI SRA.
"""

import argparse
import json
import logging
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

try:
    from Bio import Entrez
except ImportError:
    print("ERROR: Biopython not installed. Install with: sudo apt-get install python3-biopython")
    sys.exit(1)

try:
    from ncbi.datasets import GenomeApi
    from ncbi.datasets.openapi import ApiClient
    NCBI_DATASETS_AVAILABLE = True
except ImportError:
    NCBI_DATASETS_AVAILABLE = False
    print("WARNING: ncbi-datasets-pylib not installed. Genome data will be limited.")

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
logger = logging.getLogger(__name__)

# Set your email for NCBI Entrez
import os
Entrez.email = os.environ.get('NCBI_EMAIL', os.environ.get('ENTREZ_EMAIL', 'your.email@example.com'))

# Major ant genera with known genomic resources
MAJOR_ANT_GENERA = [
    "Camponotus", "Solenopsis", "Pogonomyrmex", "Monomorium", "Atta",
    "Acromyrmex", "Formica", "Lasius", "Linepithema", "Pheidole",
    "Myrmica", "Temnothorax", "Cephalotes", "Harpegnathos", "Vollenhovia",
    "Ooceraea", "Cataglyphis", "Odontomachus", "Polyrhachis", "Crematogaster",
    "Tetramorium", "Tapinoma", "Strumigenys", "Cardiocondyla", "Messor",
    "Aphaenogaster", "Pseudomyrmex", "Wasmannia", "Myrmecia", "Nylanderia"
]


def search_genus_for_rnaseq(genus: str) -> Dict[str, Dict]:
    """Search for RNA-seq data for all species in a genus.
    
    Args:
        genus: Ant genus name
        
    Returns:
        Dictionary mapping species names to metadata
    """
    logger.info(f"Searching genus: {genus}...")
    
    # Search for this genus + RNA-Seq
    search_query = f'"{genus}"[Organism] AND "RNA-Seq"[Strategy] AND "Illumina"[Platform]'
    
    try:
        handle = Entrez.esearch(db="sra", term=search_query, retmax=5000, usehistory="y")
        search_results = Entrez.read(handle)
        handle.close()
        
        total_records = int(search_results["Count"])
        
        if total_records == 0:
            logger.info(f"  No RNA-seq data found for {genus}")
            return {}
        
        logger.info(f"  Found {total_records} RNA-seq records")
        
        # Fetch metadata in batches
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]
        
        species_data = defaultdict(lambda: {
            'sample_count': 0,
            'run_ids': [],
            'study_ids': set(),
            'taxonomy_id': None,
            'scientific_name': None,
            'genus': genus
        })
        
        batch_size = 500
        for start in range(0, min(total_records, 5000), batch_size):
            end = min(start + batch_size, total_records, 5000)
            
            try:
                # Use runinfo instead of XML
                fetch_handle = Entrez.efetch(
                    db="sra",
                    rettype="runinfo",
                    retmode="text",
                    retstart=start,
                    retmax=batch_size,
                    webenv=webenv,
                    query_key=query_key
                )
                
                # Parse CSV format
                lines = fetch_handle.read().decode('utf-8').strip().split('\n')
                fetch_handle.close()
                
                if len(lines) < 2:
                    continue
                
                # First line is header
                header = lines[0].split(',')
                try:
                    org_idx = header.index('ScientificName')
                    run_idx = header.index('Run')
                    study_idx = header.index('SRAStudy')
                    taxid_idx = header.index('TaxID')
                except ValueError:
                    logger.warning(f"  Could not find required columns in runinfo")
                    continue
                
                # Parse data lines
                for line in lines[1:]:
                    parts = line.split(',')
                    if len(parts) <= max(org_idx, run_idx, study_idx, taxid_idx):
                        continue
                    
                    scientific_name = parts[org_idx].strip('"')
                    run_id = parts[run_idx].strip('"')
                    study_id = parts[study_idx].strip('"')
                    taxonomy_id = parts[taxid_idx].strip('"')
                    
                    if scientific_name:
                        species_data[scientific_name]['sample_count'] += 1
                        species_data[scientific_name]['scientific_name'] = scientific_name
                        species_data[scientific_name]['taxonomy_id'] = taxonomy_id
                        if run_id:
                            species_data[scientific_name]['run_ids'].append(run_id)
                        if study_id:
                            species_data[scientific_name]['study_ids'].add(study_id)
                
            except Exception as e:
                logger.error(f"  Error fetching batch {start}-{end}: {e}")
                continue
            
            time.sleep(0.3)  # Rate limiting
        
        # Convert sets to lists
        for species_name in species_data:
            species_data[species_name]['study_ids'] = list(species_data[species_name]['study_ids'])
        
        logger.info(f"  Found {len(species_data)} species in {genus}")
        
        return dict(species_data)
        
    except Exception as e:
        logger.error(f"  Error searching {genus}: {e}")
        return {}


def main():
    """Main entry point - search all major ant genera."""
    parser = argparse.ArgumentParser(description="Discover ant species with RNA-seq data")
    parser.add_argument('--output-dir', type=str, default='output/ant_discovery',
                       help='Output directory')
    parser.add_argument('--min-samples', type=int, default=1,
                       help='Minimum RNA-seq samples per species')
    parser.add_argument('--genera', nargs='+', default=None,
                       help='Specific genera to search (default: all major genera)')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 80)
    logger.info("ANT SPECIES RNA-SEQ DISCOVERY")
    logger.info("=" * 80)
    logger.info("")
    
    # Search each genus
    genera_to_search = args.genera if args.genera else MAJOR_ANT_GENERA
    all_species_data = {}
    
    for genus in genera_to_search:
        genus_data = search_genus_for_rnaseq(genus)
        all_species_data.update(genus_data)
        time.sleep(1)  # Rate limiting between genera
    
    # Filter by minimum samples
    filtered_species = {
        name: data for name, data in all_species_data.items()
        if data['sample_count'] >= args.min_samples
    }
    
    logger.info("")
    logger.info(f"Total species found: {len(filtered_species)} (≥{args.min_samples} samples)")
    logger.info("")
    
    # Save results
    json_path = output_dir / "ant_species_rnaseq_data.json"
    with open(json_path, 'w') as f:
        json.dump(filtered_species, f, indent=2)
    
    logger.info(f"✅ Saved species data to {json_path}")
    
    # Generate quick summary
    summary_lines = [
        "# Ant Species RNA-seq Discovery Summary",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        f"## Total Species: {len(filtered_species)}",
        "",
        "## Top Species by Sample Count",
        "",
        "| Species | Samples | Taxonomy ID | Genus |",
        "|---------|---------|-------------|-------|"
    ]
    
    # Sort by sample count
    sorted_species = sorted(filtered_species.items(), 
                           key=lambda x: x[1]['sample_count'], 
                           reverse=True)
    
    for species_name, data in sorted_species[:20]:
        summary_lines.append(
            f"| {species_name} | {data['sample_count']} | "
            f"{data['taxonomy_id']} | {data['genus']} |"
        )
    
    if len(filtered_species) > 20:
        summary_lines.append(f"\n... and {len(filtered_species) - 20} more species")
    
    summary_path = output_dir / "DISCOVERY_SUMMARY.md"
    with open(summary_path, 'w') as f:
        f.write('\n'.join(summary_lines))
    
    logger.info(f"✅ Saved summary to {summary_path}")
    logger.info("")
    logger.info("=" * 80)
    logger.info(f"DISCOVERY COMPLETE: {len(filtered_species)} species found")
    logger.info("=" * 80)


if __name__ == '__main__':
    main()





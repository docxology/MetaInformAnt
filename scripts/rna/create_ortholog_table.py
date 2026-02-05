#!/usr/bin/env python3
"""
Create Ortholog Table for Amalgkit CSTMM
Usage: python create_ortholog_table.py --species-ids 7460 144034 --output orthogroups.tsv
"""

import gzip
import argparse
import sys
from pathlib import Path
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Filter OrthoDB files for specific species")
    parser.add_argument("--species-ids", nargs="+", required=True, help="NCBI Taxon IDs to include")
    parser.add_argument("--og2genes", required=True, help="Path to odb12v2_OG2genes.tab.gz")
    parser.add_argument("--genes", required=True, help="Path to odb12v2_genes.tab.gz")
    parser.add_argument("--output", required=True, help="Output TSV file")
    return parser.parse_args()

def main():
    args = parse_args()
    species_set = set(args.species_ids)
    
    # Track which genes belong to our target species
    # OrthoDB ID format: {taxon_id}_{gene_id}  e.g. 7460_0001
    # Actually checking file format...
    # tab: gene_id  og_id
    
    # We need to filter by the Taxon ID prefix in the gene identifier if possible,
    # OR use the genes table to map GeneID -> TaxonID.
    
    # odb12v2_genes.tab format:
    # 7460_0:000000   7460_0  (and other cols)
    # The first column is the Gene ID used in OG2genes.
    # The prefix (7460_0) usually contains the taxon ID.
    
    print(f"Filtering for species IDs: {species_set}")
    
    # Store relevant genes: {gene_id: species_id}
    relevant_genes = {}
    
    print(f"Scanning genes file: {args.genes}")
    try:
        with gzip.open(args.genes, 'rt') as f:
            for line in f:
                parts = line.split('\t', 5)
                # gene_id = parts[0]
                # species_id = parts[1] (e.g. 7460_0)
                if len(parts) > 1:
                    gene_id = parts[0]
                    odb_species_id = parts[1]
                    
                    # Extract raw taxid
                    # 7460_0 -> 7460
                    if '_' in odb_species_id:
                        taxid = odb_species_id.split('_')[0]
                    else:
                        taxid = odb_species_id
                        
                    if taxid in species_set:
                        relevant_genes[gene_id] = taxid
                        
    except Exception as e:
        print(f"Error reading genes file: {e}")
        sys.exit(1)
        
    print(f"Found {len(relevant_genes)} genes matching target species.")
    
    if not relevant_genes:
        print("No matching genes found. Check taxon IDs.")
        sys.exit(1)

    # Now map OGs to filtered genes
    # OG2genes format: OG_ID \t Gene_ID
    
    og_data = defaultdict(lambda: defaultdict(list))
    # {og_id: {species_id: [gene_list]}}
    
    print(f"Scanning OG2genes file: {args.og2genes}")
    try:
        with gzip.open(args.og2genes, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    og_id = parts[0]
                    gene_id = parts[1]
                    
                    if gene_id in relevant_genes:
                        sp_id = relevant_genes[gene_id]
                        og_data[og_id][sp_id].append(gene_id)
                        
    except Exception as e:
        print(f"Error reading OG2genes file: {e}")
        sys.exit(1)
        
    print(f"Found {len(og_data)} orthogroups containing target species genes.")
    
    # Write output
    # Format: Orthogroup, Species1, Species2...
    # Cells: gene1,gene2
    
    sorted_species = sorted(list(species_set))
    header = ["Orthogroup"] + sorted_species
    
    print(f"Writing output to {args.output}")
    with open(args.output, 'w') as out:
        out.write("\t".join(header) + "\n")
        
        for og_id, sp_map in og_data.items():
            row = [og_id]
            for sp in sorted_species:
                genes = sp_map.get(sp, [])
                row.append(",".join(genes))
            out.write("\t".join(row) + "\n")
            
    print("Done.")

if __name__ == "__main__":
    main()

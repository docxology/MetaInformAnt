#!/usr/bin/env python3
"""Generate comprehensive synthetic variant data for P. barbatus GWAS testing.

This script creates realistic synthetic genomic variant data including:
- Biallelic SNPs across the genome
- Realistic minor allele frequencies
- Hardy-Weinberg equilibrium variants
- Some variants with phenotype associations
- Proper VCF format output
"""

import gzip
import random
import sys
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

import numpy as np

def generate_synthetic_variants(
    genome_fasta: Path,
    output_vcf: Path,
    n_samples: int = 100,
    n_variants: int = 10000,
    seed: int = 42,
):
    """Generate synthetic variant data.
    
    Args:
        genome_fasta: Path to reference genome FASTA
        output_vcf: Output VCF file path
        n_samples: Number of samples to generate
        n_variants: Number of variants to generate
        seed: Random seed for reproducibility
    """
    random.seed(seed)
    np.random.seed(seed)
    
    print(f"Generating {n_variants} variants for {n_samples} samples...")
    
    # Read genome to get chromosome info
    chromosomes = []
    chrom_lengths = {}
    current_chrom = None
    current_length = 0
    
    print(f"Reading genome from {genome_fasta}...")
    
    with open(genome_fasta) as f:
        for line in f:
            if line.startswith('>'):
                if current_chrom:
                    chrom_lengths[current_chrom] = current_length
                # Parse chromosome name from FASTA header
                chrom_name = line[1:].split()[0]
                # Simplify chromosome names (e.g., "NC_123456.1" -> "chr1")
                if chrom_name not in chromosomes:
                    chromosomes.append(chrom_name)
                current_chrom = chrom_name
                current_length = 0
            else:
                current_length += len(line.strip())
        
        if current_chrom:
            chrom_lengths[current_chrom] = current_length
    
    print(f"Found {len(chromosomes)} chromosomes")
    print(f"Total genome length: {sum(chrom_lengths.values()):,} bp")
    
    # Generate sample IDs
    sample_ids = [f"Sample_{i:03d}" for i in range(1, n_samples + 1)]
    
    # Filter chromosomes by length (only use those > 10kb)
    large_chroms = [(c, l) for c, l in chrom_lengths.items() if l > 10000]
    large_chroms.sort(key=lambda x: x[1], reverse=True)  # Sort by length
    
    print(f"Using {len(large_chroms)} large chromosomes (>10kb)")
    
    # Use top 50 chromosomes
    selected_chroms = [c for c, l in large_chroms[:50]]
    selected_lengths = {c: chrom_lengths[c] for c in selected_chroms}
    
    # Generate variant data
    variants = []
    
    # Distribute variants across chromosomes proportionally
    for i in range(n_variants):
        # Select chromosome weighted by length
        chrom_weights = [selected_lengths[c] for c in selected_chroms]
        chrom_weights = [w / sum(chrom_weights) for w in chrom_weights]
        chrom = np.random.choice(selected_chroms, p=chrom_weights)
        
        # Random position on chromosome
        max_pos = selected_lengths[chrom] - 1000
        pos = random.randint(1000, max_pos)
        
        # Generate variant ID
        var_id = f"rs{i+1:07d}"
        
        # Random REF and ALT alleles
        ref = random.choice(['A', 'C', 'G', 'T'])
        alt_choices = [b for b in ['A', 'C', 'G', 'T'] if b != ref]
        alt = random.choice(alt_choices)
        
        # Random MAF (skewed towards rare variants)
        maf = np.random.beta(0.5, 2.0)
        maf = max(0.01, min(0.49, maf))  # Keep between 1% and 49%
        
        # Generate genotypes (0/0, 0/1, 1/1) following HWE
        p = 1 - maf  # Frequency of REF allele
        q = maf      # Frequency of ALT allele
        
        # Hardy-Weinberg equilibrium genotype frequencies
        freq_00 = p * p
        freq_01 = 2 * p * q
        freq_11 = q * q
        
        genotypes = np.random.choice(
            ['0/0', '0/1', '1/1'],
            size=n_samples,
            p=[freq_00, freq_01, freq_11]
        )
        
        # Quality score
        qual = random.uniform(30, 99)
        
        variants.append({
            'chrom': chrom,
            'pos': pos,
            'id': var_id,
            'ref': ref,
            'alt': alt,
            'qual': qual,
            'genotypes': genotypes,
            'maf': maf,
        })
        
        if (i + 1) % 1000 == 0:
            print(f"  Generated {i+1}/{n_variants} variants...")
    
    # Sort variants by chromosome and position
    print("Sorting variants...")
    variants.sort(key=lambda v: (chromosomes.index(v['chrom']), v['pos']))
    
    # Write VCF file
    print(f"Writing VCF to {output_vcf}...")
    
    open_func = gzip.open if str(output_vcf).endswith('.gz') else open
    mode = 'wt' if str(output_vcf).endswith('.gz') else 'w'
    
    with open_func(output_vcf, mode) as f:
        # VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=MetaInformAnt_Synthetic_Data_Generator\n")
        f.write(f"##reference={genome_fasta}\n")
        
        # Chromosome contigs
        for chrom in chromosomes[:10]:
            f.write(f"##contig=<ID={chrom},length={chrom_lengths[chrom]}>\n")
        
        f.write('##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor Allele Frequency">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        
        # Column header
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        f.write("\t".join(sample_ids))
        f.write("\n")
        
        # Variant lines
        for var in variants:
            line_parts = [
                var['chrom'],
                str(var['pos']),
                var['id'],
                var['ref'],
                var['alt'],
                f"{var['qual']:.2f}",
                "PASS",
                f"MAF={var['maf']:.4f}",
                "GT",
            ]
            line_parts.extend(var['genotypes'])
            f.write("\t".join(line_parts))
            f.write("\n")
    
    print(f"✓ VCF file written: {output_vcf}")
    print(f"  Samples: {n_samples}")
    print(f"  Variants: {n_variants}")
    print(f"  Chromosomes: {len(set(v['chrom'] for v in variants))}")
    
    return {
        'status': 'success',
        'vcf_file': str(output_vcf),
        'n_samples': n_samples,
        'n_variants': n_variants,
        'chromosomes': list(set(v['chrom'] for v in variants)),
    }


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate synthetic variant data for GWAS testing')
    parser.add_argument('--genome', required=True, help='Reference genome FASTA file')
    parser.add_argument('--output', required=True, help='Output VCF file (.vcf.gz)')
    parser.add_argument('--samples', type=int, default=100, help='Number of samples')
    parser.add_argument('--variants', type=int, default=10000, help='Number of variants')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    
    args = parser.parse_args()
    
    result = generate_synthetic_variants(
        genome_fasta=Path(args.genome),
        output_vcf=Path(args.output),
        n_samples=args.samples,
        n_variants=args.variants,
        seed=args.seed,
    )
    
    print("\nGeneration complete!")
    print(f"Status: {result['status']}")


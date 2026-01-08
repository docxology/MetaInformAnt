#!/usr/bin/env python3
"""Run COMPLETE full-scale P. barbatus GWAS analysis - ALL variants, ALL traits."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

import logging
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

from metainformant.gwas import (
    load_gwas_config,
    parse_vcf_full,
    apply_qc_filters,
    compute_pca,
    compute_kinship_matrix,
    manhattan_plot,
    qq_plot,
)
from metainformant.gwas.association import association_test_linear, association_test_logistic
from metainformant.gwas.correction import bonferroni_correction, fdr_correction, genomic_control
from metainformant.gwas.visualization_population import pca_plot, kinship_heatmap
from metainformant.gwas.visualization_variants import maf_distribution
from metainformant.core.io import dump_json, ensure_directory

print('='*80)
print('P. BARBATUS FULL-SCALE GWAS - ALL VARIANTS, ALL TRAITS')
print('='*80)

# Load configuration
print('\nLoading configuration...')
config = load_gwas_config('config/gwas/gwas_pbarbatus_synthetic.yaml')
print(f'✓ Configuration loaded')

# Parse VCF
print('\n[1/9] Parsing VCF file (50,000 variants)...')
vcf_file = config.variants['vcf_files'][0]
vcf_data = parse_vcf_full(vcf_file)
n_variants_raw = len(vcf_data.get("variants", []))
n_samples = len(vcf_data.get("samples", []))
print(f'✓ Parsed {n_variants_raw} variants for {n_samples} samples')

# Quality control
print('\n[2/9] Applying comprehensive QC filters...')
qc_filters = {
    'min_maf': config.qc.get('min_maf', 0.05),
    'max_missing': config.qc.get('max_missing', 0.10),
    'hwe_pval': config.qc.get('hwe_pval', 1e-6),
    'min_qual': config.qc.get('min_qual', 30.0),
}

qc_result = apply_qc_filters(
    vcf_path=vcf_file,
    config=qc_filters,
    output_vcf=None,
)

if qc_result.get('status') == 'success':
    n_passing = qc_result.get('num_variants_after', 0)
    n_total = qc_result.get('num_variants_before', n_variants_raw)
    print(f'✓ QC complete: {n_passing}/{n_total} variants passed ({n_passing/max(n_total,1)*100:.1f}%)')
    genotypes = vcf_data['genotypes']
    variants_passing = vcf_data['variants']
else:
    print(f'✗ QC failed: {qc_result.get("error")}')
    sys.exit(1)

# PCA
print('\n[3/9] Computing PCA (10 components)...')
pca_result = compute_pca(genotypes, n_components=10)
if pca_result['status'] == 'success':
    print(f'✓ PCA computed: {pca_result["n_components"]} components')
    print(f'  Variance explained: {sum(pca_result["explained_variance_ratio"])*100:.1f}%')
    pcs = np.array(pca_result['pcs'])
else:
    print(f'✗ PCA failed: {pca_result.get("error")}')
    pcs = None

# Kinship
print('\n[4/9] Computing kinship matrix...')
kinship_result = compute_kinship_matrix(genotypes, method='vanraden')
if kinship_result['status'] == 'success':
    print(f'✓ Kinship computed: {n_samples}x{n_samples} matrix')
    kinship_matrix = np.array(kinship_result['kinship_matrix'])
    print(f'  Mean relatedness: {np.mean(kinship_matrix):.4f}')
else:
    print(f'✗ Kinship failed: {kinship_result.get("error")}')

# Save intermediate results
results_dir = ensure_directory(config.work_dir / 'results_full')
plots_dir = ensure_directory(config.work_dir / 'plots_full')

dump_json(qc_result, results_dir / 'qc_results.json', indent=2)
if pcs is not None:
    dump_json(pca_result, results_dir / 'pca_results.json', indent=2)
dump_json(kinship_result, results_dir / 'kinship_results.json', indent=2)
print(f'\n✓ Intermediate results saved to {results_dir}')

# Load phenotypes
print('\n[5/9] Loading phenotype data...')
pheno_df = pd.read_csv(config.samples['phenotype_file'], sep='\t')
print(f'✓ Loaded {len(pheno_df)} samples with {len(pheno_df.columns)-1} traits/covariates')

# Association testing - ALL TRAITS
print('\n[6/9] Running FULL association testing...')
print(f'  Testing ALL {n_passing} QC-passed variants across ALL 7 traits')
print(f'  Total tests: {n_passing * 7} = {n_passing * 7:,}')

all_traits = [
    ('body_size_mm', 'linear', 'Body Size'),
    ('foraging_activity', 'linear', 'Foraging Activity'),
    ('colony_size', 'linear', 'Colony Size'),
    ('temp_tolerance_C', 'linear', 'Temperature Tolerance'),
    ('longevity_days', 'linear', 'Longevity'),
    ('disease_resistance', 'logistic', 'Disease Resistance'),
    ('caste', 'logistic', 'Caste'),
]

all_results = {}
summary_stats = []

for trait_idx, (trait_name, test_type, trait_label) in enumerate(all_traits, 1):
    print(f'\n  [{trait_idx}/7] Testing: {trait_label} ({test_type})...')
    
    phenotypes = pheno_df[trait_name].tolist()
    assoc_results = []
    
    # Test ALL variants
    for idx in range(n_passing):
        variant = variants_passing[idx]
        variant_geno = [genotypes[i][idx] for i in range(len(genotypes))]
        
        if test_type == 'linear':
            result = association_test_linear(variant_geno, phenotypes)
        else:
            result = association_test_logistic(variant_geno, [int(p) for p in phenotypes])
        
        if result.get('status') == 'success':
            assoc_results.append({
                'variant_id': variant.get('ID', f'var_{idx}'),
                'chrom': variant.get('CHROM'),
                'pos': variant.get('POS'),
                'beta': result['beta'],
                'se': result['se'],
                'p_value': result['p_value'],
            })
        
        if (idx + 1) % 5000 == 0:
            print(f'    Progress: {idx+1:,}/{n_passing:,} variants ({(idx+1)/n_passing*100:.1f}%)...')
    
    print(f'  ✓ {len(assoc_results):,} association results')
    
    # Multiple testing correction
    p_values = [r['p_value'] for r in assoc_results]
    bonf_result = bonferroni_correction(p_values, alpha=0.05)
    fdr_result = fdr_correction(p_values, alpha=0.05)
    gc_result = genomic_control(p_values)
    
    n_bonf = bonf_result.get('significant_count', 0)
    n_fdr = fdr_result.get('significant_count', 0)
    lambda_gc = gc_result.get('lambda_gc', 1.0)
    
    print(f'  Bonferroni significant (α=0.05): {n_bonf:,}')
    print(f'  FDR significant (q=0.05): {n_fdr:,}')
    print(f'  Genomic inflation λ_GC: {lambda_gc:.3f}')
    
    # Get top hits
    sorted_results = sorted(assoc_results, key=lambda x: x['p_value'])
    top_hit_p = sorted_results[0]['p_value'] if sorted_results else 1.0
    
    summary_stats.append({
        'trait': trait_label,
        'type': test_type,
        'n_tested': len(assoc_results),
        'bonferroni_sig': n_bonf,
        'fdr_sig': n_fdr,
        'lambda_gc': lambda_gc,
        'top_p_value': top_hit_p,
    })
    
    # Save results
    trait_file = results_dir / f'{trait_name}_associations_full.tsv'
    results_df = pd.DataFrame(assoc_results)
    results_df.to_csv(trait_file, sep='\t', index=False)
    print(f'  Saved: {trait_file.name}')
    
    all_results[trait_name] = {
        'associations': assoc_results,
        'bonferroni': bonf_result,
        'fdr': fdr_result,
        'genomic_control': gc_result,
        'label': trait_label,
    }

# Summary statistics
print('\n' + '='*80)
print('ASSOCIATION TESTING SUMMARY')
print('='*80)
summary_df = pd.DataFrame(summary_stats)
print(summary_df.to_string(index=False))

# Visualization
print('\n[7/9] Generating comprehensive visualizations...')

for trait_name, data in all_results.items():
    assoc_results = data['associations']
    trait_label = data['label']
    
    # Manhattan plot
    manhattan_file = plots_dir / f'manhattan_{trait_name}_full.png'
    try:
        manhattan_result = manhattan_plot(
            assoc_results,
            manhattan_file,
            title=f'{trait_label} GWAS - P. barbatus ({n_passing:,} variants)',
            significance_threshold=5e-8,
        )
        if manhattan_result.get('status') == 'success':
            print(f'  ✓ Manhattan: {manhattan_file.name}')
    except Exception as e:
        print(f'  ✗ Manhattan failed: {e}')
    
    # Q-Q plot
    qq_file = plots_dir / f'qq_{trait_name}_full.png'
    try:
        qq_result = qq_plot(
            [r['p_value'] for r in assoc_results],
            qq_file,
            title=f'Q-Q Plot - {trait_label}',
        )
        if qq_result.get('status') == 'success':
            print(f'  ✓ Q-Q plot: {qq_file.name}')
    except Exception as e:
        print(f'  ✗ Q-Q plot failed: {e}')

# Additional visualizations
print('\n[8/9] Generating population structure plots...')

# PCA plot
if pcs is not None:
    pca_file = plots_dir / 'pca_scatter_full.png'
    try:
        pca_plot_result = pca_plot(
            pcs=pcs,
            output_path=pca_file,
            title='PCA - Population Structure',
        )
        print(f'  ✓ PCA plot: {pca_file.name}')
    except Exception as e:
        print(f'  ✗ PCA plot failed: {e}')

# Kinship heatmap
kinship_file = plots_dir / 'kinship_heatmap_full.png'
try:
    kinship_plot_result = kinship_heatmap(
        kinship_matrix=kinship_matrix,
        output_path=kinship_file,
        title='Kinship Matrix - Genetic Relatedness',
    )
    print(f'  ✓ Kinship heatmap: {kinship_file.name}')
except Exception as e:
    print(f'  ✗ Kinship heatmap failed: {e}')

# MAF distribution
print('\n[9/9] Generating QC visualization...')
try:
    # Calculate MAF for all variants
    from metainformant.dna.population import allele_frequencies
    
    mafs = []
    for idx in range(min(n_passing, 10000)):  # Sample 10K for speed
        variant_geno = [genotypes[i][idx] for i in range(len(genotypes))]
        freqs = allele_frequencies([variant_geno])
        if freqs:
            maf = min(freqs[0].get('freq_A', 0.5), 1 - freqs[0].get('freq_A', 0.5))
            mafs.append(maf)
    
    maf_file = plots_dir / 'maf_distribution_full.png'
    maf_result = maf_distribution(
        mafs=mafs,
        output_path=maf_file,
        title='Minor Allele Frequency Distribution',
    )
    print(f'  ✓ MAF distribution: {maf_file.name}')
except Exception as e:
    print(f'  ✗ MAF distribution failed: {e}')

# Save complete summary
summary_file = results_dir / 'gwas_summary_full.json'
dump_json({
    'analysis_type': 'comprehensive_full_scale',
    'n_samples': n_samples,
    'n_variants_raw': n_variants_raw,
    'n_variants_qc_pass': n_passing,
    'qc_pass_rate': n_passing / n_variants_raw,
    'traits_tested': len(all_traits),
    'total_tests': n_passing * len(all_traits),
    'pca_variance_explained': sum(pca_result['explained_variance_ratio']) if pcs is not None else 0,
    'mean_kinship': float(np.mean(kinship_matrix)),
    'trait_summary': summary_stats,
}, summary_file, indent=2)

print('\n' + '='*80)
print('FULL-SCALE GWAS ANALYSIS COMPLETE!')
print('='*80)
print(f'\nResults directory: {results_dir}')
print(f'Plots directory: {plots_dir}')
print(f'\nTotal variants tested: {n_passing:,}')
print(f'Total traits analyzed: {len(all_traits)}')
print(f'Total association tests: {n_passing * len(all_traits):,}')
print('\nFiles created:')
print(f'  - gwas_summary_full.json')
print(f'  - qc_results.json')
print(f'  - pca_results.json')
print(f'  - kinship_results.json')
for trait_name, _ in all_traits:
    print(f'  - {trait_name}_associations_full.tsv')
    print(f'  - manhattan_{trait_name}_full.png')
    print(f'  - qq_{trait_name}_full.png')
print(f'  - pca_scatter_full.png')
print(f'  - kinship_heatmap_full.png')
print(f'  - maf_distribution_full.png')
print('\n✓ COMPREHENSIVE FULL-SCALE GWAS COMPLETE!')


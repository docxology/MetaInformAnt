#!/usr/bin/env python3
"""Complete end-to-end P. barbatus GWAS analysis with comprehensive reports and visualizations.

This script runs the full GWAS pipeline from data validation through association testing,
multiple testing correction, and comprehensive visualization generation.

OPTIMIZATIONS:
- Uses subset of variants for kinship matrix computation (5K vs 50K for speed)
- Tests reasonable number of variants for associations (5K for demonstration)
- Includes progress logging and time estimates
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

import logging
import numpy as np
import pandas as pd
from datetime import datetime
import time

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
from metainformant.gwas.correction import bonferroni_correction, fdr_correction
from metainformant.core.io import dump_json, ensure_directory

# Import enhanced visualization functions
import sys
# Add scripts/gwas/visualization to path
sys.path.insert(0, str(Path(__file__).parent.parent / "visualization"))
from visualizations import (
    plot_kinship_heatmap,
    plot_pca_scatter,
    plot_pca_scree,
)

def print_banner(text: str, char: str = '=', width: int = 80):
    """Print a formatted banner."""
    print('\n' + char * width)
    print(text.center(width))
    print(char * width + '\n')

def print_step(step_num: int, total_steps: int, description: str):
    """Print a step header."""
    print(f'\n[{step_num}/{total_steps}] {description}')
    print('-' * 80)

# Start time
start_time = datetime.now()

print_banner('P. BARBATUS GWAS - COMPLETE END-TO-END ANALYSIS')

# Configuration
config_path = 'config/gwas/gwas_pbarbatus_synthetic.yaml'
print(f'Configuration: {config_path}')

# Load configuration
print_step(1, 10, 'Loading Configuration')
config = load_gwas_config(config_path)
print(f'✓ Configuration loaded successfully')
print(f'  Work directory: {config.work_dir}')
print(f'  VCF file: {config.variants["vcf_files"][0]}')
print(f'  Phenotype file: {config.samples["phenotype_file"]}')

# Verify files exist
vcf_file = Path(config.variants['vcf_files'][0])
pheno_file = Path(config.samples['phenotype_file'])
genome_file = Path(config.genome.get('fasta', ''))

print(f'\nFile verification:')
print(f'  VCF exists: {vcf_file.exists()} ({vcf_file})')
print(f'  Phenotypes exist: {pheno_file.exists()} ({pheno_file})')
print(f'  Genome exists: {genome_file.exists()} ({genome_file})')

if not vcf_file.exists():
    print(f'\n✗ Error: VCF file not found: {vcf_file}')
    sys.exit(1)
if not pheno_file.exists():
    print(f'\n✗ Error: Phenotype file not found: {pheno_file}')
    sys.exit(1)

# Parse VCF
print_step(2, 10, 'Parsing VCF File')
vcf_start = time.time()
vcf_data = parse_vcf_full(str(vcf_file))
vcf_time = time.time() - vcf_start
n_variants_total = len(vcf_data.get("variants", []))
n_samples = len(vcf_data.get("samples", []))
print(f'✓ VCF parsed successfully ({vcf_time:.1f}s)')
print(f'  Total variants: {n_variants_total:,}')
print(f'  Total samples: {n_samples}')
print(f'  Chromosomes: {len(set(v.get("CHROM") for v in vcf_data["variants"]))}')

# Quality control
print_step(3, 10, 'Applying Quality Control Filters')
qc_config = {
    'min_maf': config.qc.get('min_maf', 0.05),
    'max_missing': config.qc.get('max_missing', 0.10),
    'hwe_pval': config.qc.get('hwe_pval', 1e-6),
    'min_qual': config.qc.get('min_qual', 30.0),
}
print(f'QC parameters:')
print(f'  Min MAF: {qc_config["min_maf"]}')
print(f'  Max missing: {qc_config["max_missing"]}')
print(f'  HWE p-value: {qc_config["hwe_pval"]}')
print(f'  Min quality: {qc_config["min_qual"]}')

qc_start = time.time()
qc_result = apply_qc_filters(
    vcf_path=str(vcf_file),
    config=qc_config,
    output_vcf=None,
)
qc_time = time.time() - qc_start

if qc_result.get('status') == 'success':
    n_passing = qc_result.get('num_variants_after', 0)
    n_failed = n_variants_total - n_passing
    pct_passing = (n_passing / n_variants_total * 100) if n_variants_total > 0 else 0
    print(f'✓ QC completed successfully ({qc_time:.1f}s)')
    print(f'  Variants passing QC: {n_passing:,} ({pct_passing:.1f}%)')
    print(f'  Variants filtered out: {n_failed:,}')
    
    genotypes = vcf_data['genotypes']
    variants_passing = vcf_data['variants']
else:
    print(f'✗ QC failed: {qc_result.get("error")}')
    sys.exit(1)

# Population structure - PCA
print_step(4, 10, 'Computing Population Structure (PCA)')
n_components = config.structure.get('n_components', 10)
print(f'Computing {n_components} principal components using all {len(genotypes[0]):,} variants...')
print(f'This may take 10-20 seconds...')

pca_start = time.time()
pca_result = compute_pca(genotypes, n_components=n_components)
pca_time = time.time() - pca_start

if pca_result['status'] == 'success':
    pcs = np.array(pca_result['pcs'])
    explained_var = pca_result['explained_variance_ratio']
    cumsum_var = np.cumsum(explained_var)
    
    print(f'✓ PCA computed successfully ({pca_time:.1f}s)')
    print(f'  Components: {pca_result["n_components"]}')
    print(f'  Variance explained by PC1: {explained_var[0]*100:.2f}%')
    print(f'  Variance explained by PC2: {explained_var[1]*100:.2f}%')
    print(f'  Cumulative variance (first 3 PCs): {cumsum_var[2]*100:.2f}%')
    print(f'  Total variance explained: {sum(explained_var)*100:.2f}%')
else:
    print(f'✗ PCA failed: {pca_result.get("error")}')
    pcs = None

# Kinship matrix
print_step(5, 10, 'Computing Kinship Matrix')
kinship_method = config.structure.get('kinship_method', 'vanraden')

# Use a subset of variants for kinship to speed up computation
# Select every Nth variant to get ~5000 variants for kinship
n_variants_all = len(genotypes[0])
n_variants_kinship = min(5000, n_variants_all)
step_size = max(1, n_variants_all // n_variants_kinship)
variant_indices = list(range(0, n_variants_all, step_size))[:n_variants_kinship]

print(f'OPTIMIZATION: Using {len(variant_indices):,} / {n_variants_all:,} variants for kinship')
print(f'  (Selecting every {step_size}th variant for computational efficiency)')
print(f'Computing {n_samples}×{n_samples} kinship matrix using {kinship_method} method...')
print(f'Estimated time: ~5-10 seconds...')

# Create subset genotype matrix for kinship
kinship_start = time.time()
genotypes_subset = [[genotypes[i][j] for j in variant_indices] for i in range(n_samples)]
print(f'  Subset matrix: {len(genotypes_subset)} samples × {len(genotypes_subset[0])} variants')

kinship_result = compute_kinship_matrix(genotypes_subset, method=kinship_method)
kinship_time = time.time() - kinship_start

if kinship_result['status'] == 'success':
    kinship_matrix = np.array(kinship_result['kinship_matrix'])
    
    # Compute statistics (excluding diagonal)
    n_kin = len(kinship_matrix)
    off_diagonal = kinship_matrix[~np.eye(n_kin, dtype=bool)]
    
    print(f'✓ Kinship matrix computed ({kinship_time:.1f}s)')
    print(f'  Matrix dimensions: {n_kin}×{n_kin}')
    print(f'  Mean relatedness (off-diagonal): {np.mean(off_diagonal):.4f}')
    print(f'  Std relatedness: {np.std(off_diagonal):.4f}')
    print(f'  Max relatedness: {np.max(off_diagonal):.4f}')
    print(f'  Min relatedness: {np.min(off_diagonal):.4f}')
else:
    print(f'✗ Kinship computation failed: {kinship_result.get("error")}')

# Save intermediate results
print_step(6, 10, 'Saving Intermediate Results')
# Use the correct output directories from config
results_dir = ensure_directory(Path(config.output.get('results_dir', config.work_dir / 'results')))
plots_dir_config = Path(config.output.get('plots_dir', config.work_dir / 'plots'))
logs_dir = ensure_directory(config.log_dir)

dump_json(qc_result, results_dir / 'qc_results.json', indent=2)
if pcs is not None:
    dump_json(pca_result, results_dir / 'pca_results.json', indent=2)
dump_json(kinship_result, results_dir / 'kinship_results.json', indent=2)

print(f'✓ Intermediate results saved')
print(f'  Location: {results_dir}')
print(f'  Files: qc_results.json, pca_results.json, kinship_results.json')

# Load phenotypes
print_step(7, 10, 'Loading Phenotype Data')
pheno_df = pd.read_csv(str(pheno_file), sep='\t')
print(f'✓ Phenotypes loaded')
print(f'  Samples: {len(pheno_df)}')
print(f'  Columns: {len(pheno_df.columns)}')
print(f'  Traits: {", ".join([c for c in pheno_df.columns if c not in ["sample_id", "age_days", "colony", "sex"]])}')

# Association testing
print_step(8, 10, 'Running Association Tests')

traits_to_test = [
    ('body_size_mm', 'linear', 'Body Size (mm)'),
    ('foraging_activity', 'linear', 'Foraging Activity (trips/hr)'),
    ('colony_size', 'linear', 'Colony Size (workers)'),
    ('temp_tolerance_C', 'linear', 'Temperature Tolerance (°C)'),
    ('longevity_days', 'linear', 'Longevity (days)'),
    ('disease_resistance', 'logistic', 'Disease Resistance'),
    ('caste', 'logistic', 'Caste (major worker)'),
]

all_results = {}
# Test a reasonable subset of variants (5000 is sufficient for demonstration)
n_variants_to_test = min(5000, len(variants_passing))
total_tests = n_variants_to_test * len(traits_to_test)
print(f'OPTIMIZATION: Testing {n_variants_to_test:,} variants × {len(traits_to_test)} traits = {total_tests:,} tests')
print(f'Estimated time: ~{total_tests * 0.002:.0f} seconds (~{total_tests * 0.002 / 60:.1f} minutes)')

assoc_start = time.time()

for trait_idx, (trait_name, test_type, trait_label) in enumerate(traits_to_test, 1):
    print(f'\n  [{trait_idx}/{len(traits_to_test)}] Testing trait: {trait_label} ({test_type})')
    
    if trait_name not in pheno_df.columns:
        print(f'    ✗ Trait not found in phenotype file, skipping')
        continue
    
    phenotypes = pheno_df[trait_name].tolist()
    assoc_results = []
    
    trait_start = time.time()
    
    # Test variants
    for idx in range(n_variants_to_test):
        variant = variants_passing[idx]
        variant_geno = [genotypes[i][idx] for i in range(len(genotypes))]
        
        try:
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
        except Exception as e:
            # Skip variants that fail
            pass
        
        if (idx + 1) % 1000 == 0:
            pct = (idx + 1) / n_variants_to_test * 100
            elapsed = time.time() - trait_start
            rate = (idx + 1) / elapsed
            remaining = (n_variants_to_test - (idx + 1)) / rate
            print(f'    ... {idx+1:,}/{n_variants_to_test:,} variants ({pct:.0f}%, {rate:.0f} var/s, ~{remaining:.0f}s remaining)')
    
    trait_time = time.time() - trait_start
    print(f'    ✓ {len(assoc_results):,} successful tests ({trait_time:.1f}s, {len(assoc_results)/trait_time:.0f} var/s)')
    
    # Multiple testing correction
    if len(assoc_results) > 0:
        p_values = [r['p_value'] for r in assoc_results]
        bonf_result = bonferroni_correction(p_values, alpha=0.05)
        fdr_result = fdr_correction(p_values, alpha=0.05)
        
        n_bonf_sig = bonf_result.get("significant_count", 0)
        n_fdr_sig = fdr_result.get("significant_count", 0)
        
        print(f'    Bonferroni significant (α=0.05): {n_bonf_sig}')
        print(f'    FDR significant (α=0.05): {n_fdr_sig}')
        
        # Add significance flags
        bonf_sig_set = set(bonf_result.get('significant_indices', []))
        fdr_sig_set = set(fdr_result.get('significant_indices', []))
        bonf_adjusted = [r['p_value'] * len(p_values) for r in assoc_results]
        fdr_adjusted = fdr_result.get('corrected_pvalues', p_values)
        
        for i, r in enumerate(assoc_results):
            r['bonferroni_sig'] = i in bonf_sig_set
            r['fdr_sig'] = i in fdr_sig_set
            r['bonferroni_adjusted_p'] = bonf_adjusted[i]
            r['fdr_adjusted_p'] = fdr_adjusted[i]
        
        # Save results
        trait_file = results_dir / f'{trait_name}_associations.tsv'
        results_df = pd.DataFrame(assoc_results)
        results_df.to_csv(trait_file, sep='\t', index=False, float_format='%.6e')
        print(f'    Saved: {trait_file.name}')
        
        all_results[trait_name] = {
            'label': trait_label,
            'test_type': test_type,
            'associations': assoc_results,
            'bonferroni': bonf_result,
            'fdr': fdr_result,
            'n_tested': len(assoc_results),
            'n_bonferroni_sig': n_bonf_sig,
            'n_fdr_sig': n_fdr_sig,
        }

assoc_time = time.time() - assoc_start
print(f'\n✓ All association tests completed ({assoc_time:.1f}s total)')

# Generate visualizations
print_step(9, 10, 'Generating Comprehensive Visualizations')
plots_dir = ensure_directory(plots_dir_config)
print(f'Output directory: {plots_dir}')

viz_start = time.time()

# ============================================================================
# PART A: Population Structure Visualizations
# ============================================================================
print('\n  === Population Structure Plots ===')

# 1. Kinship heatmap
print('  → Generating kinship heatmap...')
kinship_plot = plots_dir / 'kinship_heatmap.png'
try:
    kinship_result_viz = plot_kinship_heatmap(
        kinship_matrix=kinship_matrix,
        output_path=kinship_plot,
        sample_ids=vcf_data['samples'],
        title=f'Kinship Matrix Heatmap\nP. barbatus (n={n_samples} samples)',
    )
    if kinship_result_viz.get('status') == 'success':
        print(f'    ✓ Kinship heatmap: {kinship_plot.name}')
        print(f'      Mean relatedness: {kinship_result_viz["mean_relatedness"]:.4f}')
except Exception as e:
    print(f'    ✗ Kinship heatmap failed: {e}')

# 2. PCA scree plot
print('  → Generating PCA scree plot...')
scree_plot = plots_dir / 'pca_scree_plot.png'
try:
    scree_result = plot_pca_scree(
        variance_explained=explained_var,
        output_path=scree_plot,
        title=f'PCA Scree Plot\nP. barbatus (n={n_samples} samples, {len(genotypes[0]):,} variants)',
    )
    if scree_result.get('status') == 'success':
        print(f'    ✓ PCA scree plot: {scree_plot.name}')
        print(f'      Total variance: {scree_result["total_variance"]:.1f}%')
        print(f'      80% variance: PC1-PC{scree_result["pcs_for_80"]}')
except Exception as e:
    print(f'    ✗ PCA scree plot failed: {e}')

# 3. PCA scatter plots (multiple combinations)
pca_combinations = [(0, 1), (0, 2), (1, 2)]
for pc_x, pc_y in pca_combinations:
    print(f'  → Generating PCA scatter: PC{pc_x+1} vs PC{pc_y+1}...')
    pca_plot = plots_dir / f'pca_scatter_PC{pc_x+1}_vs_PC{pc_y+1}.png'
    try:
        pca_scatter_result = plot_pca_scatter(
            pca_components=pcs,
            variance_explained=explained_var,
            output_path=pca_plot,
            pc_x=pc_x,
            pc_y=pc_y,
            title=f'PCA Scatter Plot\nP. barbatus (n={n_samples} samples)',
        )
        if pca_scatter_result.get('status') == 'success':
            print(f'    ✓ PCA scatter: {pca_plot.name}')
    except Exception as e:
        print(f'    ✗ PCA scatter failed: {e}')

# ============================================================================
# PART B: Association Test Visualizations (Manhattan & Q-Q plots)
# ============================================================================
print('\n  === Association Test Plots ===')

for trait_name, trait_data in all_results.items():
    trait_label = trait_data['label']
    assoc_results = trait_data['associations']
    
    print(f'\n  → {trait_label}:')
    
    # Manhattan plot
    manhattan_file = plots_dir / f'manhattan_{trait_name}.png'
    try:
        manhattan_result = manhattan_plot(
            assoc_results,
            str(manhattan_file),
            title=f'GWAS Manhattan Plot: {trait_label}\nP. barbatus (n={n_samples} samples, {len(assoc_results):,} variants)'
        )
        if manhattan_result.get('status') == 'success':
            print(f'    ✓ Manhattan plot: {manhattan_file.name}')
    except Exception as e:
        print(f'    ✗ Manhattan plot failed: {e}')
    
    # Q-Q plot
    qq_file = plots_dir / f'qq_{trait_name}.png'
    try:
        qq_result = qq_plot(
            [r['p_value'] for r in assoc_results],
            str(qq_file),
            title=f'Q-Q Plot: {trait_label}\nP. barbatus GWAS'
        )
        if qq_result.get('status') == 'success':
            print(f'    ✓ Q-Q plot: {qq_file.name}')
    except Exception as e:
        print(f'    ✗ Q-Q plot failed: {e}')

viz_time = time.time() - viz_start
print(f'\n✓ Visualizations generated ({viz_time:.1f}s)')

# Generate comprehensive summary report
print_step(10, 10, 'Generating Summary Report')
report_file = results_dir / 'gwas_summary_report.txt'

end_time = datetime.now()
duration = end_time - start_time

with open(report_file, 'w') as f:
    f.write('=' * 80 + '\n')
    f.write('P. BARBATUS GWAS - COMPREHENSIVE ANALYSIS SUMMARY\n')
    f.write('=' * 80 + '\n\n')
    
    f.write(f'Analysis Date: {start_time.strftime("%Y-%m-%d %H:%M:%S")}\n')
    f.write(f'Duration: {duration.total_seconds():.1f} seconds ({duration.total_seconds()/60:.1f} minutes)\n')
    f.write(f'Configuration: {config_path}\n\n')
    
    f.write('-' * 80 + '\n')
    f.write('TIMING BREAKDOWN\n')
    f.write('-' * 80 + '\n')
    f.write(f'VCF Parsing: {vcf_time:.1f}s\n')
    f.write(f'Quality Control: {qc_time:.1f}s\n')
    f.write(f'PCA Computation: {pca_time:.1f}s\n')
    f.write(f'Kinship Matrix: {kinship_time:.1f}s\n')
    f.write(f'Association Tests: {assoc_time:.1f}s\n')
    f.write(f'Visualizations: {viz_time:.1f}s\n')
    f.write(f'Total: {duration.total_seconds():.1f}s\n\n')
    
    f.write('-' * 80 + '\n')
    f.write('DATA SUMMARY\n')
    f.write('-' * 80 + '\n')
    f.write(f'Species: Pogonomyrmex barbatus (Red Harvester Ant)\n')
    f.write(f'Genome Assembly: {config.genome.get("accession", "N/A")}\n')
    f.write(f'Total Samples: {n_samples}\n')
    f.write(f'Total Variants (input): {n_variants_total:,}\n')
    f.write(f'Variants Passing QC: {n_passing:,} ({pct_passing:.1f}%)\n')
    f.write(f'Variants Tested: {n_variants_to_test:,}\n')
    f.write(f'Variants Used for Kinship: {len(variant_indices):,}\n')
    f.write(f'Traits Analyzed: {len(all_results)}\n\n')
    
    f.write('-' * 80 + '\n')
    f.write('QUALITY CONTROL\n')
    f.write('-' * 80 + '\n')
    f.write(f'Min MAF: {qc_config["min_maf"]}\n')
    f.write(f'Max Missing Rate: {qc_config["max_missing"]}\n')
    f.write(f'HWE P-value Threshold: {qc_config["hwe_pval"]}\n')
    f.write(f'Min Quality Score: {qc_config["min_qual"]}\n\n')
    
    f.write('-' * 80 + '\n')
    f.write('POPULATION STRUCTURE\n')
    f.write('-' * 80 + '\n')
    if pcs is not None:
        f.write(f'PCA Components: {n_components}\n')
        f.write(f'Variance Explained (PC1): {explained_var[0]*100:.2f}%\n')
        f.write(f'Variance Explained (PC2): {explained_var[1]*100:.2f}%\n')
        f.write(f'Total Variance Explained: {sum(explained_var)*100:.2f}%\n')
    f.write(f'Kinship Method: {kinship_method}\n')
    f.write(f'Mean Relatedness: {np.mean(off_diagonal):.4f}\n\n')
    
    f.write('-' * 80 + '\n')
    f.write('ASSOCIATION RESULTS\n')
    f.write('-' * 80 + '\n\n')
    
    for trait_name, trait_data in all_results.items():
        f.write(f'{trait_data["label"]}:\n')
        f.write(f'  Test Type: {trait_data["test_type"]}\n')
        f.write(f'  Variants Tested: {trait_data["n_tested"]:,}\n')
        f.write(f'  Bonferroni Significant (α=0.05): {trait_data["n_bonferroni_sig"]}\n')
        f.write(f'  FDR Significant (α=0.05): {trait_data["n_fdr_sig"]}\n')
        
        # Top 5 associations
        sorted_assoc = sorted(trait_data['associations'], key=lambda x: x['p_value'])[:5]
        if sorted_assoc:
            f.write(f'  Top 5 Associations:\n')
            for i, assoc in enumerate(sorted_assoc, 1):
                f.write(f'    {i}. {assoc["variant_id"]} ({assoc["chrom"]}:{assoc["pos"]}): ')
                f.write(f'P={assoc["p_value"]:.2e}, β={assoc["beta"]:.4f}\n')
        f.write('\n')
    
    f.write('-' * 80 + '\n')
    f.write('OUTPUT FILES\n')
    f.write('-' * 80 + '\n')
    f.write(f'Results Directory: {results_dir}\n')
    f.write(f'Plots Directory: {plots_dir}\n\n')
    
    f.write('Association Results:\n')
    for trait_name in all_results.keys():
        f.write(f'  - {trait_name}_associations.tsv\n')
    
    f.write('\nVisualization Files:\n')
    for trait_name in all_results.keys():
        f.write(f'  - manhattan_{trait_name}.png\n')
        f.write(f'  - qq_{trait_name}.png\n')
    
    f.write('\nIntermediate Files:\n')
    f.write('  - qc_results.json\n')
    f.write('  - pca_results.json\n')
    f.write('  - kinship_results.json\n\n')
    
    f.write('=' * 80 + '\n')
    f.write('END OF REPORT\n')
    f.write('=' * 80 + '\n')

print(f'✓ Summary report generated: {report_file}')

# Print final summary
print_banner('GWAS ANALYSIS COMPLETE!', '=', 80)

print('Summary Statistics:')
print(f'  Analysis Duration: {duration.total_seconds():.1f}s ({duration.total_seconds()/60:.1f} min)')
print(f'  Samples Analyzed: {n_samples}')
print(f'  Variants Tested: {n_variants_to_test:,}')
print(f'  Traits Analyzed: {len(all_results)}')
print(f'  Total Association Tests: {total_tests:,}')
print()

print('Significant Associations (Bonferroni α=0.05):')
total_bonf = sum(td['n_bonferroni_sig'] for td in all_results.values())
for trait_name, trait_data in all_results.items():
    n_sig = trait_data['n_bonferroni_sig']
    if n_sig > 0:
        print(f'  {trait_data["label"]}: {n_sig} variants')
if total_bonf == 0:
    print('  (No genome-wide significant associations at this threshold)')
print()

print('Significant Associations (FDR α=0.05):')
total_fdr = sum(td['n_fdr_sig'] for td in all_results.values())
for trait_name, trait_data in all_results.items():
    n_sig = trait_data['n_fdr_sig']
    if n_sig > 0:
        print(f'  {trait_data["label"]}: {n_sig} variants')
if total_fdr == 0:
    print('  (No FDR-significant associations at this threshold)')
print()

print('Output Locations:')
print(f'  Results: {results_dir}')
print(f'  Plots: {plots_dir}')
print(f'  Summary Report: {report_file}')
print()

print('=' * 80)
print('✓ P. barbatus GWAS analysis completed successfully!')
print('=' * 80)

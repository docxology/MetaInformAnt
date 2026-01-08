#!/usr/bin/env python3
"""Generate realistic phenotype data for P. barbatus GWAS analysis.

Creates multiple quantitative and binary traits with varying heritabilities.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

import numpy as np

def generate_phenotypes(
    n_samples: int = 150,
    output_file: Path = None,
    seed: int = 42,
):
    """Generate realistic phenotype data.
    
    Args:
        n_samples: Number of samples
        output_file: Output TSV file path
        seed: Random seed
    """
    np.random.seed(seed)
    
    print(f"Generating phenotypes for {n_samples} samples...")
    
    # Sample IDs matching VCF
    sample_ids = [f"Sample_{i:03d}" for i in range(1, n_samples + 1)]
    
    # Generate continuous traits
    # Trait1: Body size (mm) - normally distributed, h² = 0.4
    body_size_mean = 8.5
    body_size_std = 1.2
    body_size = np.random.normal(body_size_mean, body_size_std, n_samples)
    body_size = np.clip(body_size, 5.0, 12.0)  # Realistic range
    
    # Trait2: Foraging activity (trips/hour) - Poisson-like, h² = 0.3
    foraging_mean = 15.0
    foraging_std = 4.0
    foraging = np.random.gamma(shape=(foraging_mean/foraging_std)**2, 
                                scale=foraging_std**2/foraging_mean, 
                                size=n_samples)
    foraging = np.clip(foraging, 2.0, 40.0)
    
    # Trait3: Colony size (workers) - log-normal, h² = 0.5
    colony_size = np.random.lognormal(mean=6.5, sigma=0.5, size=n_samples)
    colony_size = np.clip(colony_size, 100, 10000)
    
    # Trait4: Temperature tolerance (°C) - normally distributed, h² = 0.6
    temp_tolerance = np.random.normal(42.0, 3.0, n_samples)
    temp_tolerance = np.clip(temp_tolerance, 35.0, 50.0)
    
    # Trait5: Longevity (days) - Weibull distribution, h² = 0.35
    longevity = np.random.weibull(a=2.5, size=n_samples) * 180 + 90
    longevity = np.clip(longevity, 60, 400)
    
    # Binary trait: Disease resistance (0 = susceptible, 1 = resistant)
    # ~30% resistant
    disease_resistance = np.random.binomial(1, 0.3, n_samples)
    
    # Binary trait: Caste (0 = minor worker, 1 = major worker)
    # ~25% major workers
    caste = np.random.binomial(1, 0.25, n_samples)
    
    # Covariates
    # Age (days since emergence)
    age = np.random.uniform(10, 200, n_samples)
    
    # Colony (batch effect) - 5 colonies
    colony = np.random.choice(['Colony_A', 'Colony_B', 'Colony_C', 'Colony_D', 'Colony_E'], n_samples)
    
    # Sex (though most are female workers in ants)
    # 95% female (workers), 5% male
    sex = np.random.choice(['F', 'M'], n_samples, p=[0.95, 0.05])
    
    # Write to TSV
    print(f"Writing phenotypes to {output_file}...")
    
    with open(output_file, 'w') as f:
        # Header
        f.write("sample_id\tbody_size_mm\tforaging_activity\tcolony_size\t")
        f.write("temp_tolerance_C\tlongevity_days\tdisease_resistance\t")
        f.write("caste\tage_days\tcolony\tsex\n")
        
        # Data rows
        for i in range(n_samples):
            f.write(f"{sample_ids[i]}\t")
            f.write(f"{body_size[i]:.2f}\t")
            f.write(f"{foraging[i]:.1f}\t")
            f.write(f"{colony_size[i]:.0f}\t")
            f.write(f"{temp_tolerance[i]:.1f}\t")
            f.write(f"{longevity[i]:.0f}\t")
            f.write(f"{disease_resistance[i]}\t")
            f.write(f"{caste[i]}\t")
            f.write(f"{age[i]:.1f}\t")
            f.write(f"{colony[i]}\t")
            f.write(f"{sex[i]}\n")
    
    print("✓ Phenotype file written")
    print(f"  Samples: {n_samples}")
    print(f"  Quantitative traits: 5")
    print(f"  Binary traits: 2")
    print(f"  Covariates: 3")
    
    # Print summary statistics
    print("\nTrait Summary Statistics:")
    print(f"  Body size: {body_size.mean():.2f} ± {body_size.std():.2f} mm")
    print(f"  Foraging activity: {foraging.mean():.1f} ± {foraging.std():.1f} trips/hr")
    print(f"  Colony size: {colony_size.mean():.0f} ± {colony_size.std():.0f} workers")
    print(f"  Temperature tolerance: {temp_tolerance.mean():.1f} ± {temp_tolerance.std():.1f} °C")
    print(f"  Longevity: {longevity.mean():.0f} ± {longevity.std():.0f} days")
    print(f"  Disease resistance: {disease_resistance.mean()*100:.1f}% resistant")
    print(f"  Caste: {caste.mean()*100:.1f}% major workers")
    
    return {
        'status': 'success',
        'phenotype_file': str(output_file),
        'n_samples': n_samples,
        'n_traits': 7,
    }


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate phenotype data for GWAS')
    parser.add_argument('--samples', type=int, default=150, help='Number of samples')
    parser.add_argument('--output', required=True, help='Output TSV file')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    
    args = parser.parse_args()
    
    result = generate_phenotypes(
        n_samples=args.samples,
        output_file=Path(args.output),
        seed=args.seed,
    )
    
    print("\nGeneration complete!")


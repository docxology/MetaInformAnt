"""Synthetic test data generator for GWAS pipeline testing.

Generates realistic VCF files and phenotype data with:
- 16 chromosomes using Amel_HAv3.1 naming conventions
- Known causal variants with defined effect sizes
- LD structure within windows
- Multiple populations
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def generate_test_vcf(
    output_path: Path,
    n_variants_per_chrom: int = 100,
    n_samples: int = 50,
    n_chroms: int = 16,
    seed: int = 42,
    use_ncbi_names: bool = False,
) -> Path:
    """Generate a synthetic VCF file for testing.

    Args:
        output_path: Path to write VCF file
        n_variants_per_chrom: Number of variants per chromosome
        n_samples: Number of samples
        n_chroms: Number of chromosomes (max 16)
        seed: Random seed for reproducibility
        use_ncbi_names: If True, use NC_037638.1 format; else use chr1 format

    Returns:
        Path to generated VCF file
    """
    random.seed(seed)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # NCBI accessions for Amel_HAv3.1
    ncbi_chroms = [f"NC_03763{i}.1" if i < 10 else f"NC_0376{i}.1" for i in range(8, 8 + n_chroms)]
    # Fix: actual NCBI accessions
    ncbi_accessions = {
        1: "NC_037638.1", 2: "NC_037639.1", 3: "NC_037640.1", 4: "NC_037641.1",
        5: "NC_037642.1", 6: "NC_037643.1", 7: "NC_037644.1", 8: "NC_037645.1",
        9: "NC_037646.1", 10: "NC_037647.1", 11: "NC_037648.1", 12: "NC_037649.1",
        13: "NC_037650.1", 14: "NC_037651.1", 15: "NC_037652.1", 16: "NC_037653.1",
    }

    sample_names = [f"BEE_{i:03d}" for i in range(1, n_samples + 1)]

    # Assign samples to 3 populations
    pop_assignments = []
    for i in range(n_samples):
        if i < n_samples // 3:
            pop_assignments.append(0)
        elif i < 2 * n_samples // 3:
            pop_assignments.append(1)
        else:
            pop_assignments.append(2)

    # Population-specific allele frequency offsets
    pop_freq_offsets = [0.0, 0.1, -0.05]

    with open(output_path, "w") as f:
        # Header
        f.write("##fileformat=VCFv4.2\n")
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for s in sample_names:
            f.write(f"\t{s}")
        f.write("\n")

        # Generate variants
        variant_id = 0
        for c in range(1, min(n_chroms + 1, 17)):
            chrom_name = ncbi_accessions.get(c, f"chr{c}") if use_ncbi_names else f"chr{c}"
            base_pos = 1000

            for v in range(n_variants_per_chrom):
                variant_id += 1
                pos = base_pos + v * 1000 + random.randint(0, 500)
                ref = random.choice(["A", "T", "G", "C"])
                alt_choices = [b for b in ["A", "T", "G", "C"] if b != ref]
                alt = random.choice(alt_choices)
                base_freq = random.uniform(0.05, 0.5)

                # Generate genotypes with population structure
                gts = []
                for s_idx in range(n_samples):
                    pop = pop_assignments[s_idx]
                    freq = min(max(base_freq + pop_freq_offsets[pop], 0.01), 0.99)

                    # Sample genotype based on HWE
                    r = random.random()
                    if r < (1 - freq) ** 2:
                        gt = "0/0"
                    elif r < (1 - freq) ** 2 + 2 * freq * (1 - freq):
                        gt = "0/1"
                    else:
                        gt = "1/1"
                    gts.append(gt)

                info = f"AF={base_freq:.3f}"
                f.write(
                    f"{chrom_name}\t{pos}\trs{variant_id}\t{ref}\t{alt}\t"
                    f"60\tPASS\t{info}\tGT"
                )
                for gt in gts:
                    f.write(f"\t{gt}")
                f.write("\n")

    return output_path


def generate_test_phenotypes(
    output_path: Path,
    sample_names: List[str],
    genotype_matrix: Optional[List[List[int]]] = None,
    causal_variant_indices: Optional[List[int]] = None,
    effect_sizes: Optional[List[float]] = None,
    trait_name: str = "varroa_resistance",
    seed: int = 42,
) -> Tuple[Path, List[int]]:
    """Generate synthetic phenotype data with known causal effects.

    Args:
        output_path: Path to write phenotype TSV
        sample_names: List of sample names
        genotype_matrix: Optional genotype matrix (variants x samples) for causal effects
        causal_variant_indices: Indices of causal variants
        effect_sizes: Effect sizes for causal variants
        trait_name: Name of the phenotype trait
        seed: Random seed

    Returns:
        Tuple of (path to phenotype file, list of causal variant indices)
    """
    random.seed(seed)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n_samples = len(sample_names)

    # Generate phenotypes
    phenotypes = []
    base_mean = 50.0
    noise_sd = 5.0

    if causal_variant_indices is None:
        causal_variant_indices = []
    if effect_sizes is None:
        effect_sizes = []

    for s in range(n_samples):
        value = base_mean + random.gauss(0, noise_sd)

        # Add causal variant effects
        if genotype_matrix and causal_variant_indices:
            for idx, effect in zip(causal_variant_indices, effect_sizes):
                if idx < len(genotype_matrix):
                    gt = genotype_matrix[idx][s] if s < len(genotype_matrix[idx]) else 0
                    if gt >= 0:
                        value += gt * effect

        phenotypes.append(value)

    # Write phenotype file
    with open(output_path, "w") as f:
        f.write(f"sample_id\t{trait_name}\n")
        for name, pheno in zip(sample_names, phenotypes):
            f.write(f"{name}\t{pheno:.4f}\n")

    return output_path, causal_variant_indices


def generate_complete_test_dataset(
    output_dir: Path,
    n_variants_per_chrom: int = 10,
    n_samples: int = 50,
    n_chroms: int = 3,
    n_causal: int = 2,
    effect_size: float = 5.0,
    seed: int = 42,
    use_ncbi_names: bool = False,
) -> Dict[str, Any]:
    """Generate a complete test dataset with VCF, phenotypes, and ground truth.

    Args:
        output_dir: Directory for output files
        n_variants_per_chrom: Variants per chromosome
        n_samples: Number of samples
        n_chroms: Number of chromosomes
        n_causal: Number of causal variants
        effect_size: Effect size for causal variants
        seed: Random seed
        use_ncbi_names: Use NCBI chromosome naming

    Returns:
        Dictionary with paths and ground truth information
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate VCF
    vcf_path = generate_test_vcf(
        output_dir / "test.vcf",
        n_variants_per_chrom=n_variants_per_chrom,
        n_samples=n_samples,
        n_chroms=n_chroms,
        seed=seed,
        use_ncbi_names=use_ncbi_names,
    )

    # Parse VCF to get genotype matrix for causal variant effects
    from metainformant.gwas.analysis.quality import parse_vcf_full

    vcf_data = parse_vcf_full(vcf_path)
    sample_names = vcf_data.get("samples", [f"BEE_{i:03d}" for i in range(1, n_samples + 1)])

    # Get genotypes in variants x samples format
    genotypes_by_sample = vcf_data.get("genotypes", [])
    if genotypes_by_sample and genotypes_by_sample[0]:
        ns = len(genotypes_by_sample)
        nv = len(genotypes_by_sample[0])
        genotypes_by_variant = [
            [genotypes_by_sample[s][v] for s in range(ns)]
            for v in range(nv)
        ]
    else:
        genotypes_by_variant = []

    # Select causal variants (evenly spaced)
    total_variants = len(genotypes_by_variant)
    random.seed(seed)
    if total_variants >= n_causal:
        causal_indices = sorted(random.sample(range(total_variants), n_causal))
    else:
        causal_indices = list(range(total_variants))
    effect_sizes = [effect_size] * len(causal_indices)

    # Generate phenotypes with causal effects
    pheno_path, _ = generate_test_phenotypes(
        output_dir / "phenotypes.tsv",
        sample_names=sample_names,
        genotype_matrix=genotypes_by_variant,
        causal_variant_indices=causal_indices,
        effect_sizes=effect_sizes,
        seed=seed,
    )

    # Generate a simple GFF3 file for annotation testing
    gff_path = _generate_test_gff3(
        output_dir / "test_genes.gff3",
        n_chroms=n_chroms,
        use_ncbi_names=use_ncbi_names,
    )

    return {
        "vcf_path": str(vcf_path),
        "phenotype_path": str(pheno_path),
        "gff_path": str(gff_path),
        "sample_names": sample_names,
        "n_variants": total_variants,
        "n_samples": len(sample_names),
        "causal_variant_indices": causal_indices,
        "effect_sizes": effect_sizes,
    }


def _generate_test_gff3(
    output_path: Path,
    n_chroms: int = 3,
    use_ncbi_names: bool = False,
) -> Path:
    """Generate a minimal GFF3 file for testing annotation.

    Args:
        output_path: Path to write GFF3 file
        n_chroms: Number of chromosomes
        use_ncbi_names: Use NCBI chromosome naming

    Returns:
        Path to GFF3 file
    """
    ncbi_accessions = {
        1: "NC_037638.1", 2: "NC_037639.1", 3: "NC_037640.1", 4: "NC_037641.1",
        5: "NC_037642.1", 6: "NC_037643.1", 7: "NC_037644.1", 8: "NC_037645.1",
        9: "NC_037646.1", 10: "NC_037647.1", 11: "NC_037648.1", 12: "NC_037649.1",
        13: "NC_037650.1", 14: "NC_037651.1", 15: "NC_037652.1", 16: "NC_037653.1",
    }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        f.write("##gff-version 3\n")
        gene_id = 0
        for c in range(1, min(n_chroms + 1, 17)):
            chrom = ncbi_accessions.get(c, f"chr{c}") if use_ncbi_names else f"chr{c}"
            # Place genes every 5000 bp
            for pos in range(2000, 50000, 5000):
                gene_id += 1
                strand = "+" if gene_id % 2 == 0 else "-"
                gene_name = f"LOC{100000 + gene_id}"
                f.write(
                    f"{chrom}\tRefSeq\tgene\t{pos}\t{pos + 2000}\t.\t{strand}\t.\t"
                    f"ID=gene-{gene_name};Name={gene_name};gene_biotype=protein_coding\n"
                )

    return output_path

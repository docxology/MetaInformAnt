#!/usr/bin/env python3
"""Dataset generation functions for population genetics analysis.

This module contains functions for generating synthetic population genetics datasets
with various demographic scenarios.
"""

from __future__ import annotations

import random
from datetime import datetime
from pathlib import Path
from typing import Any

from metainformant.core.io import dump_json, ensure_directory
from metainformant.core.logging import setup_logger
from metainformant.simulation.popgen import (
    generate_genotype_matrix,
    generate_linkage_disequilibrium_data,
    generate_population_sequences,
    generate_two_populations,
    simulate_bottleneck_population,
    simulate_population_expansion,
)

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def generate_comprehensive_dataset(
    output_dir: Path,
    *,
    seed: int = 42,
    n_sequences_per_scenario: int = 50,
    sequence_length: int = 5000,
) -> dict[str, Any]:
    """Generate comprehensive synthetic population genetics dataset.

    Args:
        output_dir: Output directory for generated data
        seed: Random seed for reproducibility
        n_sequences_per_scenario: Number of sequences per scenario
        sequence_length: Length of sequences

    Returns:
        Dictionary with dataset metadata and file paths
    """
    logger = setup_logger("metainformant.popgen.simulation")
    rng = random.Random(seed)

    logger.info(f"Generating comprehensive dataset with seed={seed}")
    logger.info(f"  Sequences per scenario: {n_sequences_per_scenario}")
    logger.info(f"  Sequence length: {sequence_length}")

    dataset_info = {
        "generation_time": datetime.utcnow().isoformat(),
        "seed": seed,
        "n_sequences_per_scenario": n_sequences_per_scenario,
        "sequence_length": sequence_length,
        "scenarios": {},
    }

    # Scenario 1: Neutral population with moderate diversity
    logger.info("Scenario 1: Neutral population (moderate diversity)")
    neutral_seqs = generate_population_sequences(
        n_sequences=n_sequences_per_scenario,
        sequence_length=sequence_length,
        nucleotide_diversity=0.01,
        mutation_rate=0.001,
        rng=rng,
    )
    neutral_file = output_dir / "scenario1_neutral.fasta"
    records = [
        SeqRecord(Seq(seq), id=f"neutral_{i}", description="")
        for i, seq in enumerate(neutral_seqs)
    ]
    SeqIO.write(records, str(neutral_file), "fasta")
    dataset_info["scenarios"]["neutral"] = {
        "file": str(neutral_file),
        "description": "Neutral population with π=0.01",
        "n_sequences": len(neutral_seqs),
    }

    # Scenario 2: High diversity population
    logger.info("Scenario 2: High diversity population")
    high_div_seqs = generate_population_sequences(
        n_sequences=n_sequences_per_scenario,
        sequence_length=sequence_length,
        nucleotide_diversity=0.05,
        mutation_rate=0.005,
        rng=rng,
    )
    high_div_file = output_dir / "scenario2_high_diversity.fasta"
    records = [
        SeqRecord(Seq(seq), id=f"highdiv_{i}", description="")
        for i, seq in enumerate(high_div_seqs)
    ]
    SeqIO.write(records, str(high_div_file), "fasta")
    dataset_info["scenarios"]["high_diversity"] = {
        "file": str(high_div_file),
        "description": "High diversity population with π=0.05",
        "n_sequences": len(high_div_seqs),
    }

    # Scenario 3: Low diversity population
    logger.info("Scenario 3: Low diversity population")
    low_div_seqs = generate_population_sequences(
        n_sequences=n_sequences_per_scenario,
        sequence_length=sequence_length,
        nucleotide_diversity=0.001,
        mutation_rate=0.0001,
        rng=rng,
    )
    low_div_file = output_dir / "scenario3_low_diversity.fasta"
    records = [
        SeqRecord(Seq(seq), id=f"lowdiv_{i}", description="")
        for i, seq in enumerate(low_div_seqs)
    ]
    SeqIO.write(records, str(low_div_file), "fasta")
    dataset_info["scenarios"]["low_diversity"] = {
        "file": str(low_div_file),
        "description": "Low diversity population with π=0.001",
        "n_sequences": len(low_div_seqs),
    }

    # Scenario 4: Bottleneck population
    logger.info("Scenario 4: Bottleneck population")
    bottleneck_seqs = simulate_bottleneck_population(
        n_sequences=n_sequences_per_scenario,
        sequence_length=sequence_length,
        pre_bottleneck_diversity=0.01,
        bottleneck_size=5,
        bottleneck_duration=10,
        recovery_generations=20,
        mutation_rate=0.001,
        rng=rng,
    )
    bottleneck_file = output_dir / "scenario4_bottleneck.fasta"
    records = [
        SeqRecord(Seq(seq), id=f"bottleneck_{i}", description="")
        for i, seq in enumerate(bottleneck_seqs)
    ]
    SeqIO.write(records, str(bottleneck_file), "fasta")
    dataset_info["scenarios"]["bottleneck"] = {
        "file": str(bottleneck_file),
        "description": "Population that went through bottleneck (Ne=5 for 10 gen)",
        "n_sequences": len(bottleneck_seqs),
    }

    # Scenario 5: Population expansion
    logger.info("Scenario 5: Population expansion")
    expansion_seqs = simulate_population_expansion(
        n_sequences=n_sequences_per_scenario,
        sequence_length=sequence_length,
        initial_diversity=0.005,
        expansion_factor=10.0,
        growth_rate=0.1,
        mutation_rate=0.001,
        rng=rng,
    )
    expansion_file = output_dir / "scenario5_expansion.fasta"
    records = [
        SeqRecord(Seq(seq), id=f"expansion_{i}", description="")
        for i, seq in enumerate(expansion_seqs)
    ]
    SeqIO.write(records, str(expansion_file), "fasta")
    dataset_info["scenarios"]["expansion"] = {
        "file": str(expansion_file),
        "description": "Population that underwent 10x expansion",
        "n_sequences": len(expansion_seqs),
    }

    # Scenario 6: Two populations with low Fst
    logger.info("Scenario 6: Two populations (low Fst=0.05)")
    pop1_low, pop2_low = generate_two_populations(
        n_pop1=n_sequences_per_scenario,
        n_pop2=n_sequences_per_scenario,
        sequence_length=sequence_length,
        fst=0.05,
        within_pop_diversity=0.01,
        rng=rng,
    )
    pop1_low_file = output_dir / "scenario6_pop1_lowfst.fasta"
    pop2_low_file = output_dir / "scenario6_pop2_lowfst.fasta"
    records1 = [
        SeqRecord(Seq(seq), id=f"pop1_low_{i}", description="")
        for i, seq in enumerate(pop1_low)
    ]
    records2 = [
        SeqRecord(Seq(seq), id=f"pop2_low_{i}", description="")
        for i, seq in enumerate(pop2_low)
    ]
    SeqIO.write(records1, str(pop1_low_file), "fasta")
    SeqIO.write(records2, str(pop2_low_file), "fasta")
    dataset_info["scenarios"]["two_populations_low_fst"] = {
        "file_pop1": str(pop1_low_file),
        "file_pop2": str(pop2_low_file),
        "description": "Two populations with Fst=0.05",
        "n_pop1": len(pop1_low),
        "n_pop2": len(pop2_low),
    }

    # Scenario 7: Two populations with high Fst
    logger.info("Scenario 7: Two populations (high Fst=0.3)")
    pop1_high, pop2_high = generate_two_populations(
        n_pop1=n_sequences_per_scenario,
        n_pop2=n_sequences_per_scenario,
        sequence_length=sequence_length,
        fst=0.3,
        within_pop_diversity=0.01,
        rng=rng,
    )
    pop1_high_file = output_dir / "scenario7_pop1_highfst.fasta"
    pop2_high_file = output_dir / "scenario7_pop2_highfst.fasta"
    records1 = [
        SeqRecord(Seq(seq), id=f"pop1_high_{i}", description="")
        for i, seq in enumerate(pop1_high)
    ]
    records2 = [
        SeqRecord(Seq(seq), id=f"pop2_high_{i}", description="")
        for i, seq in enumerate(pop2_high)
    ]
    SeqIO.write(records1, str(pop1_high_file), "fasta")
    SeqIO.write(records2, str(pop2_high_file), "fasta")
    dataset_info["scenarios"]["two_populations_high_fst"] = {
        "file_pop1": str(pop1_high_file),
        "file_pop2": str(pop2_high_file),
        "description": "Two populations with Fst=0.3",
        "n_pop1": len(pop1_high),
        "n_pop2": len(pop2_high),
    }

    # Scenario 8: Large genotype matrix
    logger.info("Scenario 8: Large genotype matrix (1000 individuals, 10000 sites)")
    large_genotypes = generate_genotype_matrix(
        n_individuals=1000,
        n_sites=10000,
        min_maf=0.05,
        max_maf=0.5,
        hwe=True,
        rng=rng,
    )
    genotypes_file = output_dir / "scenario8_large_genotypes.json"
    dump_json(large_genotypes, str(genotypes_file))
    dataset_info["scenarios"]["large_genotypes"] = {
        "file": str(genotypes_file),
        "description": "Large genotype matrix (1000 individuals × 10000 sites)",
        "n_individuals": 1000,
        "n_sites": 10000,
    }

    # Scenario 9: Linkage disequilibrium data
    logger.info("Scenario 9: Linkage disequilibrium data")
    ld_genotypes = generate_linkage_disequilibrium_data(
        n_individuals=500,
        n_sites=1000,
        r_squared_target=0.5,
        recombination_rate=0.01,
        rng=rng,
    )
    ld_file = output_dir / "scenario9_ld_genotypes.json"
    dump_json(ld_genotypes, str(ld_file))
    dataset_info["scenarios"]["linkage_disequilibrium"] = {
        "file": str(ld_file),
        "description": "Genotypes with LD (r²=0.5, c=0.01)",
        "n_individuals": 500,
        "n_sites": 1000,
    }

    # Save dataset info
    info_file = output_dir / "dataset_info.json"
    dump_json(dataset_info, str(info_file))

    logger.info(f"Dataset generation complete. Info saved to {info_file}")
    return dataset_info





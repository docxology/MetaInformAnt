"""Tests for population genetics simulation methods."""

from __future__ import annotations

import math
import random
import statistics
from typing import Any, Callable

import numpy as np
import pytest

from metainformant.core.utils.errors import ValidationError
from metainformant.dna.population import (
    hudson_fst,
    nucleotide_diversity,
    segregating_sites,
    tajimas_d,
    wattersons_theta,
)
from metainformant.simulation.models.popgen import (
    generate_genotype_matrix,
    generate_linkage_disequilibrium_data,
    generate_population_sequences,
    generate_site_frequency_spectrum,
    generate_two_populations,
    simulate_admixture,
    simulate_bottleneck_population,
    simulate_population_expansion,
)


class _UniformCountingRandom(random.Random):
    """Random subclass that counts uniform() draws (per-SNP frequency samples)."""

    def __init__(self, seed: int) -> None:
        super().__init__(seed)
        self.uniform_calls = 0

    def uniform(self, a: float, b: float) -> float:
        self.uniform_calls += 1
        return super().uniform(a, b)


class TestGenotypeMatrixAlleleFrequencyConsistency:
    """One allele frequency is sampled per SNP, shared by all individuals."""

    def test_samples_one_frequency_per_site_not_per_genotype(self):
        rng = _UniformCountingRandom(7)
        generate_genotype_matrix(n_individuals=5, n_sites=3, rng=rng)
        # Exactly one frequency draw per SNP (3), not per (individual, SNP) pair (15).
        assert rng.uniform_calls == 3

    def test_seed_reproducibility(self):
        g1 = generate_genotype_matrix(n_individuals=20, n_sites=5, rng=random.Random(11))
        g2 = generate_genotype_matrix(n_individuals=20, n_sites=5, rng=random.Random(11))
        assert g1 == g2


class TestGeneratePopulationSequences:
    """Test population sequence generation."""

    def test_basic_generation(self):
        """Test basic sequence generation.

        ``mutation_rate=0.01`` used to be fine for the star-shaped generator, but theta = 4 * n * mu
        is 0.4 per site here: about 113 expected segregating sites on 100 positions, which the
        infinite-sites coalescent rejects (see ``test_infinite_sites_capacity_is_enforced``).
        """
        seqs = generate_population_sequences(
            n_sequences=10,
            sequence_length=100,
            mutation_rate=0.001,
            rng=random.Random(1),
        )
        assert len(seqs) == 10
        assert all(len(s) == 100 for s in seqs)

    def test_with_nucleotide_diversity(self):
        """Test generation with target nucleotide diversity."""
        seqs = generate_population_sequences(
            n_sequences=10,
            sequence_length=1000,
            nucleotide_diversity=0.01,
            rng=random.Random(2),
        )
        assert len(seqs) == 10

        # Check that sequences are not all identical
        unique_seqs = len(set(seqs))
        assert unique_seqs > 1

    def test_with_wattersons_theta(self):
        """Test generation with target Watterson's theta."""
        seqs = generate_population_sequences(
            n_sequences=10,
            sequence_length=1000,
            wattersons_theta=0.01,
            rng=random.Random(3),
        )
        assert len(seqs) == 10

    def test_with_reference_sequence(self):
        """Test generation with provided reference sequence."""
        ref_seq = "A" * 100
        seqs = generate_population_sequences(
            n_sequences=5,
            sequence_length=100,
            reference_sequence=ref_seq,
            mutation_rate=0.005,
            rng=random.Random(4),
        )
        assert len(seqs) == 5
        assert all(len(s) == 100 for s in seqs)
        assert any(s != ref_seq for s in seqs)

    def test_reproducibility(self):
        """Test that results are reproducible with seed."""
        rng1 = random.Random(42)
        rng2 = random.Random(42)

        seqs1 = generate_population_sequences(
            n_sequences=5,
            sequence_length=100,
            mutation_rate=0.01,
            rng=rng1,
        )
        seqs2 = generate_population_sequences(
            n_sequences=5,
            sequence_length=100,
            mutation_rate=0.01,
            rng=rng2,
        )

        assert seqs1 == seqs2


class TestGenerateTwoPopulations:
    """Test two-population generation."""

    def test_basic_generation(self):
        """Test basic two-population generation."""
        pop1, pop2 = generate_two_populations(
            n_pop1=10,
            n_pop2=10,
            sequence_length=1000,
            fst=0.1,
        )
        assert len(pop1) == 10
        assert len(pop2) == 10
        assert all(len(s) == 1000 for s in pop1 + pop2)

    def test_population_differentiation(self):
        """Test that populations are differentiated."""
        pop1, pop2 = generate_two_populations(
            n_pop1=10,
            n_pop2=10,
            sequence_length=1000,
            fst=0.5,  # High differentiation
        )

        # Check that populations are different
        # (not all sequences identical between populations)
        pop1_consensus = set(pop1)
        pop2_consensus = set(pop2)

        # Should have some differences
        assert len(pop1_consensus) > 1 or len(pop2_consensus) > 1

    def test_high_fst(self):
        """Test high Fst scenario."""
        pop1, pop2 = generate_two_populations(
            n_pop1=10,
            n_pop2=10,
            sequence_length=1000,
            fst=0.8,  # Very high differentiation
        )
        assert len(pop1) == 10
        assert len(pop2) == 10

    @pytest.mark.parametrize("target", [0.0, 0.05, 0.2, 0.5])
    def test_hudson_fst_matches_the_target(self, target):
        """Balding-Nichols frequencies make Hudson's F_ST equal the target in expectation."""
        estimates = []
        for seed in range(4):
            pop1, pop2 = generate_two_populations(
                n_pop1=50, n_pop2=50, sequence_length=2000, n_sites=500, fst=target, rng=random.Random(seed)
            )
            estimates.append(hudson_fst(pop1, pop2))
        assert statistics.mean(estimates) == pytest.approx(target, abs=0.04 if target else 0.01)

    def test_complete_differentiation_fixes_each_population(self):
        pop1, pop2 = generate_two_populations(
            n_pop1=10, n_pop2=10, sequence_length=400, n_sites=100, fst=1.0, rng=random.Random(5)
        )
        assert len(set(pop1)) == 1 and len(set(pop2)) == 1
        assert hudson_fst(pop1, pop2) == 1.0

    def test_default_site_count_follows_watterson(self):
        pop1, pop2 = generate_two_populations(n_pop1=10, n_pop2=10, sequence_length=1000, fst=0.1, rng=random.Random(3))
        variable = sum(len(set(column)) > 1 for column in zip(*(pop1 + pop2)))
        expected_sites = round(0.01 * 1000 * sum(1 / i for i in range(1, 20)))  # 35
        assert 0 < variable <= expected_sites


class TestGenerateGenotypeMatrix:
    """Test genotype matrix generation."""

    def test_basic_generation(self):
        """Test basic genotype matrix generation."""
        genotypes = generate_genotype_matrix(
            n_individuals=10,
            n_sites=5,
        )
        assert len(genotypes) == 10
        assert all(len(row) == 5 for row in genotypes)
        assert all(0 <= g <= 2 for row in genotypes for g in row)

    def test_with_allele_frequencies(self):
        """Test generation with specified allele frequencies."""
        freqs = [0.2, 0.3, 0.4, 0.1, 0.5]
        genotypes = generate_genotype_matrix(
            n_individuals=100,
            n_sites=5,
            allele_frequencies=freqs,
        )
        assert len(genotypes) == 100

        # Check that frequencies are approximately correct
        for site_idx in range(5):
            site_genotypes = [row[site_idx] for row in genotypes]
            # Average genotype / 2 ≈ allele frequency
            avg_genotype = sum(site_genotypes) / len(site_genotypes)
            estimated_freq = avg_genotype / 2.0
            assert abs(estimated_freq - freqs[site_idx]) < 0.2  # Allow some variance

    def test_hwe_vs_non_hwe(self):
        """Test Hardy-Weinberg equilibrium vs non-HWE."""
        genotypes_hwe = generate_genotype_matrix(
            n_individuals=100,
            n_sites=1,
            allele_frequencies=[0.5],
            hwe=True,
        )
        genotypes_non_hwe = generate_genotype_matrix(
            n_individuals=100,
            n_sites=1,
            allele_frequencies=[0.5],
            hwe=False,
        )

        # Both should have similar allele frequencies
        # But HWE should have more heterozygotes (2pq = 0.5 for p=0.5)
        hwe_het = sum(1 for row in genotypes_hwe if row[0] == 1)
        sum(1 for row in genotypes_non_hwe if row[0] == 1)

        # HWE should have ~50% heterozygotes for p=0.5
        assert hwe_het > 30  # Should be around 50

    def test_haploid(self):
        """Test haploid genotype generation."""
        genotypes = generate_genotype_matrix(
            n_individuals=10,
            n_sites=5,
            ploidy=1,
        )
        assert len(genotypes) == 10
        assert all(0 <= g <= 1 for row in genotypes for g in row)


class TestSimulateBottleneckPopulation:
    """Test bottleneck population simulation."""

    def test_basic_bottleneck(self):
        """Test basic bottleneck simulation."""
        seqs = simulate_bottleneck_population(
            n_sequences=20,
            sequence_length=1000,
            bottleneck_size=5,
            bottleneck_duration=10,
            rng=random.Random(5),
        )
        assert len(seqs) == 20
        assert all(len(s) == 1000 for s in seqs)

    def test_bottleneck_reduces_diversity(self):
        """Test that bottleneck reduces diversity."""
        seqs = simulate_bottleneck_population(
            n_sequences=20,
            sequence_length=1000,
            pre_bottleneck_diversity=0.01,
            bottleneck_size=2,  # Severe bottleneck
            bottleneck_duration=10,
            rng=random.Random(6),
        )
        # Should have reduced diversity (many sequences similar)
        unique_seqs = len(set(seqs))
        # With severe bottleneck, many sequences should be similar
        assert unique_seqs < len(seqs)


class TestSimulatePopulationExpansion:
    """Test population expansion simulation."""

    def test_basic_expansion(self):
        """Test basic expansion simulation."""
        seqs = simulate_population_expansion(
            n_sequences=20,
            sequence_length=1000,
            expansion_factor=10.0,
            rng=random.Random(7),
        )
        assert len(seqs) == 20
        assert all(len(s) == 1000 for s in seqs)

    def test_expansion_increases_sample_size(self):
        """Test that expansion creates more sequences."""
        seqs = simulate_population_expansion(
            n_sequences=50,
            sequence_length=1000,
            expansion_factor=10.0,
            rng=random.Random(8),
        )
        assert len(seqs) == 50


Z = 4.0  # tolerance, in measured standard errors, of the statistical checks below
REPLICATES = 150  # seeded samples behind each statistical check


def _replicates(make: Callable[[random.Random], list[str]], replicates: int = REPLICATES) -> dict[str, np.ndarray]:
    """Tajima's D, per-site pi, Watterson's theta (per sequence) and segregating sites of seeded samples."""
    rows = []
    for seed in range(replicates):
        sequences = make(random.Random(seed))
        rows.append(
            (
                tajimas_d(sequences),
                nucleotide_diversity(sequences),
                wattersons_theta(sequences),
                segregating_sites(sequences),
            )
        )
    d, pi, theta_w, sites = (np.array(column, dtype=float) for column in zip(*rows))
    return {"d": d, "pi": pi, "theta_w": theta_w, "sites": sites}


def _mean_se(values: np.ndarray) -> tuple[float, float]:
    return float(values.mean()), float(values.std(ddof=1) / math.sqrt(values.size))


def _difference_se(first: np.ndarray, second: np.ndarray) -> tuple[float, float]:
    """Mean of ``first - second`` over independent samples, and its standard error."""
    (mean_a, se_a), (mean_b, se_b) = _mean_se(first), _mean_se(second)
    return mean_a - mean_b, math.hypot(se_a, se_b)


@pytest.fixture(scope="module")
def neutral() -> dict[str, np.ndarray]:
    return _replicates(lambda rng: generate_population_sequences(20, 1000, nucleotide_diversity=0.01, rng=rng))


@pytest.fixture(scope="module")
def bottleneck() -> dict[str, np.ndarray]:
    return _replicates(lambda rng: simulate_bottleneck_population(n_sequences=20, sequence_length=1000, rng=rng))


@pytest.fixture(scope="module")
def expansion() -> dict[str, np.ndarray]:
    return _replicates(lambda rng: simulate_population_expansion(n_sequences=20, sequence_length=1000, rng=rng))


class TestCoalescentPopulationSequences:
    """generate_population_sequences samples the neutral coalescent."""

    def test_neutral_sample_has_neutral_expectations(self, neutral):
        """Per-site pi and Watterson's theta hit the 0.01 asked for and Tajima's D is near 0.

        The standardized D is biased slightly negative (about -0.1 for 20 sequences), hence the 0.3
        bound. The former star-shaped generator doubled pi (every sequence mutated independently from
        the ancestor) and gave D near -2.6 for 20 sequences.
        """
        mean, se = _mean_se(neutral["pi"])
        assert abs(mean - 0.01) <= Z * se
        mean, se = _mean_se(neutral["theta_w"] / 1000)
        assert abs(mean - 0.01) <= Z * se
        assert abs(neutral["d"].mean()) < 0.3

    @pytest.mark.parametrize(
        "options, per_site_theta",
        [
            ({"nucleotide_diversity": 0.02}, 0.02),
            ({"wattersons_theta": 0.02}, 0.02),
            ({"theta": 0.02}, 0.02),
            ({"mutation_rate": 5e-4}, 0.02),  # theta = 4 * n_sequences * mu with n_sequences = 10
            ({}, 0.01),
        ],
    )
    def test_every_mutation_parameter_maps_to_a_per_site_theta(self, options, per_site_theta):
        stats = _replicates(lambda rng: generate_population_sequences(10, 500, rng=rng, **options), 120)
        mean, se = _mean_se(stats["pi"])
        assert abs(mean - per_site_theta) <= Z * se

    def test_the_first_parameter_given_decides(self):
        """Order: nucleotide_diversity, wattersons_theta, theta, mutation_rate; the others are ignored."""

        def make(**options: Any) -> list[str]:
            return generate_population_sequences(10, 300, rng=random.Random(1), **options)

        reference = make(nucleotide_diversity=0.01)
        assert reference == make(nucleotide_diversity=0.01, wattersons_theta=0.05, theta=0.04, mutation_rate=0.001)
        assert reference == make(wattersons_theta=0.01, theta=0.04, mutation_rate=0.001)
        assert reference == make(theta=0.01, mutation_rate=0.001)
        assert reference == make(mutation_rate=0.00025)  # 4 * 10 * 0.00025 = 0.01
        assert reference == make()  # the default is 0.01 per site
        assert reference != make(nucleotide_diversity=0.02)

    def test_sequences_share_one_genealogy(self):
        """No recombination and no recurrent mutation: every pair of variable sites shows at most three gametes."""
        for seed in range(20):
            sequences = generate_population_sequences(20, 500, nucleotide_diversity=0.01, rng=random.Random(seed))
            columns = [column for column in zip(*sequences) if len(set(column)) > 1]
            assert columns
            assert all(len(set(column)) == 2 for column in columns)
            for index, first in enumerate(columns):
                for second in columns[index + 1 :]:
                    assert len(set(zip(first, second))) <= 3

    def test_reference_sequence_is_the_ancestor(self):
        reference = ("ACGTN" * 40) + "acgtn" * 20  # 300 positions, 60 of them ambiguous
        sequences = generate_population_sequences(10, reference_sequence=reference, theta=0.05, rng=random.Random(9))
        assert {len(s) for s in sequences} == {300}
        columns = list(zip(*sequences))
        variable = {i for i, column in enumerate(columns) if len(set(column)) > 1}
        assert variable and all(reference[i] not in "Nn" for i in variable)  # only A/C/G/T positions mutate
        assert all(set(column) == {reference[i]} for i, column in enumerate(columns) if i not in variable)

    def test_n_sites_is_accepted_but_no_longer_limits_the_variation(self):
        """n_sites used to cap the candidate positions; the number of variable sites is now an outcome of theta."""
        kwargs: dict[str, Any] = dict(n_sequences=20, length=1000, nucleotide_diversity=0.01)
        sequences = generate_population_sequences(**kwargs, n_sites=5, rng=random.Random(1))
        assert segregating_sites(sequences) > 5
        assert sequences == generate_population_sequences(**kwargs, rng=random.Random(1))

    def test_infinite_sites_capacity_is_enforced(self):
        # theta = 4 * 10 * 0.01 = 0.4 per site: about 113 expected segregating sites on 100 positions.
        with pytest.raises(ValueError, match="infinite-sites"):
            generate_population_sequences(10, 100, mutation_rate=0.01, rng=random.Random(1))
        # Below that bound a draw can still need more sites than there are (here 8 expected on 10 positions).
        outcomes = []
        for seed in range(60):
            try:
                generate_population_sequences(2, 10, theta=0.8, rng=random.Random(seed))
                outcomes.append(True)
            except ValueError as error:
                assert "do not fit" in str(error)
                outcomes.append(False)
        assert any(outcomes) and not all(outcomes)

    def test_small_theta_can_leave_the_sample_monomorphic(self):
        """No mutation is forced into a sample: with theta = 1, P(S = 0) = prod j / (j + 1) = 1/5 for 5 sequences."""
        replicates = 600
        identical = sum(
            len(set(generate_population_sequences(5, 100, theta=0.01, rng=random.Random(seed)))) == 1
            for seed in range(replicates)
        )
        expected = math.prod(j / (j + 1.0) for j in range(1, 5))
        assert abs(identical / replicates - expected) <= Z * math.sqrt(expected * (1 - expected) / replicates)
        assert len(set(generate_population_sequences(10, 100, theta=0.0, rng=random.Random(1)))) == 1

    def test_random_seed_alias_matches_rng(self):
        assert generate_population_sequences(10, 200, random_seed=5) == generate_population_sequences(
            10, 200, rng=random.Random(5)
        )

    def test_invalid_parameters_are_rejected(self):
        with pytest.raises(ValidationError):
            generate_population_sequences(1, 100)
        with pytest.raises(ValidationError):
            generate_population_sequences(10, 100, theta=-0.1)
        with pytest.raises(ValidationError):
            generate_population_sequences(10, 100, n_sites=0)
        with pytest.raises(ValueError, match="length"):
            generate_population_sequences(10)
        with pytest.raises(ValueError, match="finite"):
            generate_population_sequences(10, 100, theta=float("nan"))


class TestBottleneckSequences:
    """simulate_bottleneck_population samples a coalescent population that shrank for the bottleneck duration."""

    def test_defaults_lower_diversity_and_give_positive_d(self, neutral, bottleneck):
        """Measured over 150 samples: D = +0.9 against -0.04 and pi = 0.0038 against 0.0100."""
        mean_d, se_d = _mean_se(bottleneck["d"])
        assert mean_d - Z * se_d > 0.3
        difference, se = _difference_se(neutral["pi"], bottleneck["pi"])
        assert difference - Z * se > 0.003

    def test_zero_duration_is_the_constant_size_sample(self):
        kwargs: dict[str, Any] = dict(n_sequences=20, sequence_length=300)
        assert simulate_bottleneck_population(**kwargs, bottleneck_duration=0, rng=random.Random(1)) == (
            generate_population_sequences(20, 300, theta=0.01, rng=random.Random(1))
        )

    def test_bottleneck_duration_is_honoured(self, neutral, bottleneck):
        """Longer bottlenecks coalesce more of the sample, so diversity falls steadily with the duration."""
        long = _replicates(
            lambda rng: simulate_bottleneck_population(
                n_sequences=20, sequence_length=1000, bottleneck_duration=40, rng=rng
            )
        )
        for shorter, longer in ((neutral, bottleneck), (bottleneck, long)):
            difference, se = _difference_se(shorter["pi"], longer["pi"])
            assert difference > Z * se

    def test_pre_bottleneck_diversity_scales_the_diversity(self, bottleneck):
        double = _replicates(
            lambda rng: simulate_bottleneck_population(
                n_sequences=20, sequence_length=1000, pre_bottleneck_diversity=0.02, rng=rng
            )
        )
        difference, se = _difference_se(double["pi"], 2 * bottleneck["pi"])
        assert abs(difference) <= Z * se

    def test_initial_size_sets_the_time_scale(self, bottleneck):
        """In a small ancestral population the same bottleneck spans less coalescent time, so it bites less."""
        small = _replicates(
            lambda rng: simulate_bottleneck_population(n_sequences=20, sequence_length=1000, initial_size=10, rng=rng)
        )
        difference, se = _difference_se(small["pi"], bottleneck["pi"])
        assert difference > Z * se
        difference, se = _difference_se(bottleneck["d"], small["d"])
        assert difference > Z * se

    def test_parameters_without_a_role_in_the_sequence_interface_change_nothing(self):
        """mutation_rate, final_size and generations belong to the statistics interface (see the docstring)."""
        kwargs: dict[str, Any] = dict(n_sequences=12, sequence_length=300)
        base = simulate_bottleneck_population(**kwargs, rng=random.Random(3))
        unused = dict(mutation_rate=0.5, final_size=77, generations=5)
        assert base == simulate_bottleneck_population(**kwargs, **unused, rng=random.Random(3))

    def test_invalid_parameters_are_rejected(self):
        with pytest.raises(ValueError, match="sequence_length"):
            simulate_bottleneck_population(n_sequences=10)
        with pytest.raises(ValidationError, match="initial_size"):
            simulate_bottleneck_population(n_sequences=10, sequence_length=100, bottleneck_size=50, initial_size=20)
        with pytest.raises(ValidationError, match="bottleneck_duration"):
            simulate_bottleneck_population(n_sequences=10, sequence_length=100, bottleneck_duration=-1)
        with pytest.raises(ValidationError, match="theta"):
            simulate_bottleneck_population(n_sequences=10, sequence_length=100, pre_bottleneck_diversity=-0.1)
        with pytest.raises(ValidationError, match="bottleneck_size"):
            simulate_bottleneck_population(n_sequences=10, sequence_length=100, bottleneck_size=1)


class TestExpansionSequences:
    """simulate_population_expansion samples a coalescent population that grew exponentially up to the sample."""

    def test_defaults_give_negative_d_and_raise_diversity(self, neutral, expansion):
        """Measured over 150 samples: D = -0.9 against -0.04 and pi = 0.0130 against the ancestral 0.0100."""
        mean_d, se_d = _mean_se(expansion["d"])
        assert mean_d + Z * se_d < -0.5
        difference, se = _difference_se(expansion["pi"], neutral["pi"])
        assert difference - Z * se > 0.0

    def test_factor_one_is_the_constant_size_sample(self):
        kwargs: dict[str, Any] = dict(n_sequences=20, sequence_length=300)
        assert simulate_population_expansion(**kwargs, expansion_factor=1.0, rng=random.Random(1)) == (
            generate_population_sequences(20, 300, theta=0.01, rng=random.Random(1))
        )

    def test_expansion_factor_is_honoured(self, neutral, expansion):
        """A bigger expansion leaves a stronger excess of rare variants: D falls from -0.04 to -0.9 to -1.3."""
        big = _replicates(
            lambda rng: simulate_population_expansion(
                n_sequences=20, sequence_length=1000, expansion_factor=100.0, rng=rng
            )
        )
        for weaker, stronger in ((neutral, expansion), (expansion, big)):
            difference, se = _difference_se(weaker["d"], stronger["d"])
            assert difference > Z * se

    def test_expansion_time_is_honoured(self, expansion):
        """Growth over 5 generations is too brief to matter; over 100 it leaves a clear signature."""
        brief = _replicates(
            lambda rng: simulate_population_expansion(n_sequences=20, sequence_length=1000, expansion_time=5, rng=rng)
        )
        difference, se = _difference_se(brief["d"], expansion["d"])
        assert difference > Z * se

    def test_final_size_is_another_way_to_set_the_factor(self):
        kwargs: dict[str, Any] = dict(n_sequences=20, sequence_length=300)
        by_sizes = simulate_population_expansion(**kwargs, initial_size=100, final_size=1000, rng=random.Random(2))
        assert by_sizes == simulate_population_expansion(**kwargs, expansion_factor=10.0, rng=random.Random(2))
        assert by_sizes == simulate_population_expansion(
            **kwargs, initial_size=100, final_size=1000, expansion_factor=3.0, rng=random.Random(2)
        )  # final_size wins

    def test_parameters_without_a_role_in_the_sequence_interface_change_nothing(self):
        """mutation_rate belongs to the statistics interface: the ancestral theta per site is fixed at 0.01."""
        kwargs: dict[str, Any] = dict(n_sequences=12, sequence_length=300)
        base = simulate_population_expansion(**kwargs, rng=random.Random(3))
        assert base == simulate_population_expansion(**kwargs, mutation_rate=0.5, rng=random.Random(3))

    def test_invalid_parameters_are_rejected(self):
        with pytest.raises(ValueError, match="sequence_length"):
            simulate_population_expansion(n_sequences=10)
        with pytest.raises(ValidationError, match="expansion_factor"):
            simulate_population_expansion(n_sequences=10, sequence_length=100, expansion_factor=0.5)
        with pytest.raises(ValidationError, match="final_size"):
            simulate_population_expansion(n_sequences=10, sequence_length=100, initial_size=100, final_size=50)
        with pytest.raises(ValidationError, match="expansion_time"):
            simulate_population_expansion(n_sequences=10, sequence_length=100, expansion_time=0)
        with pytest.raises(ValidationError, match="initial_size"):
            simulate_population_expansion(n_sequences=10, sequence_length=100, initial_size=0)


class TestGenerateSiteFrequencySpectrum:
    """Test site frequency spectrum generation."""

    def test_folded_sfs(self):
        """Test folded SFS generation."""
        sfs = generate_site_frequency_spectrum(
            sample_size=10,
            n_sites=100,
            folded=True,
        )
        assert len(sfs) == 5  # n//2 for n=10
        assert sum(sfs) == 100

    def test_unfolded_sfs(self):
        """Test unfolded SFS generation."""
        sfs = generate_site_frequency_spectrum(
            sample_size=10,
            n_sites=100,
            folded=False,
        )
        assert len(sfs) == 9  # n-1 for n=10
        assert sum(sfs) == 100

    def test_sfs_rare_alleles(self):
        """Test that SFS has more rare alleles."""
        sfs = generate_site_frequency_spectrum(
            sample_size=10,
            n_sites=100,
            folded=True,
        )
        # Under neutral model, rare alleles should be more common
        # (first bin should have more sites)
        assert sfs[0] > 0  # Should have some rare alleles


class TestGenerateLinkageDisequilibriumData:
    """Test linkage disequilibrium data generation."""

    def test_basic_ld_generation(self):
        """Test basic LD data generation."""
        genotypes = generate_linkage_disequilibrium_data(
            n_individuals=100,
            n_sites=10,
            r_squared_target=0.5,
        )
        assert len(genotypes) == 100
        assert all(len(row) == 10 for row in genotypes)
        assert all(0 <= g <= 2 for row in genotypes for g in row)

    def test_with_allele_frequencies(self):
        """Test LD generation with specified frequencies."""
        freqs = [0.3, 0.4, 0.2, 0.5, 0.1]
        genotypes = generate_linkage_disequilibrium_data(
            n_individuals=100,
            n_sites=5,
            allele_frequencies=freqs,
        )
        assert len(genotypes) == 100


class TestBottleneckStatsInterface:
    """Test the statistics interface of simulate_bottleneck_population.

    Regression tests: this interface previously called ``rng.poisson`` on a
    ``random.Random`` and raised AttributeError at runtime.
    """

    def test_stats_interface_returns_trajectory(self):
        """Test that the default interface returns simulation statistics."""
        result = simulate_bottleneck_population(
            initial_size=100,
            bottleneck_size=10,
            final_size=100,
            generations=6,
            rng=random.Random(42),
        )
        assert isinstance(result, dict)
        assert len(result["population_sizes"]) == 6
        assert result["population_sizes"][:2] == [100, 100]  # pre-bottleneck
        assert result["bottleneck_start"] == 2
        assert result["bottleneck_end"] == 4
        assert result["population_sizes"][3] == 10  # during bottleneck
        assert result["total_mutations"] >= 0
        assert 0.0 < result["final_diversity"] <= 1.0

    def test_stats_interface_is_deterministic(self):
        """Test reproducibility with a seeded RNG."""
        kwargs = dict(initial_size=50, bottleneck_size=5, final_size=50, generations=6)
        result_a = simulate_bottleneck_population(**kwargs, rng=random.Random(7))
        result_b = simulate_bottleneck_population(**kwargs, rng=random.Random(7))
        assert result_a == result_b


class TestExpansionStatsInterface:
    """Test the statistics interface of simulate_population_expansion.

    Regression tests: this interface previously called ``rng.poisson`` on a
    ``random.Random`` and raised AttributeError at runtime.
    """

    def test_stats_interface_returns_trajectory(self):
        """Test that the default interface returns simulation statistics."""
        result = simulate_population_expansion(
            initial_size=100,
            final_size=1000,
            expansion_time=5,
            rng=random.Random(42),
        )
        assert isinstance(result, dict)
        assert len(result["population_sizes"]) == 6  # expansion_time + 1
        assert result["population_sizes"][0] == 100
        assert result["population_sizes"][-1] == pytest.approx(1000, abs=1)
        assert result["growth_rate"] > 0
        assert result["total_mutations"] >= 0
        assert result["final_diversity"] <= 1.0


class TestSFSExpansionModel:
    """Test the expansion demographic model of generate_site_frequency_spectrum.

    Regression tests: this model previously called ``rng.zipf`` on a
    ``random.Random`` and raised AttributeError at runtime.
    """

    def test_expansion_model_default_alpha(self):
        """Test that the default alpha (1.0) produces a valid SFS."""
        sfs = generate_site_frequency_spectrum(
            n_samples=10,
            n_sites=100,
            demographic_model="expansion",
            rng=random.Random(42),
        )
        assert sum(sfs) == 100
        assert all(count >= 0 for count in sfs)

    def test_expansion_model_favors_rare_alleles(self):
        """Test that zipf-distributed frequencies favor low-frequency variants."""
        sfs = generate_site_frequency_spectrum(
            n_samples=50,
            n_sites=500,
            demographic_model="expansion",
            parameters={"alpha": 2.0},
            rng=random.Random(42),
        )
        assert sfs[0] > sfs[-1]  # singletons more common than high-frequency variants

    def test_expansion_model_custom_alpha_changes_spectrum(self):
        """Test that a provided alpha parameter is honored."""
        kwargs = dict(n_samples=10, n_sites=100, demographic_model="expansion")
        sfs_default = generate_site_frequency_spectrum(**kwargs, rng=random.Random(42))
        sfs_alpha2 = generate_site_frequency_spectrum(**kwargs, parameters={"alpha": 2.0}, rng=random.Random(42))
        assert sfs_default != sfs_alpha2


class TestSimulateAdmixture:
    """Test population admixture simulation.

    Regression tests: the drift step previously called ``rng.normal`` on a
    ``random.Random`` and raised AttributeError at runtime.
    """

    def test_basic_admixture(self):
        """Test that admixture runs and produces a frequency trajectory."""
        proportions = np.array([[0.9, 0.1], [0.1, 0.9]])
        result = simulate_admixture(2, [10, 10], proportions, 3, rng=random.Random(42))
        assert len(result["frequency_trajectory"]) == 4  # generations + 1
        assert len(result["ancestral_frequencies"]) == 2
        assert result["generations"] == 3
        assert result["migration_rate"] == 0.01
        for freqs in result["frequency_trajectory"]:
            assert all(0.0 <= f <= 1.0 for f in freqs)

    def test_admixture_is_deterministic(self):
        """Test reproducibility with a seeded RNG."""
        proportions = np.array([[0.9, 0.1], [0.1, 0.9]])
        result_a = simulate_admixture(2, [10, 10], proportions, 3, rng=random.Random(7))
        result_b = simulate_admixture(2, [10, 10], proportions, 3, rng=random.Random(7))
        assert all(
            np.array_equal(a, b) for a, b in zip(result_a["frequency_trajectory"], result_b["frequency_trajectory"])
        )

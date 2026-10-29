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
    PopulationEpoch,
    generate_genotype_matrix,
    generate_linkage_disequilibrium_data,
    generate_population_sequences,
    generate_site_frequency_spectrum,
    generate_two_populations,
    simulate_admixture,
    simulate_bottleneck_population,
    simulate_coalescent_spectra,
    simulate_population_expansion,
    simulate_selection,
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


def _tajima_numerator(spectrum: list[int], n: int) -> float:
    """pi - S / a_n of a pooled spectrum: negative when rare variants are in excess, positive when they are lacking."""
    i = np.arange(1, n)
    counts = np.asarray(spectrum, dtype=float)
    return float(counts @ (i * (n - i) / (n * (n - 1) / 2)) - counts.sum() / np.sum(1.0 / i))


def _singleton_share(spectrum: list[int]) -> float:
    return spectrum[0] / sum(spectrum)


class TestSiteFrequencySpectrumModels:
    """generate_site_frequency_spectrum draws each site's class from the coalescent of the demographic model.

    Expansion and bottleneck used to be heuristics (a Zipf draw, and a coin between two uniform ranges). They are now
    the expected branch-length spectrum of seeded coalescent genealogies, and the constant model is its closed form.
    """

    N = 20  # chromosomes
    SITES = 200_000  # sites per spectrum: multinomial error of a class share is at most 0.0011

    @staticmethod
    def _shares_within_error(spectrum: list[int], expected: np.ndarray, extra_se: np.ndarray, label: str) -> None:
        counts, total = np.asarray(spectrum, dtype=float), float(sum(spectrum))
        share = counts / total
        se = np.sqrt(expected * (1 - expected) / total + extra_se**2)
        worst = np.abs(share - expected) / se
        assert (
            worst < 4.0
        ).all(), f"{label}: classes {np.flatnonzero(worst >= 4.0) + 1} deviate (z = {worst.max():.1f})"

    def test_the_constant_model_has_the_closed_form_one_over_i(self) -> None:
        sfs = generate_site_frequency_spectrum(self.N, self.SITES, rng=random.Random(1))
        weights = 1.0 / np.arange(1, self.N)
        self._shares_within_error(sfs, weights / weights.sum(), np.zeros(self.N - 1), "constant")
        assert len(sfs) == self.N - 1 and sum(sfs) == self.SITES

    @pytest.mark.parametrize(
        "model, parameters",
        [
            ("expansion", {"expansion_factor": 10.0, "expansion_time": 0.5}),
            ("expansion", {"expansion_factor": 100.0, "expansion_time": 0.2}),
            ("bottleneck", {"bottleneck_size": 0.1, "bottleneck_duration": 0.05}),
            ("bottleneck", {"bottleneck_size": 0.05, "bottleneck_duration": 0.03, "recovery_time": 0.04}),
        ],
    )
    def test_the_demographic_models_follow_an_independent_simulation_of_their_history(
        self, model: str, parameters: dict
    ) -> None:
        """E[xi_i] / E[S] from Poisson mutations on separately simulated genealogies is the same ratio E[L_i] / E[L]."""
        if model == "expansion":
            duration, factor = parameters["expansion_time"], parameters["expansion_factor"]
            history = [PopulationEpoch(0.0, factor, math.log(factor) / duration), PopulationEpoch(duration, 1.0)]
        else:
            recovery = parameters.get("recovery_time", 0.0)
            size, duration = parameters["bottleneck_size"], parameters["bottleneck_duration"]
            history = (
                [PopulationEpoch(0.0, 1.0), PopulationEpoch(recovery, size), PopulationEpoch(recovery + duration, 1.0)]
                if recovery
                else [PopulationEpoch(0.0, size), PopulationEpoch(duration, 1.0)]
            )
        replicates = 100_000
        reference = simulate_coalescent_spectra(self.N, 60_000, theta=2.0, epochs=history, rng=3).astype(float)
        mean = reference.mean(axis=0)
        expected = mean / mean.sum()
        se_reference = reference.std(axis=0, ddof=1) / math.sqrt(reference.shape[0]) / mean.sum()
        # The function's own genealogies are few too: the class shares carry an error below sd(xi_i) / sqrt(R)
        se_function = reference.std(axis=0, ddof=1) / math.sqrt(replicates) / mean.sum()
        sfs = generate_site_frequency_spectrum(
            self.N,
            self.SITES,
            demographic_model=model,
            parameters={**parameters, "n_replicates": replicates},
            rng=random.Random(2),
        )
        self._shares_within_error(sfs, expected, np.hypot(se_reference, se_function), f"{model} {parameters}")

    def test_expansion_has_an_excess_of_rare_variants_and_a_negative_tajima_numerator(self) -> None:
        neutral = generate_site_frequency_spectrum(self.N, self.SITES, rng=random.Random(4))
        grown = generate_site_frequency_spectrum(
            self.N, self.SITES, demographic_model="expansion", rng=random.Random(4)
        )
        assert _singleton_share(grown) > _singleton_share(neutral) + 0.15  # measured 0.49 against 0.28
        assert sum(grown[self.N // 2 :]) < sum(neutral[self.N // 2 :])  # fewer common variants
        # The numerator per site is 0.000 for the neutral spectrum and -0.061 after the expansion
        assert abs(_tajima_numerator(neutral, self.N)) < 0.005 * self.SITES
        assert _tajima_numerator(grown, self.N) < -0.04 * self.SITES

    def test_a_recent_bottleneck_lacks_rare_variants_and_has_a_positive_tajima_numerator(self) -> None:
        neutral = generate_site_frequency_spectrum(self.N, self.SITES, rng=random.Random(5))
        squeezed = generate_site_frequency_spectrum(
            self.N, self.SITES, demographic_model="bottleneck", rng=random.Random(5)
        )
        assert _singleton_share(squeezed) < _singleton_share(neutral) - 0.1  # measured 0.12 against 0.28
        assert _tajima_numerator(squeezed, self.N) > 0.04 * self.SITES  # measured +0.066 per site

    def test_recovery_restores_the_rare_variants_a_bottleneck_removed(self) -> None:
        shares = [
            _singleton_share(
                generate_site_frequency_spectrum(
                    self.N,
                    self.SITES,
                    demographic_model="bottleneck",
                    parameters={"recovery_time": recovery},
                    rng=random.Random(6),
                )
            )
            for recovery in (0.0, 0.02, 0.2)
        ]
        assert shares[0] < shares[1] < shares[2]  # measured 0.12, 0.19, 0.31
        assert shares[2] > 0.282  # past the neutral share once new mutations have filled the long branches

    def test_stronger_expansion_leaves_more_rare_variants(self) -> None:
        shares = [
            _singleton_share(
                generate_site_frequency_spectrum(
                    self.N,
                    self.SITES,
                    demographic_model="expansion",
                    parameters={"expansion_factor": factor},
                    rng=random.Random(7),
                )
            )
            for factor in (1.0, 10.0, 100.0)
        ]
        assert shares[0] + 0.1 < shares[1] < shares[2] - 0.05  # measured 0.28, 0.49, 0.60

    def test_models_without_demographic_change_are_the_constant_model(self) -> None:
        neutral = 1.0 / np.arange(1, self.N)
        for model, parameters in (
            ("expansion", {"expansion_factor": 1.0}),
            ("bottleneck", {"bottleneck_size": 1.0}),
            ("bottleneck", {"bottleneck_duration": 0.0}),
        ):
            sfs = generate_site_frequency_spectrum(
                self.N, self.SITES, demographic_model=model, parameters=parameters, rng=random.Random(8)
            )
            self._shares_within_error(sfs, neutral / neutral.sum(), np.zeros(self.N - 1), model)
            # the closed form is used: the same draw as the constant model's
            assert sfs == generate_site_frequency_spectrum(self.N, self.SITES, rng=random.Random(8))

    def test_folding_adds_complementary_classes_and_keeps_every_site(self) -> None:
        for n_samples in (10, 11):
            kwargs = dict(n_sites=5000, demographic_model="expansion", parameters={"n_replicates": 500})
            unfolded = generate_site_frequency_spectrum(n_samples, rng=random.Random(9), **kwargs)
            folded = generate_site_frequency_spectrum(n_samples, folded=True, rng=random.Random(9), **kwargs)
            expected = [unfolded[i] + unfolded[n_samples - 2 - i] for i in range(n_samples // 2)]
            if n_samples % 2 == 0:
                expected[-1] = unfolded[n_samples // 2 - 1]  # the middle class is its own complement
            assert folded == expected and sum(folded) == 5000

    def test_the_spectrum_is_seeded_and_does_not_touch_the_global_numpy_state(self) -> None:
        kwargs = dict(demographic_model="bottleneck", parameters={"n_replicates": 400})
        np.random.seed(123)
        before = np.random.get_state()[1].copy()
        first = generate_site_frequency_spectrum(12, 400, rng=random.Random(3), **kwargs)
        assert np.array_equal(np.random.get_state()[1], before)
        assert first == generate_site_frequency_spectrum(12, 400, rng=random.Random(3), **kwargs)
        assert first != generate_site_frequency_spectrum(12, 400, rng=random.Random(4), **kwargs)

    def test_two_chromosomes_have_one_class(self) -> None:
        assert generate_site_frequency_spectrum(2, 50, demographic_model="expansion", rng=random.Random(1)) == [50]

    @pytest.mark.parametrize(
        "model, parameters",
        [
            ("expansion", {"alpha": 2.0}),
            ("bottleneck", {"bottleneck_strength": 0.3}),
            ("constant", {"expansion_factor": 3.0}),
            ("bottleneck", {"expansion_factor": 3.0}),
        ],
    )
    def test_names_a_model_does_not_take_are_rejected_not_ignored(self, model: str, parameters: dict) -> None:
        with pytest.raises(ValidationError, match="unknown parameter"):
            generate_site_frequency_spectrum(10, 100, demographic_model=model, parameters=parameters)

    @pytest.mark.parametrize(
        "model, parameters, message",
        [
            ("expansion", {"expansion_factor": 0.5}, "expansion_factor"),
            ("expansion", {"expansion_time": -1.0}, "expansion_time"),
            ("expansion", {"expansion_time": 0.0}, "expansion_time"),
            ("expansion", {"expansion_factor": float("nan")}, "finite"),
            ("expansion", {"n_replicates": 0}, "n_replicates"),
            ("expansion", {"n_replicates": 10.5}, "n_replicates"),
            ("bottleneck", {"bottleneck_size": 0.0}, "bottleneck_size"),
            ("bottleneck", {"bottleneck_size": 1.5}, "bottleneck_size"),
            ("bottleneck", {"bottleneck_duration": -0.1}, "bottleneck_duration"),
            ("bottleneck", {"recovery_time": -0.1}, "recovery_time"),
            ("bottleneck", {"bottleneck_duration": float("inf")}, "finite"),
        ],
    )
    def test_invalid_model_parameters_are_rejected(self, model: str, parameters: dict, message: str) -> None:
        with pytest.raises(ValidationError, match=message):
            generate_site_frequency_spectrum(10, 100, demographic_model=model, parameters=parameters)

    def test_unknown_models_and_missing_sizes_are_rejected(self) -> None:
        with pytest.raises(ValidationError, match="Unknown demographic model"):
            generate_site_frequency_spectrum(10, 100, demographic_model="migration")
        with pytest.raises(ValueError, match="n_samples"):
            generate_site_frequency_spectrum(n_sites=10)
        with pytest.raises(ValueError, match="n_sites"):
            generate_site_frequency_spectrum(10)


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


def _simulated_mean_pair_time(epochs: list[PopulationEpoch]) -> tuple[float, float]:
    """Mean and standard error of the coalescence time of two chromosomes (E[S] for theta = 1) under ``epochs``."""
    times = simulate_coalescent_spectra(2, 200_000, theta=1.0, epochs=epochs, rng=3)[:, 0].astype(float)
    return float(times.mean()), float(times.std(ddof=1) / math.sqrt(times.size))


class TestStatisticsInterfacesFromTheCoalescent:
    """The statistics interfaces report diversity, mutation counts and a spectrum that follow from the coalescent.

    Their diversity trajectory was ``size / initial_size`` floored at 0.01 (and pinned at 1.0 for expansions), their
    mutation counts used a hard-coded locus of 1000 bases and the global numpy generator.
    """

    BOTTLENECK = dict(initial_size=1000, bottleneck_size=100, final_size=1000, generations=300, mutation_rate=1e-6)

    def test_diversity_starts_at_equilibrium_and_follows_the_drift_recursion(self) -> None:
        result = simulate_bottleneck_population(**self.BOTTLENECK, rng=random.Random(1))
        sizes, mu = result["population_sizes"], 1e-6
        assert sizes[:100] == [1000] * 100 and sizes[100:200] == [100] * 100 and sizes[200:] == [1000] * 100
        pi, expected = 4 * 1000 * mu, []
        for size in sizes:  # pi' = (1 - 1 / 2N) pi + 2 mu, started at the equilibrium 4 N mu
            pi = (1 - 1 / (2 * size)) * pi + 2 * mu
            expected.append(pi)
        assert result["diversity_trajectory"] == pytest.approx(expected, rel=1e-12)
        assert result["diversity_trajectory"][0] == pytest.approx(4 * 1000 * mu, rel=1e-12)  # drift-mutation balance
        assert result["final_diversity"] == result["diversity_trajectory"][-1]
        trajectory = result["diversity_trajectory"]
        assert (np.diff(trajectory[99:200]) < 0).all() and trajectory[199] < 0.7 * trajectory[0]  # measured 0.645
        assert (np.diff(trajectory[199:]) > 0).all()  # ... and recovers slowly

    def test_final_diversity_is_the_coalescent_expectation_of_the_pairwise_time(self) -> None:
        """E[pi] = 4 N_ref mu E[T] per site, with E[T] simulated under the history the size trajectory describes."""
        result = simulate_bottleneck_population(**self.BOTTLENECK, rng=random.Random(1))
        history = [PopulationEpoch(0.0, 1.0), PopulationEpoch(0.05, 0.1), PopulationEpoch(0.10, 1.0)]
        mean_time, se = _simulated_mean_pair_time(history)
        predicted = 4 * 1000 * 1e-6 * mean_time
        # Discrete generations against the continuous approximation: a generation of timing at each bottleneck edge
        assert abs(result["final_diversity"] - predicted) <= Z * 4 * 1000 * 1e-6 * se + 0.005 * predicted

    def test_expansion_diversity_grows_toward_the_new_equilibrium(self) -> None:
        result = simulate_population_expansion(
            initial_size=100, final_size=10_000, expansion_time=2000, mutation_rate=1e-6, rng=random.Random(2)
        )
        trajectory = result["diversity_trajectory"]
        assert trajectory[0] == pytest.approx(4 * 100 * 1e-6, rel=1e-12)
        assert (np.diff(trajectory) >= 0).all()  # more and more diverse as the population grows
        # E[pi] = 4 N_ref mu E[T] for the history of this exponential growth: the old code capped it at 1.0
        duration, factor = 2000 / (2 * 100), 100.0
        mean_time, se = _simulated_mean_pair_time(
            [PopulationEpoch(0.0, factor, math.log(factor) / duration), PopulationEpoch(duration, 1.0)]
        )
        predicted = 4 * 100 * 1e-6 * mean_time
        assert abs(result["final_diversity"] - predicted) <= Z * 4 * 100 * 1e-6 * se + 0.01 * predicted
        assert result["final_diversity"] > 5 * trajectory[0]

    def test_mutation_counts_are_poisson_with_the_documented_locus_and_population(self) -> None:
        kwargs = dict(self.BOTTLENECK, mutation_rate=1e-4)
        for locus, seed in ((500, 4), (1000, 5)):
            result = simulate_bottleneck_population(**kwargs, sequence_length=locus, rng=random.Random(seed))
            means = 2 * np.asarray(result["population_sizes"]) * 1e-4 * locus  # 2 N mu L new mutations per generation
            trajectory = result["mutation_trajectory"]
            assert (np.diff(trajectory) >= 0).all() and trajectory[-1] == result["total_mutations"]
            assert all(type(count) is int for count in trajectory)
            assert abs(result["total_mutations"] - means.sum()) <= Z * math.sqrt(means.sum())
        default = simulate_bottleneck_population(**kwargs, rng=random.Random(5))  # the old hard-coded locus: 1000
        assert default == simulate_bottleneck_population(**kwargs, sequence_length=1000, rng=random.Random(5))

    def test_a_bottleneck_sampled_soon_after_has_few_rare_variants_and_a_positive_d(self) -> None:
        result = simulate_bottleneck_population(
            initial_size=1000, bottleneck_size=20, final_size=1000, generations=60, rng=random.Random(1)
        )
        assert result["tajimas_d"] - Z * result["tajimas_d_se"] > 0.4  # measured +0.75
        assert result["singleton_fraction"] < result["neutral_singleton_fraction"] - 0.1  # measured 0.14 vs 0.28
        expected = np.asarray(result["expected_sfs"])
        assert expected[0] / expected.sum() < result["neutral_singleton_fraction"]
        assert result["tajimas_d_se"] > 0 and result["n_replicates"] == 2000 and result["sample_size"] == 20

    def test_a_long_recovery_leaves_the_star_like_signature_of_expansion_instead(self) -> None:
        result = simulate_bottleneck_population(
            initial_size=1000, bottleneck_size=20, final_size=1000, generations=300, rng=random.Random(1)
        )
        assert result["tajimas_d"] + Z * result["tajimas_d_se"] < -0.5  # measured -1.64
        assert result["singleton_fraction"] > result["neutral_singleton_fraction"] + 0.1

    def test_an_expansion_has_an_excess_of_rare_variants_and_a_negative_d(self) -> None:
        result = simulate_population_expansion(rng=random.Random(2))  # 100 -> 10000 over 100 generations
        assert result["tajimas_d"] + Z * result["tajimas_d_se"] < -0.8  # measured -1.29
        assert result["singleton_fraction"] > result["neutral_singleton_fraction"] + 0.2  # measured 0.60 vs 0.28
        expected, neutral = np.asarray(result["expected_sfs"]), np.asarray(result["neutral_sfs"])
        assert expected[0] / expected.sum() > neutral[0] / neutral.sum()
        half = len(expected) // 2  # and a smaller share of the common variants
        assert expected[half:].sum() / expected.sum() < 0.7 * neutral[half:].sum() / neutral.sum()

    def test_the_expected_spectrum_follows_an_independent_simulation_of_the_history(self) -> None:
        sample, replicates = 15, 20_000
        theta_locus = 4 * 100 * 1e-4 * 250  # 4 N_ref mu L = 10
        result = simulate_population_expansion(
            initial_size=100,
            final_size=10_000,
            expansion_time=100,
            mutation_rate=1e-4,
            sequence_length=250,
            sample_size=sample,
            n_replicates=replicates,
            rng=random.Random(6),
        )
        duration = 100 / (2 * 100)
        history = [PopulationEpoch(0.0, 100.0, math.log(100.0) / duration), PopulationEpoch(duration, 1.0)]
        reference = simulate_coalescent_spectra(sample, 60_000, theta=theta_locus, epochs=history, rng=7).astype(float)
        sd = reference.std(axis=0, ddof=1)
        se = np.sqrt(sd**2 / reference.shape[0] + sd**2 / replicates)  # the function's error is below sd / sqrt(R)
        difference = np.abs(np.asarray(result["expected_sfs"]) - reference.mean(axis=0))
        assert (difference <= Z * se).all(), f"largest deviation {(difference / se).max():.1f} standard errors"
        assert result["neutral_sfs"] == pytest.approx([theta_locus / i for i in range(1, sample)])
        assert result["sample_size"] == sample and result["n_replicates"] == replicates

    def test_the_bottleneck_and_expansion_results_are_seeded_and_leave_the_global_numpy_state_alone(self) -> None:
        np.random.seed(321)
        before = np.random.get_state()[1].copy()
        for simulate, kwargs in (
            (simulate_bottleneck_population, dict(generations=30, n_replicates=100)),
            (simulate_population_expansion, dict(expansion_time=30, n_replicates=100)),
        ):
            first = simulate(rng=random.Random(8), mutation_rate=1e-4, **kwargs)
            assert first == simulate(rng=random.Random(8), mutation_rate=1e-4, **kwargs)
            assert first != simulate(rng=random.Random(9), mutation_rate=1e-4, **kwargs)
        assert np.array_equal(np.random.get_state()[1], before)

    def test_the_sample_options_are_validated(self) -> None:
        for simulate in (simulate_bottleneck_population, simulate_population_expansion):
            with pytest.raises(ValidationError, match="sample_size"):
                simulate(sample_size=3)
            with pytest.raises(ValidationError, match="n_replicates"):
                simulate(n_replicates=1)
            with pytest.raises(ValidationError, match="sequence_length"):
                simulate(sequence_length=0)


class TestAdmixtureDriftAndMigration:
    """Drift is binomial sampling of the population's chromosomes, and the migration rate takes effect."""

    PROPORTIONS = np.array([[0.7, 0.3], [0.2, 0.8]])

    def test_in_huge_populations_frequencies_and_ancestry_follow_the_mixing_matrix(self) -> None:
        rate = 0.3
        mixing = (1 - rate) * np.eye(2) + rate * self.PROPORTIONS
        result = simulate_admixture(2, [10**9, 10**9], self.PROPORTIONS, 6, migration_rate=rate, rng=random.Random(1))
        start = np.asarray(result["ancestral_frequencies"])
        for generation, (frequencies, ancestry) in enumerate(
            zip(result["frequency_trajectory"], result["admixture_trajectory"])
        ):
            power = np.linalg.matrix_power(mixing, generation)
            assert frequencies == pytest.approx(power @ start, abs=5e-4)  # binomial sd about 1e-5 per generation
            assert ancestry == pytest.approx(power, abs=1e-12)
            assert ancestry.sum(axis=1) == pytest.approx([1.0, 1.0])  # ancestry shares sum to one
        assert result["final_frequencies"] == pytest.approx(result["frequency_trajectory"][-1])

    def test_without_migration_the_populations_do_not_mix(self) -> None:
        result = simulate_admixture(2, [10**9, 10**9], self.PROPORTIONS, 5, migration_rate=0.0, rng=random.Random(2))
        assert result["frequency_trajectory"][-1] == pytest.approx(result["ancestral_frequencies"], abs=1e-3)
        assert result["admixture_trajectory"][-1] == pytest.approx(np.eye(2))
        mixed = simulate_admixture(2, [10**9, 10**9], self.PROPORTIONS, 5, migration_rate=1.0, rng=random.Random(2))
        assert mixed["admixture_trajectory"][-1] == pytest.approx(np.linalg.matrix_power(self.PROPORTIONS, 5))

    @pytest.mark.parametrize("size", [25, 400])
    def test_one_generation_of_drift_has_the_wright_fisher_variance(self, size: int) -> None:
        """Var(p' - p) = p (1 - p) / (2 N), whatever the frequency: the squared standardized change averages 1 / 2N."""
        scores = []
        for seed in range(1500):
            result = simulate_admixture(2, [size, size], np.eye(2), 1, migration_rate=0.0, rng=random.Random(seed))
            start, end = np.asarray(result["ancestral_frequencies"]), result["frequency_trajectory"][1]
            scores.extend((end - start) ** 2 / (start * (1 - start)))
        scores = np.asarray(scores)
        mean, se = scores.mean(), scores.std(ddof=1) / math.sqrt(scores.size)
        assert abs(mean - 1 / (2 * size)) <= Z * se, f"{mean:.5f} vs {1 / (2 * size):.5f} (standard error {se:.5f})"

    def test_frequencies_stay_proportions_and_the_global_numpy_state_is_untouched(self) -> None:
        np.random.seed(7)
        before = np.random.get_state()[1].copy()
        result = simulate_admixture(2, [5, 5], self.PROPORTIONS, 40, migration_rate=0.05, rng=random.Random(3))
        assert np.array_equal(np.random.get_state()[1], before)
        for frequencies in result["frequency_trajectory"]:
            assert ((frequencies >= 0) & (frequencies <= 1)).all()
        assert len(result["admixture_trajectory"]) == 41

    @pytest.mark.parametrize(
        "proportions, sizes, message",
        [
            (np.array([[0.7, 0.4], [0.2, 0.8]]), [10, 10], "every row summing to 1"),
            (np.array([[1.2, -0.2], [0.2, 0.8]]), [10, 10], "non-negative"),
            (np.array([[0.5, 0.5, 0.0], [0.2, 0.8, 0.0]]), [10, 10], "shape"),
            (np.eye(2), [10, 0], "population_sizes"),
            (np.eye(2), [10], "length"),
        ],
    )
    def test_invalid_inputs_are_rejected_not_clipped(self, proportions: np.ndarray, sizes: list, message: str) -> None:
        with pytest.raises(ValidationError, match=message):
            simulate_admixture(2, sizes, proportions, 3, rng=random.Random(1))


class TestNoSilentTruncationOrOverride:
    """Requests the generators cannot honour are errors, not quietly different requests."""

    def test_more_polymorphic_sites_than_positions_is_an_error_not_a_truncation(self) -> None:
        with pytest.raises(ValidationError, match="do not fit"):
            generate_two_populations(5, 5, 50, n_sites=80, fst=0.1, rng=random.Random(1))
        # the Watterson default for theta = 1 on 10 positions is about 35 sites
        with pytest.raises(ValidationError, match="do not fit"):
            generate_two_populations(10, 10, 10, theta=1.0, fst=0.1, rng=random.Random(1))
        pop1, pop2 = generate_two_populations(5, 5, 50, n_sites=50, fst=0.1, rng=random.Random(1))
        assert len(pop1[0]) == 50

    def test_an_explicit_hwe_deviation_is_honoured_when_hwe_is_off(self) -> None:
        """hwe=False overwrote every deviation with 0.3; a deviation of 1 leaves homozygotes only."""
        forced = generate_genotype_matrix(300, 4, hwe=False, hwe_deviation=1.0, rng=random.Random(2))
        assert {genotype for row in forced for genotype in row} == {0, 2}
        default = generate_genotype_matrix(300, 4, hwe=False, rng=random.Random(2))
        assert 1 in {genotype for row in default for genotype in row}  # hwe=False alone still means some deviation
        exact = generate_genotype_matrix(300, 4, hwe=True, hwe_deviation=0.0, rng=random.Random(2))
        assert 1 in {genotype for row in exact for genotype in row}

    def test_selection_with_a_negative_fitness_is_rejected_not_sampled_with_negative_weights(self) -> None:
        genotypes = np.array([[0, 1], [1, 2], [2, 0], [1, 1]])
        fitness = np.array([[1.0, 1.0, 1.0], [1.0, -20.0, 1.0]])  # 1 + 0.1 * -20 = -1
        with pytest.raises(ValidationError, match="fitness"):
            simulate_selection(genotypes, fitness, 3, selection_strength=0.1, rng=random.Random(1))
        result = simulate_selection(genotypes, np.ones((2, 3)), 3, selection_strength=0.1, rng=random.Random(1))
        assert result["final_genotypes"].shape == (4, 2)

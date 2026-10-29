"""Validation of the coalescent with recombination (Hudson's ancestral recombination graph).

Zero-mocks. The graph is checked structurally (every marginal tree is a valid binary tree, the trees tile the
sequence), against the non-recombining engine (rho = 0, and the per-site spectra at any rho), against Hudson's
expected number of recombination events, against the exact moments of the two-sample ancestral recombination
graph (a Markov chain solved here from first principles) and against the Ohta-Kimura decay of linkage
disequilibrium. Every tolerance is four standard errors measured from the replicates themselves, and the
simulations are seeded.
"""

from __future__ import annotations

import itertools
import math
import random

import numpy as np
import pytest
from scipy import stats

from metainformant.core.utils.errors import ValidationError
from metainformant.simulation.models import popgen
from metainformant.simulation.models.popgen import (
    PopulationEpoch,
    RecombinantGenealogy,
    RecombinantSample,
    _merge_lineages,
    _SizeHistory,
    _split_lineage,
    generate_linkage_disequilibrium_data,
    simulate_coalescent_spectra,
    simulate_recombinant_genealogy,
    simulate_recombinant_haplotypes,
)

Z = 4.0  # tolerance in measured standard errors


def _a(n: int) -> float:
    return sum(1 / i for i in range(1, n))


def _mean_se(values: np.ndarray) -> tuple[float, float]:
    return float(values.mean()), float(values.std(ddof=1) / math.sqrt(values.size))


def _assert_mean(values: np.ndarray, expected: float, label: str, se_extra: float = 0.0) -> None:
    mean, se = _mean_se(values)
    se = math.hypot(se, se_extra)
    assert abs(mean - expected) <= Z * se, f"{label}: mean {mean:.4f} vs {expected:.4f} (standard error {se:.4f})"


def _spectrum(haplotypes: np.ndarray) -> np.ndarray:
    """Mutations carried by 1..n-1 of the n sequences."""
    n = haplotypes.shape[0]
    return np.bincount(haplotypes.sum(axis=0, dtype=np.int64), minlength=n + 1)[1:n]


def _pairwise(spectra: np.ndarray, n: int) -> np.ndarray:
    i = np.arange(1, n)
    return spectra @ (i * (n - i) / (n * (n - 1) / 2))


BOTTLENECK = [PopulationEpoch(0.0, 1.0), PopulationEpoch(0.1, 0.05), PopulationEpoch(0.4, 1.0)]
GROWTH = [PopulationEpoch(0.0, 20.0, math.log(20.0) / 0.5), PopulationEpoch(0.5, 1.0)]


# --- the graph is a valid sequence of marginal trees -------------------------------------------------------


@pytest.mark.parametrize("n, length, rho", [(2, 50.0, 8.0), (6, 200.0, 12.0), (15, 1000.0, 30.0)])
def test_the_marginal_trees_are_binary_trees_that_tile_the_sequence(n: int, length: float, rho: float) -> None:
    for seed in range(6):
        genealogy = simulate_recombinant_genealogy(n, length, rho=rho, rng=seed)
        times = genealogy.node_times
        assert genealogy.breakpoints[0] == 0.0 and genealogy.breakpoints[-1] == length
        assert (np.diff(genealogy.breakpoints) > 0).all()
        assert genealogy.n_trees == genealogy.tree_lengths.size == genealogy.breakpoints.size - 1
        assert (times[:n] == 0).all() and (np.diff(times[n - 1 :]) > 0).all()  # nodes are numbered by time
        assert (times[genealogy.edge_parent] > times[genealogy.edge_child]).all()

        trees = list(genealogy.trees())
        assert len(trees) == genealogy.n_trees
        for index, tree in enumerate(trees):
            assert (tree.left, tree.right) == (genealogy.breakpoints[index], genealogy.breakpoints[index + 1])
            parent = tree.parent
            n_children = np.bincount(parent[parent >= 0], minlength=parent.size)
            roots = np.flatnonzero((parent < 0) & (n_children > 0))
            internal = np.flatnonzero(n_children > 0)
            assert roots.size == 1 and internal.size == n - 1
            assert (n_children[internal] == 2).all()  # binary: a node is in a tree with both its children or not at all
            assert (parent[:n] >= 0).all() and (n_children[:n] == 0).all()
            for leaf in range(n):  # every sample reaches the root
                node, steps = leaf, 0
                while parent[node] >= 0:
                    assert times[parent[node]] > times[node]
                    node, steps = int(parent[node]), steps + 1
                assert node == roots[0] and steps < parent.size
            lengths = times[parent[parent >= 0]] - times[np.flatnonzero(parent >= 0)]
            assert tree.total_length == pytest.approx(lengths.sum(), rel=1e-9)
            assert genealogy.tree_lengths[index] == pytest.approx(tree.total_length, rel=1e-12)
            assert tree.tmrca == times[roots[0]] == genealogy.tree_tmrca[index]
        for before, after in zip(trees, trees[1:]):  # a breakpoint changes the tree
            assert not np.array_equal(before.parent, after.parent)

        positions = [0.0, length / 3, float(np.nextafter(length, 0.0))]
        if genealogy.n_trees > 1:
            positions.append(float(genealogy.breakpoints[1]))  # a boundary belongs to the tree on its right
        for position in positions:
            tree = genealogy.tree_at(position)
            assert tree.left <= position < tree.right
        with pytest.raises(ValueError, match="position"):
            genealogy.tree_at(length)


def test_a_recombining_sample_has_many_trees_and_a_non_recombining_one_has_one() -> None:
    flat = simulate_recombinant_genealogy(12, 500.0, rng=1)
    assert (flat.n_trees, flat.n_recombinations, flat.n_gap_recombinations, flat.n_coalescences) == (1, 0, 0, 11)
    assert flat.rho == 0.0 and flat.breakpoints.tolist() == [0.0, 500.0]
    busy = simulate_recombinant_genealogy(12, 500.0, rho=60.0, rng=1)
    assert busy.n_trees > 20 and busy.n_ancestral_recombinations >= busy.n_trees - 1
    assert busy.n_recombinations == busy.n_ancestral_recombinations + busy.n_gap_recombinations
    assert busy.n_gap_recombinations > 0 and busy.n_coalescences >= 11


def test_seeds_make_the_graph_reproducible_for_every_kind_of_source() -> None:
    for make in (lambda: 5, lambda: random.Random(5), lambda: np.random.default_rng(5)):
        first = simulate_recombinant_haplotypes(9, 300.0, theta=6.0, rho=10.0, rng=make())
        second = simulate_recombinant_haplotypes(9, 300.0, theta=6.0, rho=10.0, rng=make())
        assert np.array_equal(first.haplotypes, second.haplotypes) and np.array_equal(first.positions, second.positions)
        assert np.array_equal(first.genealogy.edge_child, second.genealogy.edge_child)
    other = simulate_recombinant_haplotypes(9, 300.0, theta=6.0, rho=10.0, rng=6)
    assert not np.array_equal(other.positions, first.positions)


# --- merging and splitting lineages ------------------------------------------------------------------------


def test_merging_lineages_coalesces_overlaps_and_passes_the_rest_on() -> None:
    node_times = [0.0] * 4
    edges: list = []
    merged = _merge_lineages([(0.0, 6.0, 0, 1)], [(4.0, 10.0, 1, 1)], 4, 0.7, node_times, edges)
    assert merged == [(0.0, 4.0, 0, 1), (4.0, 6.0, 4, 2), (6.0, 10.0, 1, 1)]
    assert node_times == [0.0] * 4 + [0.7] and edges == [(4.0, 6.0, 4, 0), (4.0, 6.0, 4, 1)]

    # One event makes one new node, however many stretches overlap
    node_times, edges = [0.0] * 4, []
    merged = _merge_lineages([(0.0, 2.0, 0, 1), (4.0, 6.0, 0, 1)], [(1.0, 5.0, 1, 1)], 4, 1.5, node_times, edges)
    assert merged == [(0.0, 1.0, 0, 1), (1.0, 2.0, 4, 2), (2.0, 4.0, 1, 1), (4.0, 5.0, 4, 2), (5.0, 6.0, 0, 1)]
    assert len(node_times) == 5 and {parent for _, _, parent, _ in edges} == {4} and len(edges) == 4


def test_merging_lineages_drops_what_has_reached_its_common_ancestor() -> None:
    node_times, edges = [0.0] * 2, []
    assert _merge_lineages([(0.0, 5.0, 0, 1)], [(0.0, 5.0, 1, 1)], 2, 0.3, node_times, edges) == []
    assert edges == [(0.0, 5.0, 2, 0), (0.0, 5.0, 2, 1)]  # the root coalescence is still an edge

    # n = 3: where one sample meets two the stretch is complete and goes, where one meets one it stays
    node_times, edges = [0.0] * 3, []
    merged = _merge_lineages([(0.0, 5.0, 0, 1)], [(2.0, 8.0, 1, 2)], 3, 0.3, node_times, edges)
    assert merged == [(0.0, 2.0, 0, 1), (5.0, 8.0, 1, 2)]
    assert edges == [(2.0, 5.0, 3, 0), (2.0, 5.0, 3, 1)] and node_times == [0.0, 0.0, 0.0, 0.3]
    node_times, edges = [0.0] * 4, []
    merged = _merge_lineages([(0.0, 5.0, 0, 1)], [(2.0, 8.0, 1, 2)], 4, 0.3, node_times, edges)
    assert merged == [(0.0, 2.0, 0, 1), (2.0, 5.0, 4, 3), (5.0, 8.0, 1, 2)]


def test_merging_lineages_without_common_material_creates_no_node() -> None:
    node_times, edges = [0.0] * 4, []
    merged = _merge_lineages([(0.0, 3.0, 0, 1)], [(5.0, 8.0, 1, 1)], 4, 0.9, node_times, edges)
    assert merged == [(0.0, 3.0, 0, 1), (5.0, 8.0, 1, 1)] and node_times == [0.0] * 4 and edges == []
    # Neighbouring stretches of one node with equal sample counts are one stretch
    assert _merge_lineages([(0.0, 3.0, 2, 2)], [(3.0, 6.0, 2, 2)], 5, 0.9, node_times, edges) == [(0.0, 6.0, 2, 2)]
    assert _merge_lineages([(0.0, 3.0, 2, 2)], [(3.0, 6.0, 2, 3)], 5, 0.9, node_times, edges) == [
        (0.0, 3.0, 2, 2),
        (3.0, 6.0, 2, 3),
    ]


def test_splitting_a_lineage_tells_ancestral_cuts_from_cuts_in_a_gap() -> None:
    lineage = [(0.0, 4.0, 7, 1), (6.0, 9.0, 7, 1)]
    assert _split_lineage(lineage, 2.0) == ([(0.0, 2.0, 7, 1)], [(2.0, 4.0, 7, 1), (6.0, 9.0, 7, 1)], True)
    assert _split_lineage(lineage, 5.0) == ([(0.0, 4.0, 7, 1)], [(6.0, 9.0, 7, 1)], False)
    assert _split_lineage(lineage, 7.5) == ([(0.0, 4.0, 7, 1), (6.0, 7.5, 7, 1)], [(7.5, 9.0, 7, 1)], True)


@pytest.mark.parametrize("name", ["bottleneck", "growth"])
def test_the_scalar_hazard_maps_agree_with_the_vectorized_one(name: str) -> None:
    history = _SizeHistory(BOTTLENECK if name == "bottleneck" else GROWTH)
    hazards = np.array([0.0, 1e-4, 0.05, 0.7, 2.0, 5.0, 40.0, 300.0])
    times = history.time_at_hazard(hazards)
    for hazard, time in zip(hazards.tolist(), times.tolist()):
        assert history.time_at(hazard) == pytest.approx(time, rel=1e-12, abs=1e-15)
        assert history.hazard_at(time) == pytest.approx(hazard, rel=1e-9, abs=1e-12)


# --- without recombination: the non-recombining engine -----------------------------------------------------


def test_without_recombination_the_graph_is_one_tree_distributed_as_the_non_recombining_engine() -> None:
    n, theta, reps = 10, 8.0, 4000
    generator = np.random.default_rng(1)
    spectra, tmrca = np.zeros((reps, n - 1)), np.zeros(reps)
    for r in range(reps):
        sample = simulate_recombinant_haplotypes(n, 500.0, theta=theta, rho=0.0, rng=generator)
        assert sample.genealogy.n_trees == 1 and sample.genealogy.n_recombinations == 0
        spectra[r], tmrca[r] = _spectrum(sample.haplotypes), sample.genealogy.tree_tmrca[0]
    reference = simulate_coalescent_spectra(n, 40_000, theta=theta, rng=2)

    _assert_mean(tmrca, 2 * (1 - 1 / n), "E[TMRCA] = 2 (1 - 1/n)")
    for i in range(n - 1):
        se = math.hypot(_mean_se(reference[:, i])[1], 0.0)
        _assert_mean(spectra[:, i], reference[:, i].mean(), f"E[xi_{i + 1}]", se)
        _assert_mean(spectra[:, i], theta / (i + 1), f"E[xi_{i + 1}] = theta / {i + 1}")
    sites, ref_sites = spectra.sum(axis=1), reference.sum(axis=1)
    assert stats.ks_2samp(sites, ref_sites).pvalue > 1e-3
    assert stats.ks_2samp(_pairwise(spectra, n), _pairwise(reference, n)).pvalue > 1e-3
    a = _a(n)
    squares = (sites - sites.mean()) ** 2
    assert abs(squares.mean() - (theta * a + theta**2 * sum(1 / i**2 for i in range(1, n)))) <= (
        Z * squares.std(ddof=1) / math.sqrt(reps)
    )


def test_every_pair_of_sites_of_a_non_recombining_sample_passes_the_four_gamete_test() -> None:
    for seed in range(15):
        sample = simulate_recombinant_haplotypes(14, 100.0, theta=9.0, rho=0.0, rng=seed)
        h = sample.haplotypes
        for a, b in itertools.combinations(range(h.shape[1]), 2):
            assert len(set(2 * h[:, a] + h[:, b])) <= 3


# --- Hudson's expected number of recombination events ------------------------------------------------------


@pytest.mark.parametrize(
    "n, rho, length, reps, seed",
    [(2, 10.0, 1e6, 20_000, 1), (2, 10.0, 100.0, 20_000, 2), (5, 10.0, 1e6, 6000, 3), (10, 20.0, 1e6, 3000, 4)],
)
def test_ancestral_breakpoints_have_hudsons_expectation(
    n: int, rho: float, length: float, reps: int, seed: int
) -> None:
    """E[R] = rho * sum_{i<n} 1/i for the recombination events whose breakpoint fell in ancestral material."""
    generator = np.random.default_rng(seed)
    graphs = [simulate_recombinant_genealogy(n, length, rho=rho, rng=generator) for _ in range(reps)]
    ancestral = np.array([g.n_ancestral_recombinations for g in graphs], dtype=float)
    _assert_mean(ancestral, rho * _a(n), f"n={n}, rho={rho}: ancestral recombination events")
    total = np.array([g.n_recombinations for g in graphs], dtype=float)
    assert (total >= ancestral).all() and total.mean() > ancestral.mean()  # events in gaps come on top


def test_high_recombination_adds_events_in_stretches_that_are_not_ancestral() -> None:
    """The same expectation holds for ancestral breakpoints, while the total grows faster than linearly in rho."""
    n, rho, reps = 10, 100.0, 300
    generator = np.random.default_rng(5)
    graphs = [simulate_recombinant_genealogy(n, 1e6, rho=rho, rng=generator) for _ in range(reps)]
    ancestral = np.array([g.n_ancestral_recombinations for g in graphs], dtype=float)
    total = np.array([g.n_recombinations for g in graphs], dtype=float)
    _assert_mean(ancestral, rho * _a(n), "ancestral recombination events at rho = 100")
    assert total.mean() > 1.8 * ancestral.mean()


@pytest.mark.parametrize("epochs", [BOTTLENECK, GROWTH], ids=["bottleneck", "growth"])
def test_ancestral_breakpoints_follow_the_tree_length_under_a_demography(epochs: list[PopulationEpoch]) -> None:
    """E[R] = (rho / 2) E[tree length] = rho E[S] / theta, with E[S] from the non-recombining engine."""
    n, rho, reps = 8, 15.0, 3000
    generator = np.random.default_rng(6)
    counts = np.array(
        [
            simulate_recombinant_genealogy(n, 1e5, rho=rho, epochs=epochs, rng=generator).n_ancestral_recombinations
            for _ in range(reps)
        ],
        dtype=float,
    )
    segregating = simulate_coalescent_spectra(n, 40_000, theta=1.0, epochs=epochs, rng=7).sum(axis=1)
    mean_s, se_s = _mean_se(segregating)
    _assert_mean(counts, rho * mean_s, "ancestral recombination events", se_extra=rho * se_s)


def test_per_site_and_total_rates_describe_the_same_process() -> None:
    by_site = simulate_recombinant_genealogy(8, 400.0, rho_per_site=0.05, rng=3)
    total = simulate_recombinant_genealogy(8, 400.0, rho=20.0, rng=3)
    assert by_site.rho == pytest.approx(20.0) and total.rho == 20.0
    assert np.array_equal(by_site.edge_left, total.edge_left) and np.array_equal(by_site.node_times, total.node_times)


# --- exact moments of the two-sample ancestral recombination graph -----------------------------------------


def _exact_pair_moments(link_rates: list[float]) -> tuple[np.ndarray, np.ndarray]:
    """Exact E[T_i T_j] and E[T_i] of the pairwise coalescence times at consecutive loci, for two samples.

    The loci are joined by links with population recombination rate ``link_rates[j]`` (rate ``rate / 2`` per
    lineage that carries material on both sides of the link). A state is the multiset of lineages, each the set
    (bit mask) of loci it carries ancestry at; the chain is solved by first-step analysis on the remaining
    times to the two loci, with no use of the simulator.
    """
    m = len(link_rates) + 1
    start = (2**m - 1,) * 2
    transitions: dict[tuple[int, ...], list[tuple[float, tuple[int, ...]]]] = {}
    stack = [start]
    while stack:
        state = stack.pop()
        if state in transitions:
            continue
        moves: list[tuple[float, tuple[int, ...]]] = []
        for a, b in itertools.combinations(range(len(state)), 2):
            rest = [state[i] for i in range(len(state)) if i not in (a, b)]
            merged = state[a] ^ state[b]  # loci carried by exactly one of the two are passed on
            moves.append((1.0, tuple(sorted(rest + ([merged] if merged else [])))))
        for a, lineage in enumerate(state):
            low, high = (lineage & -lineage).bit_length() - 1, lineage.bit_length() - 1
            for link in range(low, high):
                left, right = lineage & (2 ** (link + 1) - 1), lineage >> (link + 1) << (link + 1)
                if left and right:
                    rest = [state[i] for i in range(len(state)) if i != a]
                    moves.append((link_rates[link] / 2.0, tuple(sorted(rest + [left, right]))))
        transitions[state] = moves
        stack.extend(next_state for _, next_state in moves if next_state and next_state not in transitions)
    index = {state: i for i, state in enumerate(transitions)}

    def carried(state: tuple[int, ...]) -> int:
        return int(np.bitwise_or.reduce(state)) if state else 0

    size = len(index)
    first = np.zeros((size, m))
    for locus in range(m):
        system, rhs = np.eye(size), np.zeros(size)
        for state, moves in transitions.items():
            if not (carried(state) >> locus) & 1:
                continue
            total = sum(rate for rate, _ in moves)
            rhs[index[state]] = 1 / total
            for rate, nxt in moves:
                if nxt and (carried(nxt) >> locus) & 1:
                    system[index[state], index[nxt]] -= rate / total
        first[:, locus] = np.linalg.solve(system, rhs)
    second = np.zeros((m, m))
    for i, j in itertools.product(range(m), repeat=2):
        system, rhs = np.eye(size), np.zeros(size)
        for state, moves in transitions.items():
            if not ((carried(state) >> i) & 1 and (carried(state) >> j) & 1):
                continue
            total = sum(rate for rate, _ in moves)
            row = index[state]
            rhs[row] = 2 / total**2
            for rate, nxt in moves:
                weight = rate / total
                for locus in (i, j):
                    rhs[row] += (
                        weight / total * (first[index[nxt], locus] if nxt and (carried(nxt) >> locus) & 1 else 0)
                    )
                if nxt and (carried(nxt) >> i) & 1 and (carried(nxt) >> j) & 1:
                    system[row, index[nxt]] -= weight
        second[i, j] = np.linalg.solve(system, rhs)[index[start]]
    return second, first[index[start]]


@pytest.mark.parametrize("rho", [0.0, 1.0, 5.0, 20.0])
def test_the_reference_chain_reproduces_hudsons_two_locus_covariance(rho: float) -> None:
    """Cov(T_1, T_2) = (rho + 18) / (rho^2 + 13 rho + 18) in units of 2N generations (Hudson 1983)."""
    second, first = _exact_pair_moments([rho])
    assert first.tolist() == pytest.approx([1.0, 1.0])
    assert second[0, 1] - first[0] * first[1] == pytest.approx((rho + 18) / (rho**2 + 13 * rho + 18), rel=1e-9)


@pytest.mark.parametrize(
    "probes, rho_per_site, length, reps, seed",
    [([10.0, 50.0, 110.0], 0.05, 120.0, 60_000, 11), ([10.0, 30.0, 50.0, 90.0], 0.1, 100.0, 40_000, 12)],
)
def test_two_sample_coalescence_times_follow_the_exact_moments(
    probes: list[float], rho_per_site: float, length: float, reps: int, seed: int
) -> None:
    """The times at several loci have the moments of the chain, including loci whose ancestry shares a lineage."""
    links = [rho_per_site * (b - a) for a, b in zip(probes, probes[1:])]
    second, first = _exact_pair_moments(links)
    generator = np.random.default_rng(seed)
    times = np.empty((reps, len(probes)))
    for r in range(reps):
        genealogy = simulate_recombinant_genealogy(2, length, rho_per_site=rho_per_site, rng=generator)
        times[r] = genealogy.tree_tmrca[np.searchsorted(genealogy.breakpoints, probes, side="right") - 1]
    for i in range(len(probes)):
        _assert_mean(times[:, i], first[i], f"E[T_{i}]")
        for j in range(i + 1, len(probes)):
            _assert_mean(times[:, i] * times[:, j], second[i, j], f"E[T_{i} T_{j}]")


# --- mutations placed on the marginal trees ----------------------------------------------------------------


def _replicate_spectra(
    n: int, theta: float, rho_per_site: float, length: float, reps: int, seed: int, epochs: list | None = None
) -> np.ndarray:
    """Spectra of ``reps`` recombining samples."""
    generator = np.random.default_rng(seed)
    spectra = np.zeros((reps, n - 1))
    for r in range(reps):
        sample = simulate_recombinant_haplotypes(
            n, length, theta=theta, rho_per_site=rho_per_site, epochs=epochs, rng=generator
        )
        spectra[r] = _spectrum(sample.haplotypes)
    return spectra


@pytest.mark.parametrize(
    "rho_per_site, length, reps, epochs, seed",
    [(0.02, 1000.0, 3000, None, 21), (0.5, 100.0, 1500, None, 22), (0.02, 1000.0, 3000, BOTTLENECK, 23)]
    + [(0.02, 1000.0, 3000, GROWTH, 24)],
    ids=["constant", "constant-high-rho", "bottleneck", "growth"],
)
def test_recombination_leaves_the_per_site_spectrum_diversity_and_segregating_sites_unchanged(
    rho_per_site: float, length: float, reps: int, epochs: list | None, seed: int
) -> None:
    """Each site has the Kingman genealogy, so E[xi_i], E[S] and E[pi] are those of the non-recombining engine."""
    n, theta = 10, 8.0
    spectra = _replicate_spectra(n, theta, rho_per_site, length, reps, seed, epochs)
    reference = simulate_coalescent_spectra(n, 40_000, theta=theta, epochs=epochs, rng=seed + 1)
    for i in range(n - 1):
        _assert_mean(spectra[:, i], reference[:, i].mean(), f"E[xi_{i + 1}]", _mean_se(reference[:, i])[1])
    for label, arg, ref in (
        ("S", spectra.sum(axis=1), reference.sum(axis=1)),
        ("pi", _pairwise(spectra, n), _pairwise(reference, n)),
    ):
        _assert_mean(arg, ref.mean(), f"E[{label}]", _mean_se(ref)[1])
    if epochs is None:
        _assert_mean(spectra.sum(axis=1), theta * _a(n), "E[S] = theta a_n")
        _assert_mean(_pairwise(spectra, n), theta, "E[pi] = theta")


def test_recombination_makes_the_sites_less_correlated_so_the_variance_of_s_falls() -> None:
    n, theta = 10, 8.0
    linked = _replicate_spectra(n, theta, 0.0, 400.0, 2000, 31)
    free = _replicate_spectra(n, theta, 0.05, 400.0, 2000, 32)  # rho = 20
    # Var(S) = theta a_n + theta^2 times the average covariance of the tree lengths at two sites
    assert free.sum(axis=1).var(ddof=1) < 0.6 * linked.sum(axis=1).var(ddof=1)
    assert free.sum(axis=1).var(ddof=1) > theta * _a(n)  # the Poisson variance of mutation is a floor


def test_fixed_numbers_of_mutations_are_placed_at_sorted_distinct_polymorphic_positions() -> None:
    sample = simulate_recombinant_haplotypes(12, 5000.0, segregating_sites=300, rho=40.0, rng=3)
    h, positions = sample.haplotypes, sample.positions
    assert h.shape == (12, 300) and h.dtype == np.uint8 and set(np.unique(h)) <= {0, 1}
    assert positions.shape == (300,) and (np.diff(positions) > 0).all()
    assert positions[0] >= 0.0 and positions[-1] < 5000.0
    carriers = h.sum(axis=0)
    assert ((carriers >= 1) & (carriers <= 11)).all()  # every column is polymorphic
    assert np.unique(h, axis=1).shape[1] > 20  # many different mutations, not copies of one


def test_mutations_fall_uniformly_along_the_sequence_whatever_the_recombination_rate() -> None:
    """Position of one mutation: its density is proportional to E[tree length], which is the same everywhere."""
    for rho, seed in ((0.0, 41), (10.0, 42)):
        generator = np.random.default_rng(seed)
        positions = np.array(
            [
                simulate_recombinant_haplotypes(8, 10.0, segregating_sites=1, rho=rho, rng=generator).positions[0]
                for _ in range(2500)
            ]
        )
        assert stats.kstest(positions / 10.0, "uniform").pvalue > 1e-3


def test_mutations_follow_branch_length_between_trees() -> None:
    """With S = 1 the tree that holds the mutation is chosen in proportion to width * total length."""
    genealogy = simulate_recombinant_genealogy(8, 100.0, rho=40.0, rng=7)
    masses = genealogy.tree_lengths * np.diff(genealogy.breakpoints)
    assert genealogy.n_trees > 10
    generator = np.random.default_rng(8)
    reps = 20_000
    where = np.array(
        [
            np.searchsorted(genealogy.breakpoints, genealogy.mutate(segregating_sites=1, rng=generator)[0][0], "right")
            - 1
            for _ in range(reps)
        ]
    )
    observed = np.bincount(where, minlength=genealogy.n_trees)
    expected = reps * masses / masses.sum()
    keep = expected >= 5
    chi2 = float((((observed[keep] - expected[keep]) ** 2) / expected[keep]).sum())
    assert chi2 < stats.chi2.ppf(1 - 1e-3, keep.sum() - 1)


def test_only_sites_of_one_tree_are_guaranteed_the_four_gamete_test() -> None:
    sample = simulate_recombinant_haplotypes(16, 1000.0, segregating_sites=60, rho=60.0, rng=9)
    h, genealogy = sample.haplotypes, sample.genealogy
    tree_of = np.searchsorted(genealogy.breakpoints, sample.positions, side="right") - 1
    within = between = 0
    for a, b in itertools.combinations(range(h.shape[1]), 2):
        four = len(set(2 * h[:, a] + h[:, b])) == 4
        if tree_of[a] == tree_of[b]:
            assert not four
            within += 1
        else:
            between += four
    assert within > 0 and between > 0  # recombination does break the infinite-sites compatibility between trees


def test_empty_mutation_requests_give_empty_matrices() -> None:
    genealogy = simulate_recombinant_genealogy(6, 100.0, rho=5.0, rng=1)
    for kwargs in ({"theta": 0.0}, {"segregating_sites": 0}):
        positions, haplotypes = genealogy.mutate(rng=1, **kwargs)
        assert positions.shape == (0,) and haplotypes.shape == (6, 0)
    assert isinstance(simulate_recombinant_haplotypes(6, 100.0, theta=0.0, rng=1), RecombinantSample)


# --- linkage disequilibrium ---------------------------------------------------------------------------------


def _ohta_kimura(c: np.ndarray) -> np.ndarray:
    """Ohta and Kimura's (1971) E[D^2] / E[p (1 - p) q (1 - q)] at scaled recombination distance c = 4 N r d."""
    return (10 + c) / ((2 + c) * (11 + c))


def _falling(x: np.ndarray, k: int) -> np.ndarray:
    result = np.ones_like(x)
    for i in range(k):
        result = result * (x - i)
    return result


def _unbiased_ld_terms(haplotypes: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Index pairs, and unbiased estimates of D^2 and p_A p_a p_B p_b, for every pair of sites.

    The estimators are the U-statistics that average over the ordered draws of distinct sequences, so their
    expectations do not depend on the sample size: E[D^2] = E[(p_AB p_ab - p_Ab p_aB)^2] and
    E[p (1 - p) q (1 - q)] are those of the population.
    """
    n = haplotypes.shape[0]
    h = haplotypes.astype(float)
    totals = h.sum(axis=0)
    rows, cols = np.triu_indices(h.shape[1], 1)
    a, b, both = totals[rows], totals[cols], (h.T @ h)[rows, cols]
    counts = {"AB": both, "Ab": a - both, "aB": b - both, "ab": n - a - b + both}
    ordered = n * (n - 1) * (n - 2) * (n - 3)
    d2 = (
        _falling(counts["AB"], 2) * _falling(counts["ab"], 2)
        - 2 * counts["AB"] * counts["ab"] * counts["Ab"] * counts["aB"]
        + _falling(counts["Ab"], 2) * _falling(counts["aB"], 2)
    ) / ordered
    pqpq = np.zeros_like(both)
    for picks in itertools.product(("AB", "Ab"), ("aB", "ab"), ("AB", "aB"), ("Ab", "ab")):
        term = np.ones_like(both)
        for haplotype in set(picks):
            term = term * _falling(counts[haplotype], picks.count(haplotype))
        pqpq += term
    return rows, cols, d2, pqpq / ordered


LD_BINS = np.array([0, 0.5, 1, 2, 4, 7, 12, 20, 35, 60.0])  # scaled distances rho_site * d


@pytest.mark.parametrize("n, reps, seed", [(20, 700, 51), (60, 500, 52)])
def test_linkage_disequilibrium_decays_as_ohta_and_kimura_expect(n: int, reps: int, seed: int) -> None:
    """E[D^2] / E[p (1 - p) q (1 - q)] at scaled distance rho_d is (10 + rho_d) / ((2 + rho_d) (11 + rho_d)).

    The unbiased estimators make the check independent of the sample size. For every distance bin the excess of
    the summed D^2 over the Ohta-Kimura prediction (rho_d, pair by pair, times the summed p q p' q') has a mean of
    zero, judged by its standard error over independent replicates. The prediction is itself a closure
    approximation (measured over 3000 replicates at n = 20 and 50 and 1500 at n = 100: within 0.003 of the
    simulated ratio in every bin, against standard errors of 0.0005 to 0.002).
    """
    rho_per_site, length, theta = 0.05, 1000.0, 40.0
    generator = np.random.default_rng(seed)
    excess = np.zeros((reps, LD_BINS.size - 1))  # summed D^2 - OK * (p q p' q') per replicate and bin
    mass = np.zeros_like(excess)
    for r in range(reps):
        sample = simulate_recombinant_haplotypes(n, length, theta=theta, rho_per_site=rho_per_site, rng=generator)
        rows, cols, d2, pqpq = _unbiased_ld_terms(sample.haplotypes)
        distance = rho_per_site * (sample.positions[cols] - sample.positions[rows])
        which = np.searchsorted(LD_BINS, distance, side="right") - 1
        for k in range(LD_BINS.size - 1):
            chosen = which == k
            excess[r, k] = (d2[chosen] - _ohta_kimura(distance[chosen]) * pqpq[chosen]).sum()
            mass[r, k] = pqpq[chosen].sum()
    for k in range(LD_BINS.size - 1):
        ratio_error = excess[:, k].sum() / mass[:, k].sum()
        se = excess[:, k].std(ddof=1) * math.sqrt(reps) / mass[:, k].sum()
        assert abs(ratio_error) <= Z * se, (
            f"rho_d in [{LD_BINS[k]}, {LD_BINS[k + 1]}): D^2 exceeds the prediction by {ratio_error:+.4f} "
            f"of the p q p' q' scale (standard error {se:.4f})"
        )


def test_without_recombination_linkage_disequilibrium_is_five_elevenths_at_every_distance() -> None:
    """Ohta and Kimura's value at rho_d = 0 is exactly 5/11; the distance of two sites then means nothing."""
    n, reps = 20, 2000
    generator = np.random.default_rng(53)
    numerator = np.zeros(reps)
    denominator = np.zeros(reps)
    for r in range(reps):
        sample = simulate_recombinant_haplotypes(n, 100.0, theta=20.0, rho=0.0, rng=generator)
        if sample.haplotypes.shape[1] < 2:
            continue
        _, _, d2, pqpq = _unbiased_ld_terms(sample.haplotypes)
        numerator[r], denominator[r] = d2.sum(), pqpq.sum()
    ratio = numerator.sum() / denominator.sum()
    se = (numerator - ratio * denominator).std(ddof=1) * math.sqrt(reps) / denominator.sum()
    assert abs(ratio - 5 / 11) <= Z * se, f"ratio {ratio:.4f} vs 5/11 (standard error {se:.4f})"


def test_mean_r_squared_falls_with_distance_when_the_sequence_recombines_and_not_otherwise() -> None:
    def mean_r2_by_distance(rho_per_site: float, seed: int) -> np.ndarray:
        generator = np.random.default_rng(seed)
        edges = np.array([0, 40, 120, 300, 1000.0])  # sites
        total, count = np.zeros(4), np.zeros(4)
        for _ in range(300):
            sample = simulate_recombinant_haplotypes(30, 1000.0, theta=30.0, rho_per_site=rho_per_site, rng=generator)
            h = sample.haplotypes.astype(float)
            keep = (h.mean(axis=0) >= 0.1) & (h.mean(axis=0) <= 0.9)
            h, pos = h[:, keep], sample.positions[keep]
            if h.shape[1] < 2:
                continue
            r2 = np.corrcoef(h.T)[np.triu_indices(h.shape[1], 1)] ** 2
            gap = (pos[None, :] - pos[:, None])[np.triu_indices(h.shape[1], 1)]
            which = np.searchsorted(edges, gap, side="right") - 1
            for k in range(4):
                total[k] += r2[which == k].sum()
                count[k] += (which == k).sum()
        return total / count

    decaying = mean_r2_by_distance(0.05, 61)
    assert (np.diff(decaying) < 0).all() and decaying[0] > 3 * decaying[-1]  # measured: 0.36, 0.21, 0.12, 0.08
    flat = mean_r2_by_distance(0.0, 62)  # one tree: r^2 does not depend on how far apart two sites are
    assert flat.max() - flat.min() < 0.1 * flat.mean()


# --- scale --------------------------------------------------------------------------------------------------


def test_a_hundred_sequences_over_ten_thousand_sites_with_realistic_recombination_are_cheap() -> None:
    sample = simulate_recombinant_haplotypes(100, 10_000.0, segregating_sites=10_000, rho=400.0, rng=1)
    assert sample.haplotypes.shape == (100, 10_000)
    assert sample.genealogy.n_trees > 500
    assert sample.genealogy.n_ancestral_recombinations > 0.5 * 400 * _a(100)


# --- input validation ---------------------------------------------------------------------------------------


def test_invalid_requests_are_rejected() -> None:
    with pytest.raises(ValidationError, match="n_samples"):
        simulate_recombinant_genealogy(1, 100.0, rho=1.0)
    for length in (0.0, -5.0, float("inf"), float("nan"), "100"):
        with pytest.raises(ValueError, match="sequence_length"):
            simulate_recombinant_genealogy(5, length, rho=1.0)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="at most one"):
        simulate_recombinant_genealogy(5, 100.0, rho=1.0, rho_per_site=0.01)
    for name in ("rho", "rho_per_site"):
        with pytest.raises(ValueError, match=name):
            simulate_recombinant_genealogy(5, 100.0, **{name: -1.0})
        with pytest.raises(ValueError, match="finite"):
            simulate_recombinant_genealogy(5, 100.0, **{name: float("nan")})
    with pytest.raises(ValueError, match="overflows"):
        simulate_recombinant_genealogy(5, 1e300, rho_per_site=1e300)
    with pytest.raises(ValueError, match="start at time 0"):
        simulate_recombinant_genealogy(5, 100.0, rho=1.0, epochs=[PopulationEpoch(0.1, 1.0)])
    genealogy = simulate_recombinant_genealogy(5, 100.0, rho=1.0, rng=1)
    with pytest.raises(ValueError, match="exactly one"):
        genealogy.mutate()
    with pytest.raises(ValueError, match="exactly one"):
        genealogy.mutate(theta=1.0, segregating_sites=2)
    with pytest.raises(ValueError, match="non-negative integer"):
        genealogy.mutate(segregating_sites=-1)
    with pytest.raises(ValueError, match="non-negative integer"):
        genealogy.mutate(segregating_sites=2.5)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="finite"):
        genealogy.mutate(theta=float("inf"))
    with pytest.raises(ValueError, match="exactly one"):
        simulate_recombinant_haplotypes(5, 100.0)
    assert isinstance(genealogy, RecombinantGenealogy)


# --- generate_linkage_disequilibrium_data(model="coalescent") -----------------------------------------------


def _ld_sample(rate: float, seed: int, n: int = 40, snps: int = 80, length: int = 8000, **kwargs) -> np.ndarray:
    matrix = generate_linkage_disequilibrium_data(
        n, snps, model="coalescent", recombination_rate=rate, sequence_length=length, rng=random.Random(seed), **kwargs
    )
    return np.array(matrix, dtype=float)


def _mean_r_squared(rate: float, seed: int, reps: int = 40) -> tuple[float, float, float]:
    """Mean over replicates of the mean r^2 of SNP pairs with a minor allele frequency of at least 0.1."""
    values, frequencies = [], []
    for r in range(reps):
        matrix = _ld_sample(rate, seed + r)
        frequency = matrix.mean(axis=0)
        frequencies.extend(frequency.tolist())
        common = matrix[:, (frequency >= 0.1) & (frequency <= 0.9)]
        correlation = np.corrcoef(common.T)[np.triu_indices(common.shape[1], 1)]
        values.append(float((correlation**2).mean()))
    mean, se = _mean_se(np.array(values))
    return mean, se, float(np.mean(frequencies))


class TestCoalescentLinkageDisequilibriumData:
    """The coalescent model of the generator: recombination_rate sets how fast LD decays."""

    def test_rows_are_individuals_and_columns_are_polymorphic_derived_alleles(self) -> None:
        matrix = generate_linkage_disequilibrium_data(25, 60, model="coalescent", rng=random.Random(1))
        assert len(matrix) == 25 and all(len(row) == 60 for row in matrix)
        assert all(type(allele) is int and allele in (0, 1) for row in matrix for allele in row)
        carriers = np.array(matrix).sum(axis=0)
        assert ((carriers >= 1) & (carriers <= 24)).all()  # exactly n_snps segregating sites
        assert matrix == generate_linkage_disequilibrium_data(25, 60, model="coalescent", rng=random.Random(1))
        assert matrix != generate_linkage_disequilibrium_data(25, 60, model="coalescent", rng=random.Random(2))

    def test_the_alias_for_the_number_of_sites_works(self) -> None:
        matrix = generate_linkage_disequilibrium_data(10, n_sites=7, model="coalescent", rng=random.Random(1))
        assert len(matrix[0]) == 7

    def test_allele_frequencies_follow_the_neutral_spectrum(self) -> None:
        """Mean derived allele frequency (n - 1) / (n a_n): rare derived alleles dominate, unlike the 0.5 default."""
        _, _, frequency = _mean_r_squared(1e-7, 300, reps=60)
        assert abs(frequency - 39 / (40 * _a(40))) < 0.012  # measured 0.229 against 0.229

    def test_recombination_rate_sets_the_decay_of_linkage_disequilibrium(self) -> None:
        """rho = 4 N_e r L = 0, 32 and 320 over 8 kb: the more it recombines, the lower the r^2 between SNPs."""
        (none, se_none, _), (some, se_some, _), (much, se_much, _) = (
            _mean_r_squared(rate, 100) for rate in (0.0, 1e-7, 1e-6)
        )
        assert none - some > Z * math.hypot(se_none, se_some)  # measured 0.43, 0.16, 0.054
        assert some - much > Z * math.hypot(se_some, se_much)

    def test_without_recombination_every_pair_of_snps_passes_the_four_gamete_test(self) -> None:
        for seed in range(6):
            matrix = _ld_sample(0.0, seed, n=20, snps=30).astype(int)
            for a, b in itertools.combinations(range(30), 2):
                assert len(set(2 * matrix[:, a] + matrix[:, b])) <= 3
        violating = 0
        for seed in range(6):
            matrix = _ld_sample(1e-6, seed, n=20, snps=30).astype(int)
            violating += sum(
                len(set(2 * matrix[:, a] + matrix[:, b])) == 4 for a, b in itertools.combinations(range(30), 2)
            )
        assert violating > 100

    def test_only_the_product_of_effective_size_rate_and_length_matters(self) -> None:
        """The population recombination parameter is 4 * effective_size * recombination_rate * sequence_length."""
        common = dict(model="coalescent", sequence_length=4000)
        first = generate_linkage_disequilibrium_data(
            12, 20, effective_size=8192.0, recombination_rate=2.0**-27, rng=random.Random(5), **common
        )
        second = generate_linkage_disequilibrium_data(
            12, 20, effective_size=4096.0, recombination_rate=2.0**-26, rng=random.Random(5), **common
        )
        third = generate_linkage_disequilibrium_data(
            12, 20, effective_size=4096.0, recombination_rate=2.0**-23, rng=random.Random(5), **common
        )
        assert first == second and first != third

    def test_ten_thousand_snps_of_a_hundred_sequences_use_the_default_region(self) -> None:
        matrix = generate_linkage_disequilibrium_data(100, 10_000, model="coalescent", rng=random.Random(3))
        assert len(matrix) == 100 and all(len(row) == 10_000 for row in matrix)

    def test_the_markov_model_stays_the_default_and_says_when_recombination_rate_is_ignored(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        plain = generate_linkage_disequilibrium_data(30, 8, r_squared_target=0.4, rng=random.Random(21))
        assert plain == generate_linkage_disequilibrium_data(
            30, 8, r_squared_target=0.4, model="markov", rng=random.Random(21)
        )
        with caplog.at_level("WARNING", logger=popgen.logger.name):
            generate_linkage_disequilibrium_data(30, 8, r_squared_target=0.4, rng=random.Random(21))
            assert not any("recombination_rate" in record.getMessage() for record in caplog.records)
            ignored = generate_linkage_disequilibrium_data(
                30, 8, r_squared_target=0.4, recombination_rate=0.2, rng=random.Random(21)
            )
        assert ignored == plain
        assert any("no effect" in record.getMessage() for record in caplog.records)

    @pytest.mark.parametrize(
        "options, message",
        [
            ({"r_squared_target": 0.3}, "r_squared_target"),
            ({"selection_coefficient": 0.1}, "selection_coefficient"),
            ({"allele_frequencies": [0.4, 0.6]}, "allele_frequencies"),
            ({"recombination_rate": -1e-8}, "recombination_rate"),
            ({"recombination_rate": float("nan")}, "recombination_rate"),
            ({"effective_size": 0.0}, "effective_size"),
            ({"effective_size": float("inf")}, "effective_size"),
            ({"sequence_length": 4}, "sequence_length"),
            ({"sequence_length": float("inf")}, "sequence_length"),
        ],
    )
    def test_the_coalescent_model_rejects_inputs_it_cannot_honour(self, options: dict, message: str) -> None:
        with pytest.raises(ValidationError, match=message):
            generate_linkage_disequilibrium_data(10, 5, model="coalescent", **options)

    def test_the_markov_model_rejects_coalescent_only_inputs_and_unknown_models(self) -> None:
        for options in ({"effective_size": 100.0}, {"sequence_length": 500}):
            with pytest.raises(ValidationError, match="coalescent"):
                generate_linkage_disequilibrium_data(10, 5, **options)
        with pytest.raises(ValidationError, match="model"):
            generate_linkage_disequilibrium_data(10, 5, model="ms")
        with pytest.raises(ValidationError, match="recombination_rate"):
            generate_linkage_disequilibrium_data(10, 5, recombination_rate=-1.0)

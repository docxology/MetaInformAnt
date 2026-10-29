"""Population genetics simulation utilities for generating synthetic genomic data.

This module provides functions for simulating various population genetic scenarios
including population bottlenecks, expansions, migration, and natural selection.
All functions support reproducible results through random seed control.

Sequence samples come from a seeded Kingman coalescent with infinite-sites mutation
(:func:`simulate_coalescent_haplotypes`, :func:`simulate_coalescent_spectra`), so they carry
the genealogical structure real samples have: one non-recombining genealogy, every variable
site a single mutation on one of its branches, and the neutral site frequency spectrum,
Tajima's D and diversity (or the bottleneck and expansion departures from them). Time runs
backwards from the sample in units of ``2 * N_ref`` generations, a pair of lineages coalesces at
rate ``1 / size(t)`` (``size`` relative to ``N_ref``) and mutations fall on every branch at rate
``theta / 2`` per unit time, where ``theta = 4 * N_ref * mu`` is the population mutation
parameter of the whole sequence. Size changes are described by :class:`PopulationEpoch` histories,
which rescale the coalescent rates through time exactly.

Recombination is added by Hudson's ancestral recombination graph
(:func:`simulate_recombinant_genealogy`, :func:`simulate_recombinant_haplotypes`): a lineage splits
in two at rate ``rho / 2`` per unit of the sequence it spans (``rho = 4 * N_ref * r``), every site keeps
the Kingman genealogy and the marginal trees along the sequence are correlated through the branches
they share, so linkage disequilibrium decays with distance. With ``rho = 0`` it is the non-recombining
engine above. :func:`generate_linkage_disequilibrium_data` offers it as ``model="coalescent"``.
"""

from __future__ import annotations

import math
import random
from bisect import bisect_right
from dataclasses import dataclass
from itertools import accumulate
from typing import Any, Dict, Iterator, List, Sequence, Tuple, Union

import numpy as np

from metainformant.core.data import validation
from metainformant.core.utils import errors, logging
from metainformant.simulation.rng import coerce_and_seed, coerce_rng

logger = logging.get_logger(__name__)

#: Anything the coalescent functions accept as a source of randomness.
RandomSource = Union[random.Random, np.random.Generator, int, None]

_BASES = "ACGT"
_DEFAULT_THETA_PER_SITE = 0.01
#: Replicates times samples simulated per vectorised block (bounds memory for large requests).
_BLOCK_CELLS = 2_000_000


@dataclass(frozen=True)
class PopulationEpoch:
    """One epoch of a piecewise-exponential demographic history, read backwards from the sample.

    Time ``t`` runs backwards from the sampling time (``t = 0``) in units of ``2 * N_ref``
    generations, and sizes are relative to the reference size ``N_ref`` that defines
    ``theta = 4 * N_ref * mu``. The epoch begins at ``start`` and lasts until the next epoch
    begins (for ever, in the last one). Inside it the relative size at time ``t`` is
    ``size * exp(-growth_rate * (t - start))``, so ``growth_rate`` is the exponential growth rate
    forwards in time: 0 keeps the size constant, a positive rate describes a population that was
    smaller further back (an expansion) and a negative rate one that was larger (a decline).

    A bottleneck of relative size ``b`` that ended ``t_end`` and began ``t_start`` before the
    sample is ``[PopulationEpoch(0, 1), PopulationEpoch(t_end, b), PopulationEpoch(t_start, 1)]``;
    an instantaneous expansion by a factor ``f`` that happened ``t`` ago is
    ``[PopulationEpoch(0, f), PopulationEpoch(t, 1)]``; and an exponential expansion by ``f``
    that took time ``T`` and ended at the sample is
    ``[PopulationEpoch(0, f, math.log(f) / T), PopulationEpoch(T, 1)]``.

    Attributes:
        start: Backward time at which the epoch begins; the first epoch starts at 0
        size: Population size at ``start`` relative to ``N_ref`` (positive)
        growth_rate: Forward-in-time exponential growth rate per unit time (default 0)
    """

    start: float
    size: float
    growth_rate: float = 0.0


class _SizeHistory:
    """Integrated coalescent hazard of a piecewise-exponential size history, and its inverse.

    A pair of lineages coalesces at rate ``1 / size(t)``, so ``k`` lineages next coalesce once the
    integrated hazard ``Lambda(t)``, the integral of ``1 / size``, has grown by an
    ``Exp(1) / C(k, 2)`` variate. Event times of the constant-size coalescent, mapped through
    ``Lambda^-1``, therefore rescale the coalescent rates through time exactly, whatever the
    epoch layout.
    """

    def __init__(self, epochs: Sequence[PopulationEpoch] | None) -> None:
        history = list(epochs) if epochs else [PopulationEpoch(0.0, 1.0)]
        starts = np.array([epoch.start for epoch in history], dtype=float)
        sizes = np.array([epoch.size for epoch in history], dtype=float)
        rates = np.array([epoch.growth_rate for epoch in history], dtype=float)
        if not (np.isfinite(starts).all() and np.isfinite(sizes).all() and np.isfinite(rates).all()):
            raise ValueError("epoch start, size and growth_rate must be finite numbers")
        if starts[0] != 0.0:
            raise ValueError("the first epoch must start at time 0")
        if (np.diff(starts) <= 0.0).any():
            raise ValueError("epoch start times must be strictly increasing")
        if (sizes <= 0.0).any():
            raise ValueError("epoch sizes must be positive")
        if rates[-1] < 0.0:
            raise ValueError("the last epoch needs growth_rate >= 0, or lineages may never coalesce")
        self._starts, self._sizes, self._rates = starts, sizes, rates
        hazards = [0.0]
        for index in range(len(history) - 1):
            hazards.append(hazards[-1] + self._epoch_hazard(index, float(starts[index + 1] - starts[index])))
        if not math.isfinite(hazards[-1]):
            raise ValueError("epoch growth rates are too large: the integrated coalescent hazard overflows")
        self._hazard_starts = np.array(hazards)
        # Plain-float copies for the event-by-event simulations, which map one time or hazard at a time.
        self._start_list = starts.tolist()
        self._hazard_list = hazards
        self._size_list = sizes.tolist()
        self._rate_list = rates.tolist()

    def _epoch_hazard(self, index: int, elapsed: float) -> float:
        """Hazard accumulated over ``elapsed`` time inside epoch ``index``."""
        size, rate = float(self._sizes[index]), float(self._rates[index])
        if rate == 0.0:
            return elapsed / size
        try:
            return math.expm1(rate * elapsed) / (rate * size)
        except OverflowError:
            return math.inf

    def time_at_hazard(self, hazard: np.ndarray) -> np.ndarray:
        """Backward time at which the integrated hazard reaches ``hazard`` (elementwise)."""
        index = np.searchsorted(self._hazard_starts, hazard, side="right") - 1
        remaining = hazard - self._hazard_starts[index]
        size, rate = self._sizes[index], self._rates[index]
        constant = rate == 0.0
        safe_rate = np.where(constant, 1.0, rate)
        elapsed = np.where(constant, remaining * size, np.log1p(safe_rate * size * remaining) / safe_rate)
        return np.asarray(self._starts[index] + elapsed)

    def hazard_at(self, time: float) -> float:
        """Integrated coalescent hazard accumulated by backward time ``time`` (the inverse of :meth:`time_at`)."""
        index = bisect_right(self._start_list, time) - 1
        return float(self._hazard_list[index] + self._epoch_hazard(index, time - self._start_list[index]))

    def time_at(self, hazard: float) -> float:
        """Backward time at which the integrated hazard reaches ``hazard``; one value, without array overhead."""
        index = bisect_right(self._hazard_list, hazard) - 1
        remaining = hazard - self._hazard_list[index]
        size, rate = self._size_list[index], self._rate_list[index]
        if rate == 0.0:
            return float(self._start_list[index] + remaining * size)
        return float(self._start_list[index] + math.log1p(rate * size * remaining) / rate)


def _numpy_generator(rng: RandomSource) -> np.random.Generator:
    """Return a numpy Generator for ``rng``: as is, seeded from a ``random.Random``, or from an int/None seed."""
    if isinstance(rng, np.random.Generator):
        return rng
    if isinstance(rng, random.Random):
        return np.random.default_rng(rng.getrandbits(64))
    return np.random.default_rng(rng)


def _python_random(rng: RandomSource) -> random.Random:
    """Return a ``random.Random`` for ``rng``: as is, seeded from a Generator, or from an int/None seed."""
    if isinstance(rng, random.Random):
        return rng
    if isinstance(rng, np.random.Generator):
        return random.Random(int(rng.integers(2**62)))
    return random.Random(rng)


def _check_theta(theta: float) -> None:
    """Require a finite, non-negative mutation parameter."""
    if not math.isfinite(theta):
        raise ValueError(f"theta must be a finite number, got {theta}")
    validation.validate_range(theta, min_val=0.0, name="theta")


def _coalescent_event_times(
    n_samples: int, n_replicates: int, history: _SizeHistory, generator: np.random.Generator
) -> np.ndarray:
    """Backward times of the ``n_samples - 1`` coalescent events of every replicate (one row each).

    While ``k`` lineages remain the integrated hazard grows at rate ``C(k, 2)``, so the hazard at each
    event is a cumulative sum of ``Exp(1) / C(k, 2)`` variates, mapped to time by ``history``.
    """
    lineages = np.arange(n_samples, 1, -1, dtype=float)
    pairs = lineages * (lineages - 1.0) / 2.0
    hazard = np.cumsum(generator.standard_exponential((n_replicates, n_samples - 1)) / pairs, axis=1)
    return history.time_at_hazard(hazard)


def _branch_length_spectra(times: np.ndarray, n_samples: int, generator: np.random.Generator) -> np.ndarray:
    """Total branch length carried by ``i`` of the samples, ``i = 1..n_samples-1``, for every row of ``times``.

    The pair that merges at each event is uniform among the lineages present, independently of the
    event times, so the topology is tracked only through the number of samples below each lineage.
    """
    replicates = times.shape[0]
    rows = np.arange(replicates)
    width = n_samples + 1
    clade_sizes = np.ones((replicates, n_samples), dtype=np.int64)  # samples below each current lineage
    lengths = np.zeros((replicates, width))
    previous = np.zeros(replicates)
    for step, lineages in enumerate(range(n_samples, 1, -1)):
        elapsed = times[:, step] - previous
        previous = times[:, step]
        slots = (rows[:, None] * width + clade_sizes[:, :lineages]).ravel()
        weights = np.repeat(elapsed, lineages)
        lengths += np.bincount(slots, weights=weights, minlength=replicates * width).reshape(replicates, width)
        first = generator.integers(0, lineages, size=replicates)
        second = generator.integers(0, lineages - 1, size=replicates)
        second += second >= first
        clade_sizes[rows, first] += clade_sizes[rows, second]
        clade_sizes[rows, second] = clade_sizes[:, lineages - 1]
    return lengths[:, 1:n_samples]


def simulate_coalescent_spectra(
    n_samples: int,
    n_replicates: int,
    *,
    theta: float | None = None,
    segregating_sites: int | None = None,
    epochs: Sequence[PopulationEpoch] | None = None,
    rng: RandomSource = None,
) -> np.ndarray:
    """Simulate the site frequency spectra of many neutral coalescent samples.

    Each replicate draws a genealogy of ``n_samples`` sequences from the Kingman coalescent (constant
    size, or the demographic history ``epochs``) and places infinite-sites mutations on its branches.
    Exactly one of ``theta`` and ``segregating_sites`` fixes the mutation load:

    * ``theta``: mutations are Poisson with mean ``theta / 2`` per unit of branch length, so under
      constant size ``E[S] = theta * a_n`` and ``E[xi_i] = theta / i``, with ``a_n`` the sum of
      ``1 / i`` for ``i < n_samples``.
    * ``segregating_sites``: every replicate gets exactly that many mutations, each on a branch
      chosen with probability proportional to its length, on a genealogy drawn without weighting (the
      approach of Hudson's ``ms -s``). This is the distribution of the spectrum given the number of
      segregating sites when ``theta`` is integrated out with the scale-invariant prior ``1 / theta``,
      so no mutation rate has to be estimated or assumed. For a fixed ``theta`` the genealogy given
      ``S`` would be weighted by ``T**S * exp(-theta * T / 2)`` (``T`` the total branch length); the
      difference is small: p-values built on this null are calibrated to within sampling error from a
      few tens of mutations upward and conservative below that.

    Args:
        n_samples: Number of sampled sequences (at least 2)
        n_replicates: Number of independent samples to simulate (at least 1)
        theta: Population mutation parameter of the whole sequence, ``4 * N_ref * mu``
        segregating_sites: Number of mutations in every replicate
        epochs: Demographic history (see :class:`PopulationEpoch`); constant size when omitted
        rng: ``random.Random``, ``numpy.random.Generator`` or integer seed; ``None`` draws fresh entropy

    Returns:
        Integer array of shape ``(n_replicates, n_samples - 1)`` whose entry ``[r, i - 1]`` counts the
        mutations of replicate ``r`` carried by exactly ``i`` of the sequences (ancestral state known)

    Raises:
        ValueError: If neither or both of ``theta`` and ``segregating_sites`` are given, or the
            epochs or the number of segregating sites are invalid
    """
    validation.validate_range(n_samples, min_val=2, name="n_samples")
    validation.validate_range(n_replicates, min_val=1, name="n_replicates")
    if (theta is None) == (segregating_sites is None):
        raise ValueError("give exactly one of theta and segregating_sites")
    if theta is not None:
        _check_theta(theta)
    if segregating_sites is not None:
        if not isinstance(segregating_sites, (int, np.integer)) or segregating_sites < 0:
            raise ValueError(f"segregating_sites must be a non-negative integer, got {segregating_sites!r}")
    history = _SizeHistory(epochs)
    generator = _numpy_generator(rng)

    block = max(1, _BLOCK_CELLS // n_samples)
    blocks = []
    for begin in range(0, n_replicates, block):
        times = _coalescent_event_times(n_samples, min(block, n_replicates - begin), history, generator)
        lengths = _branch_length_spectra(times, n_samples, generator)
        if segregating_sites is not None:
            blocks.append(generator.multinomial(segregating_sites, lengths / lengths.sum(axis=1, keepdims=True)))
        else:
            blocks.append(generator.poisson((theta or 0.0) / 2.0 * lengths))
    return np.concatenate(blocks)


def simulate_coalescent_haplotypes(
    n_samples: int,
    theta: float,
    *,
    epochs: Sequence[PopulationEpoch] | None = None,
    rng: RandomSource = None,
) -> np.ndarray:
    """Simulate one neutral coalescent sample with infinite-sites mutation as a haplotype matrix.

    A genealogy of ``n_samples`` sequences is drawn from the Kingman coalescent (constant size, or the
    demographic history ``epochs``) and every branch of length ``l`` receives a Poisson number of
    mutations with mean ``theta * l / 2``. A mutation is carried by exactly the sequences below its
    branch, and no site mutates twice. There is no recombination, so all sites share the genealogy:
    any two columns pass the four-gamete test.

    Args:
        n_samples: Number of sampled sequences (at least 2)
        theta: Population mutation parameter of the whole sequence, ``4 * N_ref * mu`` with ``N_ref``
            the reference size that sets the time unit of ``epochs``
        epochs: Demographic history (see :class:`PopulationEpoch`); constant size when omitted
        rng: ``random.Random``, ``numpy.random.Generator`` or integer seed; ``None`` draws fresh entropy

    Returns:
        ``uint8`` array of shape ``(n_samples, S)``: entry ``[i, j]`` is 1 when sequence ``i`` carries
        the derived allele of segregating site ``j`` and 0 for the ancestral allele. Columns are in
        random order and none is invariant; ``S`` is random (mean ``theta * a_n`` under constant size).

    Raises:
        ValueError: If ``theta`` is not finite or the epochs are invalid
    """
    validation.validate_range(n_samples, min_val=2, name="n_samples")
    _check_theta(theta)
    history = _SizeHistory(epochs)
    generator = _numpy_generator(rng)
    times = _coalescent_event_times(n_samples, 1, history, generator)[0]

    n_nodes = 2 * n_samples - 1
    below = np.zeros((n_nodes, n_samples), dtype=bool)  # which samples sit below each node
    below[np.arange(n_samples), np.arange(n_samples)] = True
    born = np.zeros(n_nodes)
    branch_length = np.zeros(n_nodes - 1)  # branch above every node but the root
    lineages = list(range(n_samples))
    for step in range(n_samples - 1):
        first, second = generator.choice(len(lineages), size=2, replace=False)
        low, high = min(first, second), max(first, second)
        child_a, child_b, parent = lineages[low], lineages[high], n_samples + step
        branch_length[child_a] = times[step] - born[child_a]
        branch_length[child_b] = times[step] - born[child_b]
        below[parent] = below[child_a] | below[child_b]
        born[parent] = times[step]
        lineages[low] = parent
        lineages[high] = lineages[-1]
        lineages.pop()

    mutations = generator.poisson(theta / 2.0 * branch_length)
    columns = np.repeat(below[: n_nodes - 1], mutations, axis=0)
    columns = columns[generator.permutation(columns.shape[0])]
    return np.ascontiguousarray(columns.T, dtype=np.uint8)


def haplotypes_to_sequences(
    haplotypes: np.ndarray | Sequence[Sequence[int]],
    length: int | None = None,
    *,
    reference_sequence: str | None = None,
    rng: RandomSource = None,
) -> List[str]:
    """Write an infinite-sites haplotype matrix as aligned DNA sequences.

    Every column becomes one segregating site at its own position, drawn uniformly without replacement
    from the A/C/G/T positions of the ancestral sequence; sequences with a 0 there keep the ancestral
    base and sequences with a 1 carry a derived base drawn uniformly from the three other bases. The
    ancestral sequence is ``reference_sequence`` or, without one, ``length`` bases drawn uniformly from
    A/C/G/T. Positions and bases are seeded through ``rng``.

    Args:
        haplotypes: Matrix of shape ``(n_sequences, n_sites)`` of 0 (ancestral) and 1 (derived) alleles,
            as returned by :func:`simulate_coalescent_haplotypes`
        length: Sequence length; required without ``reference_sequence``, which sets it otherwise
        reference_sequence: Ancestral sequence; positions holding anything but A/C/G/T (either case)
            never mutate, and a derived base keeps the case of the base it replaces
        rng: ``random.Random``, ``numpy.random.Generator`` or integer seed; ``None`` draws fresh entropy

    Returns:
        One string per row of ``haplotypes``

    Raises:
        ValueError: If the matrix is malformed, neither ``length`` nor ``reference_sequence`` is
            given (or they disagree), or there are more segregating sites than A/C/G/T positions
    """
    matrix = np.asarray(haplotypes)
    if matrix.ndim != 2 or matrix.shape[0] < 1:
        raise ValueError("haplotypes must be a 2-D array with one row per sequence")
    if not np.isin(matrix, (0, 1)).all():
        raise ValueError("haplotypes must hold only 0 (ancestral) and 1 (derived) alleles")
    source = _python_random(rng)
    if reference_sequence is not None:
        if length is not None and length != len(reference_sequence):
            raise ValueError(f"length {length} does not match the reference sequence length {len(reference_sequence)}")
        ancestral = reference_sequence
    elif length is None:
        raise ValueError("length (or reference_sequence) is required")
    else:
        validation.validate_range(length, min_val=1, name="length")
        ancestral = "".join(source.choices(_BASES, k=length))
    try:
        codes = np.frombuffer(ancestral.encode("ascii"), dtype=np.uint8)
    except UnicodeEncodeError as exc:
        raise ValueError("the ancestral sequence must be ASCII") from exc

    eligible = [position for position, base in enumerate(ancestral) if base in "ACGTacgt"]
    n_sites = matrix.shape[1]
    if n_sites > len(eligible):
        raise ValueError(
            f"{n_sites} segregating sites do not fit in {len(eligible)} A/C/G/T positions: the infinite-sites "
            "model needs a lower theta or longer sequences"
        )
    positions = np.array(source.sample(eligible, n_sites), dtype=np.intp)
    derived = bytearray()
    for position in positions:
        base = ancestral[position]
        choice = source.choice([other for other in _BASES if other != base.upper()])
        derived.append(ord(choice.lower() if base.islower() else choice))
    derived_codes = np.frombuffer(bytes(derived), dtype=np.uint8)

    sequences = []
    for carries in matrix.astype(bool):
        row = codes.copy()
        row[positions[carries]] = derived_codes[carries]
        sequences.append(row.tobytes().decode("ascii"))
    return sequences


# --- recombination: Hudson's ancestral recombination graph -------------------------------------------------

#: Most coalescent plus recombination events one genealogy may take: a guard against runaway requests.
_MAX_ARG_EVENTS = 2_000_000

#: Ancestral material of a lineage: ``(left, right, node, samples below)`` for the stretch ``[left, right)``.
_Segment = Tuple[float, float, int, int]
#: An edge of the genealogy: ``(left, right, parent, child)``, the parent of ``child`` over ``[left, right)``.
_Edge = Tuple[float, float, int, int]


@dataclass(frozen=True, eq=False)
class MarginalTree:
    """The genealogy of the sample over one stretch ``[left, right)`` of the sequence.

    Attributes:
        left: First position of the stretch
        right: End of the stretch (exclusive)
        parent: ``parent[v]`` is the parent of node ``v`` in this tree, and ``-1`` for the root and for the
            nodes that do not belong to it; nodes ``0 .. n_samples - 1`` are the sampled sequences
        total_length: Sum of the branch lengths of the tree
        tmrca: Time of its most recent common ancestor
    """

    left: float
    right: float
    parent: np.ndarray
    total_length: float
    tmrca: float


@dataclass(frozen=True, eq=False)
class RecombinantGenealogy:
    """The ancestral recombination graph of a sample, as the marginal trees along the sequence.

    The sequence is the interval ``[0, sequence_length)``. Between two breakpoints the sample has one
    genealogy, a binary tree; across a breakpoint some branch of it moves. Nodes ``0 .. n_samples - 1`` are the
    samples at time 0, and the internal nodes follow in the order (and at the times) of their coalescence, so
    a node's identifier increases with its time. Time is in the units of :class:`PopulationEpoch`. A node is
    the parent of its children on a stretch of the sequence only (an edge); the edge arrays hold, for every
    such stretch ``[edge_left, edge_right)``, the ``edge_parent`` of a node ``edge_child``, with contiguous
    stretches of one edge merged.

    Attributes:
        n_samples: Number of sampled sequences
        sequence_length: Length of the sequence, in sites
        rho: Total population recombination parameter ``4 * N_ref * r`` of the whole sequence
        node_times: Time of every node
        edge_left: Start of every edge's stretch
        edge_right: End of every edge's stretch
        edge_parent: Parent node of every edge
        edge_child: Child node of every edge
        breakpoints: Boundaries of the marginal trees, from 0 to ``sequence_length`` (one more than trees)
        tree_lengths: Total branch length of every marginal tree
        tree_tmrca: Time of the most recent common ancestor of every marginal tree
        n_coalescences: Number of common-ancestor events of the graph, including those between lineages
            with no ancestral material in common
        n_recombinations: Number of recombination events of the graph: every event that split the ancestral
            material of a lineage between two lineages
        n_gap_recombinations: Number of those events whose breakpoint fell in a stretch of the lineage that
            is not ancestral to the sample (between two ancestral stretches). They move no branch of any
            marginal tree at the moment they happen, but they separate the ancestry of the stretches on
            either side. The others, :attr:`n_ancestral_recombinations`, are the breakpoints Hudson's
            expectation ``rho * sum(1 / i for i < n_samples)`` counts
    """

    n_samples: int
    sequence_length: float
    rho: float
    node_times: np.ndarray
    edge_left: np.ndarray
    edge_right: np.ndarray
    edge_parent: np.ndarray
    edge_child: np.ndarray
    breakpoints: np.ndarray
    tree_lengths: np.ndarray
    tree_tmrca: np.ndarray
    n_coalescences: int
    n_recombinations: int
    n_gap_recombinations: int

    @property
    def n_trees(self) -> int:
        """Number of distinct marginal trees along the sequence."""
        return int(self.tree_lengths.size)

    @property
    def n_ancestral_recombinations(self) -> int:
        """Recombination events whose breakpoint fell in ancestral material of the lineage it split."""
        return self.n_recombinations - self.n_gap_recombinations

    def _edges(self) -> List[_Edge]:
        return list(
            zip(
                self.edge_left.tolist(),
                self.edge_right.tolist(),
                self.edge_parent.tolist(),
                self.edge_child.tolist(),
            )
        )

    def trees(self) -> Iterator[MarginalTree]:
        """Yield the marginal trees from the left end of the sequence to the right end."""
        for left, right, parent, _, total, tmrca in _sweep_trees(
            self.n_samples, self.node_times.tolist(), self._edges(), carriers=False
        ):
            yield MarginalTree(left, right, np.array(parent, dtype=np.int64), total, tmrca)

    def tree_at(self, position: float) -> MarginalTree:
        """The marginal tree over the stretch that holds ``position``."""
        if not 0.0 <= position < self.sequence_length:
            raise ValueError(f"position must lie in [0, {self.sequence_length}), got {position}")
        index = int(np.searchsorted(self.breakpoints, position, side="right")) - 1
        for tree_index, tree in enumerate(self.trees()):
            if tree_index == index:
                return tree
        raise RuntimeError("internal error: the marginal trees do not cover the sequence")  # pragma: no cover

    def mutate(
        self,
        *,
        theta: float | None = None,
        segregating_sites: int | None = None,
        rng: RandomSource = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Place infinite-sites mutations on the marginal trees and return where they are and who carries them.

        Exactly one of ``theta`` and ``segregating_sites`` fixes the mutation load, as in
        :func:`simulate_coalescent_spectra`:

        * ``theta``: mutations fall on a stretch of the sequence of width ``w`` with a tree of total length
          ``T`` as a Poisson number with mean ``theta * w * T / (2 * sequence_length)``, so under constant size
          the expected number of segregating sites is ``theta * sum(1 / i for i < n_samples)`` however much
          the sequence recombines.
        * ``segregating_sites``: exactly that many mutations, spread over the stretches in proportion to
          ``w * T`` (the stretch's share of the total branch length of the graph) on a genealogy drawn
          without weighting, the approach of Hudson's ``ms -s``.

        Within a stretch a mutation takes a uniformly distributed position and a branch chosen in proportion to
        its length, and is carried by the sequences below that branch. Every mutation has its own position,
        so none happens twice.

        Args:
            theta: Population mutation parameter of the whole sequence, ``4 * N_ref * mu * sequence_length``
            segregating_sites: Number of mutations to place
            rng: ``random.Random``, ``numpy.random.Generator`` or integer seed; ``None`` draws fresh entropy

        Returns:
            ``(positions, haplotypes)``: the sorted positions of the mutations, a float array of shape
            ``(S,)``, and the ``uint8`` matrix of shape ``(n_samples, S)`` whose entry ``[i, j]`` is 1 when
            sequence ``i`` carries the derived allele of mutation ``j``

        Raises:
            ValueError: If neither or both of the arguments are given, or either is invalid
        """
        if (theta is None) == (segregating_sites is None):
            raise ValueError("give exactly one of theta and segregating_sites")
        if theta is not None:
            _check_theta(theta)
        if segregating_sites is not None:
            if not isinstance(segregating_sites, (int, np.integer)) or segregating_sites < 0:
                raise ValueError(f"segregating_sites must be a non-negative integer, got {segregating_sites!r}")
        generator = _numpy_generator(rng)

        masses = self.tree_lengths * np.diff(self.breakpoints)  # a stretch's share of the branch length
        if segregating_sites is not None:
            counts = generator.multinomial(int(segregating_sites), masses / masses.sum())
        else:
            counts = generator.poisson((theta or 0.0) / (2.0 * self.sequence_length) * masses)
        total = int(counts.sum())
        if total == 0:
            return np.zeros(0), np.zeros((self.n_samples, 0), dtype=np.uint8)

        positions = np.empty(total)
        carriers: List[int] = []
        times = self.node_times
        filled = 0
        for tree_index, (left, right, parent, mask, _, _) in enumerate(
            _sweep_trees(self.n_samples, times.tolist(), self._edges(), carriers=True)
        ):
            count = int(counts[tree_index])
            if count == 0:
                continue
            parents = np.asarray(parent)
            children = np.flatnonzero(parents >= 0)
            cumulative = np.cumsum(times[parents[children]] - times[children])
            chosen = np.searchsorted(cumulative, generator.random(count) * cumulative[-1], side="right")
            chosen = np.minimum(chosen, children.size - 1)
            stretch = generator.random(count)
            positions[filled : filled + count] = np.minimum(left + (right - left) * stretch, np.nextafter(right, left))
            carriers.extend(mask[node] for node in children[chosen].tolist())
            filled += count

        order = np.argsort(positions, kind="stable")
        width = (self.n_samples + 7) // 8
        packed = np.frombuffer(b"".join(carriers[i].to_bytes(width, "little") for i in order), dtype=np.uint8)
        bits = np.unpackbits(packed.reshape(total, width), axis=1, bitorder="little")[:, : self.n_samples]
        return positions[order], np.ascontiguousarray(bits.T)


@dataclass(frozen=True, eq=False)
class RecombinantSample:
    """A sample from the coalescent with recombination: its genealogy and the mutations placed on it.

    Attributes:
        genealogy: The marginal trees along the sequence
        positions: Sorted positions of the segregating sites in ``[0, sequence_length)``
        haplotypes: ``uint8`` matrix of shape ``(n_samples, S)`` of 0 (ancestral) and 1 (derived) alleles, one
            column per segregating site in the order of ``positions``
    """

    genealogy: RecombinantGenealogy
    positions: np.ndarray
    haplotypes: np.ndarray


def _sweep_trees(
    n_samples: int, node_times: Sequence[float], edges: Sequence[_Edge], *, carriers: bool
) -> Iterator[Tuple[float, float, List[int], List[int], float, float]]:
    """Walk the marginal trees from left to right, yielding ``(left, right, parent, carrier, length, tmrca)``.

    The ``parent`` list (and the ``carrier`` list when ``carriers`` is set) is updated in place between the
    yields, edge by edge, so a consumer must copy what it keeps. ``carrier[v]`` is a bit mask of the samples
    below node ``v``. Nodes are identified by increasing time, so ``parent[v] > v``.
    """
    n_nodes = len(node_times)
    parent = [-1] * n_nodes
    n_children = [0] * n_nodes
    carrier = [1 << i for i in range(n_samples)] + [0] * (n_nodes - n_samples) if carriers else []
    inserted: Dict[float, List[int]] = {}
    removed: Dict[float, List[int]] = {}
    for index, (left, right, _, _) in enumerate(edges):
        inserted.setdefault(left, []).append(index)
        removed.setdefault(right, []).append(index)
    boundaries = sorted(inserted.keys() | removed.keys())
    roots: set[int] = set()
    length = 0.0
    for left, right in zip(boundaries, boundaries[1:]):
        for index in removed.get(left, ()):
            _, _, above, node = edges[index]
            parent[node] = -1
            length -= node_times[above] - node_times[node]
            n_children[above] -= 1
            if n_children[above] == 0:
                roots.discard(above)
            if n_children[node] > 0:
                roots.add(node)
            if carriers:
                below = carrier[node]
                while above != -1:
                    carrier[above] ^= below
                    above = parent[above]
        for index in inserted.get(left, ()):
            _, _, above, node = edges[index]
            parent[node] = above
            length += node_times[above] - node_times[node]
            n_children[above] += 1
            roots.discard(node)
            if parent[above] == -1:
                roots.add(above)
            if carriers:
                below = carrier[node]
                while above != -1:
                    carrier[above] |= below
                    above = parent[above]
        if len(roots) != 1:  # pragma: no cover - a broken genealogy
            raise RuntimeError(f"internal error: {len(roots)} roots over [{left}, {right})")
        yield left, right, parent, carrier, length, node_times[next(iter(roots))]


def _merge_lineages(
    first: Sequence[_Segment],
    second: Sequence[_Segment],
    n_samples: int,
    time: float,
    node_times: List[float],
    edges: List[_Edge],
) -> List[_Segment]:
    """Merge two lineages at ``time``: the common-ancestor event of Hudson's algorithm.

    Where the two carry ancestral material of the same stretch their ancestries coalesce into a new node
    (one per event, appended to ``node_times``) with an edge to each of the two nodes; stretches carried by
    one lineage only pass to the merged lineage unchanged. A stretch below which all ``n_samples`` samples
    have been gathered has reached its most recent common ancestor, is no longer ancestral to anything, and
    is dropped. The result may be empty.
    """
    out: List[_Segment] = []
    n_first, n_second = len(first), len(second)
    i = j = 0
    head_first, head_second = first[0], second[0]
    new_node = -1
    while i < n_first and j < n_second:
        left_a, right_a, node_a, count_a = head_first
        left_b, right_b, node_b, count_b = head_second
        if right_a <= left_b:
            out.append(head_first)
            i += 1
            if i < n_first:
                head_first = first[i]
        elif right_b <= left_a:
            out.append(head_second)
            j += 1
            if j < n_second:
                head_second = second[j]
        else:
            if left_a < left_b:
                out.append((left_a, left_b, node_a, count_a))
            elif left_b < left_a:
                out.append((left_b, left_a, node_b, count_b))
            left = max(left_a, left_b)
            right = min(right_a, right_b)
            if new_node < 0:
                new_node = len(node_times)
                node_times.append(time)
            edges.append((left, right, new_node, node_a))
            edges.append((left, right, new_node, node_b))
            if count_a + count_b < n_samples:
                out.append((left, right, new_node, count_a + count_b))
            if right_a == right:
                i += 1
                if i < n_first:
                    head_first = first[i]
            else:
                head_first = (right, right_a, node_a, count_a)
            if right_b == right:
                j += 1
                if j < n_second:
                    head_second = second[j]
            else:
                head_second = (right, right_b, node_b, count_b)
    if i < n_first:
        out.append(head_first)
        out.extend(first[i + 1 :])
    if j < n_second:
        out.append(head_second)
        out.extend(second[j + 1 :])

    joined: List[_Segment] = []
    for segment in out:  # adjacent stretches under one node with equal sample counts are one stretch
        if joined and joined[-1][1] == segment[0] and joined[-1][2:] == segment[2:]:
            joined[-1] = (joined[-1][0], segment[1], segment[2], segment[3])
        else:
            joined.append(segment)
    return joined


def _split_lineage(lineage: Sequence[_Segment], cut: float) -> Tuple[List[_Segment], List[_Segment], bool]:
    """Split a lineage at the position ``cut``; also report whether the cut fell in ancestral material."""
    index = 0
    while lineage[index][1] <= cut:
        index += 1
    left, right, node, count = lineage[index]
    if left < cut:
        return (
            list(lineage[:index]) + [(left, cut, node, count)],
            [(cut, right, node, count)] + list(lineage[index + 1 :]),
            True,
        )
    return list(lineage[:index]), list(lineage[index:]), False


def _simulate_arg(
    n_samples: int, length: float, rho: float, history: _SizeHistory, generator: np.random.Generator
) -> Tuple[List[float], List[_Edge], int, int, int]:
    """Hudson's algorithm: the coalescent with recombination, backwards in time, over ``[0, length)``.

    Every lineage carries the stretches of the sequence that are still ancestral to the sample, each with
    the node it stands for and the number of samples below it. With ``k`` lineages, any of the ``k (k - 1) / 2``
    pairs merges at rate ``1 / size(t)`` (the coalescent clock, mapped through the integrated hazard of
    ``history``), and a lineage whose material spans ``s`` units of sequence splits in two at rate
    ``rho * s / (2 * length)`` with a breakpoint uniform over that span. The earlier of the two independent
    clocks decides the next event. The process ends when every stretch of the sequence has reached its most
    recent common ancestor.

    Returns:
        ``(node_times, edges, n_coalescences, n_recombinations, n_gap_recombinations)`` with the edges
        merged along contiguous stretches
    """
    rate_per_span = rho / (2.0 * length)
    lineages: List[List[_Segment]] = [[(0.0, length, node, 1)] for node in range(n_samples)]
    spans = [length] * n_samples  # span of every lineage: from its first to its last ancestral position
    node_times: List[float] = [0.0] * n_samples
    edges: List[_Edge] = []
    exponential, uniform = generator.standard_exponential, generator.random
    time = hazard = 0.0
    n_coalescences = n_recombinations = n_gap_recombinations = 0

    def drop(index: int) -> None:
        last, last_span = lineages.pop(), spans.pop()
        if index < len(lineages):
            lineages[index], spans[index] = last, last_span

    while lineages:
        if n_coalescences + n_recombinations >= _MAX_ARG_EVENTS:
            raise ValueError(
                f"the ancestral recombination graph needs more than {_MAX_ARG_EVENTS} events: lower rho "
                f"({rho:g}), the sample size ({n_samples}) or the sequence length"
            )
        k = len(lineages)
        next_recombination = time + exponential() / (rate_per_span * sum(spans)) if rho > 0.0 else math.inf
        step = exponential() / (k * (k - 1) / 2.0)
        next_coalescence = history.time_at(hazard + step)
        if next_coalescence <= next_recombination:
            time, hazard = next_coalescence, hazard + step
            first = int(uniform() * k)
            second = int(uniform() * (k - 1))
            second += second >= first
            merged = _merge_lineages(lineages[first], lineages[second], n_samples, time, node_times, edges)
            n_coalescences += 1
            if merged:
                lineages[first], spans[first] = merged, merged[-1][1] - merged[0][0]
                drop(second)
            else:
                drop(max(first, second))
                drop(min(first, second))
        else:
            time, hazard = next_recombination, history.hazard_at(next_recombination)
            cumulative = list(accumulate(spans))
            chosen = min(bisect_right(cumulative, uniform() * cumulative[-1]), k - 1)
            lineage = lineages[chosen]
            start, end = lineage[0][0], lineage[-1][1]
            cut = start + (end - start) * uniform()
            if not start < cut < end:  # the draw sits on an end of the span: no event at all
                continue
            left_lineage, right_lineage, ancestral = _split_lineage(lineage, cut)
            lineages[chosen], spans[chosen] = left_lineage, left_lineage[-1][1] - left_lineage[0][0]
            lineages.append(right_lineage)
            spans.append(right_lineage[-1][1] - right_lineage[0][0])
            n_recombinations += 1
            n_gap_recombinations += not ancestral

    edges.sort(key=lambda edge: (edge[2], edge[3], edge[0]))
    squashed: List[_Edge] = []
    for edge in edges:  # contiguous stretches of one (parent, child) pair are one edge
        if squashed and squashed[-1][2:] == edge[2:] and squashed[-1][1] == edge[0]:
            squashed[-1] = (squashed[-1][0], edge[1], edge[2], edge[3])
        else:
            squashed.append(edge)
    return node_times, squashed, n_coalescences, n_recombinations, n_gap_recombinations


def _resolve_rho(rho: float | None, rho_per_site: float | None, sequence_length: float) -> float:
    """The total recombination parameter of a request: ``rho``, or ``rho_per_site`` times the length, or 0."""
    if rho is not None and rho_per_site is not None:
        raise ValueError("give at most one of rho (the total) and rho_per_site")
    value = 0.0
    for name, rate in (("rho", rho), ("rho_per_site", rho_per_site)):
        if rate is not None:
            if not math.isfinite(rate):
                raise ValueError(f"{name} must be a finite number, got {rate}")
            validation.validate_range(rate, min_val=0.0, name=name)
            value = rate * sequence_length if name == "rho_per_site" else rate
    if not math.isfinite(value):
        raise ValueError("the total recombination parameter overflows: lower rho_per_site or the sequence length")
    return float(value)


def simulate_recombinant_genealogy(
    n_samples: int,
    sequence_length: float,
    *,
    rho: float | None = None,
    rho_per_site: float | None = None,
    epochs: Sequence[PopulationEpoch] | None = None,
    rng: RandomSource = None,
) -> RecombinantGenealogy:
    """Simulate the genealogy of a sample over a recombining sequence: Hudson's ancestral recombination graph.

    The sequence is the interval ``[0, sequence_length)`` of ``sequence_length`` sites and recombination can
    break it anywhere (at a uniformly distributed position, so every site boundary is equally likely). Going
    back in time from the sample, each lineage carries the stretches of the sequence that are ancestral to the
    sample. A pair of lineages merges at rate ``1 / size(t)`` as in :func:`simulate_coalescent_spectra`, and a
    lineage splits in two, its material divided at a breakpoint, at rate ``rho_per_site / 2`` per site of the
    span between its first and last ancestral site (``rho / 2`` for material covering the whole sequence).
    With this clock every site, taken on its own, has the Kingman genealogy (constant size or ``epochs``),
    and the sites are correlated by the branches they share. ``rho`` is the total ``4 * N_ref * r`` of the
    sequence, with ``r`` the probability per generation of a recombination breakpoint within it and ``N_ref`` the
    reference size that sets the time unit of ``epochs``: ``rho_per_site`` is the same rate per site. With
    ``rho = 0`` the graph is one tree, drawn as :func:`simulate_coalescent_haplotypes` draws it. Under
    constant size the graph holds ``rho * sum(1 / i for i < n_samples)`` recombination events in expectation
    whose breakpoint fell in ancestral material (:attr:`RecombinantGenealogy.n_ancestral_recombinations`),
    and some further ones in stretches that are not ancestral, which grow in number with ``rho``. The cost
    of the simulation grows with the number of events, a little faster than ``rho ** 1.8``: 100 samples
    take about 0.1 s at ``rho = 400``, 1 s at ``rho = 2000`` and 5 s at ``rho = 4000``.

    Args:
        n_samples: Number of sampled sequences (at least 2)
        sequence_length: Length of the sequence in sites (positive and finite)
        rho: Total population recombination parameter of the whole sequence
        rho_per_site: Population recombination parameter per site (``rho = rho_per_site * sequence_length``).
            At most one of the two may be given; with neither, there is no recombination
        epochs: Demographic history (see :class:`PopulationEpoch`); constant size when omitted
        rng: ``random.Random``, ``numpy.random.Generator`` or integer seed; ``None`` draws fresh entropy

    Returns:
        The graph, with the marginal trees recoverable from :meth:`RecombinantGenealogy.trees`; place
        mutations on it with :meth:`RecombinantGenealogy.mutate`

    Raises:
        ValueError: If ``n_samples``, ``sequence_length`` or the recombination parameters are invalid, both
            ``rho`` and ``rho_per_site`` are given, the epochs are invalid, or the graph would need more
            than two million events
    """
    validation.validate_range(n_samples, min_val=2, name="n_samples")
    if not isinstance(sequence_length, (int, float, np.integer, np.floating)) or not 0.0 < sequence_length < math.inf:
        raise ValueError(f"sequence_length must be a positive finite number, got {sequence_length!r}")
    length = float(sequence_length)
    total_rho = _resolve_rho(rho, rho_per_site, length)
    history = _SizeHistory(epochs)
    generator = _numpy_generator(rng)

    node_times, edges, n_coalescences, n_recombinations, n_gaps = _simulate_arg(
        n_samples, length, total_rho, history, generator
    )
    lefts, rights, parents, children = (np.array(column) for column in zip(*edges))
    tree_lengths, tree_tmrca, boundaries = [], [], [0.0]
    for left, right, _, _, total, tmrca in _sweep_trees(n_samples, node_times, edges, carriers=False):
        tree_lengths.append(total)
        tree_tmrca.append(tmrca)
        boundaries.append(right)
    return RecombinantGenealogy(
        n_samples=n_samples,
        sequence_length=length,
        rho=total_rho,
        node_times=np.array(node_times),
        edge_left=lefts.astype(float),
        edge_right=rights.astype(float),
        edge_parent=parents.astype(np.int64),
        edge_child=children.astype(np.int64),
        breakpoints=np.array(boundaries),
        tree_lengths=np.array(tree_lengths),
        tree_tmrca=np.array(tree_tmrca),
        n_coalescences=n_coalescences,
        n_recombinations=n_recombinations,
        n_gap_recombinations=n_gaps,
    )


def simulate_recombinant_haplotypes(
    n_samples: int,
    sequence_length: float,
    *,
    theta: float | None = None,
    segregating_sites: int | None = None,
    rho: float | None = None,
    rho_per_site: float | None = None,
    epochs: Sequence[PopulationEpoch] | None = None,
    rng: RandomSource = None,
) -> RecombinantSample:
    """Simulate one neutral sample over a recombining sequence, with infinite-sites mutation.

    The genealogy is Hudson's ancestral recombination graph (:func:`simulate_recombinant_genealogy`) and
    the mutations are placed on its marginal trees at uniformly distributed positions
    (:meth:`RecombinantGenealogy.mutate`): exactly one of ``theta`` and ``segregating_sites`` fixes their
    number. Every site has the neutral genealogy of :func:`simulate_coalescent_haplotypes`, so the expected
    site frequency spectrum, the expected number of segregating sites and the expected diversity are those of
    the non-recombining sample (``E[S] = theta * a_n``, ``E[xi_i] = theta / i``, ``E[pi] = theta`` under
    constant size, for ``theta`` the total over the sequence); recombination changes their variances, and it
    makes pairs of sites at distance ``d`` less correlated as ``rho_per_site * d`` grows: linkage
    disequilibrium decays, and two sites can show all four gametes. Sites of one marginal tree never do.

    Args:
        n_samples: Number of sampled sequences (at least 2)
        sequence_length: Length of the sequence in sites
        theta: Population mutation parameter of the whole sequence, ``4 * N_ref * mu * sequence_length``
        segregating_sites: Number of mutations (alternative to ``theta``)
        rho: Total population recombination parameter ``4 * N_ref * r`` of the whole sequence
        rho_per_site: The same rate per site (at most one of ``rho`` and ``rho_per_site``)
        epochs: Demographic history (see :class:`PopulationEpoch`); constant size when omitted
        rng: ``random.Random``, ``numpy.random.Generator`` or integer seed; ``None`` draws fresh entropy

    Returns:
        The genealogy, the sorted mutation positions and the haplotype matrix

    Raises:
        ValueError: As :func:`simulate_recombinant_genealogy` and :meth:`RecombinantGenealogy.mutate`
    """
    if (theta is None) == (segregating_sites is None):
        raise ValueError("give exactly one of theta and segregating_sites")
    generator = _numpy_generator(rng)
    genealogy = simulate_recombinant_genealogy(
        n_samples, sequence_length, rho=rho, rho_per_site=rho_per_site, epochs=epochs, rng=generator
    )
    positions, haplotypes = genealogy.mutate(theta=theta, segregating_sites=segregating_sites, rng=generator)
    return RecombinantSample(genealogy=genealogy, positions=positions, haplotypes=haplotypes)


def _coalescent_sequences(
    n_sequences: int,
    length: int,
    theta_per_site: float,
    *,
    rng: random.Random,
    epochs: Sequence[PopulationEpoch] | None = None,
    reference_sequence: str | None = None,
) -> List[str]:
    """Sample aligned sequences from the coalescent; ``theta_per_site`` applies to every A/C/G/T site."""
    ancestral = reference_sequence if reference_sequence is not None else "".join(rng.choices(_BASES, k=length))
    n_sites = sum(base in "ACGTacgt" for base in ancestral)
    theta = theta_per_site * n_sites
    if epochs is None:
        expected_sites = theta * sum(1.0 / i for i in range(1, n_sequences))
        if expected_sites > n_sites:
            raise ValueError(
                f"theta of {theta_per_site:g} per site gives about {expected_sites:.0f} expected segregating sites, "
                f"more than the {n_sites} A/C/G/T positions available; the infinite-sites model needs "
                "theta * a_n well below the sequence length (lower theta or lengthen the sequences)"
            )
    haplotypes = simulate_coalescent_haplotypes(n_sequences, theta, epochs=epochs, rng=rng)
    return haplotypes_to_sequences(haplotypes, reference_sequence=ancestral, rng=rng)


def generate_population_sequences(
    n_sequences: int,
    length: int | None = None,
    *,
    theta: float | None = None,
    n_sites: int = 1000,
    rng: random.Random | None = None,
    # Parameter aliases for backward compatibility
    sequence_length: int | None = None,
    random_seed: int | None = None,
    mutation_rate: float | None = None,
    nucleotide_diversity: float | None = None,
    wattersons_theta: float | None = None,
    reference_sequence: str | None = None,
) -> List[str]:
    """Generate aligned DNA sequences as one neutral coalescent sample.

    The sequences are a sample from the standard neutral coalescent (constant population size, no
    recombination, infinite-sites mutation; see :func:`simulate_coalescent_haplotypes`). They share a
    single genealogy, so every variable site is one mutation on one branch of it, and the site
    frequency spectrum, Tajima's D and the diversity follow neutral expectations: with ``theta`` the
    parameter of the whole sequence (defined below), ``E[S] = theta * a_n`` segregating sites,
    ``E[pi] = theta`` pairwise differences, ``E[xi_i] = theta / i`` mutations carried by ``i``
    sequences and a mean Tajima's D just below 0 (about -0.1 for 20 sequences; the standardized
    statistic is slightly biased). One locus is one draw from this distribution, so its statistics
    scatter widely around these means.

    Mutation parameter: every way of asking for variation is mapped to ``theta_site``, the population
    mutation parameter per site (``4 * N * mu`` with ``mu`` the per-site mutation rate), and the
    sample is drawn with ``theta = theta_site * L`` for the ``L`` A/C/G/T sites of the ancestral
    sequence. When several are given the first of these decides and the others are ignored:

    1. ``nucleotide_diversity``: ``theta_site`` equals the per-site pi asked for, because
       ``E[pi] = theta_site``
    2. ``wattersons_theta``: ``theta_site`` equals the per-site Watterson's theta asked for, because
       ``E[theta_W] = theta_site``
    3. ``theta``: used as ``theta_site`` (per site, as in :func:`generate_two_populations`)
    4. ``mutation_rate``: per-site rate ``mu`` per generation. The population size is taken to be
       ``n_sequences``, as the former ``theta = 4 * n_sequences * mutation_rate`` did, so
       ``theta_site = 4 * n_sequences * mutation_rate``
    5. none of them: ``theta_site = 0.01``

    The diversity targets are expectations over genealogies; the pi and Watterson's theta of a
    realised sample scatter around them.

    Args:
        n_sequences: Number of sequences to generate (at least 2)
        length: Length of each sequence
        theta: Per-site population mutation parameter ``4*N*mu`` (see above)
        n_sites: Accepted and validated for backward compatibility, but not used. The former
            star-shaped generator mutated a pool of this many candidate sites; under the
            infinite-sites coalescent the number of variable sites is an outcome of the genealogy
            and ``theta`` (mean ``theta * a_n``), not an input, and the variable sites are spread
            uniformly over the sequence whatever the size of a candidate pool
        rng: Random number generator for reproducibility
        sequence_length: Alias for length (backward compatibility)
        random_seed: Random seed for reproducibility (creates rng if provided)
        mutation_rate: Per-site mutation rate per generation (alternative to theta; see above)
        nucleotide_diversity: Target per-site nucleotide diversity pi (see above)
        wattersons_theta: Target per-site Watterson's theta (see above)
        reference_sequence: Ancestral sequence, whose length replaces ``length``. Positions that do not
            segregate carry its bases, and only its A/C/G/T positions can mutate

    Returns:
        List of ``n_sequences`` equal-length DNA sequences. When ``theta * a_n`` is small a sample can be
        monomorphic: the probability is the product of ``j / (j + theta)`` over ``j < n_sequences``.

    Raises:
        ValueError: If ``length`` is missing, or the infinite-sites model cannot hold the variation
            (more expected or drawn segregating sites than A/C/G/T positions)
        ValidationError: If a numeric parameter is out of range
    """
    # Handle parameter aliases
    if sequence_length is not None and length is None:
        length = sequence_length
    if reference_sequence is not None:
        length = len(reference_sequence)
    if length is None:
        raise ValueError("length (or sequence_length) is required")
    if random_seed is not None and rng is None:
        rng = random.Random(random_seed)

    # Per-site population mutation parameter: the first parameter given decides
    if nucleotide_diversity is not None:
        theta_per_site = nucleotide_diversity
    elif wattersons_theta is not None:
        theta_per_site = wattersons_theta
    elif theta is not None:
        theta_per_site = theta
    elif mutation_rate is not None:
        theta_per_site = 4 * n_sequences * mutation_rate
    else:
        theta_per_site = _DEFAULT_THETA_PER_SITE

    validation.validate_range(n_sequences, min_val=2, name="n_sequences")
    validation.validate_range(length, min_val=1, name="length")
    _check_theta(theta_per_site)
    validation.validate_range(n_sites, min_val=1, name="n_sites")

    return _coalescent_sequences(
        n_sequences, length, theta_per_site, rng=coerce_rng(rng), reference_sequence=reference_sequence
    )


def generate_two_populations(
    n_pop1: int,
    n_pop2: int,
    length: int | None = None,
    *,
    f_st: float | None = None,
    fst: float | None = None,
    theta: float = 0.01,
    n_sites: int | None = None,
    rng: random.Random | None = None,
    sequence_length: int | None = None,
) -> Tuple[List[str], List[str]]:
    """Generate two populations differentiated by F_ST under the Balding-Nichols model.

    Each polymorphic site gets an ancestral derived-allele frequency p drawn
    uniformly from [0.05, 0.95]. Each population then has frequency
    p_k ~ Beta(p(1 - F)/F, (1 - p)(1 - F)/F), so E[p_k] = p and
    Var(p_k) = F p(1 - p), which makes Hudson's F_ST equal F in expectation.
    Every sequence carries the derived allele at each site independently with
    its population's frequency. F = 0 gives identical frequencies and F = 1
    fixes each population for one allele.

    Args:
        n_pop1: Size of population 1
        n_pop2: Size of population 2
        length: Sequence length
        f_st: Desired F_ST value (0-1)
        fst: Alias for f_st (lowercase)
        theta: Per-site population mutation rate; sets the default ``n_sites``
        n_sites: Number of polymorphic sites, at most ``length`` (a larger number is an error, not a
            truncation). Defaults to Watterson's expectation ``theta * length * a_n`` for the
            ``n_pop1 + n_pop2`` sequences (at least 1), which must also fit in ``length``.
        rng: Random number generator for reproducibility
        sequence_length: Alias for length

    Returns:
        Tuple of (population1_sequences, population2_sequences)

    Raises:
        ValueError: If parameters are invalid, or more polymorphic sites are asked for (or follow from
            ``theta``) than the sequences have positions
    """
    # Handle parameter aliases
    if sequence_length is not None and length is None:
        length = sequence_length
    if fst is not None and f_st is None:
        f_st = fst
    if f_st is None:
        f_st = 0.1
    if length is None:
        raise ValueError("length (or sequence_length) is required")

    validation.validate_range(n_pop1, min_val=2, name="n_pop1")
    validation.validate_range(n_pop2, min_val=2, name="n_pop2")
    validation.validate_range(length, min_val=1, name="length")
    validation.validate_range(f_st, min_val=0.0, max_val=1.0, name="f_st")
    validation.validate_range(theta, min_val=0.0, name="theta")
    if n_sites is None:
        a_n = sum(1.0 / i for i in range(1, n_pop1 + n_pop2))
        n_sites = max(1, round(theta * length * a_n))
    validation.validate_range(n_sites, min_val=1, name="n_sites")
    if n_sites > length:
        raise errors.ValidationError(
            f"{n_sites} polymorphic sites do not fit in sequences of length {length}: lower theta or n_sites, "
            "or lengthen the sequences"
        )

    rng = coerce_and_seed(rng)

    ancestral = [rng.choice("ACGT") for _ in range(length)]
    populations = [[list(ancestral) for _ in range(n_pop1)], [list(ancestral) for _ in range(n_pop2)]]
    for position in sorted(rng.sample(range(length), n_sites)):
        derived = rng.choice([base for base in "ACGT" if base != ancestral[position]])
        p = rng.uniform(0.05, 0.95)
        for population in populations:
            if f_st == 0:
                frequency = p
            elif f_st == 1:
                frequency = 1.0 if rng.random() < p else 0.0
            else:
                frequency = rng.betavariate(p * (1 - f_st) / f_st, (1 - p) * (1 - f_st) / f_st)
            for sequence in population:
                if rng.random() < frequency:
                    sequence[position] = derived

    pop1 = ["".join(sequence) for sequence in populations[0]]
    pop2 = ["".join(sequence) for sequence in populations[1]]
    return pop1, pop2


def generate_genotype_matrix(
    n_individuals: int,
    n_snps: int | None = None,
    *,
    n_sites: int | None = None,
    maf_min: float = 0.05,
    maf_max: float = 0.5,
    hwe_deviation: float = 0.0,
    hwe: bool = True,
    allele_frequencies: List[float] | None = None,
    ploidy: int = 2,
    rng: random.Random | None = None,
) -> List[List[int]]:
    """Generate a genotype matrix for population genetics analysis.

    Args:
        n_individuals: Number of individuals
        n_snps: Number of SNPs (alias: n_sites)
        n_sites: Number of sites (alias for n_snps)
        maf_min: Minimum minor allele frequency
        maf_max: Maximum minor allele frequency
        hwe_deviation: Deviation from Hardy-Weinberg equilibrium (0 = perfect HWE)
        hwe: Whether to generate genotypes under HWE (True) or not (False). ``hwe=False`` with no
            ``hwe_deviation`` uses a deviation of 0.3; an explicit ``hwe_deviation`` is always honoured
        allele_frequencies: List of allele frequencies for each site
        ploidy: Ploidy level (1=haploid, 2=diploid)
        rng: Random number generator for reproducibility

    Returns:
        Genotype matrix as list of lists (n_individuals x n_sites) with values 0-ploidy

    Raises:
        ValueError: If parameters are invalid
    """
    # Handle parameter aliases
    if n_sites is not None and n_snps is None:
        n_snps = n_sites
    if n_snps is None:
        raise ValueError("n_snps (or n_sites) is required")

    # Use allele_frequencies to infer n_snps if provided
    if allele_frequencies is not None:
        n_snps = len(allele_frequencies)

    validation.validate_range(n_individuals, min_val=1, name="n_individuals")
    validation.validate_range(n_snps, min_val=1, name="n_snps")
    validation.validate_range(maf_min, min_val=0.0, max_val=0.5, name="maf_min")
    validation.validate_range(maf_max, min_val=maf_min, max_val=0.5, name="maf_max")
    validation.validate_range(hwe_deviation, min_val=0.0, max_val=1.0, name="hwe_deviation")
    validation.validate_range(ploidy, min_val=1, max_val=2, name="ploidy")

    # hwe=False without a deviation of its own means "some deviation from HWE"
    if not hwe and hwe_deviation == 0:
        hwe_deviation = 0.3

    rng = coerce_and_seed(rng)

    # Sample one allele frequency per SNP up front so every individual at a
    # site shares the same population allele frequency.
    if allele_frequencies is not None:
        snp_mafs = list(allele_frequencies)
    else:
        snp_mafs = [rng.uniform(maf_min, maf_max) for _ in range(n_snps)]

    genotype_matrix = []

    for ind_idx in range(n_individuals):
        ind_genotypes = []
        for snp_idx in range(n_snps):
            maf = snp_mafs[snp_idx]

            # Generate genotype based on ploidy
            if ploidy == 1:
                # Haploid: 0 or 1
                genotype = 1 if rng.random() < maf else 0
            else:
                # Diploid: 0, 1, or 2
                if hwe_deviation == 0:
                    # Perfect HWE
                    p_aa = (1 - maf) ** 2
                    p_ab = 2 * maf * (1 - maf)
                    p_bb = maf**2

                    genotype = rng.choices([0, 1, 2], weights=[p_aa, p_ab, p_bb])[0]
                else:
                    # Deviate from HWE
                    if rng.random() < hwe_deviation:
                        genotype = rng.choice([0, 2])  # Favor homozygotes
                    else:
                        p_aa = (1 - maf) ** 2
                        p_ab = 2 * maf * (1 - maf)
                        p_bb = maf**2
                        genotype = rng.choices([0, 1, 2], weights=[p_aa, p_ab, p_bb])[0]

            ind_genotypes.append(genotype)
        genotype_matrix.append(ind_genotypes)

    return genotype_matrix


#: Mutations in every simulated locus behind the Tajima's D reported by the statistics interfaces.
_SIGNATURE_SEGREGATING_SITES = 50


def _check_summary_sample(sample_size: int, n_replicates: int) -> None:
    """Validate the size of the simulated sample behind the statistics interfaces."""
    validation.validate_range(sample_size, min_val=4, name="sample_size")
    validation.validate_range(n_replicates, min_val=2, name="n_replicates")


def _expansion_epochs(factor: float, duration: float) -> List[PopulationEpoch]:
    """Exponential growth by ``factor`` over the last ``duration`` time units before the sample."""
    if factor == 1.0:
        return [PopulationEpoch(0.0, 1.0)]
    return [PopulationEpoch(0.0, factor, math.log(factor) / duration), PopulationEpoch(duration, 1.0)]


def _bottleneck_epochs(size: float, duration: float, recovery: float = 0.0) -> List[PopulationEpoch]:
    """A bottleneck of relative ``size`` and ``duration`` that ended ``recovery`` time units before the sample."""
    epochs = (
        [PopulationEpoch(0.0, 1.0), PopulationEpoch(recovery, size)] if recovery > 0.0 else [PopulationEpoch(0.0, size)]
    )
    return epochs + [PopulationEpoch(recovery + duration, 1.0)]


def _epochs_from_sizes(sizes: Sequence[int], reference_size: float) -> List[PopulationEpoch]:
    """The history of a per-generation size trajectory read backwards from its last generation.

    Runs of equal sizes become one epoch (time in units of ``2 * reference_size`` generations), and the
    size of the first generation continues into the past for ever.
    """
    runs: List[List[int]] = []  # [size, generations], latest first
    for size in reversed(sizes):
        if runs and runs[-1][0] == size:
            runs[-1][1] += 1
        else:
            runs.append([size, 1])
    epochs, start = [], 0
    for size, length in runs:
        epochs.append(PopulationEpoch(start / (2.0 * reference_size), size / reference_size))
        start += length
    return epochs


def _expected_branch_spectrum(
    n_samples: int, epochs: Sequence[PopulationEpoch] | None, n_replicates: int, generator: np.random.Generator
) -> np.ndarray:
    """Mean total branch length carried by ``i`` of the samples, ``i = 1..n_samples-1``, over seeded genealogies.

    The expected site frequency spectrum of a neutral sample is proportional to it (``E[xi_i] = theta E[L_i] / 2``);
    its relative sampling error is about ``1 / sqrt(n_replicates)`` per class.
    """
    history = _SizeHistory(epochs)
    total = np.zeros(n_samples - 1)
    block = max(1, _BLOCK_CELLS // n_samples)
    for begin in range(0, n_replicates, block):
        times = _coalescent_event_times(n_samples, min(block, n_replicates - begin), history, generator)
        total += _branch_length_spectra(times, n_samples, generator).sum(axis=0)
    return total / n_replicates


def _tajimas_d(spectra: np.ndarray, n_samples: int) -> np.ndarray:
    """Tajima's (1989) D of every row of ``spectra`` (mutations carried by 1..n-1 samples); 0 without mutations."""
    i = np.arange(1, n_samples, dtype=float)
    a1, a2 = float(np.sum(1.0 / i)), float(np.sum(1.0 / i**2))
    b1 = (n_samples + 1) / (3.0 * (n_samples - 1))
    b2 = 2.0 * (n_samples**2 + n_samples + 3) / (9.0 * n_samples * (n_samples - 1))
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n_samples + 2) / (a1 * n_samples) + a2 / a1**2
    segregating = spectra.sum(axis=1).astype(float)
    pairwise = spectra @ (i * (n_samples - i) / (n_samples * (n_samples - 1) / 2.0))
    variance = c1 / a1 * segregating + c2 / (a1**2 + a2) * segregating * (segregating - 1.0)
    statistic: np.ndarray = np.divide(
        pairwise - segregating / a1, np.sqrt(variance), out=np.zeros_like(pairwise), where=variance > 0
    )
    return statistic


def _trajectory_statistics(
    population_sizes: Sequence[int],
    epochs: Sequence[PopulationEpoch],
    *,
    initial_size: int,
    mutation_rate: float,
    locus_length: int,
    sample_size: int,
    n_replicates: int,
    generator: np.random.Generator,
) -> Dict[str, Any]:
    """Statistics of a size trajectory that follow from the coalescent.

    ``diversity_trajectory`` is the expected per-site nucleotide diversity after every generation, from the
    recursion ``pi' = (1 - 1 / (2 N)) pi + 2 mu`` (two lineages share a parent with probability ``1 / (2 N)``
    and each accrues ``mu`` mutations per site), started at the mutation-drift equilibrium ``4 N mu`` of the
    initial size. ``mutation_trajectory`` counts the new mutations of the whole population, Poisson with mean
    ``2 N mu * locus_length`` per generation. The rest describes a sample of ``sample_size`` chromosomes drawn
    after the last generation, from ``n_replicates`` genealogies of the history ``epochs``.
    """
    pi = 4.0 * initial_size * mutation_rate
    diversity = []
    for size in population_sizes:
        pi = (1.0 - 1.0 / (2.0 * size)) * pi + 2.0 * mutation_rate
        diversity.append(pi)
    new_mutations = generator.poisson(2.0 * np.asarray(population_sizes, dtype=float) * mutation_rate * locus_length)
    mutation_trajectory = np.cumsum(new_mutations)

    theta_locus = 4.0 * initial_size * mutation_rate * locus_length
    branch = _expected_branch_spectrum(sample_size, epochs, n_replicates, generator)
    classes = np.arange(1, sample_size, dtype=float)
    spectra = simulate_coalescent_spectra(
        sample_size, n_replicates, segregating_sites=_SIGNATURE_SEGREGATING_SITES, epochs=epochs, rng=generator
    )
    d = _tajimas_d(spectra, sample_size)
    return {
        "diversity_trajectory": diversity,
        "mutation_trajectory": mutation_trajectory.tolist(),
        "total_mutations": int(mutation_trajectory[-1]),
        "final_diversity": diversity[-1],
        "sample_size": sample_size,
        "n_replicates": n_replicates,
        "expected_sfs": (theta_locus / 2.0 * branch).tolist(),
        "neutral_sfs": (theta_locus / classes).tolist(),
        "singleton_fraction": float(branch[0] / branch.sum()),
        "neutral_singleton_fraction": float(1.0 / np.sum(1.0 / classes)),
        "tajimas_d": float(d.mean()),
        "tajimas_d_se": float(d.std(ddof=1) / math.sqrt(n_replicates)),
    }


def simulate_bottleneck_population(
    initial_size: int | None = None,
    bottleneck_size: int | None = None,
    final_size: int | None = None,
    generations: int | None = None,
    *,
    mutation_rate: float = 1e-8,
    rng: random.Random | None = None,
    # Alternative interface for sequence generation
    n_sequences: int | None = None,
    sequence_length: int | None = None,
    bottleneck_duration: int | None = None,
    pre_bottleneck_diversity: float | None = None,
    # Statistics interface: the sample behind the reported site frequency spectrum and Tajima's D
    sample_size: int = 20,
    n_replicates: int = 2000,
) -> List[str] | Dict[str, Any]:
    """Simulate a population bottleneck followed by recovery.

    Returns either aligned sequences (the sequence interface, selected by passing ``n_sequences``) or a
    statistics dictionary (the default interface).

    Sequence interface: the sequences are a neutral coalescent sample (see
    :func:`simulate_coalescent_haplotypes`; one non-recombining genealogy) from a population of
    ``initial_size`` individuals that shrank to ``bottleneck_size`` for the last ``bottleneck_duration``
    generations. The sample is drawn at the end of the bottleneck, before any recovery. In units of
    ``2 * initial_size`` generations the bottleneck has relative size ``bottleneck_size / initial_size``
    and lasts ``bottleneck_duration / (2 * initial_size)``; its strength is set mainly by
    ``bottleneck_duration / (2 * bottleneck_size)``, the coalescent time it spans. Lineages coalesce
    quickly inside it and leave a few long ancestral branches, so diversity falls, rare variants are
    lost and Tajima's D is positive on average, the textbook signature of a bottleneck observed shortly
    afterwards (a severe bottleneck can leave the sample monomorphic). The ancestral population has
    per-site population mutation parameter ``pre_bottleneck_diversity`` (``4 * N * mu`` per site, its
    expected nucleotide diversity).

    Statistics interface: the deterministic size trajectory (``initial_size`` for the first third of
    ``generations``, ``bottleneck_size`` for the middle third, ``final_size`` for the last) and what the
    coalescent says about it, with the sample drawn after the last generation (see
    :func:`_trajectory_statistics`): ``diversity_trajectory`` is the expected per-site nucleotide diversity
    after every generation (an exact recursion, started at the equilibrium ``4 * initial_size * mutation_rate``
    of the initial population), ``mutation_trajectory`` the cumulative count of new mutations in the whole
    population (Poisson, ``2 * size * mutation_rate * sequence_length`` per generation), and ``expected_sfs``,
    ``singleton_fraction`` and ``tajimas_d`` describe a sample of ``sample_size`` chromosomes simulated from
    the coalescent of that history: ``expected_sfs[i - 1]`` is the expected number of mutations in the locus
    carried by ``i`` chromosomes, averaged over ``n_replicates`` seeded genealogies (``neutral_sfs`` is the
    constant-size expectation ``theta / i`` for comparison), and ``tajimas_d`` is the mean Tajima's D of
    ``n_replicates`` loci that carry 50 mutations each (``tajimas_d_se`` is its standard error). A bottleneck
    sampled while its signature lasts shows a deficit of rare variants (a ``singleton_fraction`` below
    ``neutral_singleton_fraction``) and a positive D; recovery for long enough turns both around.

    Args:
        initial_size: Population size before the bottleneck. Sequence interface: default
            ``max(1000, 10 * bottleneck_size)``, at least ``bottleneck_size``. Statistics interface:
            default 1000
        bottleneck_size: Population size during the bottleneck (at least 2). Default
            ``max(2, n_sequences // 4)`` for sequences and 100 for statistics
        final_size: Population size after recovery (statistics interface only, default
            ``initial_size``). The sequence sample is drawn at the end of the bottleneck, before the
            recovery, so ``final_size`` cannot affect it and is not used
        generations: Total generations to simulate (statistics interface only, default 100); not used
            for sequences, where ``bottleneck_duration`` sets the length of the bottleneck
        mutation_rate: Per-base mutation rate for the statistics interface. Not used for sequences:
            ``pre_bottleneck_diversity`` already fixes their ``4 * N * mu``
        rng: Random number generator for reproducibility
        n_sequences: Number of sequences to generate (selects the sequence interface)
        sequence_length: Length of each sequence (required with ``n_sequences``). Statistics interface:
            length of the locus in base pairs, which scales the mutation counts and the expected
            spectrum (default 1000)
        bottleneck_duration: Duration of the bottleneck in generations (sequence interface, default 10;
            0 means no bottleneck)
        pre_bottleneck_diversity: Per-site population mutation parameter, the expected nucleotide
            diversity, before the bottleneck (sequence interface, default 0.01)
        sample_size: Chromosomes in the sample behind the spectrum and Tajima's D (statistics interface,
            default 20, at least 4)
        n_replicates: Genealogies behind the spectrum and Tajima's D (statistics interface, default 2000,
            at least 2); their sampling error is about ``1 / sqrt(n_replicates)``

    Returns:
        List of sequences if n_sequences provided, otherwise dictionary with ``population_sizes``,
        ``diversity_trajectory``, ``mutation_trajectory``, ``bottleneck_start``, ``bottleneck_end``,
        ``total_mutations``, ``final_diversity`` and the coalescent summary keys above

    Raises:
        ValueError: If ``sequence_length`` is missing with ``n_sequences``, or the infinite-sites model
            cannot hold the variation (see :func:`generate_population_sequences`)
        ValidationError: If a numeric parameter is out of range
    """
    rng = coerce_rng(rng)

    # Alternative interface: return sequences
    if n_sequences is not None:
        if sequence_length is None:
            raise ValueError("sequence_length required when using n_sequences interface")
        if bottleneck_size is None:
            bottleneck_size = max(2, n_sequences // 4)
        if bottleneck_duration is None:
            bottleneck_duration = 10
        if initial_size is None:
            initial_size = max(1000, 10 * bottleneck_size)
        theta_per_site = _DEFAULT_THETA_PER_SITE if pre_bottleneck_diversity is None else pre_bottleneck_diversity

        validation.validate_range(n_sequences, min_val=2, name="n_sequences")
        validation.validate_range(sequence_length, min_val=1, name="sequence_length")
        validation.validate_range(bottleneck_size, min_val=2, name="bottleneck_size")
        validation.validate_range(initial_size, min_val=bottleneck_size, name="initial_size")
        validation.validate_range(bottleneck_duration, min_val=0, name="bottleneck_duration")
        _check_theta(theta_per_site)

        # Time in units of 2 * initial_size generations: the last bottleneck_duration generations
        # before the sample have relative size bottleneck_size / initial_size, then the size is 1.
        epochs = None
        if bottleneck_duration > 0:
            epochs = [
                PopulationEpoch(0.0, bottleneck_size / initial_size),
                PopulationEpoch(bottleneck_duration / (2.0 * initial_size), 1.0),
            ]
        return _coalescent_sequences(n_sequences, sequence_length, theta_per_site, rng=rng, epochs=epochs)

    # Original interface: return statistics
    if initial_size is None:
        initial_size = 1000
    if bottleneck_size is None:
        bottleneck_size = 100
    if final_size is None:
        final_size = initial_size
    if generations is None:
        generations = 100
    locus_length = 1000 if sequence_length is None else sequence_length

    validation.validate_range(initial_size, min_val=2, name="initial_size")
    validation.validate_range(bottleneck_size, min_val=2, name="bottleneck_size")
    validation.validate_range(final_size, min_val=2, name="final_size")
    validation.validate_range(generations, min_val=1, name="generations")
    validation.validate_range(mutation_rate, min_val=0.0, name="mutation_rate")
    validation.validate_range(locus_length, min_val=1, name="sequence_length")
    _check_summary_sample(sample_size, n_replicates)

    bottleneck_start = generations // 3
    bottleneck_end = 2 * generations // 3
    population_sizes = [
        initial_size if gen < bottleneck_start else bottleneck_size if gen < bottleneck_end else final_size
        for gen in range(generations)
    ]
    statistics = _trajectory_statistics(
        population_sizes,
        _epochs_from_sizes(population_sizes, initial_size),
        initial_size=initial_size,
        mutation_rate=mutation_rate,
        locus_length=locus_length,
        sample_size=sample_size,
        n_replicates=n_replicates,
        generator=_numpy_generator(rng),
    )
    return {
        "population_sizes": population_sizes,
        "bottleneck_start": bottleneck_start,
        "bottleneck_end": bottleneck_end,
        **statistics,
    }


def simulate_population_expansion(
    initial_size: int | None = None,
    final_size: int | None = None,
    expansion_time: int | None = None,
    *,
    mutation_rate: float = 1e-8,
    rng: random.Random | None = None,
    # Alternative interface for sequence generation
    n_sequences: int | None = None,
    sequence_length: int | None = None,
    expansion_factor: float | None = None,
    # Statistics interface: the sample behind the reported site frequency spectrum and Tajima's D
    sample_size: int = 20,
    n_replicates: int = 2000,
) -> List[str] | Dict[str, Any]:
    """Simulate exponential population expansion.

    Returns either aligned sequences (the sequence interface, selected by passing ``n_sequences``) or a
    statistics dictionary (the default interface).

    Sequence interface: the sequences are a neutral coalescent sample (see
    :func:`simulate_coalescent_haplotypes`; one non-recombining genealogy) from a population of
    ``initial_size`` individuals that grew exponentially by a factor of ``expansion_factor`` over
    the last ``expansion_time`` generations; the sample is drawn when the growth ends. In units of
    ``2 * initial_size`` generations the expansion lasts ``expansion_time / (2 * initial_size)`` and
    its growth rate is ``log(expansion_factor)`` divided by that. The large recent population barely
    coalesces, so the genealogy is star-like with long external branches: rare variants are in
    excess and Tajima's D is negative on average (about -0.9 for 20 sequences with the defaults),
    the textbook signature of a recent expansion. The ancestral population has the fixed per-site
    population mutation parameter 0.01 (``4 * N * mu`` per site), as the founders of the former
    generator did.

    Statistics interface: the deterministic size trajectory (exponential growth from ``initial_size`` to
    ``final_size`` over ``expansion_time`` generations, sizes truncated to whole individuals) and what the
    coalescent says about it, with the sample drawn when the growth ends (see
    :func:`_trajectory_statistics`): ``diversity_trajectory`` is the expected per-site nucleotide diversity
    after every generation (an exact recursion, started at the equilibrium ``4 * initial_size * mutation_rate``
    of the initial population), ``mutation_trajectory`` the cumulative count of new mutations in the whole
    population (Poisson, ``2 * size * mutation_rate * sequence_length`` per generation), and ``expected_sfs``,
    ``singleton_fraction`` and ``tajimas_d`` describe a sample of ``sample_size`` chromosomes simulated from
    the coalescent of that history: ``expected_sfs[i - 1]`` is the expected number of mutations in the locus
    carried by ``i`` chromosomes, averaged over ``n_replicates`` seeded genealogies (``neutral_sfs`` is the
    constant-size expectation ``theta / i`` for comparison), and ``tajimas_d`` is the mean Tajima's D of
    ``n_replicates`` loci that carry 50 mutations each (``tajimas_d_se`` is its standard error). Growth
    leaves an excess of rare variants (a ``singleton_fraction`` above ``neutral_singleton_fraction``) and a
    negative D.

    Args:
        initial_size: Population size before the expansion (default 100, both interfaces)
        final_size: Population size when the expansion ends (statistics interface default 10000). For
            sequences it is another way to set the factor, ``final_size / initial_size``, and takes
            precedence over ``expansion_factor``; it must be at least ``initial_size``
        expansion_time: Generations over which the expansion occurs (default 100, both interfaces)
        mutation_rate: Per-base mutation rate for the statistics interface. Not used for sequences,
            whose ancestral ``4 * N * mu`` per site is fixed at 0.01
        rng: Random number generator for reproducibility
        n_sequences: Number of sequences to generate (selects the sequence interface)
        sequence_length: Length of each sequence (required with ``n_sequences``). Statistics interface:
            length of the locus in base pairs, which scales the mutation counts and the expected
            spectrum (default 1000)
        expansion_factor: Factor by which the population expanded (sequence interface, default 10.0,
            at least 1; 1 means no expansion)
        sample_size: Chromosomes in the sample behind the spectrum and Tajima's D (statistics interface,
            default 20, at least 4)
        n_replicates: Genealogies behind the spectrum and Tajima's D (statistics interface, default 2000,
            at least 2); their sampling error is about ``1 / sqrt(n_replicates)``

    Returns:
        List of sequences if n_sequences provided, otherwise dictionary with ``population_sizes``,
        ``diversity_trajectory``, ``mutation_trajectory``, ``growth_rate``, ``expansion_time``,
        ``total_mutations``, ``final_diversity`` and the coalescent summary keys above

    Raises:
        ValueError: If ``sequence_length`` is missing with ``n_sequences``, or the infinite-sites model
            cannot hold the variation (see :func:`generate_population_sequences`)
        ValidationError: If a numeric parameter is out of range
    """
    rng = coerce_rng(rng)

    # Alternative interface: return sequences
    if n_sequences is not None:
        if sequence_length is None:
            raise ValueError("sequence_length required when using n_sequences interface")
        if initial_size is None:
            initial_size = 100
        if expansion_time is None:
            expansion_time = 100
        validation.validate_range(initial_size, min_val=2, name="initial_size")
        if final_size is not None:
            validation.validate_range(final_size, min_val=initial_size, name="final_size")
            expansion_factor = final_size / initial_size
        elif expansion_factor is None:
            expansion_factor = 10.0

        validation.validate_range(n_sequences, min_val=2, name="n_sequences")
        validation.validate_range(sequence_length, min_val=1, name="sequence_length")
        validation.validate_range(expansion_time, min_val=1, name="expansion_time")
        validation.validate_range(expansion_factor, min_val=1.0, name="expansion_factor")

        # Time in units of 2 * initial_size generations: the size falls from expansion_factor at the
        # sample to 1 at expansion_time ago, exponentially, and stays 1 before that.
        epochs = None
        if expansion_factor > 1.0:
            duration = expansion_time / (2.0 * initial_size)
            epochs = [
                PopulationEpoch(0.0, expansion_factor, math.log(expansion_factor) / duration),
                PopulationEpoch(duration, 1.0),
            ]
        return _coalescent_sequences(n_sequences, sequence_length, _DEFAULT_THETA_PER_SITE, rng=rng, epochs=epochs)

    # Original interface: return statistics
    if initial_size is None:
        initial_size = 100
    if final_size is None:
        final_size = 10000
    if expansion_time is None:
        expansion_time = 100
    locus_length = 1000 if sequence_length is None else sequence_length

    validation.validate_range(initial_size, min_val=2, name="initial_size")
    validation.validate_range(final_size, min_val=initial_size, name="final_size")
    validation.validate_range(expansion_time, min_val=1, name="expansion_time")
    validation.validate_range(mutation_rate, min_val=0.0, name="mutation_rate")
    validation.validate_range(locus_length, min_val=1, name="sequence_length")
    _check_summary_sample(sample_size, n_replicates)

    # Exponential growth model
    growth_rate = (np.log(final_size) - np.log(initial_size)) / expansion_time
    population_sizes = [
        min(int(initial_size * np.exp(growth_rate * gen)), final_size) for gen in range(expansion_time + 1)
    ]
    # In units of 2 * initial_size generations the growth lasts expansion_time / (2 * initial_size)
    epochs = _expansion_epochs(final_size / initial_size, expansion_time / (2.0 * initial_size))
    statistics = _trajectory_statistics(
        population_sizes,
        epochs,
        initial_size=initial_size,
        mutation_rate=mutation_rate,
        locus_length=locus_length,
        sample_size=sample_size,
        n_replicates=n_replicates,
        generator=_numpy_generator(rng),
    )
    return {
        "population_sizes": population_sizes,
        "growth_rate": growth_rate,
        "expansion_time": expansion_time,
        **statistics,
    }


#: Parameters of each demographic model of :func:`generate_site_frequency_spectrum`, with their defaults.
_SFS_MODEL_PARAMETERS: Dict[str, Dict[str, float]] = {
    "constant": {},
    "expansion": {"expansion_factor": 10.0, "expansion_time": 0.5, "n_replicates": 10_000},
    "bottleneck": {"bottleneck_size": 0.1, "bottleneck_duration": 0.05, "recovery_time": 0.0, "n_replicates": 10_000},
}


def generate_site_frequency_spectrum(
    n_samples: int | None = None,
    n_sites: int | None = None,
    *,
    sample_size: int | None = None,
    folded: bool = False,
    demographic_model: str = "constant",
    parameters: Dict[str, float] | None = None,
    rng: random.Random | None = None,
) -> List[int]:
    """Generate a site frequency spectrum under different demographic models.

    Every one of the ``n_sites`` polymorphic sites has the genealogy of its own (unlinked sites) from the
    neutral coalescent of the model, carries one mutation, and falls in the class of the number of the
    ``n_samples`` chromosomes that carry it. A site in class ``i`` is therefore drawn with the probability
    ``E[L_i] / E[L]``, the expected share of the branch length of a genealogy that supports ``i`` chromosomes
    (the mutation is more likely on a longer branch), and the counts are one multinomial draw:

    * ``"constant"``: constant size. The probabilities are the closed form ``(1 / i) / a_n``, with ``a_n`` the
      sum of ``1 / j`` for ``j < n_samples``.
    * ``"expansion"``: exponential growth by ``expansion_factor`` over the last ``expansion_time`` before the
      sample (constant size before). Growth leaves long external branches: an excess of rare variants and a
      negative Tajima's D.
    * ``"bottleneck"``: size shrinks to ``bottleneck_size`` of the ancestral size for ``bottleneck_duration``
      and ends ``recovery_time`` before the sample. Lineages coalesce inside it, which removes rare variants
      and leaves a positive Tajima's D until recovery has added new mutations to the long external branches.

    For the last two the probabilities are not available in closed form: they are the mean branch-length
    spectrum of ``n_replicates`` genealogies simulated from the history with the seeded coalescent
    (:func:`simulate_coalescent_spectra` and its time change), 10000 by default, so each probability carries
    a relative sampling error of about ``1 / sqrt(n_replicates)`` on top of the multinomial draw.

    Args:
        n_samples: Number of sampled chromosomes (alias: sample_size)
        n_sites: Number of polymorphic sites
        sample_size: Alias for n_samples
        folded: Whether to return folded SFS (default False)
        demographic_model: Demographic model ("constant", "expansion", "bottleneck")
        parameters: Model-specific parameters, with times in units of ``2 * N_ref`` generations and sizes
            relative to the ancestral size ``N_ref``; unknown names are rejected. ``"expansion"``:
            ``expansion_factor`` (final over ancestral size, at least 1, default 10), ``expansion_time``
            (default 0.5), ``n_replicates`` (default 10000). ``"bottleneck"``: ``bottleneck_size`` (in
            (0, 1], default 0.1), ``bottleneck_duration`` (default 0.05), ``recovery_time`` (default 0, the
            sample is taken as the bottleneck ends), ``n_replicates``. ``"constant"`` has none
        rng: Random number generator for reproducibility

    Returns:
        Site frequency spectrum as list of counts: ``n_samples - 1`` entries (entry ``i - 1`` counts the
        sites carried by ``i`` chromosomes), or the ``n_samples // 2`` folded classes

    Raises:
        ValueError: If parameters are invalid
        ValidationError: If the model is unknown or ``parameters`` holds a name the model does not take
    """
    # Handle parameter aliases
    if sample_size is not None and n_samples is None:
        n_samples = sample_size
    if n_samples is None:
        raise ValueError("n_samples (or sample_size) is required")
    if n_sites is None:
        raise ValueError("n_sites is required")

    validation.validate_range(n_samples, min_val=2, name="n_samples")
    validation.validate_range(n_sites, min_val=1, name="n_sites")

    if demographic_model not in _SFS_MODEL_PARAMETERS:
        raise errors.ValidationError(f"Unknown demographic model: {demographic_model}")
    defaults = _SFS_MODEL_PARAMETERS[demographic_model]
    unknown = sorted(set(parameters or {}) - set(defaults))
    if unknown:
        raise errors.ValidationError(
            f"unknown parameter(s) {unknown} for the {demographic_model!r} model; it takes {sorted(defaults) or 'none'}"
        )
    values = {**defaults, **(parameters or {})}
    generator = _numpy_generator(coerce_rng(rng))

    epochs: List[PopulationEpoch] | None = None  # None: constant size, whose spectrum has a closed form
    replicates = 1
    if demographic_model != "constant":
        for name, value in values.items():
            if not math.isfinite(value):
                raise errors.ValidationError(f"{name} must be a finite number, got {value}")
        replicates = int(values["n_replicates"])
        if replicates != values["n_replicates"] or replicates < 1:
            raise errors.ValidationError(f"n_replicates must be a positive integer, got {values['n_replicates']}")
        if demographic_model == "expansion":
            validation.validate_range(values["expansion_factor"], min_val=1.0, name="expansion_factor")
            validation.validate_range(values["expansion_time"], min_val=0.0, name="expansion_time")
            if values["expansion_factor"] > 1.0:
                if values["expansion_time"] <= 0.0:
                    raise errors.ValidationError("expansion_time must be positive when expansion_factor exceeds 1")
                epochs = _expansion_epochs(values["expansion_factor"], values["expansion_time"])
        else:
            if not 0.0 < values["bottleneck_size"] <= 1.0:
                raise errors.ValidationError(f"bottleneck_size must be in (0, 1], got {values['bottleneck_size']}")
            validation.validate_range(values["bottleneck_duration"], min_val=0.0, name="bottleneck_duration")
            validation.validate_range(values["recovery_time"], min_val=0.0, name="recovery_time")
            if values["bottleneck_duration"] > 0.0 and values["bottleneck_size"] < 1.0:
                epochs = _bottleneck_epochs(
                    values["bottleneck_size"], values["bottleneck_duration"], values["recovery_time"]
                )
    if epochs is None:
        weights = 1.0 / np.arange(1, n_samples, dtype=float)
    else:
        weights = _expected_branch_spectrum(n_samples, epochs, replicates, generator)

    sfs: List[int] = generator.multinomial(n_sites, weights / weights.sum()).tolist()

    # Fold SFS if requested
    if folded:
        folded_sfs = []
        mid = n_samples // 2
        for i in range(mid):
            # Combine frequency class i and n-1-i (complement frequencies)
            j = n_samples - 2 - i
            if i < len(sfs) and j < len(sfs) and i != j:
                folded_sfs.append(sfs[i] + sfs[j])
            elif i < len(sfs):
                # Middle bin (only when i == j) or edge case
                folded_sfs.append(sfs[i])
        return folded_sfs

    return sfs


def _adjacent_ld_conditionals(p_prev: float, p_cur: float, correlation: float) -> Tuple[float, float, bool]:
    """Allele-1 probabilities of a site given the allele at the previous site.

    The pair gets correlation ``correlation`` (so r-squared is ``correlation`` squared) while the
    marginal frequency of the current site stays exactly ``p_cur``. With the disequilibrium
    D = r * sqrt(p_prev (1 - p_prev) p_cur (1 - p_cur)):

        P(cur=1 | prev=1) = p_cur + D / p_prev
        P(cur=1 | prev=0) = p_cur - D / (1 - p_prev)

    A positive D cannot exceed min(p_prev (1 - p_cur), (1 - p_prev) p_cur), because both
    probabilities must stay within [0, 1]; a larger request is capped there.

    Returns:
        (P(cur=1 | prev=1), P(cur=1 | prev=0), whether the request had to be capped)
    """
    if p_prev <= 0.0 or p_prev >= 1.0 or p_cur <= 0.0 or p_cur >= 1.0:
        # A monomorphic site has no variance to correlate with.
        return p_cur, p_cur, False
    d = correlation * (p_prev * (1.0 - p_prev) * p_cur * (1.0 - p_cur)) ** 0.5
    d_max = min(p_prev * (1.0 - p_cur), (1.0 - p_prev) * p_cur)
    capped = d > d_max * (1.0 + 1e-12)
    d = min(d, d_max)
    return min(1.0, p_cur + d / p_prev), max(0.0, p_cur - d / (1.0 - p_prev)), capped


#: The models of :func:`generate_linkage_disequilibrium_data`.
_LD_MODELS = ("markov", "coalescent")
#: Region length per SNP, in base pairs, of the coalescent model when no ``sequence_length`` is given.
_DEFAULT_BP_PER_SNP = 100
#: Per-generation recombination probability per base pair of the coalescent model by default (a human average).
_DEFAULT_RECOMBINATION_RATE = 1e-8
#: Effective population size of the coalescent model by default.
_DEFAULT_EFFECTIVE_SIZE = 10_000.0


def generate_linkage_disequilibrium_data(
    n_individuals: int,
    n_snps: int | None = None,
    *,
    n_sites: int | None = None,
    recombination_rate: float | None = None,
    selection_coefficient: float = 0.0,
    r_squared_target: float | None = None,
    allele_frequencies: List[float] | None = None,
    rng: random.Random | None = None,
    model: str = "markov",
    effective_size: float | None = None,
    sequence_length: int | None = None,
) -> List[List[int]]:
    """Generate 0/1 allele data with linkage disequilibrium between sites.

    ``model`` chooses how the linkage disequilibrium arises. Both models return the same kind of matrix, one
    row per individual (a haploid sequence) and one column per site, holding 0/1 alleles.

    ``model="markov"`` (the default): sites are generated left to right as a first-order Markov chain.
    Every site keeps its own allele frequency exactly, and each adjacent pair is given the
    disequilibrium D = r * sqrt(p1 (1 - p1) p2 (1 - p2)) with r = sqrt(``r_squared_target``), so
    the expected r-squared between adjacent sites equals ``r_squared_target`` (a sample r-squared
    scatters around it by about 2 r (1 - r-squared) / sqrt(``n_individuals``)). For sites with
    equal frequencies the correlation falls to r ** k between sites k apart. The model has no recombination
    and no physical distance, so ``recombination_rate`` has no effect on it.

    ``model="coalescent"``: the sample is drawn from the neutral coalescent with recombination
    (:func:`simulate_recombinant_haplotypes`, Hudson's ancestral recombination graph) over a region of
    ``sequence_length`` base pairs in a population of ``effective_size`` diploid individuals, with a
    recombination probability of ``recombination_rate`` per base pair per generation, so the population
    recombination parameter is ``rho = 4 * effective_size * recombination_rate * sequence_length``. Exactly
    ``n_snps`` mutations are placed on the genealogy (the segregating-sites conditioning of Hudson's
    ``ms -s``; the mutation rate is therefore not an input), at positions along the region that follow the
    marginal trees; the columns are those segregating sites in order of position, and allele 1 is the derived
    allele. Allele frequencies and r-squared are outcomes of the genealogy, not inputs: they follow the
    neutral site frequency spectrum (rare derived alleles dominate), and r-squared decays with the distance
    between sites as recombination breaks the shared genealogy (the quantity ``E[D^2] / E[p (1 - p) q (1 - q)]``
    is close to ``(10 + rho_d) / ((2 + rho_d) (11 + rho_d))`` for two sites ``rho_d`` apart, with ``rho_d`` the
    recombination parameter between them). With ``recombination_rate = 0`` all columns come from one tree and
    any two of them show at most three of the four gametes. To get the SNP positions, or a mutation rate
    instead of a fixed number of sites, use :func:`simulate_recombinant_haplotypes`.

    Args:
        n_individuals: Number of individuals (rows), at least 2
        n_snps: Number of sites (columns), at least 2 (alias: n_sites). Taken from the length of
            ``allele_frequencies`` when that is given. Under the coalescent model, the exact number of
            segregating sites simulated
        n_sites: Alias for n_snps
        recombination_rate: Probability of a recombination breakpoint per base pair per generation,
            default 1e-8. Used by ``model="coalescent"`` only; ``model="markov"`` accepts and validates it
            (it must not be negative) but it does not change the output, and says so in the log
        selection_coefficient: Selection coefficient (0 = neutral). When non-zero, its absolute
            value replaces sqrt(``r_squared_target``) as the correlation r between adjacent
            sites (capped at 1); a heuristic coupling, not a model of selection. Markov model only
        r_squared_target: Target r-squared between adjacent sites, within [0, 1]. The default,
            0.25, is the r-squared this generator produced before a target was honoured. A pair
            whose allele frequencies cannot reach the target (a rare allele next to a common one
            limits r-squared) is generated at the largest attainable r-squared and a warning is
            logged. Markov model only
        allele_frequencies: Frequency of allele 1 at each site, each within [0, 1] (default 0.5
            at every site). Markov model only
        rng: Random number generator for reproducibility
        model: ``"markov"`` or ``"coalescent"``
        effective_size: Effective population size ``N_e`` of the coalescent model (default 10000), which
            sets how many generations of recombination one unit of coalescent time stands for. Coalescent
            model only
        sequence_length: Length of the simulated region in base pairs, at least ``n_snps`` (default
            100 base pairs per SNP). Coalescent model only

    Returns:
        List of lists with one row per individual and one column per site, holding 0/1 alleles

    Raises:
        ValueError: If neither n_snps nor n_sites is given
        ValidationError: If a parameter is out of range, r_squared_target, an allele frequency
            or selection_coefficient is not a valid number, or a parameter is given that the chosen
            model does not use (``effective_size`` and ``sequence_length`` under ``"markov"``,
            ``r_squared_target``, ``selection_coefficient`` and ``allele_frequencies`` under ``"coalescent"``)
    """
    if model not in _LD_MODELS:
        raise errors.ValidationError(f"model must be one of {_LD_MODELS}, got {model!r}")

    # Handle parameter aliases
    if n_sites is not None and n_snps is None:
        n_snps = n_sites
    if allele_frequencies is not None:
        n_snps = len(allele_frequencies)
    if n_snps is None:
        raise ValueError("n_snps (or n_sites) is required")

    validation.validate_range(n_individuals, min_val=2, name="n_individuals")
    validation.validate_range(n_snps, min_val=2, name="n_snps")
    if recombination_rate is not None:
        if not math.isfinite(recombination_rate):
            raise errors.ValidationError(f"recombination_rate must be finite, got {recombination_rate}")
        validation.validate_range(recombination_rate, min_val=0.0, name="recombination_rate")
    if not np.isfinite(selection_coefficient):
        raise errors.ValidationError(f"selection_coefficient must be finite, got {selection_coefficient}")

    if model == "coalescent":
        return _coalescent_ld_data(
            n_individuals,
            n_snps,
            recombination_rate=_DEFAULT_RECOMBINATION_RATE if recombination_rate is None else recombination_rate,
            effective_size=effective_size,
            sequence_length=sequence_length,
            unused={
                "r_squared_target": r_squared_target is not None,
                "selection_coefficient": selection_coefficient != 0,
                "allele_frequencies": allele_frequencies is not None,
            },
            rng=rng,
        )

    for name, value in (("effective_size", effective_size), ("sequence_length", sequence_length)):
        if value is not None:
            raise errors.ValidationError(f"{name} is used by model='coalescent' only, not by model='markov'")
    if recombination_rate is not None:
        logger.warning(
            "recombination_rate=%g has no effect on model='markov', which has no recombination step; "
            "use model='coalescent' for linkage disequilibrium that depends on it",
            recombination_rate,
        )
    if r_squared_target is not None and not 0.0 <= r_squared_target <= 1.0:
        raise errors.ValidationError(f"r_squared_target must be within [0, 1], got {r_squared_target}")
    if allele_frequencies is not None and not all(0.0 <= freq <= 1.0 for freq in allele_frequencies):
        raise errors.ValidationError("allele_frequencies must all be within [0, 1]")

    rng = coerce_and_seed(rng)

    frequencies = [0.5] * n_snps if allele_frequencies is None else [float(freq) for freq in allele_frequencies]

    # Correlation r between adjacent sites. The target is an r-squared, so r is its square root;
    # a non-zero selection coefficient replaces it and |s| is used directly as r.
    if selection_coefficient != 0:
        correlation = min(abs(selection_coefficient), 1.0)
    else:
        correlation = (0.25 if r_squared_target is None else r_squared_target) ** 0.5

    # P(allele 1 | previous allele 1) and P(allele 1 | previous allele 0) at each site; they depend
    # only on the allele frequencies and the target, so every individual shares them.
    given_one = [0.0] * n_snps
    given_zero = [0.0] * n_snps
    n_capped = 0
    for pos in range(1, n_snps):
        given_one[pos], given_zero[pos], capped = _adjacent_ld_conditionals(
            frequencies[pos - 1], frequencies[pos], correlation
        )
        n_capped += capped
    if n_capped:
        logger.warning(
            "r-squared %.4g is not attainable at the given allele frequencies for %d of %d adjacent site pairs; "
            "those pairs use the largest attainable r-squared",
            correlation**2,
            n_capped,
            n_snps - 1,
        )

    genotype_matrix: List[List[int]] = []
    for _ in range(n_individuals):
        allele = 1 if rng.random() < frequencies[0] else 0
        row = [allele]
        for pos in range(1, n_snps):
            allele = 1 if rng.random() < (given_one[pos] if allele else given_zero[pos]) else 0
            row.append(allele)
        genotype_matrix.append(row)

    return genotype_matrix


def _coalescent_ld_data(
    n_individuals: int,
    n_snps: int,
    *,
    recombination_rate: float,
    effective_size: float | None,
    sequence_length: int | None,
    unused: Dict[str, bool],
    rng: random.Random | None,
) -> List[List[int]]:
    """The ``model="coalescent"`` branch of :func:`generate_linkage_disequilibrium_data`."""
    for name, given in unused.items():
        if given:
            raise errors.ValidationError(
                f"{name} is an input of model='markov' only: under model='coalescent' it is an outcome of the "
                "genealogy"
            )
    size = _DEFAULT_EFFECTIVE_SIZE if effective_size is None else effective_size
    if not math.isfinite(size) or size <= 0:
        raise errors.ValidationError(f"effective_size must be a positive finite number, got {size}")
    length = _DEFAULT_BP_PER_SNP * n_snps if sequence_length is None else sequence_length
    validation.validate_range(length, min_val=n_snps, name="sequence_length")
    if not math.isfinite(length):
        raise errors.ValidationError(f"sequence_length must be finite, got {length}")

    rho = 4.0 * size * recombination_rate * length
    sample = simulate_recombinant_haplotypes(
        n_individuals, length, segregating_sites=n_snps, rho=rho, rng=_numpy_generator(rng)
    )
    matrix: List[List[int]] = sample.haplotypes.tolist()
    return matrix


def simulate_admixture(
    n_populations: int,
    population_sizes: List[int],
    admixture_proportions: np.ndarray,
    generations: int,
    *,
    migration_rate: float = 0.01,
    rng: random.Random | None = None,
) -> Dict[str, Any]:
    """Simulate population admixture: migration among populations and genetic drift at one locus.

    Each population starts with its own allele frequency, drawn uniformly from [0.1, 0.9]. Every generation
    a fraction ``migration_rate`` of each population is replaced by migrants whose allele frequency is the
    mixture ``admixture_proportions @ p`` of the populations' frequencies ``p`` (row ``i`` holds the
    proportions in which population ``i`` receives its migrants from every population), so the frequencies
    become ``(1 - m) p + m * admixture_proportions @ p``; ``migration_rate = 1`` replaces everyone. Genetic
    drift then resamples the ``2 * population_sizes[i]`` chromosomes of population ``i`` binomially, so a
    small population drifts more (the variance of its frequency change is ``p (1 - p) / (2 N)``) and a very
    large one hardly at all.

    Args:
        n_populations: Number of ancestral populations
        population_sizes: Number of diploid individuals in each ancestral population (positive integers)
        admixture_proportions: Mixing proportions (n_populations x n_populations): non-negative, every
            row summing to 1
        generations: Generations since admixture
        migration_rate: Fraction of every population replaced by migrants each generation, within [0, 1]
        rng: Random number generator for reproducibility

    Returns:
        Dictionary with ``ancestral_frequencies``, the ``frequency_trajectory`` (one array of the
        populations' frequencies per generation, ``generations + 1`` of them), ``final_frequencies``,
        ``admixture_trajectory`` (per generation, a matrix whose entry ``[i, j]`` is the expected share of
        population ``i``'s ancestry that comes from the original population ``j``: the identity at the
        start, ``((1 - m) I + m * admixture_proportions) ** t`` after ``t`` generations, drift aside),
        ``generations`` and ``migration_rate``

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_populations, min_val=2, name="n_populations")
    validation.validate_type(population_sizes, list, "population_sizes")
    validation.validate_range(generations, min_val=0, name="generations")
    validation.validate_range(migration_rate, min_val=0.0, max_val=1.0, name="migration_rate")

    if len(population_sizes) != n_populations:
        raise errors.ValidationError(f"population_sizes must have length {n_populations}")
    for size in population_sizes:
        validation.validate_range(size, min_val=1, name="population_sizes entries")

    proportions = np.asarray(admixture_proportions, dtype=float)
    if proportions.shape != (n_populations, n_populations):
        raise errors.ValidationError(f"admixture_proportions must have shape {(n_populations, n_populations)}")
    if not np.isfinite(proportions).all() or (proportions < 0).any() or not np.allclose(proportions.sum(axis=1), 1.0):
        raise errors.ValidationError("admixture_proportions must be non-negative with every row summing to 1")

    rng = coerce_rng(rng)
    generator = _numpy_generator(rng)
    chromosomes = 2 * np.asarray(population_sizes)

    # Initialize ancestral populations
    ancestral_frequencies = []
    for i in range(n_populations):
        freq = rng.uniform(0.1, 0.9)  # Allele frequency in ancestral population
        ancestral_frequencies.append(freq)

    # Migration mixes the populations' frequencies; drift resamples every population's chromosomes
    mixing = (1.0 - migration_rate) * np.eye(n_populations) + migration_rate * proportions
    current_frequencies = np.array(ancestral_frequencies)
    frequency_trajectory = [current_frequencies.copy()]
    ancestry = np.eye(n_populations)
    admixture_trajectory = [ancestry.copy()]

    for gen in range(generations):
        mixed = np.clip(mixing @ current_frequencies, 0.0, 1.0)
        current_frequencies = generator.binomial(chromosomes, mixed) / chromosomes
        frequency_trajectory.append(current_frequencies.copy())
        ancestry = mixing @ ancestry
        admixture_trajectory.append(ancestry.copy())

    return {
        "ancestral_frequencies": ancestral_frequencies,
        "frequency_trajectory": frequency_trajectory,
        "admixture_trajectory": admixture_trajectory,
        "final_frequencies": current_frequencies,
        "generations": generations,
        "migration_rate": migration_rate,
    }


def simulate_selection(
    genotype_matrix: np.ndarray,
    fitness_effects: np.ndarray,
    generations: int,
    *,
    selection_strength: float = 0.1,
    rng: random.Random | None = None,
) -> Dict[str, Any]:
    """Simulate natural selection on genotypes.

    Args:
        genotype_matrix: Initial genotype matrix (individuals x SNPs)
        fitness_effects: Fitness effect for each genotype (3 values per SNP: AA, AB, BB)
        generations: Generations to simulate
        selection_strength: Overall strength of selection
        rng: Random number generator for reproducibility

    Returns:
        Dictionary with selection simulation results

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_type(genotype_matrix, np.ndarray, "genotype_matrix")
    validation.validate_type(fitness_effects, np.ndarray, "fitness_effects")
    validation.validate_range(generations, min_val=1, name="generations")
    validation.validate_range(selection_strength, min_val=0.0, name="selection_strength")

    n_individuals, n_snps = genotype_matrix.shape

    if fitness_effects.shape != (n_snps, 3):
        raise errors.ValidationError(f"fitness_effects must have shape ({n_snps}, 3)")
    if not (1 + selection_strength * fitness_effects > 0).all():
        raise errors.ValidationError(
            "1 + selection_strength * fitness_effects must be positive for every genotype: a fitness cannot be negative"
        )

    rng = coerce_and_seed(rng)

    # Track allele frequencies over time
    allele_frequencies = []
    current_genotypes = genotype_matrix.copy()

    for gen in range(generations + 1):
        # Calculate current allele frequencies
        freqs = []
        for snp in range(n_snps):
            genotypes = current_genotypes[:, snp]
            alt_allele_count = np.sum(genotypes)  # 0=AA, 1=AB, 2=BB
            freq = alt_allele_count / (2 * n_individuals)
            freqs.append(freq)
        allele_frequencies.append(freqs)

        if gen == generations:
            break

        # Calculate fitness for each individual
        fitness_values = np.ones(n_individuals)

        for snp in range(n_snps):
            genotypes = current_genotypes[:, snp]
            fitness_snp = np.zeros(n_individuals)

            for i, genotype in enumerate(genotypes):
                fitness_snp[i] = fitness_effects[snp, genotype]

            fitness_values *= 1 + selection_strength * fitness_snp

        # Normalize fitness (Wright-Fisher model)
        fitness_values /= np.mean(fitness_values)

        # Generate next generation
        next_genotypes = np.zeros_like(current_genotypes)

        for snp in range(n_snps):
            # Sample genotypes proportional to fitness
            genotypes = current_genotypes[:, snp]

            # Create weighted sampling
            weights = fitness_values.tolist()
            selected_indices = rng.choices(range(n_individuals), weights=weights, k=n_individuals)

            next_genotypes[:, snp] = genotypes[selected_indices]

        current_genotypes = next_genotypes

    return {
        "initial_genotypes": genotype_matrix,
        "final_genotypes": current_genotypes,
        "allele_frequencies": allele_frequencies,
        "fitness_effects": fitness_effects,
        "generations": generations,
        "selection_strength": selection_strength,
    }

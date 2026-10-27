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
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Any, Dict, List, Sequence, Tuple, Union

import numpy as np

from metainformant.core.data import validation
from metainformant.core.utils import errors, logging
from metainformant.simulation.rng import coerce_and_seed, coerce_rng, seed_numpy_from_rng

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
        n_sites: Number of polymorphic sites, at most ``length``. Defaults to
            Watterson's expectation ``theta * length * a_n`` for the
            ``n_pop1 + n_pop2`` sequences (at least 1).
        rng: Random number generator for reproducibility
        sequence_length: Alias for length

    Returns:
        Tuple of (population1_sequences, population2_sequences)

    Raises:
        ValueError: If parameters are invalid
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

    rng = coerce_and_seed(rng)

    ancestral = [rng.choice("ACGT") for _ in range(length)]
    populations = [[list(ancestral) for _ in range(n_pop1)], [list(ancestral) for _ in range(n_pop2)]]
    for position in sorted(rng.sample(range(length), min(n_sites, length))):
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
        hwe: Whether to generate genotypes under HWE (True) or not (False)
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

    # Convert hwe boolean to hwe_deviation
    if not hwe:
        hwe_deviation = 0.3  # Some deviation from HWE

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

    Statistics interface: a deterministic size trajectory with heuristic diversity and Poisson
    mutation counts. It is not a coalescent simulation. The bottleneck occupies the middle third of
    ``generations``.

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
        sequence_length: Length of each sequence (required with ``n_sequences``)
        bottleneck_duration: Duration of the bottleneck in generations (sequence interface, default 10;
            0 means no bottleneck)
        pre_bottleneck_diversity: Per-site population mutation parameter, the expected nucleotide
            diversity, before the bottleneck (sequence interface, default 0.01)

    Returns:
        List of sequences if n_sequences provided, otherwise dictionary with simulation results

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

    validation.validate_range(initial_size, min_val=2, name="initial_size")
    validation.validate_range(bottleneck_size, min_val=2, name="bottleneck_size")
    validation.validate_range(final_size, min_val=2, name="final_size")
    validation.validate_range(generations, min_val=1, name="generations")
    validation.validate_range(mutation_rate, min_val=0.0, name="mutation_rate")

    seed_numpy_from_rng(rng)

    # Simulate population size changes
    bottleneck_start = generations // 3
    bottleneck_end = 2 * generations // 3

    population_sizes = []

    for gen in range(generations):
        if gen < bottleneck_start:
            size = initial_size
        elif gen < bottleneck_end:
            size = bottleneck_size
        else:
            size = final_size
        population_sizes.append(size)

    # Simulate genetic diversity changes (simplified)
    initial_diversity = 1.0
    diversity_trajectory = []

    for size in population_sizes:
        diversity = initial_diversity * (size / initial_size)
        diversity = max(diversity, 0.01)
        diversity_trajectory.append(diversity)

    # Simulate mutations
    total_mutations = 0
    mutation_trajectory = []

    for gen, size in enumerate(population_sizes):
        expected_mutations = size * 1000 * mutation_rate
        mutations_this_gen = np.random.poisson(expected_mutations)
        total_mutations += mutations_this_gen
        mutation_trajectory.append(total_mutations)

    return {
        "population_sizes": population_sizes,
        "diversity_trajectory": diversity_trajectory,
        "mutation_trajectory": mutation_trajectory,
        "bottleneck_start": bottleneck_start,
        "bottleneck_end": bottleneck_end,
        "total_mutations": total_mutations,
        "final_diversity": diversity_trajectory[-1],
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

    Statistics interface: a deterministic size trajectory with heuristic diversity and Poisson
    mutation counts. It is not a coalescent simulation.

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
        sequence_length: Length of each sequence (required with ``n_sequences``)
        expansion_factor: Factor by which the population expanded (sequence interface, default 10.0,
            at least 1; 1 means no expansion)

    Returns:
        List of sequences if n_sequences provided, otherwise dictionary with simulation results

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

    validation.validate_range(initial_size, min_val=2, name="initial_size")
    validation.validate_range(final_size, min_val=initial_size, name="final_size")
    validation.validate_range(expansion_time, min_val=1, name="expansion_time")
    validation.validate_range(mutation_rate, min_val=0.0, name="mutation_rate")

    seed_numpy_from_rng(rng)

    # Exponential growth model
    growth_rate = (np.log(final_size) - np.log(initial_size)) / expansion_time

    population_sizes = []
    diversity_trajectory = []
    mutation_trajectory = []

    total_mutations = 0
    current_diversity = 1.0

    for gen in range(expansion_time + 1):
        size = int(initial_size * np.exp(growth_rate * gen))
        size = min(size, final_size)
        population_sizes.append(size)

        diversity_increase = (size - initial_size) / (final_size - initial_size) * 0.1
        current_diversity = min(1.0, current_diversity + diversity_increase)
        diversity_trajectory.append(current_diversity)

        expected_mutations = size * 1000 * mutation_rate
        mutations_this_gen = np.random.poisson(expected_mutations)
        total_mutations += mutations_this_gen
        mutation_trajectory.append(total_mutations)

    return {
        "population_sizes": population_sizes,
        "diversity_trajectory": diversity_trajectory,
        "mutation_trajectory": mutation_trajectory,
        "growth_rate": growth_rate,
        "expansion_time": expansion_time,
        "total_mutations": total_mutations,
        "final_diversity": diversity_trajectory[-1],
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

    Args:
        n_samples: Number of sampled chromosomes (alias: sample_size)
        n_sites: Number of polymorphic sites
        sample_size: Alias for n_samples
        folded: Whether to return folded SFS (default False)
        demographic_model: Demographic model ("constant", "expansion", "bottleneck")
        parameters: Model-specific parameters
        rng: Random number generator for reproducibility

    Returns:
        Site frequency spectrum as list of counts

    Raises:
        ValueError: If parameters are invalid
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

    if demographic_model not in ["constant", "expansion", "bottleneck"]:
        raise errors.ValidationError(f"Unknown demographic model: {demographic_model}")

    rng = coerce_and_seed(rng)

    if parameters is None:
        parameters = {}

    # Unfolded SFS has n-1 bins
    sfs = [0] * (n_samples - 1)

    for site in range(n_sites):
        if demographic_model == "constant":
            # Neutral coalescent - 1/i distribution
            weights = [1.0 / i for i in range(1, n_samples)]
            total = sum(weights)
            weights = [w / total for w in weights]
            freq = rng.choices(range(1, n_samples), weights=weights)[0]
        elif demographic_model == "expansion":
            # Favor low frequency variants (zipf requires exponent > 1; for
            # alpha <= 1 fall back to the neutral 1/i weighting)
            alpha = parameters.get("alpha", 1.0)
            if alpha > 1.0:
                freq = min(int(np.random.zipf(alpha)), n_samples - 1)
            else:
                weights = [1.0 / i for i in range(1, n_samples)]
                total = sum(weights)
                weights = [w / total for w in weights]
                freq = rng.choices(range(1, n_samples), weights=weights)[0]
        elif demographic_model == "bottleneck":
            # Favor high frequency variants
            bottleneck_strength = parameters.get("bottleneck_strength", 0.1)
            if rng.random() < bottleneck_strength:
                freq = rng.randint(n_samples // 2, n_samples - 1)
            else:
                freq = rng.randint(1, max(2, n_samples // 2))

        # Ensure valid frequency
        freq = max(1, min(freq, n_samples - 1))

        sfs[freq - 1] += 1

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


def generate_linkage_disequilibrium_data(
    n_individuals: int,
    n_snps: int | None = None,
    *,
    n_sites: int | None = None,
    recombination_rate: float = 1e-8,
    selection_coefficient: float = 0.0,
    r_squared_target: float | None = None,
    allele_frequencies: List[float] | None = None,
    rng: random.Random | None = None,
) -> List[List[int]]:
    """Generate 0/1 allele data with linkage disequilibrium between adjacent sites.

    Sites are generated left to right as a first-order Markov chain, one row per individual. Every
    site keeps its own allele frequency exactly, and each adjacent pair is given the
    disequilibrium D = r * sqrt(p1 (1 - p1) p2 (1 - p2)) with r = sqrt(``r_squared_target``), so
    the expected r-squared between adjacent sites equals ``r_squared_target`` (a sample r-squared
    scatters around it by about 2 r (1 - r-squared) / sqrt(``n_individuals``)). For sites with
    equal frequencies the correlation falls to r ** k between sites k apart.

    Args:
        n_individuals: Number of individuals (rows), at least 2
        n_snps: Number of sites (columns), at least 2 (alias: n_sites). Taken from the length of
            ``allele_frequencies`` when that is given
        n_sites: Alias for n_snps
        recombination_rate: Accepted and validated (it must not be negative) but it does not change
            the output. The model has no recombination step: LD between adjacent sites is set by
            ``r_squared_target`` alone and decays geometrically with distance
        selection_coefficient: Selection coefficient (0 = neutral). When non-zero, its absolute
            value replaces sqrt(``r_squared_target``) as the correlation r between adjacent
            sites (capped at 1); a heuristic coupling, not a model of selection
        r_squared_target: Target r-squared between adjacent sites, within [0, 1]. The default,
            0.25, is the r-squared this generator produced before a target was honoured. A pair
            whose allele frequencies cannot reach the target (a rare allele next to a common one
            limits r-squared) is generated at the largest attainable r-squared and a warning is
            logged
        allele_frequencies: Frequency of allele 1 at each site, each within [0, 1] (default 0.5
            at every site)
        rng: Random number generator for reproducibility

    Returns:
        List of lists with one row per individual and one column per site, holding 0/1 alleles

    Raises:
        ValueError: If neither n_snps nor n_sites is given
        ValidationError: If a parameter is out of range, or r_squared_target, an allele frequency
            or selection_coefficient is not a valid number
    """
    # Handle parameter aliases
    if n_sites is not None and n_snps is None:
        n_snps = n_sites
    if allele_frequencies is not None:
        n_snps = len(allele_frequencies)
    if n_snps is None:
        raise ValueError("n_snps (or n_sites) is required")

    validation.validate_range(n_individuals, min_val=2, name="n_individuals")
    validation.validate_range(n_snps, min_val=2, name="n_snps")
    validation.validate_range(recombination_rate, min_val=0.0, name="recombination_rate")
    if r_squared_target is not None and not 0.0 <= r_squared_target <= 1.0:
        raise errors.ValidationError(f"r_squared_target must be within [0, 1], got {r_squared_target}")
    if allele_frequencies is not None and not all(0.0 <= freq <= 1.0 for freq in allele_frequencies):
        raise errors.ValidationError("allele_frequencies must all be within [0, 1]")
    if not np.isfinite(selection_coefficient):
        raise errors.ValidationError(f"selection_coefficient must be finite, got {selection_coefficient}")

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


def simulate_admixture(
    n_populations: int,
    population_sizes: List[int],
    admixture_proportions: np.ndarray,
    generations: int,
    *,
    migration_rate: float = 0.01,
    rng: random.Random | None = None,
) -> Dict[str, Any]:
    """Simulate population admixture.

    Args:
        n_populations: Number of ancestral populations
        population_sizes: Size of each ancestral population
        admixture_proportions: Mixing proportions (n_populations x n_populations)
        generations: Generations since admixture
        migration_rate: Migration rate between populations
        rng: Random number generator for reproducibility

    Returns:
        Dictionary with admixture simulation results

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_populations, min_val=2, name="n_populations")
    validation.validate_type(population_sizes, list, "population_sizes")
    validation.validate_range(generations, min_val=0, name="generations")
    validation.validate_range(migration_rate, min_val=0.0, max_val=1.0, name="migration_rate")

    if len(population_sizes) != n_populations:
        raise errors.ValidationError(f"population_sizes must have length {n_populations}")

    if admixture_proportions.shape != (n_populations, n_populations):
        raise errors.ValidationError(f"admixture_proportions must have shape {(n_populations, n_populations)}")

    rng = coerce_and_seed(rng)

    # Initialize ancestral populations
    ancestral_frequencies = []
    for i in range(n_populations):
        freq = rng.uniform(0.1, 0.9)  # Allele frequency in ancestral population
        ancestral_frequencies.append(freq)

    # Simulate admixture over generations
    current_frequencies = np.array(ancestral_frequencies)
    frequency_trajectory = [current_frequencies.copy()]

    for gen in range(generations):
        # Apply migration/admixture
        new_frequencies = admixture_proportions @ current_frequencies

        # Add some drift
        for i in range(n_populations):
            drift_effect = np.random.normal(0, 0.01)
            new_frequencies[i] += drift_effect
            new_frequencies[i] = np.clip(new_frequencies[i], 0.0, 1.0)

        current_frequencies = new_frequencies
        frequency_trajectory.append(current_frequencies.copy())

    # Calculate admixture proportions over time
    admixture_trajectory = []
    for freqs in frequency_trajectory:
        admixture_trajectory.append(freqs / np.sum(freqs))

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

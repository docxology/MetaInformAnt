"""Simulated null distributions for the neutrality tests (DNA-simulation integration).

The neutrality tests of :mod:`metainformant.dna.population.analysis` (Tajima's D, Fu and Li's D* and F*,
Fay and Wu's H) have their exact null distribution only through the neutral coalescent, whose simulator belongs
to the simulation domain (:func:`metainformant.simulation.models.popgen.simulate_coalescent_spectra`). Package
boundaries let a domain depend on another one only from an ``integration`` module, so this adapter is where the
two meet: ``analysis`` imports it lazily, when a caller asks for ``p_value_method="simulation"``, and the
coalescent kernel itself stays in ``simulation``.
"""

from __future__ import annotations

import random
from typing import Dict, List, Sequence, Tuple

import numpy as np

from metainformant.dna.population import core
from metainformant.simulation.models.popgen import simulate_coalescent_spectra

__all__ = ["monte_carlo_two_sided_p", "simulated_neutrality_tests"]

#: Tests whose weights are symmetric in ``i`` and ``n - i``: they score the folded spectrum and need no outgroup.
FOLDED_TESTS = ("tajima_d", "fu_li_d_star", "fu_li_f_star")
#: Tests that need the outgroup-polarized (unfolded) spectrum.
POLARIZED_TESTS = ("fay_wu_h",)


def monte_carlo_two_sided_p(observed: float, simulated: np.ndarray) -> float:
    """Two-sided Monte Carlo p-value of ``observed`` against ``simulated`` null statistics.

    Each one-sided tail counts the simulated values at least as extreme as the observed one (ties
    included, up to rounding) and gets the +1 correction, so the observed statistic counts as one more
    draw from the null and the p-value is never 0: ``(1 + #{sim <= obs}) / (R + 1)`` below and
    ``(1 + #{sim >= obs}) / (R + 1)`` above. The two-sided p-value is twice the smaller tail, at most 1:
    the equal-tailed version of the normal approximation's ``2 * min(Phi(z), 1 - Phi(z))``, which stays
    valid when the null is not symmetric.

    Args:
        observed: The observed statistic
        simulated: Statistics of the null replicates

    Returns:
        The p-value, at least ``2 / (len(simulated) + 1)`` and at most 1
    """
    tolerance = 1e-9 * max(1.0, abs(observed))
    lower = (1 + np.count_nonzero(simulated <= observed + tolerance)) / (simulated.size + 1)
    upper = (1 + np.count_nonzero(simulated >= observed - tolerance)) / (simulated.size + 1)
    return float(min(1.0, 2.0 * min(lower, upper)))


def _tests_against_null(
    sequences: Sequence[str],
    tests: List[str],
    outgroup: str | None,
    n_simulations: int,
    random_state: int | random.Random | np.random.Generator | None,
) -> Dict[str, Tuple[float, float]]:
    """Score tests that share one observed spectrum against one simulated null."""
    n = len(sequences)
    observed_spectrum = core._site_frequency_spectrum(sequences, outgroup)
    mutations = int(observed_spectrum.sum())
    if mutations == 0:
        return {test: (0.0, 1.0) for test in tests}
    null_spectra = simulate_coalescent_spectra(n, n_simulations, segregating_sites=mutations, rng=random_state)
    results: Dict[str, Tuple[float, float]] = {}
    for test in tests:
        observed = float(core._neutrality_statistics_from_spectra(observed_spectrum[1:], n, test)[0])
        null = core._neutrality_statistics_from_spectra(null_spectra, n, test)
        results[test] = (observed, monte_carlo_two_sided_p(observed, null))
    return results


def simulated_neutrality_tests(
    sequences: Sequence[str],
    tests: Sequence[str],
    *,
    n_simulations: int,
    random_state: int | random.Random | np.random.Generator | None = None,
    outgroup: str | None = None,
) -> Dict[str, Tuple[float, float]]:
    """Observed statistic and Monte Carlo p-value of each requested neutrality test.

    The null is the standard neutral coalescent for the same number of sequences (constant size, infinite
    sites, no recombination), with exactly as many mutations in every replicate as the observed spectrum
    has, dropped on an unweighted genealogy in proportion to branch length
    (:func:`metainformant.simulation.models.popgen.simulate_coalescent_spectra`). The observed spectrum of
    the :data:`FOLDED_TESTS` is folded, as in :mod:`metainformant.dna.population.core`; they are symmetric in
    ``i`` and ``n - i``, so the unfolded simulated spectra give the same statistics. The
    :data:`POLARIZED_TESTS` use the spectrum polarized by ``outgroup``, with the number of mutations that
    spectrum holds. Each group of tests is scored against its own simulation. A sample without mutations
    has statistic 0 and p-value 1, as with the normal approximation.

    Args:
        sequences: Aligned DNA sequences (at least 4 for a meaningful test)
        tests: Names from :data:`FOLDED_TESTS` and :data:`POLARIZED_TESTS`
        n_simulations: Number of coalescent replicates per group of tests
        random_state: Seed of the replicates: an integer, a ``random.Random`` or a
            ``numpy.random.Generator``; ``None`` draws fresh entropy
        outgroup: Outgroup sequence aligned to ``sequences``; required for the polarized tests

    Returns:
        Mapping from test name to ``(observed statistic, two-sided p-value)``

    Raises:
        ValueError: If a test name is unknown or a polarized test is requested without an outgroup
    """
    unknown = [test for test in tests if test not in FOLDED_TESTS + POLARIZED_TESTS]
    if unknown:
        raise ValueError(f"unknown neutrality test(s) {unknown}; choose from {FOLDED_TESTS + POLARIZED_TESTS}")
    polarized = [test for test in tests if test in POLARIZED_TESTS]
    if polarized and outgroup is None:
        raise ValueError(f"{polarized} need an outgroup to polarize the spectrum")

    results: Dict[str, Tuple[float, float]] = {}
    for group, group_outgroup in ((FOLDED_TESTS, None), (POLARIZED_TESTS, outgroup)):
        wanted = [test for test in tests if test in group]
        if wanted:
            results.update(_tests_against_null(sequences, wanted, group_outgroup, n_simulations, random_state))
    return results

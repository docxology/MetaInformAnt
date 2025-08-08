from __future__ import annotations

from metainformant.math import (
    narrow_sense_heritability,
    breeders_equation_response,
    lande_equation_response,
)


def test_narrow_sense_heritability_bounds():
    assert abs(narrow_sense_heritability(0.2, 1.0) - 0.2) < 1e-12
    assert narrow_sense_heritability(1.0, 0.0) == 0.0
    assert narrow_sense_heritability(2.0, 1.0) == 1.0


def test_breeders_equation_response():
    assert abs(breeders_equation_response(1.5, 0.4) - 0.6) < 1e-12


def test_lande_equation_response():
    beta = [0.2, -0.1]
    G = [[1.0, 0.5], [0.5, 2.0]]
    resp = lande_equation_response(beta, G)
    # [1*0.2 + 0.5*(-0.1), 0.5*0.2 + 2*(-0.1)] = [0.15, -0.1]
    assert len(resp) == 2
    assert abs(resp[0] - 0.15) < 1e-12
    assert abs(resp[1] + 0.1) < 1e-12




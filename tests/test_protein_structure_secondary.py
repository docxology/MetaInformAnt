from __future__ import annotations


def test_simple_secondary_structure_prediction_signature():
    from metainformant.protein.structure.general.general.general.general.general.general.general.secondary import simple_helix_coil_propensity

    seq = "MAAAAAGGGGLLLLPPPP"
    probs = simple_helix_coil_propensity(seq)
    assert isinstance(probs, list)
    assert len(probs) == len(seq)
    assert all(0.0 <= p <= 1.0 for p in probs)

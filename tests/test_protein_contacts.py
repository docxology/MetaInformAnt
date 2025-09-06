from __future__ import annotations


def test_contact_pairs_threshold():
    from metainformant.protein.contacts import compute_ca_contact_pairs

    coords = [
        (0.0, 0.0, 0.0),  # 0
        (3.0, 0.0, 0.0),  # 1 (within 4.0)
        (10.0, 0.0, 0.0),  # 2 (far)
    ]
    pairs = compute_ca_contact_pairs(coords, threshold=4.0)
    assert (0, 1) in pairs and all(i < j for i, j in pairs)
    assert (0, 2) not in pairs and (1, 2) not in pairs

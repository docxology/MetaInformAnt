"""Tests for protein contact analysis functions."""

from __future__ import annotations

import numpy as np

from metainformant.protein.structure.contacts import compute_ca_contact_pairs


def test_contact_pairs_threshold():
    """Test CA contact pairs with distance threshold."""
    coords = [
        (0.0, 0.0, 0.0),  # 0
        (3.0, 0.0, 0.0),  # 1 (within 4.0 of 0)
        (10.0, 0.0, 0.0),  # 2 (far from 0 and 1)
    ]
    pairs = compute_ca_contact_pairs(coords, threshold=4.0)
    assert (0, 1) in pairs and all(i < j for i, j in pairs)
    assert (0, 2) not in pairs and (1, 2) not in pairs


def test_contact_pairs_all_close():
    """Test when all atoms are within threshold."""
    coords = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
    ]
    pairs = compute_ca_contact_pairs(coords, threshold=2.0)
    assert len(pairs) == 3  # All pairs within 2.0


def test_contact_pairs_none_close():
    """Test when no atoms are within threshold."""
    coords = [
        (0.0, 0.0, 0.0),
        (100.0, 0.0, 0.0),
        (0.0, 100.0, 0.0),
    ]
    pairs = compute_ca_contact_pairs(coords, threshold=1.0)
    assert len(pairs) == 0


def test_contact_pairs_empty():
    """Test with empty coordinates."""
    pairs = compute_ca_contact_pairs([], threshold=8.0)
    assert pairs == []


def test_contact_pairs_single():
    """Test with single coordinate."""
    pairs = compute_ca_contact_pairs([(0.0, 0.0, 0.0)], threshold=8.0)
    assert pairs == []


def test_residue_contacts_matrix():
    """Test residue contact map calculation."""
    from metainformant.protein.structure.contacts import calculate_residue_contacts

    coords = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],  # residue 0
        [5.0, 0.0, 0.0], [6.0, 0.0, 0.0],  # residue 1
        [50.0, 0.0, 0.0], [51.0, 0.0, 0.0],  # residue 2 (far away)
    ])
    residue_ranges = [(0, 2), (2, 4), (4, 6)]
    contacts = calculate_residue_contacts(coords, residue_ranges, threshold=8.0)

    assert contacts.shape == (3, 3)
    assert contacts[0, 1] == 1  # residue 0 and 1 in contact
    assert contacts[0, 2] == 0  # residue 0 and 2 not in contact


def test_contact_network_analysis():
    """Test contact network analysis."""
    from metainformant.protein.structure.contacts import analyze_contact_network

    contact_map = np.array([
        [0, 1, 1, 0],
        [1, 0, 1, 0],
        [1, 1, 0, 1],
        [0, 0, 1, 0],
    ])
    analysis = analyze_contact_network(contact_map)

    assert analysis["n_residues"] == 4
    assert analysis["n_contacts"] == 4  # (0,1), (0,2), (1,2), (2,3)
    assert analysis["n_components"] == 1


def test_contact_persistence():
    """Test contact persistence across structures."""
    from metainformant.protein.structure.contacts import calculate_contact_persistence

    map1 = np.array([[0, 1], [1, 0]])
    map2 = np.array([[0, 0], [0, 0]])
    persistence = calculate_contact_persistence([map1, map2])

    assert persistence.shape == (2, 2)
    assert persistence[0, 1] == 0.5

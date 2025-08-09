from __future__ import annotations

from typing import Iterable, List, Tuple


def compute_ca_contact_pairs(coords: Iterable[tuple[float, float, float]], *, threshold: float = 8.0) -> List[Tuple[int, int]]:
    """Return index pairs (i,j) with i<j for which Euclidean distance < threshold.

    coords: iterable of (x,y,z). Threshold in Angstroms.
    """
    pts = list(coords)
    pairs: List[Tuple[int, int]] = []
    thr2 = float(threshold) ** 2
    for i in range(len(pts)):
        xi, yi, zi = pts[i]
        for j in range(i + 1, len(pts)):
            xj, yj, zj = pts[j]
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            if dx * dx + dy * dy + dz * dz < thr2:
                pairs.append((i, j))
    return pairs



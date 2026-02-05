#!/usr/bin/env python3
"""
Verification script for phenotype module expansion.
Attempts to import and instantiate classes from all new submodules.
"""
import sys
from datetime import datetime

import numpy as np

try:
    print("Verifying behavior module...")
    from metainformant.phenotype.behavior import BehaviorSequence, Ethogram

    ethogram = Ethogram({"w": "walk", "r": "rest"})
    print("  - Ethogram created")

    # Mock event class for sequence since we didn't export the event class directly in init
    # checking sequence.py to see if we can use a dict or if it needs the dataclass
    # sequence.py expects BehaviorEvent objects.
    from metainformant.phenotype.behavior.sequence import BehaviorEvent

    events = [
        BehaviorEvent(timestamp=0.0, code="w", duration=1.0),
        BehaviorEvent(timestamp=1.0, code="r", duration=2.0),
    ]
    seq = BehaviorSequence(events, ethogram)
    budget = seq.calculate_time_budget()
    print(f"  - Behavior sequence validated. Budget: {budget}")

    print("Verifying chemical module...")
    from metainformant.phenotype.chemical import ChemicalProfile, Compound

    c1 = Compound(name="C25", retention_time=20.5)
    prof = ChemicalProfile("s1", {c1: 0.5})
    norm = prof.normalize()
    print("  - Chemical profile created and normalized")

    print("Verifying sonic module...")
    from metainformant.phenotype.sonic import AcousticSignal

    # Create dummy waveform
    dummy_wave = np.sin(np.linspace(0, 10, 1000))
    sig = AcousticSignal(waveform=dummy_wave, sample_rate=100)
    print(f"  - Acoustic signal created. Duration: {sig.duration}")

    print("Verifying morphological module...")
    from metainformant.phenotype.morphological import Measurement, MorphometricProfile

    m1 = Measurement("HeadWidth", 2.5, "mm")
    m2 = Measurement("HeadLength", 3.0, "mm")
    morph = MorphometricProfile("spec1", [m1, m2])
    idx = morph.calculate_index("CI", "HeadWidth", "HeadLength")
    print(f"  - Morphometric profile created. Index: {idx}")

    print("Verifying electronic module...")
    from metainformant.phenotype.electronic import TrackingPoint, Trajectory

    p1 = TrackingPoint(0, 0, 0)
    p2 = TrackingPoint(1, 1, 1)
    traj = Trajectory("ant1", [p1, p2])
    dist = traj.total_distance()
    print(f"  - Trajectory created. Distance: {dist}")

    print("\nSUCCESS: All modules verified.")

except Exception as e:
    print(f"\nFAILURE: {e}")
    import traceback

    traceback.print_exc()
    sys.exit(1)

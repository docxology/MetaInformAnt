from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from scipy import stats

from metainformant.core.utils.errors import ValidationError

from .ethogram import Ethogram


@dataclass
class BehaviorEvent:
    """A single occurrence of a behavior."""

    timestamp: float
    code: str
    duration: float = 0.0
    modifiers: Dict[str, Any] = field(default_factory=dict)


class BehaviorSequence:
    """A sequence of behavioral events with analysis capabilities.

    Provides time budgets, transition matrices, Shannon diversity,
    bout analysis, and first-order Markov chain stationarity tests.
    """

    def __init__(self, events: List[BehaviorEvent], ethogram: Ethogram):
        self.events = sorted(events, key=lambda x: x.timestamp)
        self.ethogram = ethogram
        self._validate_events()

    def _validate_events(self) -> None:
        for event in self.events:
            if not self.ethogram.validate(event.code):
                raise ValidationError(f"Behavior code '{event.code}' not found in ethogram.")

    def calculate_time_budget(self) -> Dict[str, float]:
        """Proportion of total time spent in each behavior (duration-weighted)."""
        total_duration = sum(e.duration for e in self.events)
        if total_duration == 0:
            return {}

        budget: Dict[str, float] = defaultdict(float)
        for event in self.events:
            budget[event.code] += event.duration

        return {code: dur / total_duration for code, dur in budget.items()}

    def transition_matrix(self) -> Dict[str, Dict[str, float]]:
        """First-order transition probability matrix between behaviors."""
        counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

        for i in range(len(self.events) - 1):
            counts[self.events[i].code][self.events[i + 1].code] += 1

        matrix: Dict[str, Dict[str, float]] = {}
        for start, targets in counts.items():
            total = sum(targets.values())
            matrix[start] = {t: c / total for t, c in targets.items()}
        return matrix

    def shannon_diversity(self) -> float:
        """Shannon diversity index (H') of behavior frequencies.

        Higher values indicate more even distribution of behaviors.
        Returns 0.0 for empty sequences.
        """
        if not self.events:
            return 0.0

        counts: Dict[str, int] = defaultdict(int)
        for e in self.events:
            counts[e.code] += 1

        total = len(self.events)
        h = 0.0
        for count in counts.values():
            p = count / total
            if p > 0:
                h -= p * math.log(p)
        return h

    def simpson_diversity(self) -> float:
        """Simpson's diversity index (1 - D).

        Probability that two randomly chosen events are different behaviors.
        Range [0, 1], higher = more diverse.
        """
        if len(self.events) < 2:
            return 0.0

        counts: Dict[str, int] = defaultdict(int)
        for e in self.events:
            counts[e.code] += 1

        n = len(self.events)
        d = sum(c * (c - 1) for c in counts.values()) / (n * (n - 1))
        return 1.0 - d

    def bout_analysis(self, min_bout_duration: float = 0.0) -> Dict[str, Dict[str, Any]]:
        """Analyze behavioral bouts (consecutive runs of the same behavior).

        A bout is a contiguous sequence of events with the same behavior code.

        Args:
            min_bout_duration: Minimum total duration to count as a bout.

        Returns:
            Per-behavior dict with bout_count, mean_duration, total_duration, max_duration.
        """
        if not self.events:
            return {}

        bouts: Dict[str, List[float]] = defaultdict(list)
        current_code = self.events[0].code
        current_dur = self.events[0].duration

        for event in self.events[1:]:
            if event.code == current_code:
                current_dur += event.duration
            else:
                if current_dur >= min_bout_duration:
                    bouts[current_code].append(current_dur)
                current_code = event.code
                current_dur = event.duration

        # Don't forget the last bout
        if current_dur >= min_bout_duration:
            bouts[current_code].append(current_dur)

        results: Dict[str, Dict[str, Any]] = {}
        for code, durations in bouts.items():
            results[code] = {
                "bout_count": len(durations),
                "mean_duration": sum(durations) / len(durations) if durations else 0.0,
                "total_duration": sum(durations),
                "max_duration": max(durations) if durations else 0.0,
            }
        return results

    def event_rate(self, time_window: Optional[float] = None) -> float:
        """Events per unit time.

        Args:
            time_window: Custom time window. If None, uses sequence span.

        Returns:
            The number of events divided by the window (or span). A sequence whose events all share one timestamp
            has a span of zero, over which a rate is not defined; the number of events is returned then (the
            rate over one unit of time). A non-positive ``time_window`` returns 0.0.
        """
        if not self.events:
            return 0.0

        if time_window is not None:
            return len(self.events) / time_window if time_window > 0 else 0.0

        span = self.events[-1].timestamp - self.events[0].timestamp
        if span <= 0:
            return float(len(self.events))
        return len(self.events) / span

    def markov_stationarity_chi2(self) -> Dict[str, Any]:
        """Chi-squared test of first-order Markov stationarity (homogeneity of the transitions over time).

        Splits the sequence in half and compares the transition counts of the two halves state by state: for each
        behavior observed in both halves, the 2 x K table (half by next behavior) is tested for homogeneity with
        Pearson's chi-squared, its expected counts coming from the table's own margins
        (``row total * column total / state total``), so halves that contribute different numbers of transitions
        from a state are not mistaken for a change in the transition probabilities. The statistics add up over
        states, with ``sum(K - 1)`` degrees of freedom. The chi-squared reference distribution is asymptotic and
        unreliable when expected counts are small.

        Returns:
            Dict with ``chi2``, ``df``, ``p_value`` and ``sufficient_data``. With fewer than six events, or no
            behavior followed by a transition in both halves, there is nothing to compare:
            ``sufficient_data`` is False, ``chi2`` is 0.0, ``df`` is 0 and ``p_value`` is ``nan``.
        """
        if len(self.events) < 6:
            return {"chi2": 0.0, "df": 0, "p_value": math.nan, "sufficient_data": False}

        mid = len(self.events) // 2
        first_half = self.events[:mid]
        second_half = self.events[mid:]

        def _counts(evts: List[BehaviorEvent]) -> Dict[str, Dict[str, int]]:
            c: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
            for i in range(len(evts) - 1):
                c[evts[i].code][evts[i + 1].code] += 1
            return c

        c1 = _counts(first_half)
        c2 = _counts(second_half)

        chi2 = 0.0
        df = 0
        for state in set(c1.keys()) & set(c2.keys()):
            row1, row2 = c1[state], c2[state]
            n1, n2 = sum(row1.values()), sum(row2.values())
            targets = set(row1) | set(row2)
            for target in targets:
                o1, o2 = row1.get(target, 0), row2.get(target, 0)
                column = o1 + o2
                expected1 = n1 * column / (n1 + n2)
                expected2 = n2 * column / (n1 + n2)
                chi2 += (o1 - expected1) ** 2 / expected1 + (o2 - expected2) ** 2 / expected2
            df += len(targets) - 1

        if df == 0:
            return {"chi2": 0.0, "df": 0, "p_value": math.nan, "sufficient_data": False}
        return {"chi2": chi2, "df": df, "p_value": float(stats.chi2.sf(chi2, df)), "sufficient_data": True}

    def latency_to_first(self, behavior_code: str) -> Optional[float]:
        """Time from sequence start to first occurrence of a behavior."""
        if not self.events:
            return None
        start_time = self.events[0].timestamp
        for event in self.events:
            if event.code == behavior_code:
                return event.timestamp - start_time
        return None

    def behavior_counts(self) -> Dict[str, int]:
        """Frequency count of each behavior code."""
        counts: Dict[str, int] = defaultdict(int)
        for e in self.events:
            counts[e.code] += 1
        return dict(counts)

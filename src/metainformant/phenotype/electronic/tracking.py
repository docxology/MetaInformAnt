"""Electronic tracking and movement analysis for RFID, video, and GPS data."""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple


@dataclass(frozen=True)
class TrackingPoint:
    """A single point in a spatiotemporal trajectory."""

    x: float
    y: float
    timestamp: float
    z: float = 0.0
    confidence: float = 1.0

    def distance_to(self, other: TrackingPoint) -> float:
        """Euclidean distance to another point (3D)."""
        return math.sqrt((self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2)


class Trajectory:
    """A sequence of tracking points representing movement.

    Provides path analysis, velocity profiles, home range estimation,
    activity classification, sinuosity metrics, and interaction detection.
    """

    def __init__(self, entity_id: str, points: List[TrackingPoint]) -> None:
        self.entity_id = entity_id
        self.points = sorted(points, key=lambda p: p.timestamp)

    @property
    def duration(self) -> float:
        """Total time span in seconds."""
        if len(self.points) < 2:
            return 0.0
        return self.points[-1].timestamp - self.points[0].timestamp

    @property
    def n_points(self) -> int:
        return len(self.points)

    def total_distance(self) -> float:
        """Total path length (sum of step distances)."""
        if len(self.points) < 2:
            return 0.0
        return sum(self.points[i].distance_to(self.points[i + 1]) for i in range(len(self.points) - 1))

    def net_displacement(self) -> float:
        """Straight-line distance from first to last point."""
        if len(self.points) < 2:
            return 0.0
        return self.points[0].distance_to(self.points[-1])

    def sinuosity(self) -> float:
        """Sinuosity index: total_distance / net_displacement.

        Value of 1.0 means perfectly straight; higher means more tortuous.
        Returns inf if net displacement is zero (circular path).
        """
        net = self.net_displacement()
        if net == 0:
            total = self.total_distance()
            return float("inf") if total > 0 else 1.0
        return self.total_distance() / net

    def velocity_profile(self) -> List[float]:
        """Instantaneous velocity between consecutive points (units/sec)."""
        velocities = []
        for i in range(len(self.points) - 1):
            p1, p2 = self.points[i], self.points[i + 1]
            dt = p2.timestamp - p1.timestamp
            dd = p1.distance_to(p2)
            velocities.append(dd / dt if dt > 0 else 0.0)
        return velocities

    def mean_velocity(self) -> float:
        """Mean velocity over the trajectory."""
        dur = self.duration
        if dur <= 0:
            return 0.0
        return self.total_distance() / dur

    def max_velocity(self) -> float:
        """Maximum instantaneous velocity."""
        vels = self.velocity_profile()
        return max(vels) if vels else 0.0

    def acceleration_profile(self) -> List[float]:
        """Acceleration between consecutive velocity measurements (units/sec^2)."""
        vels = self.velocity_profile()
        if len(vels) < 2:
            return []

        accels = []
        for i in range(len(vels) - 1):
            dt = self.points[i + 2].timestamp - self.points[i + 1].timestamp
            if dt > 0:
                accels.append((vels[i + 1] - vels[i]) / dt)
            else:
                accels.append(0.0)
        return accels

    def turning_angles(self) -> List[float]:
        """Turning angles between consecutive segments (radians, [0, pi])."""
        if len(self.points) < 3:
            return []

        angles = []
        for i in range(len(self.points) - 2):
            p0, p1, p2 = self.points[i], self.points[i + 1], self.points[i + 2]
            v1x, v1y = p1.x - p0.x, p1.y - p0.y
            v2x, v2y = p2.x - p1.x, p2.y - p1.y

            dot = v1x * v2x + v1y * v2y
            mag1 = math.sqrt(v1x**2 + v1y**2)
            mag2 = math.sqrt(v2x**2 + v2y**2)

            if mag1 == 0 or mag2 == 0:
                angles.append(0.0)
            else:
                cos_angle = max(-1.0, min(1.0, dot / (mag1 * mag2)))
                angles.append(math.acos(cos_angle))
        return angles

    def mean_turning_angle(self) -> float:
        """Mean turning angle in radians."""
        angles = self.turning_angles()
        return statistics.mean(angles) if angles else 0.0

    def home_range_mcp(self) -> float:
        """Minimum convex polygon (MCP) home range area estimate (2D).

        Uses the convex hull of all x,y points. Returns area in square units.
        """
        if len(self.points) < 3:
            return 0.0

        pts = [(p.x, p.y) for p in self.points]
        hull = _convex_hull(pts)
        return _polygon_area(hull)

    def bounding_box(self) -> Dict[str, float]:
        """Axis-aligned bounding box of the trajectory."""
        if not self.points:
            return {"x_min": 0, "x_max": 0, "y_min": 0, "y_max": 0}
        xs = [p.x for p in self.points]
        ys = [p.y for p in self.points]
        return {
            "x_min": min(xs),
            "x_max": max(xs),
            "y_min": min(ys),
            "y_max": max(ys),
        }

    def activity_states(self, velocity_threshold: float = 1.0) -> List[str]:
        """Classify each step as 'active' or 'resting' based on velocity threshold."""
        vels = self.velocity_profile()
        return ["active" if v >= velocity_threshold else "resting" for v in vels]

    def activity_budget(self, velocity_threshold: float = 1.0) -> Dict[str, float]:
        """Fraction of time spent active vs resting."""
        states = self.activity_states(velocity_threshold)
        if not states:
            return {"active": 0.0, "resting": 0.0}
        n = len(states)
        active_count = states.count("active")
        return {
            "active": active_count / n,
            "resting": (n - active_count) / n,
        }

    def segment_by_time(self, interval: float) -> List[Trajectory]:
        """Split trajectory into fixed-time segments."""
        if not self.points or interval <= 0:
            return []

        segments = []
        start_time = self.points[0].timestamp
        current_points: List[TrackingPoint] = []

        for pt in self.points:
            if pt.timestamp - start_time >= interval and current_points:
                segments.append(Trajectory(self.entity_id, current_points))
                start_time = pt.timestamp
                current_points = []
            current_points.append(pt)

        if current_points:
            segments.append(Trajectory(self.entity_id, current_points))

        return segments


def detect_interactions(
    trajectories: List[Trajectory],
    distance_threshold: float,
    time_threshold: float,
) -> List[Dict[str, Any]]:
    """Detect spatial interactions between trajectories.

    Two entities interact when they are within distance_threshold
    at approximately the same time (within time_threshold).

    Args:
        trajectories: List of Trajectory objects.
        distance_threshold: Max distance for interaction.
        time_threshold: Max temporal offset for co-occurrence.

    Returns:
        List of interaction events with entity ids, time, and distance.
    """
    interactions = []
    for i in range(len(trajectories)):
        for j in range(i + 1, len(trajectories)):
            t_a, t_b = trajectories[i], trajectories[j]
            for pa in t_a.points:
                for pb in t_b.points:
                    if abs(pa.timestamp - pb.timestamp) <= time_threshold:
                        d = pa.distance_to(pb)
                        if d <= distance_threshold:
                            interactions.append(
                                {
                                    "entity_a": t_a.entity_id,
                                    "entity_b": t_b.entity_id,
                                    "timestamp": (pa.timestamp + pb.timestamp) / 2,
                                    "distance": d,
                                }
                            )
    return interactions


def _convex_hull(points: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """Andrew's monotone chain convex hull algorithm. O(n log n)."""
    points = sorted(set(points))
    if len(points) <= 1:
        return list(points)

    def cross(o: Tuple[float, float], a: Tuple[float, float], b: Tuple[float, float]) -> float:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: List[Tuple[float, float]] = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper: List[Tuple[float, float]] = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    return lower[:-1] + upper[:-1]


def _polygon_area(vertices: List[Tuple[float, float]]) -> float:
    """Shoelace formula for polygon area."""
    n = len(vertices)
    if n < 3:
        return 0.0
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += vertices[i][0] * vertices[j][1]
        area -= vertices[j][0] * vertices[i][1]
    return abs(area) / 2.0

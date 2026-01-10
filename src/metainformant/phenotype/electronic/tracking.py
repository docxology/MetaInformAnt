from dataclasses import dataclass
from typing import List, Optional
import math

@dataclass(frozen=True)
class TrackingPoint:
    """
    A single point in a spatiotemporal trajectory.
    """
    x: float
    y: float
    timestamp: float
    z: float = 0.0
    confidence: float = 1.0
    
    def distance_to(self, other: 'TrackingPoint') -> float:
        """Euclidean distance to another point (3D)."""
        return math.sqrt(
            (self.x - other.x)**2 + 
            (self.y - other.y)**2 + 
            (self.z - other.z)**2
        )

class Trajectory:
    """
    A sequence of tracking points representing movement.
    """
    
    def __init__(self, entity_id: str, points: List[TrackingPoint]):
        self.entity_id = entity_id
        # Sort by timestamp
        self.points = sorted(points, key=lambda p: p.timestamp)
        
    def total_distance(self) -> float:
        """Calculate total path length."""
        if len(self.points) < 2:
            return 0.0
            
        dist = 0.0
        for i in range(len(self.points) - 1):
            dist += self.points[i].distance_to(self.points[i+1])
            
        return dist
        
    def velocity_profile(self) -> List[float]:
        """
        Calculate velocity between consecutive points (units/sec).
        """
        velocities = []
        for i in range(len(self.points) - 1):
            p1 = self.points[i]
            p2 = self.points[i+1]
            
            dt = p2.timestamp - p1.timestamp
            dd = p1.distance_to(p2)
            
            if dt > 0:
                velocities.append(dd / dt)
            else:
                velocities.append(0.0) # Instant teleporation?
                
        return velocities

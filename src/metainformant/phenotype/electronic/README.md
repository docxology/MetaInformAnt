# Electronic Phenotype Module

The `electronic` module handles data from electronic sensors, such as RFID tags, GPS tracking, and automated monitoring systems.

## Overview
This module handles:
- **Tracking Data**: Spatiotemporal points (x, y, z, t).
- **Sensor Readings**: Environmental data associated with phenotypes.
- **RFID**: Interaction events from tagged individuals.

## Components
- `TrackingPoint`: Class for a single spatiotemporal datum.
- `Trajectory`: Sequence of tracking points.

## Usage
```python
from metainformant.phenotype.electronic import TrackingPoint, Trajectory

# Define a point
p1 = TrackingPoint(x=10.5, y=20.0, timestamp=1620000000)
p2 = TrackingPoint(x=10.6, y=20.1, timestamp=1620000001)

# Create trajectory
traj = Trajectory(id="ant_001", points=[p1, p2])
print(f"Total distance: {traj.distance()} units")
```

# Electronic Module Technical Specification

## Architecture
The `electronic` module provides tools for handling high-frequency sensor data.

## Components

### `TrackingPoint`
- **Attributes**:
    - `x`, `y`, `z`: Coordinates (float).
    - `timestamp`: Unix timestamp or datetime.
    - `confidence`: Optional detection confidence.

### `Trajectory`
- **Attributes**:
    - `id`: Unique identifier.
    - `points`: List of `TrackingPoint`.
- **Methods**:
    - `distance()`: Total path length.
    - `velocity_profile()`: Speed over time.
    - `resample(rate)`: Interpolate to fixed frame rate.

## Integration
- Can be converted to `BehaviorSequence` (in `metainformant.phenotype.behavior`) via spatial discretization.

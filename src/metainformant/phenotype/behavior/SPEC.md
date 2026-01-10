# Behavior Module Technical Specification

## Architecture
The `behavior` module is designed to capture and analyze discrete and continuous behavioral data.

## Components

### `Ethogram`
- **Purpose**: Defines a controlled vocabulary of behaviors.
- **Attributes**:
    - `behaviors`: Dict connecting behavior codes to descriptions/metadata.
- **Methods**:
    - `validate(code)`: Checks if a code exists in the ethogram.
    - `load(path)`: Loads ethogram from JSON/YAML.

### `BehaviorSequence`
- **Purpose**: Represents a sequence of behavioral events.
- **Attributes**:
    - `events`: List of (timestamp, behavior_code, duration, modifiers).
    - `ethogram`: Reference `Ethogram`.
- **Methods**:
    - `time_budget()`: Calculates proportion of time in each state.
    - `transition_matrix()`: Calculates transition probabilities between states.

## Data Structures
- **Events**: Should support discrete events (point in time) and states (duration).

# Sonic Module Technical Specification

## Architecture
The `sonic` module encapsulates acoustic data handling, focusing on the properties relevant to insect communication (stridulation, drumming).

## Components

### `AcousticSignal`
- **Attributes**:
    - `waveform`: Numpy array of samples.
    - `sample_rate`: Sampling frequency (Hz).
    - `metadata`: Recording conditions.
- **Properties**:
    - `duration`: Length in seconds.
    - `dominant_frequency`: Frequency with highest energy.
    - `pulse_count`: Number of distinct pulses (requires thresholding).

### `SignalAnalyst`
- **Methods**:
    - `bandpass_filter(signal, low, high)`: Applies Butterworth filter.
    - `detect_pulses(signal, threshold)`: Identifies pulse start/stop times.

## Dependencies
- `numpy`, `scipy.signal` (for processing).
- `soundfile` or `librosa` (for I/O - optional).

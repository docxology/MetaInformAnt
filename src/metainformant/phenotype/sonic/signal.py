import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Any, Optional


@dataclass
class AcousticSignal:
    """
    Represents an acoustic signal (e.g., stridulation).
    """

    waveform: np.ndarray
    sample_rate: int
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def duration(self) -> float:
        """Duration in seconds."""
        return len(self.waveform) / self.sample_rate

    @property
    def time_axis(self) -> np.ndarray:
        """Time points for the waveform."""
        return np.arange(len(self.waveform)) / self.sample_rate

    @classmethod
    def from_file(cls, path: str) -> "AcousticSignal":
        """
        Placeholder for loading from file.
        Requires soundfile or librosa which might not be installed in all envs.
        For now, raising NotImplementedError or importing dynamically.
        """
        raise NotImplementedError("File loading requires 'soundfile' or 'librosa' dependency.")

    def dominant_frequency(self) -> float:
        """
        Calculate dominant frequency using FFT.
        """
        if len(self.waveform) == 0:
            return 0.0

        # Simple FFT
        spectrum = np.fft.rfft(self.waveform)
        frequencies = np.fft.rfftfreq(len(self.waveform), d=1 / self.sample_rate)

        peak_idx = np.argmax(np.abs(spectrum))
        return frequencies[peak_idx]

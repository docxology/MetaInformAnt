"""Acoustic signal analysis for stridulation, vibrational communication, and bioacoustics."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


@dataclass
class AcousticSignal:
    """Represents an acoustic signal with spectral and temporal analysis.

    Provides FFT-based frequency analysis, spectrogram computation,
    RMS energy, syllable detection, and temporal pattern extraction.
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

    @property
    def n_samples(self) -> int:
        return len(self.waveform)

    @classmethod
    def from_file(cls, path: str) -> AcousticSignal:
        """Load from audio file. Requires 'soundfile' package."""
        try:
            import soundfile as sf
        except ImportError:
            raise ImportError("File loading requires the 'soundfile' package: uv pip install soundfile")

        data, sr = sf.read(path, dtype="float64")
        if data.ndim > 1:
            data = data.mean(axis=1)  # Mix to mono
        return cls(waveform=data, sample_rate=sr)

    @classmethod
    def generate_tone(cls, frequency: float, duration: float, sample_rate: int = 44100) -> AcousticSignal:
        """Generate a pure sine tone for testing."""
        t = np.arange(int(sample_rate * duration)) / sample_rate
        waveform = np.sin(2 * np.pi * frequency * t)
        return cls(waveform=waveform, sample_rate=sample_rate, metadata={"type": "tone", "frequency": frequency})

    def dominant_frequency(self) -> float:
        """Dominant frequency using FFT (Hz)."""
        if len(self.waveform) == 0:
            return 0.0

        spectrum = np.fft.rfft(self.waveform)
        frequencies = np.fft.rfftfreq(len(self.waveform), d=1 / self.sample_rate)
        peak_idx = np.argmax(np.abs(spectrum))
        return float(frequencies[peak_idx])

    def power_spectrum(self) -> Tuple[np.ndarray, np.ndarray]:
        """Compute power spectrum density.

        Returns:
            (frequencies, power) arrays.
        """
        spectrum = np.fft.rfft(self.waveform)
        frequencies = np.fft.rfftfreq(len(self.waveform), d=1 / self.sample_rate)
        power = np.abs(spectrum) ** 2 / len(self.waveform)
        return frequencies, power

    def spectral_centroid(self) -> float:
        """Spectral centroid (weighted mean frequency, Hz)."""
        frequencies, power = self.power_spectrum()
        total_power = np.sum(power)
        if total_power == 0:
            return 0.0
        return float(np.sum(frequencies * power) / total_power)

    def spectral_bandwidth(self) -> float:
        """Spectral bandwidth (weighted std of frequency, Hz)."""
        frequencies, power = self.power_spectrum()
        total_power = np.sum(power)
        if total_power == 0:
            return 0.0
        centroid = np.sum(frequencies * power) / total_power
        variance = np.sum(power * (frequencies - centroid) ** 2) / total_power
        return float(np.sqrt(variance))

    def band_energy(self, low_hz: float, high_hz: float) -> float:
        """Total energy in a frequency band (fraction of total)."""
        frequencies, power = self.power_spectrum()
        total = np.sum(power)
        if total == 0:
            return 0.0
        mask = (frequencies >= low_hz) & (frequencies <= high_hz)
        return float(np.sum(power[mask]) / total)

    def rms_energy(self) -> float:
        """Root mean square energy of the waveform."""
        if len(self.waveform) == 0:
            return 0.0
        return float(np.sqrt(np.mean(self.waveform**2)))

    def peak_amplitude(self) -> float:
        """Maximum absolute amplitude."""
        if len(self.waveform) == 0:
            return 0.0
        return float(np.max(np.abs(self.waveform)))

    def zero_crossing_rate(self) -> float:
        """Zero crossing rate (crossings per second)."""
        if len(self.waveform) < 2:
            return 0.0
        signs = np.sign(self.waveform)
        crossings = np.sum(np.abs(np.diff(signs)) > 0)
        return float(crossings / self.duration) if self.duration > 0 else 0.0

    def spectrogram(self, window_size: int = 1024, hop_size: int = 512) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Short-time Fourier transform spectrogram.

        Args:
            window_size: FFT window size in samples.
            hop_size: Hop between windows in samples.

        Returns:
            (times, frequencies, magnitude) where magnitude is shape (n_freq, n_frames).
        """
        n_frames = 1 + (len(self.waveform) - window_size) // hop_size
        if n_frames <= 0:
            return np.array([]), np.array([]), np.array([[]])

        n_freq = window_size // 2 + 1
        magnitude = np.zeros((n_freq, n_frames))
        window = np.hanning(window_size)

        for i in range(n_frames):
            start = i * hop_size
            frame = self.waveform[start : start + window_size] * window
            spectrum = np.fft.rfft(frame)
            magnitude[:, i] = np.abs(spectrum)

        frequencies = np.fft.rfftfreq(window_size, d=1 / self.sample_rate)
        times = np.arange(n_frames) * hop_size / self.sample_rate

        return times, frequencies, magnitude

    def detect_syllables(
        self, threshold_factor: float = 2.0, min_gap_ms: float = 10.0, window_ms: float = 5.0
    ) -> List[Dict[str, float]]:
        """Detect syllables based on amplitude envelope thresholding.

        Args:
            threshold_factor: Multiplier over mean RMS for onset detection.
            min_gap_ms: Minimum gap between syllables (ms).
            window_ms: Envelope smoothing window (ms).

        Returns:
            List of dicts with onset, offset, duration (in seconds).
        """
        if len(self.waveform) == 0:
            return []

        # Compute smoothed amplitude envelope
        window_samples = max(1, int(self.sample_rate * window_ms / 1000))
        abs_signal = np.abs(self.waveform)
        if len(abs_signal) < window_samples:
            envelope = abs_signal
        else:
            kernel = np.ones(window_samples) / window_samples
            envelope = np.convolve(abs_signal, kernel, mode="same")

        threshold = np.mean(envelope) * threshold_factor
        min_gap_samples = int(self.sample_rate * min_gap_ms / 1000)

        syllables = []
        in_syllable = False
        onset_idx = 0

        for i, val in enumerate(envelope):
            if not in_syllable and val >= threshold:
                in_syllable = True
                onset_idx = i
            elif in_syllable and val < threshold:
                in_syllable = False
                offset_idx = i
                # Check minimum gap before starting new syllable
                if syllables:
                    prev_offset = int(syllables[-1]["offset"] * self.sample_rate)
                    if onset_idx - prev_offset < min_gap_samples:
                        # Merge with previous
                        syllables[-1]["offset"] = offset_idx / self.sample_rate
                        syllables[-1]["duration"] = syllables[-1]["offset"] - syllables[-1]["onset"]
                        continue

                syllables.append(
                    {
                        "onset": onset_idx / self.sample_rate,
                        "offset": offset_idx / self.sample_rate,
                        "duration": (offset_idx - onset_idx) / self.sample_rate,
                    }
                )

        # Handle case where signal ends during a syllable
        if in_syllable:
            syllables.append(
                {
                    "onset": onset_idx / self.sample_rate,
                    "offset": len(self.waveform) / self.sample_rate,
                    "duration": (len(self.waveform) - onset_idx) / self.sample_rate,
                }
            )

        return syllables

    def temporal_pattern(self, threshold_factor: float = 2.0) -> Dict[str, Any]:
        """Extract temporal pattern from syllable detection.

        Returns:
            Dict with syllable_count, mean_duration, mean_inter_syllable_interval,
            duty_cycle, and regularity (CV of intervals).
        """
        syllables = self.detect_syllables(threshold_factor=threshold_factor)

        if not syllables:
            return {
                "syllable_count": 0,
                "mean_duration": 0.0,
                "mean_isi": 0.0,
                "duty_cycle": 0.0,
                "regularity": 0.0,
            }

        durations = [s["duration"] for s in syllables]
        mean_dur = sum(durations) / len(durations)

        # Inter-syllable intervals
        isis = []
        for i in range(len(syllables) - 1):
            isi = syllables[i + 1]["onset"] - syllables[i]["offset"]
            isis.append(isi)

        mean_isi = sum(isis) / len(isis) if isis else 0.0

        total_active = sum(durations)
        duty_cycle = total_active / self.duration if self.duration > 0 else 0.0

        # Regularity = 1 - CV of ISIs (higher = more regular)
        if len(isis) > 1 and mean_isi > 0:
            import statistics as stats

            cv = stats.stdev(isis) / mean_isi
            regularity = max(0.0, 1.0 - cv)
        else:
            regularity = 1.0 if isis else 0.0

        return {
            "syllable_count": len(syllables),
            "mean_duration": mean_dur,
            "mean_isi": mean_isi,
            "duty_cycle": duty_cycle,
            "regularity": regularity,
        }

    def trim(self, start_sec: float = 0.0, end_sec: Optional[float] = None) -> AcousticSignal:
        """Return a trimmed copy of the signal."""
        start_idx = int(start_sec * self.sample_rate)
        end_idx = int(end_sec * self.sample_rate) if end_sec is not None else len(self.waveform)
        start_idx = max(0, min(start_idx, len(self.waveform)))
        end_idx = max(start_idx, min(end_idx, len(self.waveform)))
        return AcousticSignal(
            waveform=self.waveform[start_idx:end_idx].copy(),
            sample_rate=self.sample_rate,
            metadata=self.metadata,
        )

    def normalize_amplitude(self) -> AcousticSignal:
        """Return a copy with amplitude normalized to [-1, 1]."""
        peak = self.peak_amplitude()
        if peak == 0:
            return AcousticSignal(self.waveform.copy(), self.sample_rate, self.metadata)
        return AcousticSignal(self.waveform / peak, self.sample_rate, self.metadata)

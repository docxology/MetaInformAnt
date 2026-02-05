"""Comprehensive tests for the phenotype module.

Tests all submodules: behavior, chemical, electronic, morphological, sonic,
analysis, workflow, and integration. No mocking — all real computations.
"""

from __future__ import annotations

import math
import json
from pathlib import Path
from typing import Dict

import numpy as np
import pytest

# ============================================================
# BEHAVIOR SUBMODULE
# ============================================================


class TestEthogram:
    def test_create_from_strings(self):
        from metainformant.phenotype.behavior.ethogram import Ethogram

        eth = Ethogram({"F": "Foraging", "G": "Grooming", "R": "Resting"})
        assert len(eth) == 3
        assert eth.validate("F")
        assert not eth.validate("X")

    def test_create_from_definitions(self):
        from metainformant.phenotype.behavior.ethogram import Ethogram, BehaviorDefinition

        bd = BehaviorDefinition(code="F", name="Forage", description="Searching for food")
        eth = Ethogram({"F": bd})
        assert eth.get("F").name == "Forage"

    def test_invalid_definition_raises(self):
        from metainformant.phenotype.behavior.ethogram import Ethogram
        from metainformant.core.utils.errors import ValidationError

        with pytest.raises(ValidationError):
            Ethogram({"X": 123})


class TestBehaviorSequence:
    @pytest.fixture
    def ethogram(self):
        from metainformant.phenotype.behavior.ethogram import Ethogram

        return Ethogram({"F": "Foraging", "G": "Grooming", "R": "Resting", "T": "Trophallaxis"})

    @pytest.fixture
    def sequence(self, ethogram):
        from metainformant.phenotype.behavior.sequence import BehaviorSequence, BehaviorEvent

        events = [
            BehaviorEvent(timestamp=0, code="F", duration=10.0),
            BehaviorEvent(timestamp=10, code="F", duration=5.0),
            BehaviorEvent(timestamp=15, code="G", duration=8.0),
            BehaviorEvent(timestamp=23, code="R", duration=12.0),
            BehaviorEvent(timestamp=35, code="F", duration=7.0),
            BehaviorEvent(timestamp=42, code="T", duration=3.0),
            BehaviorEvent(timestamp=45, code="G", duration=5.0),
            BehaviorEvent(timestamp=50, code="R", duration=10.0),
        ]
        return BehaviorSequence(events, ethogram)

    def test_time_budget(self, sequence):
        budget = sequence.calculate_time_budget()
        assert len(budget) == 4
        total = sum(budget.values())
        assert abs(total - 1.0) < 1e-9

    def test_transition_matrix(self, sequence):
        matrix = sequence.transition_matrix()
        assert "F" in matrix
        for probs in matrix.values():
            total = sum(probs.values())
            assert abs(total - 1.0) < 1e-9

    def test_shannon_diversity(self, sequence):
        h = sequence.shannon_diversity()
        assert h > 0
        # Max entropy for 4 categories = ln(4) ≈ 1.386
        assert h <= math.log(4) + 0.01

    def test_simpson_diversity(self, sequence):
        d = sequence.simpson_diversity()
        assert 0 <= d <= 1.0

    def test_bout_analysis(self, sequence):
        bouts = sequence.bout_analysis()
        assert "F" in bouts
        assert bouts["F"]["bout_count"] >= 1
        assert bouts["F"]["total_duration"] > 0

    def test_event_rate(self, sequence):
        rate = sequence.event_rate()
        assert rate > 0

    def test_markov_stationarity(self, sequence):
        result = sequence.markov_stationarity_chi2()
        assert "chi2" in result
        assert "df" in result
        assert result["sufficient_data"] is True

    def test_latency_to_first(self, sequence):
        lat = sequence.latency_to_first("G")
        assert lat == 15.0
        assert sequence.latency_to_first("F") == 0.0
        assert sequence.latency_to_first("X") is None

    def test_behavior_counts(self, sequence):
        counts = sequence.behavior_counts()
        assert counts["F"] == 3
        assert counts["G"] == 2

    def test_invalid_code_raises(self, ethogram):
        from metainformant.phenotype.behavior.sequence import BehaviorSequence, BehaviorEvent
        from metainformant.core.utils.errors import ValidationError

        events = [BehaviorEvent(timestamp=0, code="INVALID", duration=1.0)]
        with pytest.raises(ValidationError):
            BehaviorSequence(events, ethogram)


# ============================================================
# CHEMICAL SUBMODULE
# ============================================================


class TestChemicalProfile:
    @pytest.fixture
    def compounds(self):
        from metainformant.phenotype.chemical.compound import Compound

        return {
            Compound("n-C25"): 30.0,
            Compound("n-C27"): 45.0,
            Compound("n-C29"): 20.0,
            Compound("3-MeC25"): 5.0,
        }

    @pytest.fixture
    def profile(self, compounds):
        from metainformant.phenotype.chemical.profile import ChemicalProfile

        return ChemicalProfile("sample_1", compounds)

    @pytest.fixture
    def profile_b(self):
        from metainformant.phenotype.chemical.compound import Compound
        from metainformant.phenotype.chemical.profile import ChemicalProfile

        compounds = {
            Compound("n-C25"): 10.0,
            Compound("n-C27"): 60.0,
            Compound("n-C29"): 25.0,
            Compound("n-C31"): 5.0,
        }
        return ChemicalProfile("sample_2", compounds)

    def test_normalize_tic(self, profile):
        normed = profile.normalize("total_ion_current")
        assert abs(sum(normed.compounds.values()) - 1.0) < 1e-9

    def test_normalize_max_peak(self, profile):
        normed = profile.normalize("max_peak")
        assert max(normed.compounds.values()) == 1.0

    def test_bray_curtis(self, profile, profile_b):
        d = profile.bray_curtis_distance(profile_b)
        assert 0 <= d <= 1.0

    def test_euclidean_distance(self, profile, profile_b):
        d = profile.euclidean_distance(profile_b)
        assert d > 0

    def test_cosine_similarity(self, profile, profile_b):
        s = profile.cosine_similarity(profile_b)
        assert -1 <= s <= 1
        self_sim = profile.cosine_similarity(profile)
        assert abs(self_sim - 1.0) < 1e-9

    def test_shannon_diversity(self, profile):
        h = profile.shannon_diversity()
        assert h > 0

    def test_richness(self, profile):
        assert profile.richness() == 4
        assert profile.richness(threshold=10.0) == 3

    def test_top_compounds(self, profile):
        top = profile.top_compounds(2)
        assert len(top) == 2
        assert top[0][1] >= top[1][1]

    def test_filter_by_abundance(self, profile):
        filtered = profile.filter_by_abundance(min_abundance=10.0)
        assert filtered.n_compounds == 3

    def test_relative_abundance(self, profile, compounds):
        from metainformant.phenotype.chemical.compound import Compound

        ra = profile.relative_abundance(Compound("n-C27"))
        assert abs(ra - 45.0 / 100.0) < 1e-9

    def test_compound_ratio(self, profile):
        from metainformant.phenotype.chemical.compound import Compound

        ratio = profile.compound_ratio(Compound("n-C27"), Compound("n-C25"))
        assert abs(ratio - 1.5) < 1e-9


class TestChemicalDistanceMatrix:
    def test_distance_matrix(self):
        from metainformant.phenotype.chemical.compound import Compound
        from metainformant.phenotype.chemical.profile import ChemicalProfile, distance_matrix

        c1, c2 = Compound("A"), Compound("B")
        profiles = [
            ChemicalProfile("s1", {c1: 10, c2: 20}),
            ChemicalProfile("s2", {c1: 15, c2: 25}),
            ChemicalProfile("s3", {c1: 50, c2: 5}),
        ]
        mat = distance_matrix(profiles, "bray_curtis")
        assert len(mat) == 3
        assert mat[0][0] == 0.0
        assert mat[0][1] == mat[1][0]

    def test_distance_matrix_euclidean(self):
        from metainformant.phenotype.chemical.compound import Compound
        from metainformant.phenotype.chemical.profile import ChemicalProfile, distance_matrix

        c1 = Compound("A")
        p1 = ChemicalProfile("s1", {c1: 0})
        p2 = ChemicalProfile("s2", {c1: 3})
        mat = distance_matrix([p1, p2], "euclidean")
        assert abs(mat[0][1] - 3.0) < 1e-9


class TestMarkerCompounds:
    def test_identify_markers(self):
        from metainformant.phenotype.chemical.compound import Compound
        from metainformant.phenotype.chemical.profile import ChemicalProfile, identify_marker_compounds

        c1, c2, c3 = Compound("A"), Compound("B"), Compound("C")
        group_a = [
            ChemicalProfile("a1", {c1: 100, c2: 10, c3: 5}),
            ChemicalProfile("a2", {c1: 90, c2: 12, c3: 6}),
        ]
        group_b = [
            ChemicalProfile("b1", {c1: 10, c2: 80, c3: 5}),
            ChemicalProfile("b2", {c1: 12, c2: 75, c3: 6}),
        ]
        markers = identify_marker_compounds(group_a, group_b, min_fold_change=2.0)
        assert len(markers) >= 2
        # c1 and c2 should be markers
        compound_names = [m["compound"].name for m in markers]
        assert "A" in compound_names
        assert "B" in compound_names


# ============================================================
# ELECTRONIC / TRACKING SUBMODULE
# ============================================================


class TestTrajectory:
    @pytest.fixture
    def trajectory(self):
        from metainformant.phenotype.electronic.tracking import TrackingPoint, Trajectory

        # Square path: (0,0) -> (10,0) -> (10,10) -> (0,10) -> (0,0)
        points = [
            TrackingPoint(x=0, y=0, timestamp=0),
            TrackingPoint(x=10, y=0, timestamp=1),
            TrackingPoint(x=10, y=10, timestamp=2),
            TrackingPoint(x=0, y=10, timestamp=3),
            TrackingPoint(x=0, y=0, timestamp=4),
        ]
        return Trajectory("ant_1", points)

    @pytest.fixture
    def straight_trajectory(self):
        from metainformant.phenotype.electronic.tracking import TrackingPoint, Trajectory

        points = [
            TrackingPoint(x=0, y=0, timestamp=0),
            TrackingPoint(x=5, y=0, timestamp=1),
            TrackingPoint(x=10, y=0, timestamp=2),
        ]
        return Trajectory("ant_2", points)

    def test_total_distance(self, trajectory):
        dist = trajectory.total_distance()
        assert abs(dist - 40.0) < 1e-9

    def test_net_displacement(self, trajectory):
        # Returns to origin
        assert trajectory.net_displacement() < 1e-9

    def test_sinuosity_circular(self, trajectory):
        sin = trajectory.sinuosity()
        assert sin == float("inf")

    def test_sinuosity_straight(self, straight_trajectory):
        sin = straight_trajectory.sinuosity()
        assert abs(sin - 1.0) < 1e-9

    def test_velocity_profile(self, trajectory):
        vels = trajectory.velocity_profile()
        assert len(vels) == 4
        for v in vels:
            assert abs(v - 10.0) < 1e-9

    def test_mean_velocity(self, trajectory):
        assert abs(trajectory.mean_velocity() - 10.0) < 1e-9

    def test_duration(self, trajectory):
        assert trajectory.duration == 4.0

    def test_turning_angles(self, trajectory):
        angles = trajectory.turning_angles()
        assert len(angles) == 3
        for a in angles:
            assert abs(a - math.pi / 2) < 0.01

    def test_home_range_mcp(self, trajectory):
        area = trajectory.home_range_mcp()
        assert abs(area - 100.0) < 1e-6

    def test_bounding_box(self, trajectory):
        bb = trajectory.bounding_box()
        assert bb["x_min"] == 0
        assert bb["x_max"] == 10
        assert bb["y_min"] == 0
        assert bb["y_max"] == 10

    def test_activity_budget(self, trajectory):
        budget = trajectory.activity_budget(velocity_threshold=5.0)
        assert budget["active"] == 1.0

    def test_segment_by_time(self, trajectory):
        segs = trajectory.segment_by_time(2.0)
        assert len(segs) >= 2


class TestDetectInteractions:
    def test_interactions(self):
        from metainformant.phenotype.electronic.tracking import (
            TrackingPoint,
            Trajectory,
            detect_interactions,
        )

        t1 = Trajectory("ant_a", [TrackingPoint(0, 0, 0), TrackingPoint(1, 0, 1)])
        t2 = Trajectory("ant_b", [TrackingPoint(0.5, 0, 0), TrackingPoint(5, 5, 1)])

        interactions = detect_interactions([t1, t2], distance_threshold=1.0, time_threshold=0.5)
        assert len(interactions) >= 1
        assert interactions[0]["entity_a"] == "ant_a"


# ============================================================
# MORPHOLOGICAL SUBMODULE
# ============================================================


class TestMorphometricProfile:
    @pytest.fixture
    def profile(self):
        from metainformant.phenotype.morphological import Measurement, MorphometricProfile

        return MorphometricProfile(
            "specimen_1",
            measurements=[
                Measurement("head_width", 1.2, "mm"),
                Measurement("head_length", 1.5, "mm"),
                Measurement("thorax_length", 2.0, "mm"),
                Measurement("femur_length", 1.8, "mm"),
            ],
        )

    def test_calculate_index(self, profile):
        ci = profile.calculate_index("cephalic_index", "head_width", "head_length")
        assert abs(ci - 80.0) < 1e-9

    def test_geometric_mean_size(self, profile):
        gm = profile.geometric_mean_size()
        assert gm > 0
        expected = math.exp((math.log(1.2) + math.log(1.5) + math.log(2.0) + math.log(1.8)) / 4)
        assert abs(gm - expected) < 1e-9

    def test_size_corrected_values(self, profile):
        corrected = profile.size_corrected_values()
        assert len(corrected) == 4
        gm = profile.geometric_mean_size()
        assert abs(corrected["head_width"] - 1.2 / gm) < 1e-9

    def test_values_dict(self, profile):
        vals = profile.values_dict()
        assert vals["head_width"] == 1.2

    def test_measurement_names(self, profile):
        names = profile.measurement_names
        assert "head_width" in names
        assert len(names) == 4


class TestMorphologicalFunctions:
    @pytest.fixture
    def profiles(self):
        from metainformant.phenotype.morphological import Measurement, MorphometricProfile

        profiles = []
        # Create specimens with known scaling
        for i, (hw, hl, tl, fl) in enumerate(
            [
                (1.0, 1.2, 1.8, 1.5),
                (1.5, 1.8, 2.7, 2.25),
                (2.0, 2.4, 3.6, 3.0),
                (2.5, 3.0, 4.5, 3.75),
                (3.0, 3.6, 5.4, 4.5),
            ]
        ):
            profiles.append(
                MorphometricProfile(
                    f"sp_{i}",
                    measurements=[
                        Measurement("head_width", hw, "mm"),
                        Measurement("head_length", hl, "mm"),
                        Measurement("thorax_length", tl, "mm"),
                        Measurement("femur_length", fl, "mm"),
                    ],
                )
            )
        return profiles

    def test_coefficient_of_variation(self, profiles):
        from metainformant.phenotype.morphological.profile import coefficient_of_variation

        cv = coefficient_of_variation(profiles, "head_width")
        assert cv is not None
        assert cv > 0

    def test_allometric_regression(self, profiles):
        from metainformant.phenotype.morphological.profile import allometric_regression

        result = allometric_regression(profiles, "head_width", "thorax_length")
        assert "slope" in result
        assert result["r_squared"] > 0.99  # Perfect scaling in test data
        assert abs(result["slope"] - 1.0) < 0.1  # Near isometric

    def test_compare_profiles(self, profiles):
        from metainformant.phenotype.morphological.profile import compare_profiles

        result = compare_profiles(profiles[0], profiles[4])
        assert result["shared_measurements"] == 4
        assert "head_width" in result["comparisons"]

    def test_summary_statistics(self, profiles):
        from metainformant.phenotype.morphological.profile import summary_statistics

        stats = summary_statistics(profiles)
        assert "head_width" in stats
        assert stats["head_width"]["n"] == 5
        assert stats["head_width"]["mean"] == 2.0


class TestMeasurementConversion:
    def test_mm_to_cm(self):
        from metainformant.phenotype.morphological import Measurement

        m = Measurement("hw", 15.0, "mm")
        converted = m.convert("cm")
        assert abs(converted.value - 1.5) < 1e-9

    def test_cm_to_mm(self):
        from metainformant.phenotype.morphological import Measurement

        m = Measurement("hw", 1.5, "cm")
        converted = m.convert("mm")
        assert abs(converted.value - 15.0) < 1e-9


# ============================================================
# SONIC SUBMODULE
# ============================================================


class TestAcousticSignal:
    @pytest.fixture
    def tone_440(self):
        from metainformant.phenotype.sonic.signal import AcousticSignal

        return AcousticSignal.generate_tone(440.0, duration=1.0, sample_rate=44100)

    @pytest.fixture
    def multi_tone(self):
        """Signal with 200Hz + 800Hz components."""
        sr = 44100
        t = np.arange(sr) / sr
        waveform = np.sin(2 * np.pi * 200 * t) + 0.5 * np.sin(2 * np.pi * 800 * t)
        from metainformant.phenotype.sonic.signal import AcousticSignal

        return AcousticSignal(waveform=waveform, sample_rate=sr)

    def test_duration(self, tone_440):
        assert abs(tone_440.duration - 1.0) < 0.001

    def test_n_samples(self, tone_440):
        assert tone_440.n_samples == 44100

    def test_dominant_frequency(self, tone_440):
        freq = tone_440.dominant_frequency()
        assert abs(freq - 440.0) < 2.0  # Within 2 Hz

    def test_power_spectrum(self, tone_440):
        freqs, power = tone_440.power_spectrum()
        assert len(freqs) == len(power)
        peak_idx = np.argmax(power)
        assert abs(freqs[peak_idx] - 440.0) < 2.0

    def test_spectral_centroid(self, tone_440):
        centroid = tone_440.spectral_centroid()
        assert abs(centroid - 440.0) < 5.0

    def test_spectral_bandwidth(self, multi_tone):
        bw = multi_tone.spectral_bandwidth()
        assert bw > 0  # Multi-tone should have bandwidth

    def test_band_energy(self, multi_tone):
        low_energy = multi_tone.band_energy(100, 400)
        high_energy = multi_tone.band_energy(600, 1000)
        assert low_energy > 0
        assert high_energy > 0
        # 200Hz is stronger than 800Hz (amplitude 1.0 vs 0.5)
        assert low_energy > high_energy

    def test_rms_energy(self, tone_440):
        rms = tone_440.rms_energy()
        # RMS of sin wave = 1/sqrt(2) ≈ 0.707
        assert abs(rms - 1.0 / math.sqrt(2)) < 0.01

    def test_peak_amplitude(self, tone_440):
        peak = tone_440.peak_amplitude()
        assert abs(peak - 1.0) < 0.01

    def test_zero_crossing_rate(self, tone_440):
        zcr = tone_440.zero_crossing_rate()
        # 440 Hz tone crosses zero ~880 times per second
        assert abs(zcr - 880) < 10

    def test_spectrogram(self, tone_440):
        times, freqs, mag = tone_440.spectrogram(window_size=1024, hop_size=512)
        assert len(times) > 0
        assert len(freqs) == 513  # 1024/2 + 1
        assert mag.shape[0] == 513

    def test_trim(self, tone_440):
        trimmed = tone_440.trim(0.1, 0.5)
        assert abs(trimmed.duration - 0.4) < 0.001

    def test_normalize_amplitude(self, multi_tone):
        normed = multi_tone.normalize_amplitude()
        assert abs(normed.peak_amplitude() - 1.0) < 1e-9

    def test_generate_tone_metadata(self, tone_440):
        assert tone_440.metadata["type"] == "tone"
        assert tone_440.metadata["frequency"] == 440.0


class TestSyllableDetection:
    def test_detect_syllables_in_bursts(self):
        """Create signal with clear on/off bursts."""
        from metainformant.phenotype.sonic.signal import AcousticSignal

        sr = 44100
        silence = np.zeros(int(sr * 0.1))  # 100ms silence
        burst = np.sin(2 * np.pi * 1000 * np.arange(int(sr * 0.05)) / sr)  # 50ms burst

        waveform = np.concatenate([silence, burst, silence, burst, silence, burst, silence])
        signal = AcousticSignal(waveform=waveform, sample_rate=sr)

        syllables = signal.detect_syllables(threshold_factor=1.5)
        assert len(syllables) >= 2  # At least 2 detected

    def test_temporal_pattern(self):
        from metainformant.phenotype.sonic.signal import AcousticSignal

        sr = 44100
        silence = np.zeros(int(sr * 0.1))
        burst = np.sin(2 * np.pi * 1000 * np.arange(int(sr * 0.05)) / sr)
        waveform = np.concatenate([silence, burst, silence, burst, silence])
        signal = AcousticSignal(waveform=waveform, sample_rate=sr)

        pattern = signal.temporal_pattern()
        assert "syllable_count" in pattern
        assert "duty_cycle" in pattern
        assert 0 <= pattern["duty_cycle"] <= 1


# ============================================================
# ANALYSIS / LIFE COURSE
# ============================================================


class TestLifeCourseAnalysis:
    @pytest.fixture
    def sequences(self):
        from metainformant.phenotype.analysis.life_course import Event, EventSequence

        seqs = []
        for i in range(5):
            events = [
                Event(timestamp=0, event_type="job_start", description="First job"),
                Event(timestamp=365, event_type="social_event", description="Social"),
                Event(timestamp=730, event_type="job_end", description="Left job"),
                Event(timestamp=750, event_type="job_start", description="New job"),
                Event(timestamp=1000, event_type="health_event", description="Checkup"),
            ]
            if i > 2:
                events.extend(
                    [
                        Event(timestamp=1100, event_type="job_end"),
                        Event(timestamp=1200, event_type="job_start"),
                        Event(timestamp=1300, event_type="job_end"),
                        Event(timestamp=1400, event_type="job_start"),
                        Event(timestamp=1500, event_type="health_event"),
                        Event(timestamp=1600, event_type="health_event"),
                        Event(timestamp=1700, event_type="health_event"),
                        Event(timestamp=1800, event_type="health_event"),
                        Event(timestamp=1900, event_type="health_event"),
                        Event(timestamp=2000, event_type="health_event"),
                    ]
                )
            seqs.append(EventSequence(person_id=f"person_{i}", events=events))
        return seqs

    def test_extract_phenotypes(self, sequences):
        from metainformant.phenotype.analysis.life_course import extract_phenotypes_from_events

        phenotypes = extract_phenotypes_from_events(sequences[0])
        assert "career_changes" in phenotypes

    def test_identify_trajectory_patterns(self, sequences):
        from metainformant.phenotype.analysis.life_course import identify_trajectory_patterns

        patterns = identify_trajectory_patterns(sequences)
        assert "patterns" in patterns
        assert "transitions" in patterns
        assert patterns["n_sequences"] == 5

    def test_identify_critical_periods(self, sequences):
        from metainformant.phenotype.analysis.life_course import identify_critical_periods

        result = identify_critical_periods(sequences, [(0, 500), (500, 1000), (1000, 2000)])
        assert "period_importance" in result
        assert len(result["period_importance"]) == 3

    def test_analyze_trajectories(self, sequences):
        from metainformant.phenotype.analysis.life_course import analyze_life_course_trajectories

        result = analyze_life_course_trajectories(sequences)
        assert result["total_sequences"] == 5
        assert "trajectory_patterns" in result

    def test_analyze_life_course(self, sequences):
        from metainformant.phenotype.analysis.life_course import analyze_life_course

        result = analyze_life_course(sequences)
        assert result["n_sequences"] == 5
        assert "trajectory_patterns" in result

    def test_create_report(self, sequences, tmp_path):
        from metainformant.phenotype.analysis.life_course import create_life_course_report

        report = create_life_course_report(sequences, output_path=tmp_path / "report.txt")
        assert "LIFE COURSE PHENOTYPE ANALYSIS REPORT" in report
        assert (tmp_path / "report.txt").exists()

    def test_predict_outcomes(self, sequences):
        from metainformant.phenotype.analysis.life_course import predict_life_course_outcomes

        result = predict_life_course_outcomes(sequences, prediction_horizon=2.0)
        assert result["prediction_horizon"] == 2.0


# ============================================================
# WORKFLOW PIPELINE
# ============================================================


class TestPipelineConfig:
    def test_create_default(self):
        from metainformant.phenotype.workflow.pipeline import PipelineConfig

        config = PipelineConfig()
        assert config.name == "phenotype_pipeline"
        assert "load" in config.steps

    def test_from_dict(self):
        from metainformant.phenotype.workflow.pipeline import PipelineConfig

        config = PipelineConfig.from_dict(
            {
                "name": "test_pipeline",
                "phenotype_types": ["behavioral", "chemical"],
                "steps": ["load", "analyze"],
            }
        )
        assert config.name == "test_pipeline"
        assert len(config.phenotype_types) == 2

    def test_validate(self):
        from metainformant.phenotype.workflow.pipeline import PipelineConfig

        config = PipelineConfig(phenotype_types=["unknown_type"])
        issues = config.validate()
        assert len(issues) > 0


class TestPhenotypePipeline:
    def test_run_with_data(self):
        from metainformant.phenotype.workflow.pipeline import PhenotypePipeline, PipelineConfig

        config = PipelineConfig(
            name="test",
            phenotype_types=["morphological"],
            steps=["load", "validate", "analyze", "summarize"],
        )
        pipeline = PhenotypePipeline(config)
        result = pipeline.run(data=[{"specimen": "ant_1", "hw": 1.2}])

        assert result.success
        assert len(result.outputs) == 4
        assert result.metrics["steps_completed"] == 4

    def test_run_empty(self):
        from metainformant.phenotype.workflow.pipeline import PhenotypePipeline, PipelineConfig

        config = PipelineConfig(steps=["load", "validate", "summarize"])
        pipeline = PhenotypePipeline(config)
        result = pipeline.run()
        assert result.success  # No data is handled gracefully

    def test_result_to_dict(self):
        from metainformant.phenotype.workflow.pipeline import PhenotypePipeline, PipelineConfig

        config = PipelineConfig(steps=["load", "summarize"])
        pipeline = PhenotypePipeline(config)
        result = pipeline.run(data={"test": True})
        d = result.to_dict()
        assert "success" in d
        assert "timestamp" in d

    def test_result_save_json(self, tmp_path):
        from metainformant.phenotype.workflow.pipeline import PhenotypePipeline, PipelineConfig

        config = PipelineConfig(steps=["load"])
        pipeline = PhenotypePipeline(config)
        result = pipeline.run(data=[])
        out_path = tmp_path / "result.json"
        result.save_json(out_path)
        assert out_path.exists()
        loaded = json.loads(out_path.read_text())
        assert "success" in loaded

    def test_export_step(self, tmp_path):
        from metainformant.phenotype.workflow.pipeline import PhenotypePipeline, PipelineConfig

        config = PipelineConfig(
            steps=["load", "validate", "analyze", "summarize", "export"],
            output_path=str(tmp_path / "export.json"),
        )
        pipeline = PhenotypePipeline(config)
        result = pipeline.run(data=[{"x": 1}])
        assert result.success
        assert (tmp_path / "export.json").exists()

    def test_custom_step(self):
        from metainformant.phenotype.workflow.pipeline import PhenotypePipeline, PipelineConfig

        config = PipelineConfig(steps=["load", "custom_check"])
        pipeline = PhenotypePipeline(config)
        pipeline.register_step("custom_check", lambda: {"status": "custom_ok"})
        result = pipeline.run(data=[1, 2, 3])
        assert result.success
        assert result.outputs["custom_check"]["status"] == "custom_ok"


# ============================================================
# INTEGRATION MODULE
# ============================================================


class TestPhenotypeGenotypeAssociation:
    def test_basic_association(self):
        from metainformant.phenotype.integration.cross_omic import phenotype_genotype_association

        phenotypes = {"s1": 10.0, "s2": 20.0, "s3": 30.0, "s4": 40.0, "s5": 50.0}
        genotypes = {
            "var_1": [0, 1, 1, 2, 2],
            "var_2": [2, 2, 1, 0, 0],
        }
        result = phenotype_genotype_association(phenotypes, genotypes)
        assert result["n_samples"] == 5
        assert "var_1" in result["associations"]
        # var_1 should have positive correlation with phenotype
        assert result["associations"]["var_1"]["beta"] > 0

    def test_empty_input(self):
        from metainformant.phenotype.integration.cross_omic import phenotype_genotype_association

        result = phenotype_genotype_association({}, {})
        assert "error" in result


class TestTraitExpressionCorrelation:
    def test_basic_correlation(self):
        from metainformant.phenotype.integration.cross_omic import trait_expression_correlation

        traits = {"s1": 1.0, "s2": 2.0, "s3": 3.0, "s4": 4.0}
        expression = {
            "gene_A": {"s1": 10, "s2": 20, "s3": 30, "s4": 40},
            "gene_B": {"s1": 40, "s2": 30, "s3": 20, "s4": 10},
        }
        result = trait_expression_correlation(traits, expression)
        assert result["n_genes_tested"] == 2
        assert result["correlations"]["gene_A"]["correlation"] > 0.9
        assert result["correlations"]["gene_B"]["correlation"] < -0.9


class TestMultiPhenotypeIntegration:
    def test_cross_phenotype(self):
        from metainformant.phenotype.integration.cross_omic import multi_phenotype_integration

        phenotype_matrices = {
            "body_size": {"s1": 1.0, "s2": 2.0, "s3": 3.0, "s4": 4.0},
            "head_width": {"s1": 0.5, "s2": 1.0, "s3": 1.5, "s4": 2.0},
            "aggression": {"s1": 4.0, "s2": 3.0, "s3": 2.0, "s4": 1.0},
        }
        result = multi_phenotype_integration(phenotype_matrices)
        assert result["n_phenotypes"] == 3
        assert result["n_common_samples"] == 4
        # Body size and head width should be strongly positive
        assert result["correlation_matrix"][0][1] > 0.9


class TestGxEInteraction:
    def test_basic_gxe(self):
        from metainformant.phenotype.integration.cross_omic import phenotype_environment_interaction

        phenotypes = {"s1": 10, "s2": 20, "s3": 30, "s4": 40, "s5": 50}
        genotypes = {"var_1": [0, 1, 1, 2, 2]}
        environment = {"s1": 1.0, "s2": 2.0, "s3": 3.0, "s4": 4.0, "s5": 5.0}

        result = phenotype_environment_interaction(phenotypes, genotypes, environment)
        assert result["n_samples"] == 5
        assert "var_1" in result["interactions"]


# ============================================================
# TOP-LEVEL IMPORTS
# ============================================================


class TestTopLevelImports:
    def test_all_submodules_importable(self):
        from metainformant import phenotype

        assert hasattr(phenotype, "behavior")
        assert hasattr(phenotype, "chemical")
        assert hasattr(phenotype, "electronic")
        assert hasattr(phenotype, "morphological")
        assert hasattr(phenotype, "sonic")
        assert hasattr(phenotype, "analysis")
        assert hasattr(phenotype, "data")
        assert hasattr(phenotype, "visualization")
        assert hasattr(phenotype, "workflow")
        assert hasattr(phenotype, "integration")

    def test_key_classes_importable(self):
        from metainformant.phenotype import (
            Ethogram,
            BehaviorSequence,
            Compound,
            ChemicalProfile,
            TrackingPoint,
            Trajectory,
            Measurement,
            MorphometricProfile,
            AcousticSignal,
            PhenotypePipeline,
            PipelineConfig,
            PipelineResult,
        )

        assert all(
            [
                Ethogram,
                BehaviorSequence,
                Compound,
                ChemicalProfile,
                TrackingPoint,
                Trajectory,
                Measurement,
                MorphometricProfile,
                AcousticSignal,
                PhenotypePipeline,
                PipelineConfig,
                PipelineResult,
            ]
        )

"""Tests for the diagnostics module."""

import os
import pytest
import tempfile
from pathlib import Path

from beast2py.config import ConfigParser
from beast2py.diagnostics import (
    ConflictDetector,
    DiagnosticEngine,
    ReportGenerator,
)

from .conftest import TEST_FASTA, TEST_MULTI_CAL_YAML


class TestConflictDetector:
    """Test the ConflictDetector class."""

    def setup_method(self, method):
        """Create a temp directory for each test."""
        self.tmpdir = Path(tempfile.mkdtemp())
        fasta_path = self.tmpdir / "test.fasta"
        fasta_path.write_text(TEST_FASTA)

    def test_detect_no_conflicts(self):
        """Test that consistent calibrations produce no conflicts."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_MULTI_CAL_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        conflicts = ConflictDetector.detect_conflicts(config.calibrations)

        # calA (5.0) is nested in calB (10.0) which is nested in rootCal (15-30)
        # These should be temporally consistent
        error_conflicts = [c for c in conflicts if c.severity == "error"]
        assert len(error_conflicts) == 0

    def test_detect_temporal_conflict(self):
        """Test detection of temporal inconsistency."""
        from beast2py.models import (
            CalibrationPoint,
            DistributionConfig,
        )

        # Parent clade younger than nested clade
        cal_parent = CalibrationPoint(
            name="parent",
            taxa=["A", "B", "C"],
            monophyletic=True,
            distribution=DistributionConfig(
                type="normal",
                parameters={"mean": 3.0, "sigma": 0.5},
            ),
        )
        cal_child = CalibrationPoint(
            name="child",
            taxa=["A", "B"],
            monophyletic=True,
            distribution=DistributionConfig(
                type="normal",
                parameters={"mean": 10.0, "sigma": 1.0},
            ),
        )

        conflicts = ConflictDetector.detect_conflicts([cal_parent, cal_child])
        error_conflicts = [c for c in conflicts if c.severity == "error"]
        assert len(error_conflicts) > 0
        assert error_conflicts[0].conflict_type == "temporal"

    def test_detect_monophyly_conflict(self):
        """Test detection of monophyly conflict."""
        from beast2py.models import (
            CalibrationPoint,
            DistributionConfig,
        )

        # Overlapping but non-nested taxon sets
        cal_a = CalibrationPoint(
            name="calA",
            taxa=["A", "B", "C"],
            monophyletic=True,
            distribution=DistributionConfig(
                type="normal",
                parameters={"mean": 5.0, "sigma": 0.5},
            ),
        )
        cal_b = CalibrationPoint(
            name="calB",
            taxa=["B", "C", "D"],
            monophyletic=True,
            distribution=DistributionConfig(
                type="normal",
                parameters={"mean": 10.0, "sigma": 1.0},
            ),
        )

        conflicts = ConflictDetector.detect_conflicts([cal_a, cal_b])
        monophyly_conflicts = [c for c in conflicts if c.conflict_type == "monophyly"]
        assert len(monophyly_conflicts) > 0

    def test_no_conflict_single_calibration(self):
        """Test that single calibration produces no conflicts."""
        from beast2py.models import (
            CalibrationPoint,
            DistributionConfig,
        )

        cal = CalibrationPoint(
            name="only",
            taxa=["A", "B"],
            monophyletic=True,
            distribution=DistributionConfig(
                type="normal",
                parameters={"mean": 5.0, "sigma": 0.5},
            ),
        )

        conflicts = ConflictDetector.detect_conflicts([cal])
        assert len(conflicts) == 0


class TestDiagnosticEngine:
    """Test the DiagnosticEngine class."""

    def setup_method(self, method):
        """Create a temp directory for each test."""
        self.tmpdir = Path(tempfile.mkdtemp())
        fasta_path = self.tmpdir / "test.fasta"
        fasta_path.write_text(TEST_FASTA)

    def test_run_diagnostics(self):
        """Test running the full diagnostic workflow."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_MULTI_CAL_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        output_dir = self.tmpdir / "diagnostics"
        report = DiagnosticEngine.run(
            config,
            output_dir=output_dir,
            run_sensitivity=False,
            run_visualization=False,
        )

        assert len(report.summaries) == 3
        assert report.config == config

    def test_run_diagnostics_with_report(self):
        """Test generating an HTML report."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_MULTI_CAL_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        output_dir = self.tmpdir / "diagnostics"
        report = DiagnosticEngine.run(
            config,
            output_dir=output_dir,
            run_sensitivity=False,
            run_visualization=False,
        )

        report_path = self.tmpdir / "report.html"
        ReportGenerator.generate_report(
            config=config,
            conflicts=report.conflicts,
            summaries=report.summaries,
            sensitivity_xmls=[],
            visualization_file=None,
            output_file=report_path,
        )

        assert os.path.exists(report_path)
        with open(report_path, "r") as f:
            content = f.read()
        assert "Beast2Py Diagnostic Report" in content
        assert "test_multi_cal" in content


class TestOverlapCriterionIsRelativeToDefault:
    """Supplement §S1.2 documents the ratio as the default criterion.

    W grows with node depth, so the absolute criterion is opt-in: two priors that
    are far apart relative to their own width must not be called redundant unless
    the caller asks for the absolute measure.
    """

    WIDTH = 2 * 1.96 * 0.05          # each prior's own 95% interval width

    def _pair(self):
        from beast2py.models import CalibrationPoint, DistributionConfig
        return [
            CalibrationPoint(name="A", taxa=["t1", "t2"],
                             distribution=DistributionConfig(
                                 type="normal", parameters={"mean": 10.0, "sigma": 0.05})),
            CalibrationPoint(name="B", taxa=["t3", "t4"],
                             distribution=DistributionConfig(
                                 type="normal", parameters={"mean": 10.5, "sigma": 0.05})),
        ]

    @staticmethod
    def _overlap(conflicts):
        return [c for c in conflicts if c.conflict_type == "distribution_overlap"]

    def test_thresholds_are_the_ones_the_supplement_quotes(self):
        assert ConflictDetector.OVERLAP_RATIO_THRESHOLD == 0.15
        assert ConflictDetector.OVERLAP_DISTANCE_THRESHOLD == 2.0

    def test_absolute_criterion_is_not_applied_by_default(self):
        # W = 0.5 is below the absolute 2.0, but R = 0.5 / 0.196 = 2.55 is far
        # above 0.15, so this pair is informative, not redundant.
        assert 0.5 < ConflictDetector.OVERLAP_DISTANCE_THRESHOLD
        out = ConflictDetector.detect_conflicts(self._pair())
        assert self._overlap(out) == []

    @pytest.mark.parametrize("measure", ["absolute", "both"])
    def test_absolute_criterion_applies_when_requested(self, measure):
        out = ConflictDetector.detect_conflicts(self._pair(),
                                                settings={"overlap_measure": measure})
        hits = self._overlap(out)
        assert len(hits) == 1
        assert hits[0].severity == "warning"

    def test_relative_requested_explicitly_matches_the_default(self):
        assert self._overlap(ConflictDetector.detect_conflicts(
            self._pair(), settings={"overlap_measure": "relative"})) == []

    def test_unknown_measure_is_rejected_not_silently_ignored(self):
        from beast2py.config import ConfigError
        with pytest.raises((ConfigError, ValueError)):
            ConflictDetector.detect_conflicts(
                self._pair(), settings={"overlap_measure": "whichever"})

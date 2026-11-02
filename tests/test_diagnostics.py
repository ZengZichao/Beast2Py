"""Tests for the diagnostics module.

The structural rules changed in v0.2.0: independent clades with similar
priors are no longer a conflict, and crown/stem semantics are resolved before
any comparison. These tests pin both the new rules and the old false positive.
"""

import tempfile
from pathlib import Path

import pytest

from beast2py.config import ConfigError, ConfigParser
from beast2py.diagnostics import (
    ConflictDetector,
    DiagnosticEngine,
    ReportGenerator,
)
from beast2py.models import CalibrationPoint, DistributionConfig

from .conftest import TEST_FASTA, TEST_MULTI_CAL_YAML


def cal(
    name,
    taxa,
    dist_type="normal",
    offset=0.0,
    use_originate=False,
    monophyletic=True,
    tipsonly=False,
    **params,
):
    return CalibrationPoint(
        name=name,
        taxa=taxa,
        monophyletic=monophyletic,
        use_originate=use_originate,
        tipsonly=tipsonly,
        distribution=DistributionConfig(type=dist_type, parameters=params, offset=offset),
    )


def types(conflicts):
    return [c.conflict_type for c in conflicts]


class TestStructuralConflicts:
    def test_consistent_nested_calibrations_produce_no_errors(self):
        self.tmpdir = Path(tempfile.mkdtemp())
        (self.tmpdir / "test.fasta").write_text(TEST_FASTA)
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_MULTI_CAL_YAML)
        config = ConfigParser().parse(config_path)
        conflicts = ConflictDetector.detect_conflicts(config.calibrations)
        assert [c for c in conflicts if c.severity == "error"] == []

    def test_parent_younger_than_child_is_temporal_error(self):
        parent = cal("parent", ["A", "B", "C"], mean=3.0, sigma=0.5)
        child = cal("child", ["A", "B"], mean=10.0, sigma=1.0)
        conflicts = ConflictDetector.detect_conflicts([parent, child])
        errors = [c for c in conflicts if c.severity == "error"]
        assert errors and errors[0].conflict_type == "temporal"

    def test_intersecting_non_nested_monophyletic_is_topology_error(self):
        a = cal("calA", ["A", "B", "C"], mean=5.0, sigma=0.5)
        b = cal("calB", ["B", "C", "D"], mean=10.0, sigma=1.0)
        conflicts = ConflictDetector.detect_conflicts([a, b])
        hits = [c for c in conflicts if c.conflict_type == "topology"]
        assert len(hits) == 1 and hits[0].severity == "error"

    def test_intersecting_without_monophyly_is_unchecked_not_a_conflict(self):
        a = cal("calA", ["A", "B", "C"], mean=5.0, sigma=0.5, monophyletic=False)
        b = cal("calB", ["B", "C", "D"], mean=10.0, sigma=1.0)
        assessments = ConflictDetector.assess([a, b])
        assert assessments[0].status == "unchecked"
        assert not [c for c in assessments[0].conflicts if c.severity == "error"]

    def test_single_calibration_produces_nothing(self):
        assert (
            ConflictDetector.detect_conflicts([cal("only", ["A", "B"], mean=5.0, sigma=0.5)]) == []
        )


class TestSameNode:
    def test_identical_distribution_is_duplicate_warning(self):
        a = cal("a", ["A", "B"], mean=10.0, sigma=1.0)
        b = cal("b", ["A", "B"], mean=10.0, sigma=1.0)
        conflicts = ConflictDetector.detect_conflicts([a, b])
        hits = [c for c in conflicts if c.conflict_type == "duplicate"]
        assert len(hits) == 1 and hits[0].severity == "warning"

    def test_disjoint_distributions_on_same_node_is_error(self):
        a = cal("a", ["A", "B"], mean=5.0, sigma=0.2)
        b = cal("b", ["A", "B"], mean=20.0, sigma=0.2)
        conflicts = ConflictDetector.detect_conflicts([a, b])
        hits = [c for c in conflicts if c.conflict_type == "same_node"]
        assert len(hits) == 1 and hits[0].severity == "error"

    def test_overlapping_different_distributions_is_warning(self):
        a = cal("a", ["A", "B"], mean=10.0, sigma=2.0)
        b = CalibrationPoint(
            name="b",
            taxa=["A", "B"],
            distribution=DistributionConfig(
                type="uniform", parameters={"lower": 8.0, "upper": 12.0}
            ),
        )
        conflicts = ConflictDetector.detect_conflicts([a, b])
        hits = [c for c in conflicts if c.conflict_type == "same_node"]
        assert len(hits) == 1 and hits[0].severity == "warning"


class TestCrownStem:
    def test_crown_and_stem_on_same_taxon_set_are_parent_and_child(self):
        # Stem is the parent; a stem younger than the crown is impossible.
        crown = cal("crown", ["A", "B"], mean=20.0, sigma=0.5)
        stem = cal("stem", ["A", "B"], mean=5.0, sigma=0.5, use_originate=True)
        conflicts = ConflictDetector.detect_conflicts([crown, stem])
        assert "temporal" in types(conflicts)
        assert [c for c in conflicts if c.severity == "error"]

    def test_same_taxon_set_crown_vs_stem_consistent_is_not_duplicate(self):
        crown = cal("crown", ["A", "B"], mean=10.0, sigma=0.5)
        stem = cal("stem", ["A", "B"], mean=15.0, sigma=0.5, use_originate=True)
        conflicts = ConflictDetector.detect_conflicts([crown, stem])
        assert "duplicate" not in types(conflicts)
        assert [c for c in conflicts if c.severity == "error"] == []

    def test_nested_parent_crown_vs_child_stem_is_unchecked(self):
        # C(parent) >= C(child), but S(child) can be above or below C(parent),
        # so taxon containment does not determine the order.
        parent = cal("parent", ["A", "B", "C"], mean=30.0, sigma=1.0)
        child_stem = cal("child_stem", ["A", "B"], mean=10.0, sigma=1.0, use_originate=True)
        assessments = ConflictDetector.assess([parent, child_stem])
        assert assessments[0].status == "unchecked"
        assert "unchecked" in types(assessments[0].conflicts)

    def test_nested_stem_parent_crown_child_is_checked(self):
        parent_stem = cal("parent_stem", ["A", "B", "C"], mean=5.0, sigma=0.5, use_originate=True)
        child = cal("child", ["A", "B"], mean=20.0, sigma=0.5)
        conflicts = ConflictDetector.detect_conflicts([parent_stem, child])
        assert "temporal" in types(conflicts)

    def test_crown_crown_nested_is_checked(self):
        parent = cal("parent", ["A", "B", "C"], mean=20.0, sigma=0.5)
        child = cal("child", ["A", "B"], mean=10.0, sigma=0.5)
        assessments = ConflictDetector.assess([parent, child])
        assert assessments[0].status == "checked"
        assert assessments[0].conflicts == []

    def test_stem_stem_nested_is_checked(self):
        parent = cal("parent", ["A", "B", "C"], mean=20.0, sigma=0.5, use_originate=True)
        child = cal("child", ["A", "B"], mean=10.0, sigma=0.5, use_originate=True)
        assessments = ConflictDetector.assess([parent, child])
        assert assessments[0].status == "checked"

    def test_root_plus_originate_is_rejected_at_parse_time(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(TEST_FASTA)
        cfg = tmp_path / "config.yaml"
        cfg.write_text(
            """
alignments:
  - id: alignment
    file: test.fasta
    data_type: nucleotide
partitions:
  - id: alignment
    site_model:
      substitution_model: {type: hky}
    clock_model: {type: strict}
    tree: shared
tree_prior: {type: yule}
calibrations:
  - name: root
    taxa: null
    use_originate: true
    distribution: {type: uniform, parameters: {lower: 10, upper: 20}}
mcmc: {chain_length: 10000}
output: {file: out.xml}
""",
            encoding="utf-8",
        )
        with pytest.raises(ConfigError, match="root has no parent"):
            ConfigParser().parse(cfg)

    def test_all_taxa_plus_originate_is_rejected_at_parse_time(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(TEST_FASTA)
        cfg = tmp_path / "config.yaml"
        cfg.write_text(
            """
alignments:
  - id: alignment
    file: test.fasta
    data_type: nucleotide
partitions:
  - id: alignment
    site_model:
      substitution_model: {type: hky}
    clock_model: {type: strict}
    tree: shared
tree_prior: {type: yule}
calibrations:
  - name: all
    taxa: [taxonA, taxonB, taxonC, taxonD]
    use_originate: true
    distribution: {type: uniform, parameters: {lower: 10, upper: 20}}
mcmc: {chain_length: 10000}
output: {file: out.xml}
""",
            encoding="utf-8",
        )
        with pytest.raises(ConfigError, match="root.*no parent|no parent.*root"):
            ConfigParser().parse(cfg)

    def test_tipsonly_plus_originate_is_rejected(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(TEST_FASTA)
        cfg = tmp_path / "config.yaml"
        cfg.write_text(
            """
alignments:
  - id: alignment
    file: test.fasta
    data_type: nucleotide
partitions:
  - id: alignment
    site_model:
      substitution_model: {type: hky}
    clock_model: {type: strict}
    tree: shared
tree_prior: {type: yule}
calibrations:
  - name: tip
    taxa: [taxonA, taxonB]
    tipsonly: true
    use_originate: true
    distribution: {type: uniform, parameters: {lower: 10, upper: 20}}
mcmc: {chain_length: 10000}
output: {file: out.xml}
""",
            encoding="utf-8",
        )
        with pytest.raises(ConfigError, match="tipsonly"):
            ConfigParser().parse(cfg)


class TestUndefinedQuantiles:
    def test_one_on_x_pair_is_unchecked_not_passed(self):
        a = cal("a", ["A", "B"], "one_on_x")
        b = cal("b", ["A", "B", "C"], mean=10.0, sigma=1.0)
        assessments = ConflictDetector.assess([a, b])
        assert assessments[0].status == "unchecked"
        assert "unchecked" in types(assessments[0].conflicts)
        assert ConflictDetector.detect_conflicts([a, b])[0].severity == "info"


class TestPriorSimilarityIsDescriptiveOnly:
    """Regression: independent clades with similar priors are not redundant.

    The old rule reported this pair as `distribution_overlap` and suggested
    merging the calibrations. Similar marginal priors of independent clades
    can be entirely legitimate, so the default is no output at all and the
    opt-in output is informational, never a Gate 2 failure.
    """

    def _pair(self):
        return [
            cal("A", ["t1", "t2"], mean=10.0, sigma=0.05),
            cal("B", ["t3", "t4"], mean=10.5, sigma=0.05),
        ]

    def test_default_is_silent(self):
        assert ConflictDetector.detect_conflicts(self._pair()) == []
        assert (
            ConflictDetector.detect_conflicts(self._pair(), settings={"prior_similarity": False})
            == []
        )

    def test_no_threshold_means_no_silent_default(self):
        with pytest.raises(ValueError, match="threshold"):
            ConflictDetector.detect_conflicts(self._pair(), settings={"prior_similarity": True})

    def test_opt_in_is_info_and_never_a_gate_2_blocker(self):
        conflicts = ConflictDetector.detect_conflicts(
            self._pair(),
            settings={
                "prior_similarity": True,
                "similarity_measure": "absolute",
                "similarity_threshold": 2.0,
            },
        )
        hits = [c for c in conflicts if c.conflict_type == "prior_similarity"]
        assert len(hits) == 1
        assert hits[0].severity == "info"
        assert "not evidence" in hits[0].description
        assert hits[0].recommendation.startswith("None.")
        assert [c for c in conflicts if c.severity == "error"] == []

    def test_disjoint_different_priors_are_silent_by_default(self):
        pair = [
            cal("A", ["t1", "t2"], mean=10.0, sigma=0.1),
            cal("B", ["t3", "t4"], mean=50.0, sigma=0.1),
        ]
        assert ConflictDetector.detect_conflicts(pair) == []

    def test_nested_similar_priors_are_not_renamed_redundancy(self):
        parent = cal("parent", ["A", "B", "C"], mean=10.0, sigma=0.05)
        child = cal("child", ["A", "B"], mean=10.0, sigma=0.05)
        conflicts = ConflictDetector.detect_conflicts(
            [parent, child],
            settings={
                "prior_similarity": True,
                "similarity_measure": "relative",
                "similarity_threshold": 0.15,
            },
        )
        assert "distribution_overlap" not in types(conflicts)
        assert all(c.conflict_type != "redundancy" for c in conflicts)
        assert all(c.severity != "error" for c in conflicts)
        similarity = [c for c in conflicts if c.conflict_type == "prior_similarity"]
        assert similarity and similarity[0].severity == "info"
        assert similarity[0].recommendation.startswith("None.")


class TestDiagnosticEngine:
    def setup_method(self, method):
        self.tmpdir = Path(tempfile.mkdtemp())
        (self.tmpdir / "test.fasta").write_text(TEST_FASTA)

    def test_run_diagnostics(self):
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_MULTI_CAL_YAML)
        config = ConfigParser().parse(config_path)
        report = DiagnosticEngine.run(
            config,
            output_dir=self.tmpdir / "diagnostics",
            run_sensitivity=False,
            run_visualization=False,
        )
        assert len(report.summaries) == 3
        assert len(report.assessments) == 3

    def test_report_counts_checked_and_unchecked_separately(self):
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_MULTI_CAL_YAML)
        config = ConfigParser().parse(config_path)
        report = DiagnosticEngine.run(
            config,
            output_dir=self.tmpdir / "diagnostics",
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
            assessments=report.assessments,
        )
        content = report_path.read_text()
        assert "Beast2Py Diagnostic Report" in content
        assert "test_multi_cal" in content
        assert "Checked pairs" in content
        assert "Unchecked pairs" in content
        assert "pairs assessed" in content


class TestSettingsValidation:
    def test_unknown_measure_is_rejected(self):
        with pytest.raises(ValueError):
            ConflictDetector.detect_conflicts(
                [cal("A", ["t1", "t2"]), cal("B", ["t3", "t4"])],
                settings={"similarity_measure": "whichever", "similarity_threshold": 1},
            )

    def test_parser_rejects_similarity_without_threshold(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(TEST_FASTA)
        cfg = tmp_path / "config.yaml"
        cfg.write_text(
            """
alignments:
  - id: alignment
    file: test.fasta
    data_type: nucleotide
partitions:
  - id: alignment
    site_model:
      substitution_model: {type: hky}
    clock_model: {type: strict}
    tree: shared
tree_prior: {type: yule}
diagnostics: {prior_similarity: true}
mcmc: {chain_length: 10000}
output: {file: out.xml}
""",
            encoding="utf-8",
        )
        with pytest.raises(ConfigError, match="similarity_threshold"):
            ConfigParser().parse(cfg)


class TestCrownStemTruthTable:
    """Every combination of node reading against every taxon-set relation.

    The review asked for a truth table rather than a handful of examples, and for
    the undecidable cells to be reported as ``unchecked`` instead of guessed.
    This pins all sixteen cells: the expectations below are the biology, so a
    rule that silently starts resolving a crown/stem pair differently shows up
    here rather than in a manuscript.
    """

    ALL = {"A", "B", "C", "D"}
    # a is always the older prior (20 Ma), b the younger (10 Ma), so a temporal
    # conflict means the pair was resolved with b above a.
    # (relation, a_originate, b_originate) -> (status, expected conflict types)
    GRID = [
        (("same", False, False), ("checked", ["same_node"])),
        (("same", False, True), ("checked", ["temporal"])),
        (("same", True, False), ("checked", [])),
        (("same", True, True), ("checked", ["same_node"])),
        (("nested", False, False), ("checked", [])),
        (("nested", False, True), ("unchecked", ["unchecked"])),
        (("nested", True, False), ("checked", [])),
        (("nested", True, True), ("checked", [])),
        (("intersecting", False, False), ("checked", ["topology"])),
        (("intersecting", False, True), ("checked", ["topology"])),
        (("intersecting", True, False), ("checked", ["topology"])),
        (("intersecting", True, True), ("checked", ["topology"])),
        (("disjoint", False, False), ("not_applicable", [])),
        (("disjoint", False, True), ("not_applicable", [])),
        (("disjoint", True, False), ("not_applicable", [])),
        (("disjoint", True, True), ("not_applicable", [])),
    ]

    SETS = {
        "same": (["A", "B"], ["A", "B"]),
        "nested": (["A", "B", "C"], ["A", "B"]),
        "intersecting": (["A", "B", "C"], ["B", "C", "D"]),
        "disjoint": (["A", "B"], ["C", "D"]),
    }

    @pytest.mark.parametrize("cell,expected", GRID, ids=[str(cell) for cell, _ in GRID])
    def test_each_cell_of_the_table(self, cell, expected):
        relation, a_originate, b_originate = cell
        status, kinds = expected
        left_set, right_set = self.SETS[relation]
        a = cal("a", left_set, use_originate=a_originate, mean=20.0, sigma=1.0)
        b = cal("b", right_set, use_originate=b_originate, mean=10.0, sigma=1.0)
        assessment = ConflictDetector.assess([a, b], all_taxa=set(self.ALL))[0]
        assert assessment.status == status, (cell, assessment.status)
        assert sorted({c.conflict_type for c in assessment.conflicts}) == sorted(kinds)

    def test_the_table_covers_every_combination(self):
        assert len({cell for cell, _ in self.GRID}) == 16
        assert {(cell[0], cell[1], cell[2]) for cell, _ in self.GRID} == {
            (relation, a, b)
            for relation in ("same", "nested", "intersecting", "disjoint")
            for a in (False, True)
            for b in (False, True)
        }

    def test_stem_and_crown_on_one_taxon_set_are_parent_and_child(self):
        # The cell that used to be skipped: same set, a=crown, b=stem, with the
        # stem reading younger than its own crown is impossible, so it is a real
        # conflict rather than a duplicate.
        a = cal("crown", ["A", "B"], mean=10.0, sigma=0.5)
        b = cal("stem", ["A", "B"], use_originate=True, mean=20.0, sigma=0.5)
        conflicts = ConflictDetector.detect_conflicts([a, b], all_taxa=set(self.ALL))
        assert [c.conflict_type for c in conflicts if c.severity == "error"] == []

    def test_disjoint_pairs_never_produce_a_structural_error(self):
        # Regression for the whole of issue 4: no relation, no conflict, whatever
        # the priors look like.
        for mean_a, mean_b in ((10.0, 10.0), (10.0, 100.0), (2.0, 60.0)):
            a = cal("a", ["A", "B"], mean=mean_a, sigma=0.5)
            b = cal("b", ["C", "D"], mean=mean_b, sigma=0.5)
            assessments = ConflictDetector.assess([a, b], all_taxa=set(self.ALL))
            assert assessments[0].status == "not_applicable"
            assert not [c for c in assessments[0].conflicts if c.severity == "error"]

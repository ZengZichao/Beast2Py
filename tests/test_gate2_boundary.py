"""Gate 2's stated boundary: the lists exist, and the outputs print them.

The review's complaint about the static gate was not that it was wrong but that
a green tick invites a reader to assume more than structure was checked.  These
tests pin the remedy: the boundary is a single source in ``diagnostics.py``, it is
rendered into the CLI output and the HTML report, and the documentation cannot
claim a detection the code does not implement.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from beast2py import diagnostics
from beast2py.diagnostics import (
    DETECTABLE_CHECKS,
    NOT_DETECTABLE_CHECKS,
    NOT_DETECTABLE_LABELS,
    PRIOR_LAYERS,
    gate2_scope_html,
    gate2_scope_lines,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DOC_FILES = sorted((REPO_ROOT / "docs").glob("*.md")) + [
    REPO_ROOT / "README.md",
    REPO_ROOT / "README_zh.md",
]


class TestBoundaryLists:
    def test_both_lists_are_non_empty_and_aligned(self):
        assert DETECTABLE_CHECKS, "the gate must state what it decides"
        assert len(NOT_DETECTABLE_CHECKS) == len(NOT_DETECTABLE_LABELS)
        # The short labels are the same items as the long sentences, so the one
        # line the CLI prints cannot drift from the list the report prints.
        for label, full in zip(NOT_DETECTABLE_LABELS, NOT_DETECTABLE_CHECKS):
            assert label.split(" ")[-1] in full or label.replace("-", " ") in full

    def test_the_four_prior_layers_are_named_in_order(self):
        names = [name for name, _ in PRIOR_LAYERS]
        assert names == [
            "specified marginal prior",
            "effective marginal prior",
            "joint effective prior",
            "posterior",
        ]

    def test_no_detectable_claim_is_a_posterior_claim(self):
        # Gate 2 is static; nothing it can detect may be phrased as a statement
        # about a posterior, or the boundary would be self-contradictory.
        for item in DETECTABLE_CHECKS:
            assert "posterior" not in item

    def test_the_list_names_the_four_things_a_static_gate_cannot_see(self):
        joined = " ".join(NOT_DETECTABLE_CHECKS).lower()
        for concept in ("effective marginal prior", "posterior dominance", "ess", "data"):
            assert concept in joined, concept


class TestRendering:
    def test_text_rendering_lists_both_halves(self):
        lines = gate2_scope_lines()
        text = "\n".join(lines)
        assert "Gate 2 can detect" in text
        assert "cannot detect" in text
        for item in DETECTABLE_CHECKS:
            assert item in text
        for item in NOT_DETECTABLE_CHECKS:
            assert item in text

    def test_html_rendering_escapes_and_includes_the_layers(self):
        html = gate2_scope_html()
        for name, gloss in PRIOR_LAYERS:
            assert name in html
        assert "<li>" in html
        assert "<script" not in html


class TestOutputsCarryTheBoundary:
    @pytest.fixture
    def config_with_calibrations(self, tmp_path):
        from .conftest import TEST_FASTA, TEST_MULTI_CAL_YAML

        (tmp_path / "test.fasta").write_text(TEST_FASTA, encoding="utf-8")
        config = tmp_path / "config.yaml"
        config.write_text(TEST_MULTI_CAL_YAML, encoding="utf-8")
        return config

    def test_generate_prints_the_boundary_with_the_passing_gate(
        self, config_with_calibrations, tmp_path, capsys
    ):
        from beast2py.main import main

        out = tmp_path / "gen.xml"
        assert main(["generate", "-c", str(config_with_calibrations), "-o", str(out)]) == 0
        printed = capsys.readouterr().out
        assert "Gate 2/3 passed" in printed
        assert "It cannot rule out:" in printed
        for label in NOT_DETECTABLE_LABELS:
            assert label in printed

    def test_diagnose_prints_the_full_lists_and_the_report_contains_them(
        self, config_with_calibrations, tmp_path, capsys
    ):
        from beast2py.main import main

        report = tmp_path / "report.html"
        ret = main(["diagnose", "-c", str(config_with_calibrations), "--report", str(report)])
        printed = capsys.readouterr().out
        assert ret == 0, printed
        assert "Gate 2 can detect (static structure only):" in printed
        for item in NOT_DETECTABLE_CHECKS:
            assert item in printed
        html = report.read_text(encoding="utf-8")
        assert "What this gate decides" in html
        for name, _ in PRIOR_LAYERS:
            assert name in html

    def test_the_boundary_lists_survive_html_escaping(self):
        # The items contain "95%" and "MCMC"; escaping must not mangle them into
        # something a reader cannot search for in the report.
        html = gate2_scope_html()
        assert "MCMC non-convergence" in html


class TestNoUncalibratedThresholdSurvives:
    def test_detector_holds_no_threshold_constant(self):
        # The old rules shipped 0.15 and 2.0 as scientific defaults.  Nothing may
        # quietly reintroduce a number that no reference set has calibrated.
        leftovers = [
            name for name in vars(diagnostics.ConflictDetector) if "THRESHOLD" in name.upper()
        ]
        assert leftovers == []
        for name in ("OVERLAP_RATIO_THRESHOLD", "OVERLAP_DISTANCE_THRESHOLD"):
            assert not hasattr(diagnostics, name), name
            assert not hasattr(diagnostics.ConflictDetector, name), name

    @pytest.mark.parametrize("path", DOC_FILES, ids=lambda p: p.name)
    def test_documentation_never_presents_a_default_threshold(self, path):
        text = path.read_text(encoding="utf-8")
        for match in re.finditer(r"similarity_threshold[^\n]*", text):
            line = match.group(0)
            assert "0.15" not in line, f"{path.name}: {line}"
            assert "default" not in line.lower() or "no calibrated default" in line.lower()


class TestDocumentationMatchesTheCode:
    """A manual that describes a detection the code lacks is a defect."""

    def test_no_manual_promises_a_mean_only_temporal_warning(self):
        # The detector reports only a wholly younger parent interval; a doc that
        # says "warning when only the means disagree" promises an alarm the tool
        # deliberately does not raise.
        offenders = []
        for path in DOC_FILES:
            text = path.read_text(encoding="utf-8")
            if "only the means disagree" in text or "只有均值冲突" in text:
                offenders.append(path.name)
        assert offenders == []

    def test_no_manual_calls_a_prior_only_run_a_posterior_comparison(self):
        offenders = []
        for path in DOC_FILES:
            text = path.read_text(encoding="utf-8")
            # "compare the tree-height posteriors" was the old wording for the
            # leave-one-out runs, which sample from the prior.
            if re.search(r"tree-height posteriors|树高的后验", text, re.I):
                offenders.append(path.name)
        assert offenders == []

    def test_no_manual_asserts_a_short_chain_is_enough(self):
        offenders = []
        for path in DOC_FILES:
            text = path.read_text(encoding="utf-8")
            # A prior-only run is not automatically cheap and not automatically
            # converged: no doc may assert a length is enough without an ESS
            # criterion.
            if re.search(
                r"(a short chain is enough|short chain suffices|短链即可"
                r"|1 million generations is enough)",
                text,
                re.I,
            ):
                offenders.append(path.name)
        assert offenders == []

    def test_both_manuals_describe_the_three_structural_rules(self):
        for path in (
            REPO_ROOT / "docs" / "user_manual_en.md",
            REPO_ROOT / "docs" / "user_manual_zh.md",
        ):
            text = path.read_text(encoding="utf-8")
            assert "same_node" in text or "同一节点" in text, path.name
            assert "topology" in text or "单系" in text, path.name

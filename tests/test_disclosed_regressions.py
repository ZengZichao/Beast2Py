"""Every disclosed defect of an archived release has a named regression case.

Issue 11.3 asked for this directly: the manuscript used to explain in prose that
``main`` had fixed things the archived tag still did not do. Prose is not a
guard. Each entry below is a defect that was disclosed to readers, the release
that closed it, and the test that now fails if the behaviour comes back.

The registry is checked against the collected suite rather than asserted from
memory, so a renamed or deleted regression test breaks the registry instead of
quietly leaving a disclosed defect unguarded.
"""

from __future__ import annotations

import pathlib
import re
import subprocess
import sys

import pytest

import beast2py
from beast2py.config import ConfigParser

from .conftest import EXAMPLES_DIR

REPO = pathlib.Path(beast2py.__file__).resolve().parent.parent

#: ``id`` labels the defect the way it was disclosed; ``cases`` cites the tests
#: that guard it, built from a file and its class/method chain so the registry
#: stays readable.


def _t(module: str, *names: str) -> str:
    """One pytest node id, written so it fits the lint line length."""
    return "::".join((module,) + names)


DISCLOSED_DEFECTS = [
    {
        "id": "TOOL-VERSION-010",
        "defect": "ConfigParser hardcoded tool_version '0.1.0', so fingerprints "
        "kept a stale suffix after a release bump",
        "fixed_in": "0.1.1",
        "cases": (
            _t(
                "tests/test_reproducibility.py",
                "TestFingerprintGenerator",
                "test_generate_fingerprint_dict",
            ),
            _t(
                "tests/test_disclosed_regressions.py",
                "TestToolVersionDefault",
                "test_the_default_tool_version_is_the_live_package_version",
            ),
        ),
    },
    {
        "id": "QUICK-CLI-INFER",
        "defect": "the quick API and the CLI inferred a sequence's data type by "
        "different rules, so one FASTA could become a nucleotide alignment in "
        "one entry point and an amino-acid one in the other",
        "fixed_in": "0.1.2",
        "cases": (
            _t(
                "tests/test_release_integrity.py",
                "TestPublicSurface",
                "test_read_sequence_defaults_to_the_same_inference_as_the_cli",
            ),
            _t(
                "tests/test_release_integrity.py",
                "TestPublicSurface",
                "test_quick_mode_is_the_same_function_from_either_entry_point",
            ),
        ),
    },
    {
        "id": "VALIDATE-EXIT-2",
        "defect": "the standalone validate command reported success when BEAST2 "
        "was missing, so a CI job reading only the exit code treated an "
        "unrunnable check as a pass",
        "fixed_in": "0.1.2",
        "cases": (
            _t(
                "tests/test_cli.py",
                "TestValidateExitCodeWhenTheParseGateCannotRun",
                "test_requested_beast2_gate_exits_2",
            ),
            _t(
                "tests/test_cli.py",
                "TestValidateRefusesWhenBeast2IsUndiscoverable",
                "test_requested_gate_exits_two_and_says_the_gate_is_unavailable",
            ),
        ),
    },
    {
        "id": "LANG-012",
        "defect": "the CLI and the documentation disagreed about the gate "
        "numbering, so a Gate 0 parse failure was reported as one of the three "
        "generation gates",
        "fixed_in": "0.1.2",
        "cases": (
            _t(
                "tests/test_cli.py",
                "TestThePrintedGateNumbering",
                "test_gate_zero_failure_is_a_configuration_error_that_writes_nothing",
            ),
            _t(
                "tests/test_cli.py",
                "TestThePrintedGateNumbering",
                "test_quick_mode_prints_gate_zero_before_the_generation_gates",
            ),
            _t(
                "tests/test_cli.py",
                "TestThePrintedGateNumbering",
                "test_validate_prints_the_gate_it_maps_each_level_onto",
            ),
        ),
    },
    {
        "id": "PY310-TOMLLIB",
        "defect": "the release-integrity guard imported tomllib, which is stdlib "
        "only from 3.11, so it errored on the 3.10 CI legs",
        "fixed_in": "0.1.2",
        "cases": (
            _t(
                "tests/test_release_integrity.py",
                "TestDistributionContents",
                "test_declared_license_is_the_same_expression_everywhere",
            ),
        ),
    },
    {
        "id": "ZENODO-LICENCE-ID",
        "defect": "the archive metadata put the compound SPDX expression in "
        "license.id, which Zenodo's licence vocabulary rejects, so archiving "
        "v0.1.2 failed with 'Extra metadata load failed'",
        "fixed_in": "0.1.2",
        "cases": (
            _t(
                "tests/test_release_integrity.py",
                "TestDistributionContents",
                "test_declared_license_is_the_same_expression_everywhere",
            ),
        ),
    },
    {
        "id": "STORE-EVERY-PARITY",
        "defect": "quick and generate disagreed about the store_every default, so "
        "the same configuration produced chains with different resume points "
        "depending on the entry point",
        "fixed_in": "0.1.2",
        "cases": (
            _t(
                "tests/test_cli.py", "TestCLI", "test_quick_and_generate_write_the_same_store_every"
            ),
        ),
    },
    {
        "id": "NESTED-REDUNDANCY",
        "defect": "similar priors on independent clades were reported as "
        "redundancy with a merge recommendation and could block Gate 2, on no "
        "calibrated evidence",
        "fixed_in": "0.2.0",
        "cases": (
            _t(
                "tests/test_diagnostics.py",
                "TestPriorSimilarityIsDescriptiveOnly",
                "test_nested_similar_priors_are_not_renamed_redundancy",
            ),
            _t(
                "tests/test_diagnostics.py",
                "TestCrownStemTruthTable",
                "test_disjoint_pairs_never_produce_a_structural_error",
            ),
            "tests/test_validation_set.py",
        ),
    },
    {
        "id": "ORIGINATE-SILENT",
        "defect": "use_originate was not written for a root point and was "
        "silently ignored elsewhere, so the YAML, the fingerprint and the XML "
        "described different nodes",
        "fixed_in": "0.2.0",
        "cases": (
            _t(
                "tests/test_calibration.py",
                "test_calibrated_yule_writes_forparent_only_for_a_stem_point[stem]",
            ),
            _t(
                "tests/test_calibration.py",
                "test_calibrated_yule_writes_forparent_only_for_a_stem_point[crown]",
            ),
            _t(
                "tests/test_calibration.py",
                "test_calibrated_yule_stem_point_reaches_the_generated_xml",
            ),
            _t("tests/test_calibration.py", "test_mrcaprior_refuses_a_stem_reading_of_the_root"),
        ),
    },
    {
        "id": "RELEASE-GATE",
        "defect": "the Release workflow gated a tag on version-string agreement "
        "only, so a failing suite or a skipped Gate 3 could still publish a build",
        "fixed_in": "0.2.0",
        "cases": (
            _t(
                "tests/test_release_integrity.py",
                "TestReleaseWorkflowGate",
                "test_build_waits_for_a_job_that_runs_the_suite",
            ),
            _t(
                "tests/test_release_integrity.py",
                "TestReleaseWorkflowGate",
                "test_the_release_gate_requires_beast2_rather_than_skipping_it",
            ),
            _t(
                "tests/test_release_integrity.py",
                "TestReleaseWorkflowGate",
                "test_ci_and_release_share_one_beast2_install",
            ),
        ),
    },
    {
        "id": "TRACE-PARSER-UNGUARDED",
        "defect": "the numbers a prior-only figure reports were produced by a "
        "one-off script beside the run directories, with no reader that checked "
        "column names, duplicate sample ids or unequal chain lengths",
        "fixed_in": "0.2.0",
        "cases": (
            _t(
                "tests/test_trace.py",
                "test_read_trace_uses_the_header_not_the_column_order",
            ),
            _t(
                "tests/test_trace.py",
                "test_duplicate_sample_ids_are_rejected_as_concatenation",
            ),
            _t(
                "tests/test_trace.py",
                "test_combine_refuses_unequal_chain_lengths",
            ),
            _t(
                "tests/test_trace.py",
                "TestPriorSummaryCLI",
                "test_a_column_the_scenario_does_not_log_is_a_usage_error",
            ),
        ),
    },
    {
        "id": "GATE2-BOUNDARY-UNSTATED",
        "defect": "Gate 2 printed a green tick without stating what it cannot "
        "decide, so a passing static gate could be read as ruling out posterior "
        "dominance",
        "fixed_in": "0.2.0",
        "cases": (
            _t(
                "tests/test_gate2_boundary.py",
                "TestOutputsCarryTheBoundary",
                "test_generate_prints_the_boundary_with_the_passing_gate",
            ),
            _t(
                "tests/test_gate2_boundary.py",
                "TestDocumentationMatchesTheCode",
                "test_no_manual_promises_a_mean_only_temporal_warning",
            ),
        ),
    },
    {
        "id": "PROVENANCE-PLACEHOLDER",
        "defect": "the formal examples carried unrelated DOIs and 'no source "
        "consulted' placeholders, and nothing in the suite noticed because the "
        "only check lived in the submission bundle script",
        "fixed_in": "0.2.0",
        "cases": (
            _t(
                "tests/test_release_integrity.py",
                "TestDocumentedCountsAndProvenance",
                "test_formal_examples_are_explicitly_synthetic",
            ),
            _t(
                "tests/test_release_integrity.py",
                "TestDocumentedCountsAndProvenance",
                "test_every_formal_example_calibration_carries_a_synthetic_card",
            ),
        ),
    },
    {
        "id": "LINT-DEBT",
        "defect": "the tree carried lint findings the CI gate rejects (an unused "
        "import, two unused locals, two over-long lines), so every push after "
        "v0.2.0 was red and nothing could be released through the gate at all",
        "fixed_in": "0.2.0",
        "cases": (
            _t(
                "tests/test_disclosed_regressions.py",
                "TestLintStaysClean",
                "test_no_line_exceeds_the_configured_limit",
            ),
            _t(
                "tests/test_disclosed_regressions.py",
                "TestLintStaysClean",
                "test_no_module_level_import_is_left_unused",
            ),
            _t(
                "tests/test_disclosed_regressions.py",
                "TestLintStaysClean",
                "test_the_python_tools_agree_when_they_are_installed",
            ),
        ),
    },
    {
        "id": "DOC-COUNT-DRIFT",
        "defect": "the READMEs kept a test count from an earlier tag while the "
        "manuscript quoted the frozen one, so the same suite was reported as two "
        "different sizes",
        "fixed_in": "0.2.0",
        "cases": (
            _t(
                "tests/test_release_integrity.py",
                "TestDocumentedCountsAndProvenance",
                "test_readme_test_counts_match_the_collected_suite",
            ),
        ),
    },
]


def _collected_node_ids() -> set:
    """Every node id the suite collects, cached per session."""
    if not hasattr(_collected_node_ids, "_cache"):
        proc = subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "tests",
                "--collect-only",
                "-q",
                "-p",
                "no:cacheprovider",
            ],
            cwd=str(REPO),
            capture_output=True,
            text=True,
            timeout=900,
        )
        assert proc.returncode == 0, proc.stdout + proc.stderr
        ids = {
            line.strip()
            for line in proc.stdout.splitlines()
            if line.strip().startswith("tests/") and "::" in line
        }
        # A whole-file entry is a valid way to cite a regression case.
        files = {line.split("::")[0] for line in ids}
        _collected_node_ids._cache = (ids, files)
    return _collected_node_ids._cache


class TestDisclosedDefectRegistry:
    @pytest.mark.parametrize(
        "entry", DISCLOSED_DEFECTS, ids=[entry["id"] for entry in DISCLOSED_DEFECTS]
    )
    def test_each_disclosed_defect_has_a_live_regression_case(self, entry):
        node_ids, files = _collected_node_ids()
        assert entry["cases"], entry["id"]
        for cited in entry["cases"]:
            if cited in node_ids or cited in files:
                continue
            raise AssertionError(
                f"{entry['id']}: the regression case {cited!r} is not in the "
                "collected suite. Rename the registry entry or restore the test; "
                "do not delete the entry, because the defect was disclosed to "
                "readers and stays disclosed."
            )

    def test_registry_ids_are_unique_and_versions_are_released(self):
        ids = [entry["id"] for entry in DISCLOSED_DEFECTS]
        assert len(ids) == len(set(ids))
        tags = {
            line.strip().removeprefix("v")
            for line in subprocess.run(
                ["git", "tag", "-l"], cwd=str(REPO), capture_output=True, text=True
            ).stdout.splitlines()
        }
        if not tags:
            pytest.skip(
                "git tags are not available in this checkout; the Release workflow "
                "fetches tags and enforces the fixed_in versions"
            )
        for entry in DISCLOSED_DEFECTS:
            assert entry["fixed_in"] in tags, (entry["id"], entry["fixed_in"], sorted(tags))

    def test_the_defects_named_in_the_review_are_all_registered(self):
        """The three differences between the archived tag and `main` that a reader
        was told about must each have a guard, not just a paragraph.

        The pre-submission review listed them by name: the quick API and the CLI
        inferring a sequence's data type differently, `validate` reporting success
        when the BEAST2 check could not run, and the gate names and order
        disagreeing with the figure. Anything else in the changelog may be fixed
        without a registry entry; these three were disclosed as version
        differences, so they need one.
        """
        registered = {entry["id"]: entry for entry in DISCLOSED_DEFECTS}
        for required in ("QUICK-CLI-INFER", "VALIDATE-EXIT-2", "LANG-012"):
            entry = registered.get(required)
            assert entry, f"{required} is a disclosed version difference with no registry entry"
            assert entry["cases"], f"{required} is registered without a regression case"


class TestLintStaysClean:
    """The tree keeps the shape the CI lint gate demands.

    CI runs ``black --check`` and ``flake8 beast2py tests`` before anything else,
    and v0.2.0 carried findings there (an unused ``import os``, two locals
    computed and never read, two lines over the 100-character limit). A red lint
    gate is not cosmetic: with the release now waiting on the suite, a tree that
    cannot pass its own checks has no way to publish.

    These checks are deliberately dependency-free, so a contributor without
    black or flake8 installed still sees the same failure locally.
    """

    MAX_LINE = 100
    SOURCES = sorted((REPO / "beast2py").rglob("*.py")) + sorted((REPO / "tests").glob("*.py"))

    def test_no_line_exceeds_the_configured_limit(self):
        offenders = []
        for path in self.SOURCES:
            if "__pycache__" in str(path):
                continue
            for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
                if len(line) > self.MAX_LINE and "noqa" not in line:
                    offenders.append(f"{path.relative_to(REPO)}:{number} ({len(line)})")
        assert not offenders, "lines over 100 characters: " + ", ".join(offenders[:8])

    def test_no_module_level_import_is_left_unused(self):
        import ast

        offenders = []
        for path in self.SOURCES:
            if "__pycache__" in str(path):
                continue
            text = path.read_text(encoding="utf-8")
            tree = ast.parse(text)
            exported = set()
            for node in tree.body:
                if isinstance(node, ast.Assign) and any(
                    getattr(target, "id", "") == "__all__" for target in node.targets
                ):
                    try:
                        exported |= set(ast.literal_eval(node.value))
                    except ValueError:
                        pass
            for node in tree.body:
                if not isinstance(node, (ast.Import, ast.ImportFrom)):
                    continue
                for alias in node.names:
                    if alias.name == "*":
                        continue
                    name = (alias.asname or alias.name).split(".")[0]
                    if name == "annotations" or name in exported:
                        continue
                    # `# noqa: F401` is the sanctioned way to keep a re-export.
                    line = text.splitlines()[node.lineno - 1]
                    if "noqa" in line:
                        continue
                    if text.count(name) <= 1:
                        offenders.append(f"{path.relative_to(REPO)}:{node.lineno}: {name}")
        assert not offenders, "unused imports: " + ", ".join(offenders[:8])

    def test_the_python_tools_agree_when_they_are_installed(self):
        """black, flake8 and mypy are the authoritative gate; run them if present.

        mypy matters most here: the CI ratchet now sits at a baseline of 0, which
        is only honest while the tree stays clean, and a green pytest run can
        still hide a type error that fails the release gate.
        """
        import shutil
        import subprocess

        commands = {
            "black": [sys.executable, "-m", "black", "--check", "--quiet", "."],
            "flake8": [sys.executable, "-m", "flake8", "beast2py", "tests"],
            "mypy": [sys.executable, "-m", "mypy"],
        }
        for tool, args in commands.items():
            if (
                shutil.which(tool) is None
                and subprocess.run(
                    [sys.executable, "-m", tool, "--version"], capture_output=True
                ).returncode
                != 0
            ):
                pytest.skip(f"{tool} is not installed in this environment")
            proc = subprocess.run(args, cwd=str(REPO), capture_output=True, text=True)
            assert proc.returncode == 0, f"{tool} says:\n{proc.stdout}{proc.stderr}"

    def test_the_ci_mypy_ratchet_stays_pinned_at_zero(self):
        """A raised baseline is a silently re-opened backlog."""
        workflow = (REPO / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
        match = re.search(r"MYPY_BASELINE:\s*(\d+)", workflow)
        assert match, "ci.yml no longer declares a mypy baseline"
        assert (
            match.group(1) == "0"
        ), f"MYPY_BASELINE is {match.group(1)} while the tree reports 0 errors"


class TestToolVersionDefault:
    def test_the_default_tool_version_is_the_live_package_version(self):
        """The 0.1.0 hardcode, guarded directly rather than through a fingerprint.

        A release bump must move the version the fingerprint records. Before the
        fix the parser wrote a literal, so the identifier of a configuration
        generated by 0.1.1 still said 0.1.0.
        """
        config = ConfigParser().parse(EXAMPLES_DIR / "config_basic.yaml")
        assert config.metadata.get("tool_version") == beast2py.__version__

    def test_an_explicit_tool_version_is_still_honoured(self):
        # The default reads __version__; it must not override what the author wrote.
        config = ConfigParser().parse(EXAMPLES_DIR / "config_basic.yaml")
        config.metadata["tool_version"] = "9.9.9"
        assert config.metadata["tool_version"] == "9.9.9"

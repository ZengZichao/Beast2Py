"""Tests for the CLI main module."""

import os
import shutil
import tempfile
import xml.etree.ElementTree as ET
import pytest
from pathlib import Path

from beast2py.main import build_parser, main

from .conftest import EXAMPLES_DIR, TEST_FASTA, TEST_CONFIG_YAML


def _store_every_values(path):
    """Return every ``storeEvery`` attribute in an emitted XML, as ints.

    Reads them all rather than the first: the MCMC run element and the state
    element must agree, and a default that only reached one of the two is the
    kind of split this test exists to catch.
    """
    values = []
    for elem in ET.parse(path).iter():
        raw = elem.get("storeEvery")
        if raw is not None:
            values.append(int(raw))
    return values


class TestCLI:
    """Test the CLI interface."""

    def setup_method(self, method):
        """Create a temp directory for each test."""
        self.tmpdir = Path(tempfile.mkdtemp())
        fasta_path = self.tmpdir / "test.fasta"
        fasta_path.write_text(TEST_FASTA)

    def test_parser_has_subcommands(self):
        """Test that the parser has all expected subcommands."""
        parser = build_parser()
        actions = [a for a in parser._actions if hasattr(a, "choices")]
        subcommand_choices = None
        for action in actions:
            if action.dest == "command":
                subcommand_choices = action.choices
                break
        assert subcommand_choices is not None
        assert "generate" in subcommand_choices
        assert "diagnose" in subcommand_choices
        assert "quick" in subcommand_choices
        assert "validate" in subcommand_choices
        assert "methods" in subcommand_choices
        assert "pipeline" in subcommand_choices
        assert "fingerprint" in subcommand_choices
        assert "list-models" in subcommand_choices

    def test_no_command_prints_help(self):
        """Test that running with no command prints help."""
        ret = main([])
        assert ret == 0

    def test_version_flag(self, capsys):
        """Test the --version flag."""
        with pytest.raises(SystemExit):
            main(["--version"])
        captured = capsys.readouterr()
        # Track the single source of truth instead of a literal that rots on
        # every release bump.
        from beast2py import __version__

        assert f"v{__version__}" in captured.out

    def test_generate_command(self):
        """Test the generate subcommand."""
        # Create config file
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        output_path = self.tmpdir / "output.xml"
        ret = main(["generate", "-c", str(config_path), "-o", str(output_path), "-v"])
        assert ret == 0
        assert os.path.exists(output_path)

        with open(output_path, "r") as f:
            xml = f.read()
        assert "<beast" in xml
        assert "</beast>" in xml

    def test_generate_with_fingerprint(self):
        """Test the generate command with --fingerprint flag."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        output_path = self.tmpdir / "output.xml"
        ret = main(
            [
                "generate",
                "-c",
                str(config_path),
                "-o",
                str(output_path),
                "--fingerprint",
                "-v",
            ]
        )
        assert ret == 0

        fp_path = str(Path(output_path).with_suffix(".fingerprint.json"))
        assert os.path.exists(fp_path)

    def test_generate_with_methods(self):
        """Test the generate command with --methods flag."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        output_path = self.tmpdir / "output.xml"
        ret = main(
            [
                "generate",
                "-c",
                str(config_path),
                "-o",
                str(output_path),
                "--methods",
                "-v",
            ]
        )
        assert ret == 0

        methods_path = str(Path(output_path).with_suffix(".methods.tex"))
        assert os.path.exists(methods_path)

    def test_validate_command(self):
        """Test the validate subcommand."""
        # First generate XML
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        output_path = self.tmpdir / "output.xml"
        main(["generate", "-c", str(config_path), "-o", str(output_path)])

        # Then validate it
        ret = main(["validate", "--xml", str(output_path)])
        assert ret == 0

    def test_methods_command(self):
        """Test the methods subcommand."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        ret = main(["methods", "-c", str(config_path)])
        assert ret == 0

    def test_methods_command_with_output(self):
        """Test the methods subcommand with output file."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        methods_path = self.tmpdir / "methods.tex"
        ret = main(["methods", "-c", str(config_path), "-o", str(methods_path)])
        assert ret == 0
        assert os.path.exists(methods_path)

    def test_fingerprint_command(self):
        """Test the fingerprint subcommand."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        ret = main(["fingerprint", "-c", str(config_path), "-v"])
        assert ret == 0

    def test_list_models_command(self):
        """Test the list-models subcommand."""
        ret = main(["list-models"])
        assert ret == 0

    def test_pipeline_snakemake_command(self):
        """Test the pipeline subcommand with snakemake."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        output_dir = self.tmpdir / "pipeline"
        ret = main(
            [
                "pipeline",
                "-c",
                str(config_path),
                "--type",
                "snakemake",
                "--output-dir",
                str(output_dir),
            ]
        )
        assert ret == 0
        assert (Path(output_dir) / "Snakefile").exists()

    def test_pipeline_threads_reaches_the_emitted_files(self):
        """CODE-008: `--threads` is a knob, not the old welded-in `-threads 4`.

        Both engines are checked, plus the rejection path, because the point of
        the knob is a 2-core laptop or a resource-restricted HPC partition.
        """
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        default_dir = self.tmpdir / "pipe_default"
        assert (
            main(
                [
                    "pipeline",
                    "-c",
                    str(config_path),
                    "--type",
                    "both",
                    "--output-dir",
                    str(default_dir),
                ]
            )
            == 0
        )
        for name in ("Snakefile", "main.nf"):
            text = (default_dir / name).read_text()
            assert "-threads 4" in text, name
            assert "-threads 1" not in text, name

        one_dir = self.tmpdir / "pipe_one_thread"
        assert (
            main(
                [
                    "pipeline",
                    "-c",
                    str(config_path),
                    "--type",
                    "both",
                    "--output-dir",
                    str(one_dir),
                    "--threads",
                    "1",
                ]
            )
            == 0
        )
        for name in ("Snakefile", "main.nf"):
            text = (one_dir / name).read_text()
            assert "-threads 1" in text, name
            assert "-threads 4" not in text, name

        bad_dir = self.tmpdir / "pipe_bad"
        assert (
            main(
                [
                    "pipeline",
                    "-c",
                    str(config_path),
                    "--output-dir",
                    str(bad_dir),
                    "--threads",
                    "0",
                ]
            )
            == 1
        )
        assert not (bad_dir / "Snakefile").exists(), "a rejected thread count wrote a file"

    def test_quick_command(self):
        """Test the quick subcommand."""
        fasta_path = self.tmpdir / "primates.fasta"
        fasta_path.write_text(TEST_FASTA)

        cal_path = self.tmpdir / "cals.yaml"
        cal_path.write_text("""
- name: rootCal
  taxa: null
  monophyletic: true
  distribution:
    type: normal
    parameters: {mean: 10.0, sigma: 1.0}
""")

        output_path = self.tmpdir / "quick_output.xml"
        ret = main(
            [
                "quick",
                "-a",
                str(fasta_path),
                "-o",
                str(output_path),
                "--tree-prior",
                "yule",
                "--subst-model",
                "hky",
                "--clock-model",
                "strict",
                "--calibration-yaml",
                str(cal_path),
                "--chain-length",
                "100000",
                "-v",
            ]
        )
        assert ret == 0
        assert os.path.exists(output_path)

    def test_quick_and_generate_write_the_same_store_every(self):
        """CODE-007: one documented ``store_every`` default for every path.

        ``quick`` used to build ``MCMCConfig`` from a dataclass default of
        10000, while ``generate`` derived it from the chain length as
        ``chain_length // 10000``. Both XMLs passed every gate, but the two
        chains checkpointed at different rates, so the entry point silently
        changed the sampling/resume behaviour of an otherwise identical run.
        """
        from beast2py.api import Beast2Py
        from beast2py.models import MCMCConfig, default_store_every

        chain_length = 1_000_000
        # The documented rule, spelled out so a change here is a deliberate one.
        expected = max(1, chain_length // 10_000)
        assert expected == default_store_every(chain_length) == 100

        fasta = self.tmpdir / "primates.fasta"
        fasta.write_text(TEST_FASTA)

        quick_xml = self.tmpdir / "quick.xml"
        assert (
            main(
                [
                    "quick",
                    "-a",
                    str(fasta),
                    "-o",
                    str(quick_xml),
                    "--tree-prior",
                    "yule",
                    "--chain-length",
                    str(chain_length),
                ]
            )
            == 0
        )

        config_path = self.tmpdir / "gen_config.yaml"
        config_path.write_text(
            TEST_CONFIG_YAML.replace("chain_length: 100000", f"chain_length: {chain_length}")
        )
        generated_xml = self.tmpdir / "generated.xml"
        assert main(["generate", "-c", str(config_path), "-o", str(generated_xml)]) == 0

        api_xml = self.tmpdir / "api_quick.xml"
        Beast2Py().quick_generate(
            alignment=str(fasta),
            output=str(api_xml),
            tree_prior="yule",
            chain_length=chain_length,
        )

        for label, path in (
            ("quick", quick_xml),
            ("generate", generated_xml),
            ("api", api_xml),
        ):
            values = _store_every_values(path)
            # Every storeEvery in the document, not just the first: a partial
            # writer that dropped the scaling rule would otherwise slip through.
            assert values, f"{label} emitted no storeEvery attribute"
            assert set(values) == {expected}, (label, values, expected)

        # The dataclass default is now the same rule, not a fixed 10000.
        assert MCMCConfig(chain_length=chain_length).store_every == expected

    def test_store_every_floor_and_explicit_value_still_win(self):
        """The documented default is a floor of 1 and never overrides a user value."""
        from beast2py.models import default_store_every

        assert default_store_every(9_999) == 1  # short chain, no storeEvery="0"
        assert default_store_every(100_000_000) == 10_000

        fasta = self.tmpdir / "primates.fasta"
        fasta.write_text(TEST_FASTA)
        config_path = self.tmpdir / "explicit.yaml"
        config_path.write_text(
            TEST_CONFIG_YAML.replace("  pre_burnin: 0\n", "  pre_burnin: 0\n  store_every: 7\n")
        )
        out = self.tmpdir / "explicit.xml"
        assert main(["generate", "-c", str(config_path), "-o", str(out)]) == 0
        assert set(_store_every_values(out)) == {7}

    def test_generate_nonexistent_config(self):
        """Test that nonexistent config file returns error."""
        ret = main(["generate", "-c", "/nonexistent/config.yaml", "-o", "out.xml"])
        assert ret == 1


class TestUserYamlErrorReporting:
    """A hand-edited YAML file is the likeliest thing a new user breaks.

    Every such mistake must reach the terminal as one ``[ERROR]`` line that
    names the file, not as a Python traceback: the traceback path in
    :func:`main` is reserved for bugs, and a stack trace on a typo teaches
    the reader nothing about their own file.
    """

    def setup_method(self, method):
        self.tmpdir = Path(tempfile.mkdtemp())
        (self.tmpdir / "test.fasta").write_text(TEST_FASTA)

    def _run(self, capsys, argv):
        ret = main(argv)
        captured = capsys.readouterr()
        return ret, captured.out + captured.err

    def _assert_clean_error(self, capsys, argv, needle):
        ret, text = self._run(capsys, argv)
        assert ret == 1, text
        assert "Traceback" not in text, text
        assert "[ERROR]" in text, text
        assert needle in text, text

    def test_malformed_config_yaml(self, capsys):
        bad = self.tmpdir / "bad.yaml"
        bad.write_text("partitions:\n  - id: alignment\n\t broken: :\n")
        self._assert_clean_error(
            capsys,
            ["generate", "-c", str(bad), "-o", str(self.tmpdir / "out.xml")],
            "Invalid YAML in",
        )

    def test_malformed_calibration_yaml(self, capsys):
        bad = self.tmpdir / "cals.yaml"
        bad.write_text("- name: rootCal\n   taxa: [oops\n")
        self._assert_clean_error(
            capsys,
            [
                "quick",
                "-a",
                str(self.tmpdir / "test.fasta"),
                "-o",
                str(self.tmpdir / "out.xml"),
                "--calibration-yaml",
                str(bad),
            ],
            "Invalid YAML in",
        )

    def test_missing_calibration_file(self, capsys):
        self._assert_clean_error(
            capsys,
            [
                "quick",
                "-a",
                str(self.tmpdir / "test.fasta"),
                "-o",
                str(self.tmpdir / "out.xml"),
                "--calibration-yaml",
                str(self.tmpdir / "nope.yaml"),
            ],
            "Calibration file not found",
        )

    def test_calibration_file_that_is_not_a_list(self, capsys):
        mapping = self.tmpdir / "cals.yaml"
        mapping.write_text("name: rootCal\n")
        self._assert_clean_error(
            capsys,
            [
                "quick",
                "-a",
                str(self.tmpdir / "test.fasta"),
                "-o",
                str(self.tmpdir / "out.xml"),
                "--calibration-yaml",
                str(mapping),
            ],
            "must be a list",
        )

    def test_calibration_that_names_an_unknown_taxon(self, capsys):
        cal = self.tmpdir / "cals.yaml"
        cal.write_text(
            "- name: ghostCal\n"
            "  taxa: [ghostA, ghostB]\n"
            "  monophyletic: true\n"
            "  distribution:\n"
            "    type: normal\n"
            "    parameters: {mean: 5.0, sigma: 0.5}\n"
        )
        self._assert_clean_error(
            capsys,
            [
                "quick",
                "-a",
                str(self.tmpdir / "test.fasta"),
                "-o",
                str(self.tmpdir / "out.xml"),
                "--calibration-yaml",
                str(cal),
            ],
            "not found in any alignment",
        )

    def test_valid_calibration_still_succeeds(self, capsys):
        """The guards must not swallow the working path."""
        cal = self.tmpdir / "cals.yaml"
        cal.write_text(
            "- name: rootCal\n"
            "  taxa: null\n"
            "  monophyletic: true\n"
            "  distribution:\n"
            "    type: normal\n"
            "    parameters: {mean: 10.0, sigma: 1.0}\n"
        )
        out = self.tmpdir / "quick.xml"
        ret, text = self._run(
            capsys,
            [
                "quick",
                "-a",
                str(self.tmpdir / "test.fasta"),
                "-o",
                str(out),
                "--calibration-yaml",
                str(cal),
            ],
        )
        assert ret == 0, text
        assert out.exists()


class TestValidateExitCodeWhenTheParseGateCannotRun:
    """``validate --beast2`` must not go green on a check that never happened.

    The archived release printed a JVM "could not find or load main class"
    message and exited 1 for a broken installation, which is exactly the code it
    uses for "BEAST2 rejected this model" -- so a missing checker and a rejected
    analysis were indistinguishable, and before that the same case warned and
    exited 0. The gate now refuses with its own code, 2, and an English
    tool-level sentence.

    ``validate`` with no ``--beast2`` still exits 0 on a structural pass, and that
    is deliberate: generate gates a write, validate reports what it checked.
    """

    def setup_method(self, method):
        """Write a structurally valid XML file for the command to read."""
        self.tmpdir = Path(tempfile.mkdtemp())
        self.xml = self.tmpdir / "sample.xml"
        shutil.copy(EXAMPLES_DIR / "output_basic.xml", self.xml)

    def _run_without_validator(self, monkeypatch, argv):
        """Run *argv* with the bundled launcher pointed at one lacking its helper."""
        import beast2py.validator as validator_module

        hollow = self.tmpdir / "beast2_validate.sh"
        hollow.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
        monkeypatch.setattr(
            validator_module.XMLValidator,
            "_bundled_validator_script",
            staticmethod(lambda: hollow),
        )
        return main(argv)

    def test_requested_beast2_gate_exits_2(self, monkeypatch, capsys):
        """Asking for BEAST2 and not getting it is a refusal, not a pass or a 1."""
        ret = self._run_without_validator(
            monkeypatch, ["validate", "--xml", str(self.xml), "--beast2"]
        )
        assert ret == 2
        err = capsys.readouterr().err
        assert "Level 2/2 not run" in err
        assert "NOT been verified by BEAST2" in err

    def test_the_failure_carries_no_jvm_crash_text(self, monkeypatch, capsys):
        """The report is the tool's own sentence, not a Java stack trace."""
        self._run_without_validator(monkeypatch, ["validate", "--xml", str(self.xml), "--beast2"])
        captured = capsys.readouterr()
        text = captured.out + captured.err
        for fragment in ("ClassNotFound", "main class", "\u9519\u8bef", "\u539f\u56e0"):
            assert fragment not in text, fragment

    def test_the_exemption_is_explicit_and_named_like_generate(self, monkeypatch, capsys):
        """--allow-unvalidated is the opt-in downgrade, matching the generate flag."""
        ret = self._run_without_validator(
            monkeypatch,
            ["validate", "--xml", str(self.xml), "--beast2", "--allow-unvalidated"],
        )
        assert ret == 0
        assert "Level 2/2 not run" in capsys.readouterr().err

    def test_structural_only_validation_still_exits_zero(self, capsys):
        """The documented asymmetry: no --beast2 means no gate was asked for."""
        assert main(["validate", "--xml", str(self.xml)]) == 0
        out = capsys.readouterr().out
        assert "Level 1/2 passed" in out
        assert "Level 2/2" not in out

    def test_structural_only_validation_still_exits_one_on_a_bad_file(self, capsys):
        """Refusing to claim a BEAST2 verdict is not the same as passing anything."""
        broken = self.tmpdir / "broken.xml"
        broken.write_text('<?xml version="1.0"?><foo/>', encoding="utf-8")
        assert main(["validate", "--xml", str(broken)]) == 1

    def test_a_real_rejection_is_still_reported_as_one(self, monkeypatch, capsys):
        """Non-vacuity: 2 means "could not run"; a rejected model must stay 1."""
        import beast2py.validator as validator_module

        monkeypatch.setattr(
            validator_module.XMLValidator,
            "beast2_validation_status",
            staticmethod(lambda *a, **k: (False, "INVALID: sample.xml  bad model", "failed")),
        )
        assert main(["validate", "--xml", str(self.xml), "--beast2"]) == 1
        assert "BEAST2 rejected the model" in capsys.readouterr().err

    def test_a_real_pass_is_still_reported_as_zero(self, monkeypatch, capsys):
        """And the happy path is untouched: a parsed model exits 0."""
        import beast2py.validator as validator_module

        monkeypatch.setattr(
            validator_module.XMLValidator,
            "beast2_validation_status",
            staticmethod(lambda *a, **k: (True, "VALID: sample.xml", "passed")),
        )
        assert main(["validate", "--xml", str(self.xml), "--beast2"]) == 0
        assert "Level 2/2 passed" in capsys.readouterr().out

"""Tests for the CLI main module."""

import os
import tempfile
import pytest
from pathlib import Path

from beast2py.main import build_parser, main

from .conftest import TEST_FASTA, TEST_CONFIG_YAML


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
        assert "v0.1.0" in captured.out

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

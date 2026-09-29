"""Tests for the reproducibility module."""

import os
import re
import shutil
import tempfile
from pathlib import Path

from beast2py.config import ConfigParser
from beast2py.reproducibility import (
    FingerprintGenerator,
    MethodsGenerator,
    PipelineGenerator,
)

from .conftest import TEST_FASTA, TEST_CONFIG_YAML


class TestFingerprintGenerator:
    """Test the FingerprintGenerator class."""

    def setup_method(self, method):
        """Create a temp directory for each test."""
        self.tmpdir = Path(tempfile.mkdtemp())
        fasta_path = self.tmpdir / "test.fasta"
        fasta_path.write_text(TEST_FASTA)

    def test_generate_fingerprint(self):
        """Test fingerprint generation."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        fingerprint = FingerprintGenerator.generate_fingerprint(config)
        assert fingerprint.startswith("B2P-")
        parts = fingerprint.split("-")
        assert len(parts) >= 3  # BNC, hash, version, date

    def test_generate_fingerprint_dict(self):
        """Test detailed fingerprint dictionary."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        fp_dict = FingerprintGenerator.generate_fingerprint_dict(config)
        assert "fingerprint" in fp_dict
        assert "full_hash" in fp_dict
        assert "tool_version" in fp_dict
        assert "beast2_version" in fp_dict
        assert "generation_time" in fp_dict

    def test_fingerprint_consistency(self):
        """Test that the same config produces the same fingerprint hash."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        fp1 = FingerprintGenerator.generate_fingerprint(config)
        fp2 = FingerprintGenerator.generate_fingerprint(config)

        # The identifier is a pure function of the analysis: no calendar date,
        # so two runs of the same configuration are equal in full, not merely
        # in their hash part ( /.
        assert fp1 == fp2
        assert len(fp1.split("-")) == 3, f"unexpected date component: {fp1}"

    def test_write_fingerprint_file(self):
        """Test writing fingerprint to a JSON file."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        fp_path = self.tmpdir / "fingerprint.json"
        FingerprintGenerator.write_fingerprint_file(config, fp_path)

        assert os.path.exists(fp_path)
        import json

        with open(fp_path, "r") as f:
            data = json.load(f)
        assert "fingerprint" in data


class TestMethodsGenerator:
    """Test the MethodsGenerator class."""

    def setup_method(self, method):
        """Create a temp directory for each test."""
        self.tmpdir = Path(tempfile.mkdtemp())
        fasta_path = self.tmpdir / "test.fasta"
        fasta_path.write_text(TEST_FASTA)

    def test_generate_methods(self):
        """Test methods description generation."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        methods = MethodsGenerator.generate_methods(config)
        assert isinstance(methods, str)
        assert "BEAST2" in methods
        assert "HKY" in methods
        assert "Yule" in methods
        assert "Analysis Fingerprint" in methods
        assert "B2P-" in methods

    def test_methods_contains_calibration_info(self):
        """Test that methods description contains calibration information."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        methods = MethodsGenerator.generate_methods(config)
        assert "calibration" in methods.lower()
        assert "Normal" in methods


class TestPipelineGenerator:
    """Test the PipelineGenerator class."""

    def setup_method(self, method):
        """Create a temp directory for each test."""
        self.tmpdir = Path(tempfile.mkdtemp())
        fasta_path = self.tmpdir / "test.fasta"
        fasta_path.write_text(TEST_FASTA)

    def test_generate_snakemake(self):
        """Test Snakemake pipeline generation."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        output_dir = self.tmpdir / "pipeline"
        snakefile = PipelineGenerator.generate_snakemake(
            config,
            output_dir=output_dir,
            config_file="config.yaml",
        )

        assert os.path.exists(snakefile)
        with open(snakefile, "r") as f:
            content = f.read()
        assert "Snakemake" in content
        assert "rule generate_xml" in content
        assert "rule run_beast" in content
        assert "beast2py" in content

    def test_snakemake_sequence_paths_are_relative_to_the_snakefile(self):
        """The generated Snakefile must not pin the pipeline to one machine.

        Snakemake resolves input paths against the workflow directory, so an
        absolute path here records the generating host and breaks on any other.
        """
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)
        config = ConfigParser().parse(config_path)

        output_dir = self.tmpdir / "pipeline"
        snakefile = PipelineGenerator.generate_snakemake(
            config, output_dir=output_dir, config_file="config.yaml"
        )
        content = Path(snakefile).read_text()

        assert "sequences = [" in content
        assert str(self.tmpdir) not in content, "absolute host path leaked into Snakefile"
        assert "../test.fasta" in content

    def test_snakemake_runs_from_a_relocated_workflow_directory(self):
        """A copied pipeline directory still points at its own sequence file."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)
        config = ConfigParser().parse(config_path)
        output_dir = self.tmpdir / "pipeline"
        PipelineGenerator.generate_snakemake(
            config, output_dir=output_dir, config_file="config.yaml"
        )

        moved = Path(tempfile.mkdtemp()) / "pipeline"
        shutil.copytree(output_dir, moved)
        shutil.copy(self.tmpdir / "test.fasta", moved.parent / "test.fasta")
        content = (moved / "Snakefile").read_text()
        referenced = re.findall(r"'([^']*\.fasta)'", content)
        assert len(referenced) == 1
        target = os.path.realpath(str(moved / referenced[0]))
        assert target == os.path.realpath(str(moved.parent / "test.fasta")), (
            "Snakefile still points outside the relocated tree: %s" % referenced[0])

    def test_generate_nextflow(self):
        """Test Nextflow pipeline generation."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        output_dir = self.tmpdir / "pipeline"
        nf_file = PipelineGenerator.generate_nextflow(
            config,
            output_dir=output_dir,
            config_file="config.yaml",
        )

        assert os.path.exists(nf_file)
        with open(nf_file, "r") as f:
            content = f.read()
        assert "Nextflow" in content
        assert "process generate_xml" in content
        assert "process run_beast" in content

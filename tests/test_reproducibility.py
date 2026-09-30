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

    def _methods_for_example(self, name):
        cfg = Path(__file__).resolve().parents[1] / "examples" / name
        return MethodsGenerator.generate_methods(ConfigParser().parse(cfg))

    def test_citations_use_the_author_year_in_text_form(self):
        """The paragraph is meant to be pasted into a methods section.

        In-text citations follow the author-year form (Jones 1970), so a
        parenthetical may not nest another parenthesis, carry journal titles, or
        spell co-authors with an ampersand.
        """
        methods = self._methods_for_example("config_basic.yaml")
        assert "&" not in methods, "ampersand is not the author-year in-text form"
        assert not re.search(r"\([^()]*\(", methods), "a citation nests parentheses"
        for paren in re.findall(r"\(([^()]*)\)", methods):
            if re.search(r"\b(?:19|20)\d{2}\b", paren):
                for journal in ("PLoS", "J Mol Evol", "Syst Biol", "BMC", "Mol Biol Evol"):
                    assert journal not in paren, "full reference inside a citation: %s" % paren
        assert "(Bouckaert et al. 2014)" in methods
        assert "(Hasegawa et al. 1985)" in methods

    def test_full_references_are_kept_after_the_paragraph(self):
        """Shortening in-text citations must not lose the bibliography entries."""
        methods = self._methods_for_example("config_basic.yaml")
        tail = methods[methods.find("\n\n"):] if "\n\n" in methods else ""
        assert "PLoS Comput Biol 10: e1003537" in tail
        assert "Hasegawa" in tail and "1985" in tail

    def test_tree_prior_is_named_in_words_not_config_keys(self):
        """Every tree prior gets prose a reviewer can read, not a YAML keyword."""
        methods = self._methods_for_example("config_bd_skyline.yaml")
        assert "bd_skyline_serial" not in methods
        assert re.search(r"birth[- ]death skyline", methods, re.I)

    def test_committed_methods_artefact_is_what_the_code_generates(self):
        """The shipped example paragraph must not age behind the generator."""
        ex = Path(__file__).resolve().parents[1] / "examples"
        generated = MethodsGenerator.generate_methods(
            ConfigParser().parse(ex / "config_basic.yaml"))
        assert (ex / "output_basic.methods.tex").read_text(encoding="utf-8") == generated


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

    def _generated(self, kind):
        """Generate one workflow file for TEST_CONFIG_YAML and return (config, text)."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)
        config = ConfigParser().parse(config_path)
        gen = (PipelineGenerator.generate_snakemake if kind == "sm"
               else PipelineGenerator.generate_nextflow)
        path = gen(config, output_dir=self.tmpdir / ("pipeline_" + kind),
                   config_file="config.yaml")
        return config, Path(path).read_text()

    def test_run_rule_declares_only_outputs_the_run_produces(self):
        """`rule run_beast` may not name log files BEAST2 never writes.

        The names used to come from metadata.analysis_name, while the run writes
        whatever the XML's own loggers say, so the declared outputs of the rule
        could never appear and the DAG could not be satisfied.
        """
        for kind in ("sm", "nf"):
            config, content = self._generated(kind)
            analysis = config.metadata.get("analysis_name", "analysis")
            assert "%s.log" % analysis not in content, kind
            assert "%s.trees" % analysis not in content, kind
            assert "beast2_run.done" in content, kind
            assert "touch" in content, kind

    def test_run_rule_records_where_beast2_writes(self):
        """The real logger names from the config are stated, not invented ones."""
        for kind in ("sm", "nf"):
            config, content = self._generated(kind)
            for logger in (config.loggers.trace_log, config.loggers.tree_log):
                assert logger["file_name"] in content, (kind, logger["file_name"])

    def test_nextflow_carries_every_documented_stage(self):
        """The manual promises both engines cover methods description too."""
        _, content = self._generated("nf")
        assert "process methods" in content
        assert re.search(r"workflow\s*\{[^}]*methods\(\)", content, re.S)

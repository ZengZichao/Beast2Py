"""Tests for the config module."""

import tempfile
from pathlib import Path
import pytest

from beast2py.config import ConfigParser, ConfigError
from beast2py.models import (
    BEASTConfig,
    CalibrationMethod,
    TreePriorType,
)

from .conftest import TEST_FASTA, TEST_CONFIG_YAML, TEST_MULTI_CAL_YAML


class TestConfigParser:
    """Test the ConfigParser class."""

    def setup_method(self, method):
        """Create a temp directory for each test."""
        self.tmpdir = Path(tempfile.mkdtemp())
        # Create test FASTA file
        self.fasta_path = self.tmpdir / "test.fasta"
        self.fasta_path.write_text(TEST_FASTA)

    def test_parse_basic_config(self):
        """Test parsing a basic configuration file."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        assert isinstance(config, BEASTConfig)
        assert config.metadata["analysis_name"] == "test_analysis"
        assert len(config.partitions) == 1
        assert config.partitions[0].id == "alignment"
        assert config.tree_prior_type == TreePriorType.YULE
        assert len(config.calibrations) == 1
        assert config.calibrations[0].name == "rootCalibration"

    def test_parse_multi_calibration(self):
        """Test parsing a config with multiple calibrations."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_MULTI_CAL_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        assert len(config.calibrations) == 3
        assert config.calibrations[0].name == "calA"
        assert config.calibrations[1].name == "calB"
        assert config.calibrations[2].name == "rootCal"

    def test_parse_invalid_file(self):
        """Test that parsing a non-existent file raises ConfigError."""
        parser = ConfigParser()
        with pytest.raises(ConfigError):
            parser.parse("/nonexistent/config.yaml")

    def test_parse_alignment_not_found(self):
        """Test that missing alignment file raises error."""
        bad_config = TEST_CONFIG_YAML.replace("test.fasta", "nonexistent.fasta")
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(bad_config)

        parser = ConfigParser()
        with pytest.raises(ConfigError):
            parser.parse(config_path)

    def test_parse_invalid_substitution_model(self):
        """Test that invalid substitution model raises error."""
        bad_config = TEST_CONFIG_YAML.replace('type: "hky"', 'type: "invalid_model"')
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(bad_config)

        parser = ConfigParser()
        with pytest.raises(ConfigError):
            parser.parse(config_path)

    def test_parse_invalid_tree_prior(self):
        """Test that invalid tree prior raises error."""
        bad_config = TEST_CONFIG_YAML.replace('type: "yule"', 'type: "invalid_prior"')
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(bad_config)

        parser = ConfigParser()
        with pytest.raises(ConfigError):
            parser.parse(config_path)

    def test_parse_calibration_method(self):
        """Test parsing calibration method."""
        config_with_method = TEST_CONFIG_YAML + "calibration_method: mrca_prior\n"
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(config_with_method)

        parser = ConfigParser()
        config = parser.parse(config_path)
        assert config.calibration_method == CalibrationMethod.MRCA_PRIOR

    def test_parse_provenance(self):
        """Test that provenance is correctly parsed."""
        config_with_prov = """
metadata:
  analysis_name: "test_prov"
  analysis_description: "Test with provenance"

alignments:
  - id: "alignment"
    file: "test.fasta"
    format: "auto"
    data_type: "nucleotide"

partitions:
  - id: "alignment"
    site_model:
      substitution_model:
        type: "hky"
    clock_model:
      type: "strict"
    tree: "shared"

tree_prior:
  type: "yule"

calibrations:
  - name: "rootCalibration"
    taxa: null
    monophyletic: true
    distribution:
      type: normal
      parameters: {mean: 10.0, sigma: 1.0}
    provenance:
      source: fossil
      reference: DOI:10.1038/xxx
      calibration_type: soft
      notes: Test fossil

mcmc:
  type: "standard"
  chain_length: 100000

output:
  file: "test_output.xml"
"""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(config_with_prov)

        parser = ConfigParser()
        config = parser.parse(config_path)
        assert config.calibrations[0].provenance is not None
        assert config.calibrations[0].provenance.source == "fossil"
        assert config.calibrations[0].provenance.calibration_type == "soft"

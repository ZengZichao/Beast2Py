"""Tests for the XML writer module."""

import tempfile
from pathlib import Path

from beast2py.config import ConfigParser
from beast2py.xml_writer import XMLWriter
from beast2py.validator import XMLValidator

from .conftest import TEST_FASTA, TEST_CONFIG_YAML, TEST_MULTI_CAL_YAML


class TestXMLWriter:
    """Test the XMLWriter class."""

    def setup_method(self, method):
        """Create a temp directory for each test."""
        self.tmpdir = Path(tempfile.mkdtemp())
        # Create test FASTA file
        fasta_path = self.tmpdir / "test.fasta"
        fasta_path.write_text(TEST_FASTA)

    def test_generate_basic_xml(self):
        """Test generating XML from a basic config."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        writer = XMLWriter(config)
        xml_str = writer.generate_xml()

        assert xml_str.startswith("<?xml")
        assert "<beast" in xml_str
        assert "</beast>" in xml_str

    def test_generate_xml_has_data(self):
        """Test that generated XML contains data elements."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        writer = XMLWriter(config)
        xml_str = writer.generate_xml()

        assert "<data" in xml_str
        assert 'id="alignment"' in xml_str
        assert "<sequence" in xml_str

    def test_generate_xml_has_mcmc(self):
        """Test that generated XML contains MCMC run element."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        writer = XMLWriter(config)
        xml_str = writer.generate_xml()

        assert "<run" in xml_str
        assert "MCMC" in xml_str
        assert "chainLength" in xml_str

    def test_generate_xml_has_calibration(self):
        """Test that generated XML contains calibration elements."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        writer = XMLWriter(config)
        xml_str = writer.generate_xml()

        assert "MRCAPrior" in xml_str
        assert "rootCalibration" in xml_str

    def test_generate_multi_cal_xml(self):
        """Test generating XML with multiple calibrations."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_MULTI_CAL_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        writer = XMLWriter(config)
        xml_str = writer.generate_xml()

        assert "calA" in xml_str
        assert "calB" in xml_str
        assert "rootCal" in xml_str

    def test_generated_xml_passes_validation(self):
        """Test that generated XML passes structural validation."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        writer = XMLWriter(config)
        xml_str = writer.generate_xml()

        result = XMLValidator.validate(xml_str)
        assert result.is_valid, f"Validation errors: {result.errors}"

    def test_xml_has_fingerprint(self):
        """Test that XML contains analysis fingerprint."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        writer = XMLWriter(config)
        xml_str = writer.generate_xml()

        assert "Analysis Fingerprint" in xml_str
        assert "B2P-" in xml_str

    def test_xml_has_loggers(self):
        """Test that XML contains logger elements."""
        config_path = self.tmpdir / "config.yaml"
        config_path.write_text(TEST_CONFIG_YAML)

        parser = ConfigParser()
        config = parser.parse(config_path)

        writer = XMLWriter(config)
        xml_str = writer.generate_xml()

        assert "tracelog" in xml_str
        assert "treelog" in xml_str
        assert "screenlog" in xml_str

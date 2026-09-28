"""Tests for the validator module."""

from beast2py.validator import XMLValidator, ValidationResult


class TestXMLValidator:
    """Test the XMLValidator class."""

    def test_validate_valid_xml(self):
        """Test validation of a valid XML structure."""
        xml = """<?xml version="1.0" encoding="UTF-8"?>
<beast version="2.0" namespace="beast.base.core:beast.base.inference">
    <data id="alignment" dataType="nucleotide">
        <sequence taxon="A">ACGT</sequence>
        <sequence taxon="B">ACGT</sequence>
    </data>
    <run id="mcmc" spec="MCMC" chainLength="1000000">
        <state id="state">
            <tree id="tree" spec="beast.base.evolution.tree.Tree"/>
        </state>
        <distribution idref="posterior"/>
        <logger id="tracelog" fileName="out.log" logEvery="1000"/>
        <logger id="treelog" fileName="out.trees" logEvery="1000" mode="tree"/>
        <logger id="screenlog" logEvery="10000"/>
    </run>
    <distribution id="posterior" spec="CompoundDistribution"/>
</beast>"""
        result = XMLValidator.validate(xml)
        assert result.is_valid

    def test_validate_invalid_xml_syntax(self):
        """Test that malformed XML is detected."""
        xml = "<beast><data></beast>"
        result = XMLValidator.validate(xml)
        assert not result.is_valid
        assert any("parse" in e.lower() for e in result.errors)

    def test_validate_wrong_root(self):
        """Test that wrong root element is detected."""
        xml = '<?xml version="1.0"?><foo/>'
        result = XMLValidator.validate(xml)
        assert not result.is_valid
        assert any("root" in e.lower() for e in result.errors)

    def test_validate_missing_data(self):
        """Test that missing data element is detected."""
        xml = """<?xml version="1.0"?>
<beast version="2.0" namespace="beast.base.core">
    <run id="mcmc" spec="MCMC" chainLength="1000000">
        <state id="state"/>
        <distribution idref="posterior"/>
        <logger id="tracelog" fileName="out.log" logEvery="1000"/>
    </run>
</beast>"""
        result = XMLValidator.validate(xml)
        assert any("data" in e.lower() for e in result.errors)

    def test_validate_missing_run(self):
        """Test that missing run element is detected."""
        xml = """<?xml version="1.0"?>
<beast version="2.0" namespace="beast.base.core">
    <data id="alignment" dataType="nucleotide">
        <sequence taxon="A">ACGT</sequence>
    </data>
</beast>"""
        result = XMLValidator.validate(xml)
        assert any("run" in e.lower() for e in result.errors)

    def test_validate_dangling_idref(self):
        """Test that dangling idref is detected."""
        xml = """<?xml version="1.0"?>
<beast version="2.0" namespace="beast.base.core">
    <data id="alignment" dataType="nucleotide">
        <sequence taxon="A">ACGT</sequence>
    </data>
    <run id="mcmc" spec="MCMC" chainLength="1000000">
        <state id="state"/>
        <distribution idref="nonexistent"/>
        <logger id="tracelog" fileName="out.log" logEvery="1000"/>
    </run>
</beast>"""
        result = XMLValidator.validate(xml)
        assert any("nonexistent" in e for e in result.errors)

    def test_validation_result_str(self):
        """Test ValidationResult string representation."""
        result = ValidationResult()
        result.add_error("Test error")
        result.add_warning("Test warning")
        s = str(result)
        assert "ERRORS" in s
        assert "WARNINGS" in s
        assert "Test error" in s

    def test_validation_result_valid_str(self):
        """Test ValidationResult string when valid."""
        result = ValidationResult()
        s = str(result)
        assert "valid" in s.lower()

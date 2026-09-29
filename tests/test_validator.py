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


class TestPriorCoverageRule:
    """Gate 1 asks for a prior on continuous parameters only.

    A flat prior on an unbounded continuous support is improper, so the gate must
    refuse it.  Rate-category indices, skyline group sizes and indicator bits are
    integer or boolean parameters whose domain is fixed and finite by construction,
    so they are exempt -- and the shipped examples rely on that exemption.
    """

    HEAD = """<?xml version="1.0" encoding="UTF-8"?>
<beast version="2.0" namespace="beast.base.core:beast.base.inference">
    <data id="alignment" dataType="nucleotide">
        <sequence taxon="A">ACGT</sequence>
        <sequence taxon="B">ACGT</sequence>
    </data>
    <run id="mcmc" spec="MCMC" chainLength="1000000">
        <state id="state">
%(nodes)s
        </state>
        <distribution idref="posterior"/>
        <logger id="tracelog" fileName="out.log" logEvery="1000"/>
        <logger id="treelog" fileName="out.trees" logEvery="1000" mode="tree"/>
        <logger id="screenlog" logEvery="10000"/>
    </run>
    <distribution id="posterior" spec="CompoundDistribution">%(priors)s</distribution>
</beast>"""

    PRIOR = """
        <distribution id="x.prior" spec="Prior" x="@x">
            <distr spec="beast.base.inference.distribution.Uniform" lower="0" upper="10"/>
        </distribution>"""

    def _xml(self, nodes, priors=""):
        return self.HEAD % {"nodes": nodes, "priors": priors}

    def test_continuous_parameter_without_a_prior_is_refused(self):
        xml = self._xml(
            '<parameter id="x" spec="beast.base.inference.parameter.RealParameter" '
            'value="1" lower="0" estimate="true"/>')
        result = XMLValidator.validate(xml)
        assert not result.is_valid
        assert any("has no <Prior>" in e for e in result.errors), result.errors

    def test_a_continuous_parameter_with_a_prior_passes_the_rule(self):
        # non-vacuity: the same file with the prior written down is accepted
        node = ('<parameter id="x" spec="beast.base.inference.parameter.RealParameter" '
                'value="1" lower="0" estimate="true"/>')
        result = XMLValidator.validate(self._xml(node, self.PRIOR))
        assert not any("has no <Prior>" in e for e in result.errors), result.errors

    def test_integer_and_boolean_parameters_are_exempt(self):
        nodes = (
            '<parameter id="rateCategories" '
            'spec="beast.base.inference.parameter.IntegerParameter" value="0" '
            'lower="0" dimension="4" estimate="true"/> '
            '<parameter id="indicators" '
            'spec="beast.base.inference.parameter.BooleanParameter" value="0 0 0 0" '
            'dimension="4" estimate="true"/>')
        result = XMLValidator.validate(self._xml(nodes))
        assert not any("has no <Prior>" in e for e in result.errors), result.errors

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
            'value="1" lower="0" estimate="true"/>'
        )
        result = XMLValidator.validate(xml)
        assert not result.is_valid
        assert any("has no <Prior>" in e for e in result.errors), result.errors

    def test_a_continuous_parameter_with_a_prior_passes_the_rule(self):
        # non-vacuity: the same file with the prior written down is accepted
        node = (
            '<parameter id="x" spec="beast.base.inference.parameter.RealParameter" '
            'value="1" lower="0" estimate="true"/>'
        )
        result = XMLValidator.validate(self._xml(node, self.PRIOR))
        assert not any("has no <Prior>" in e for e in result.errors), result.errors

    def test_integer_and_boolean_parameters_are_exempt(self):
        nodes = (
            '<parameter id="rateCategories" '
            'spec="beast.base.inference.parameter.IntegerParameter" value="0" '
            'lower="0" dimension="4" estimate="true"/> '
            '<parameter id="indicators" '
            'spec="beast.base.inference.parameter.BooleanParameter" value="0 0 0 0" '
            'dimension="4" estimate="true"/>'
        )
        result = XMLValidator.validate(self._xml(nodes))
        assert not any("has no <Prior>" in e for e in result.errors), result.errors


class TestTheParseGateRefusesToPretendItRan:
    """A missing validator is a setup error, never a verdict about the model.

    The bundled launcher assembles its class path from ``beast2py/tools``. When
    that directory is absent -- a source install that omitted it, a partial
    checkout -- java used to be reached anyway and answered with a localized
    "could not find or load main class" and exit 1, which a caller gating on
    exit codes reads as "BEAST2 rejected this XML". Both the launcher and this
    layer now stop before the JVM, and say so in English.
    """

    @staticmethod
    def _launcher_without_validator(tmp_path):
        """A real launcher file whose sibling ``tools`` directory does not exist."""
        script = tmp_path / "beast2_validate.sh"
        script.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
        return script

    def test_the_setup_problem_is_named_before_any_jvm_starts(self, tmp_path, monkeypatch):
        """The check must short-circuit: a started JVM is already the defect."""
        import subprocess

        def _no_process_may_start(*args, **kwargs):
            raise AssertionError("a JVM was reached for a validator that is not installed")

        monkeypatch.setattr(subprocess, "run", _no_process_may_start)
        monkeypatch.setattr(
            XMLValidator,
            "_bundled_validator_script",
            staticmethod(lambda: self._launcher_without_validator(tmp_path)),
        )
        monkeypatch.setattr("shutil.which", lambda name: None)

        success, output = XMLValidator.validate_with_beast2("in.xml")
        assert not success
        assert output.startswith("SETUP-ERROR")

    def test_the_diagnostic_says_which_gate_could_not_run_and_why(self, tmp_path, monkeypatch):
        """It names the parse gate and the missing helper, not a stack trace."""
        monkeypatch.setattr(
            XMLValidator,
            "_bundled_validator_script",
            staticmethod(lambda: self._launcher_without_validator(tmp_path)),
        )
        _, output = XMLValidator.validate_with_beast2("in.xml")
        assert "BEAST2 parse gate could not be run" in output
        assert "Beast2Validator" in output
        assert "not a pass and not a rejection" in output
        assert "ClassNotFound" not in output

    def test_the_status_is_unavailable_rather_than_failed(self, tmp_path, monkeypatch):
        """The status that lets a caller exit 2 instead of reporting a failure."""
        monkeypatch.setattr(
            XMLValidator,
            "_bundled_validator_script",
            staticmethod(lambda: self._launcher_without_validator(tmp_path)),
        )
        success, _, status = XMLValidator.beast2_validation_status("in.xml")
        assert not success
        assert status == "unavailable"

    def test_a_jvm_class_path_failure_also_reads_as_unavailable(self, monkeypatch):
        """The locale-independent half of the old crash is caught too.

        An installation carrying the *previous* launcher still reaches java, so
        the status layer keeps the signature: the exception name survives any
        locale, while the readable sentence around it does not.
        """
        crash = (
            "\u9519\u8bef: \u627e\u4e0d\u5230\u6216\u65e0\u6cd5\u52a0\u8f7d"
            "\u4e3b\u7c7b Beast2Validator\n"
            "\u539f\u56e0: java.lang.ClassNotFoundException: Beast2Validator"
        )
        monkeypatch.setattr(
            XMLValidator, "validate_with_beast2", staticmethod(lambda *a, **k: (False, crash))
        )
        _, _, status = XMLValidator.beast2_validation_status("in.xml")
        assert status == "unavailable"

    def test_a_launcher_that_does_ship_the_validator_is_left_alone(self, tmp_path, monkeypatch):
        """Non-vacuity: the guard fires on the missing helper, not on every path."""
        script = self._launcher_without_validator(tmp_path)
        real = XMLValidator._bundled_validator_script()
        assert (
            real is not None
        ), "the package must ship its own launcher for this test to mean anything"
        assert XMLValidator._validator_setup_problem(real, True) is None
        assert XMLValidator._validator_setup_problem(script, True) is not None
        # a script the user named is taken at face value: only its exit code decides
        assert XMLValidator._validator_setup_problem(script, False) is None

    def test_validate_file_still_reports_the_structural_verdict(self):
        """The structural level is untouched by the gate's new behaviour.

        The shipped worked example is the strictest file in the repository, so it
        is the one worth asserting level 1 still accepts on its own merits.
        """
        from .conftest import EXAMPLES_DIR

        committed = EXAMPLES_DIR / "output_basic.xml"
        assert committed.is_file(), "the worked example must ship with the package"
        assert XMLValidator.validate_file(str(committed)).is_valid

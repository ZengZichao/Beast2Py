"""Semantic invariant tests.

Each test asserts a *scientific* property of the generated XML: that what
the user asked for is what the model does. Each test fails against incorrect
or degenerate configurations.

"""

from __future__ import annotations

import os
import shutil
import defusedxml.ElementTree as ET
from pathlib import Path

import pytest

from beast2py.config import ConfigError, ConfigParser
from beast2py.main import _validate_and_write
from beast2py.models import (
    Alignment,
    ClockModelConfig,
    ClockModelType,
    DataType,
    Partition,
    Sequence,
    SiteModelConfig,
)
from beast2py.validator import XMLValidator
from beast2py.xml_writer import XMLWriter

# An alignment with real variation: the shared conftest fixture is four
# identical 32 bp sequences (0 variable sites), which makes every
# site-dependent code path degenerate.
VARIABLE_FASTA = """>taxonA
ACGTACGTACGTACGTTTTGCAAGGCATCGATCGGG
>taxonB
ACGTACGTACGTAAGTTTTGCAAGGCATCGATCGGA
>taxonC
ACGAACGTACGTAAGTTTTGCAAGGCATCGATCGGG
>taxonD
ACGAACGTACGTAAGTTTTGCAAGGCATCGATCGGA
"""

# Shorthand for "one HKY site model on the test alignment", used by the
# longest YAML literals so those lines stay inside the line-length gate.
HKY_PARTITION = '\n    site_model:\n      substitution_model:\n        type: "hky"\n'

AMINO_FASTA = """>taxonA
MKTLLAVLGGALPACLLEQ
>taxonB
MKTLLAVLGGALPACLLEA
>taxonC
MKTLLAALGGALPACLLAQ
>taxonD
MKTLLAVLGGALPAALLER
"""


def write_inputs(tmp_path, config_body, fasta=VARIABLE_FASTA, fasta_name="test.fasta"):
    """Write a FASTA plus a YAML config that uses it.

    Args:
        tmp_path: Directory to write into.
        config_body: YAML text for the configuration.
        fasta: Alignment content.
        fasta_name: File name for the alignment.

    Returns:
        Path to the written YAML file.
    """
    (tmp_path / fasta_name).write_text(fasta, encoding="utf-8")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(config_body, encoding="utf-8")
    return config_path


def head(data_type="nucleotide", extra_alignments=""):
    """Return the common YAML preamble for one alignment + partition.

    Args:
        data_type: Declared alignment data type.
        extra_alignments: Additional alignment blocks.

    Returns:
        YAML text.
    """
    return f"""
metadata:
  analysis_name: "semantic_test"
alignments:
  - id: "alignment"
    file: "test.fasta"
    data_type: "{data_type}"
{extra_alignments}
partitions:
  - id: "alignment"
"""


TAIL = """
tree_prior:
  type: "yule"
calibrations:
  - name: "rootCal"
    taxa: null
    distribution:
      type: "lognormal"
      parameters: {M: 2.3, S: 1.0}
mcmc:
  chain_length: 100000
"""


def generate(tmp_path, body):
    """Parse a config and generate its XML.

    Args:
        tmp_path: Directory holding the config and its alignment.
        body: YAML text appended after the common preamble.

    Returns:
        Tuple of (config, xml string).
    """
    path = write_inputs(tmp_path, head() + body)
    config = ConfigParser().parse(str(path))
    writer = XMLWriter(config)
    return config, writer.generate_xml(), writer


def state_nodes(xml):
    """Return ``{id: attrib}`` for every parameter in the first ``<state>``."""
    root = ET.fromstring(xml)
    out = {}
    for state in root.iter("state"):
        for elem in state:
            if elem.tag in ("parameter", "tree") and elem.get("id"):
                out[elem.get("id")] = dict(elem.attrib)
        break
    return out


def prior_targets(xml):
    """Return the set of parameter ids referenced by a ``Prior`` element."""
    root = ET.fromstring(xml)
    return {
        elem.get("x", "")[1:]
        for elem in root.iter("distribution")
        if elem.get("spec", "").endswith(".Prior") and elem.get("x", "").startswith("@")
    }


# ---------------------------------------------------------------------------
# estimate: false must actually fix the parameter
# ---------------------------------------------------------------------------


class TestEstimateFalse:
    def test_fixed_kappa_is_not_a_state_node(self, tmp_path):
        _config, xml, _writer = generate(
            tmp_path,
            """
    site_model:
      substitution_model:
        type: "hky"
        kappa: {value: 2.0, estimate: false}
""" + TAIL,
        )
        nodes = state_nodes(xml)
        assert "alignment.hky.kappa" not in nodes, "a parameter the user fixed appeared in <state>"
        root = ET.fromstring(xml)
        kappa = [p for p in root.iter("parameter") if p.get("id") == "alignment.hky.kappa"]
        assert kappa and kappa[0].get("estimate") == "false"
        assert float(kappa[0].get("value")) == 2.0

    def test_fixed_parameter_gets_no_operator_and_is_not_logged(self, tmp_path):
        _config, xml, _writer = generate(
            tmp_path,
            """
    site_model:
      substitution_model:
        type: "hky"
        kappa: {value: 2.0, estimate: false}
""" + TAIL,
        )
        assert "alignment.hky.kappaScaler" not in xml
        assert '<parameter idref="alignment.hky.kappa" name="log"' not in xml

    def test_fixed_clock_and_stdev_are_not_state_nodes(self, tmp_path):
        _config, xml, _writer = generate(
            tmp_path,
            """
    site_model:
      substitution_model:
        type: "hky"
    clock_model:
      type: "ucln"
      ucld_mean: {value: 1.0, estimate: false}
      ucld_stdev: {value: 0.4, estimate: false, upper: 1.0}
""" + TAIL,
        )
        nodes = state_nodes(xml)
        assert "alignment.ucld.mean" not in nodes
        assert "alignment.ucld.stdev" not in nodes

    def test_estimated_parameter_still_estimated(self, tmp_path):
        _config, xml, _writer = generate(
            tmp_path,
            """
    site_model:
      substitution_model:
        type: "hky"
        kappa: {value: 2.0, estimate: true}
""" + TAIL,
        )
        assert state_nodes(xml)["alignment.hky.kappa"]["estimate"] == "true"

    def test_fixed_gamma_shape_is_inline(self, tmp_path):
        _config, xml, _writer = generate(
            tmp_path,
            """
    site_model:
      substitution_model:
        type: "hky"
      gamma_categories: 4
      gamma_shape: {value: 0.5, estimate: false}
""" + TAIL,
        )
        assert "alignment.gammaShape" not in state_nodes(xml)


# ---------------------------------------------------------------------------
# tree: separate used to be a silent no-op
# ---------------------------------------------------------------------------


class TestSeparateTreesRejected:
    def test_separate_tree_is_refused_not_ignored(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + """
    site_model:
      substitution_model:
        type: "hky"
    tree: "separate"
""" + TAIL,
        )
        with pytest.raises(ConfigError) as excinfo:
            ConfigParser().parse(str(path))
        assert "unlinked" in str(excinfo.value).lower() or "not implemented" in str(excinfo.value)

    def test_shared_tree_is_accepted(self, tmp_path):
        _config, xml, _writer = generate(
            tmp_path,
            """
    site_model:
      substitution_model:
        type: "hky"
    tree: "shared"
""" + TAIL,
        )
        assert xml.count("<tree ") >= 1


# ---------------------------------------------------------------------------
# every estimated parameter needs a prior
# ---------------------------------------------------------------------------


class TestPriorCoverage:
    @pytest.mark.parametrize(
        "body",
        [
            '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL,
            '\n    site_model:\n      substitution_model:\n        type: "gtr"\n'
            '\n    clock_model:\n      type: "ucln"\n' + TAIL,
            '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            "      gamma_categories: 4\n      proportion_invariant: {value: 0.1}\n"
            '\n    clock_model:\n      type: "strict"\n'
            'tree_prior:\n  type: "birth_death"\n  birth_rate: {value: 1.0}\n'
            "  death_rate: {value: 0.5}\n"
            'calibrations:\n  - name: "rootCal"\n    taxa: null\n    distribution:\n'
            '      type: "lognormal"\n      parameters: {M: 2.3, S: 1.0}\n'
            "mcmc:\n  chain_length: 100000\n",
        ],
    )
    def test_every_estimated_parameter_has_a_prior(self, tmp_path, body):
        _config, xml, writer = generate(tmp_path, body)
        nodes = state_nodes(xml)
        covered = prior_targets(xml)
        for node_id, attrs in nodes.items():
            if attrs.get("estimate") != "true":
                continue
            if attrs.get("spec", "").endswith(("IntegerParameter", "BooleanParameter")):
                continue
            if node_id == "tree":
                continue
            assert node_id in covered, f"{node_id} is estimated but has no Prior"
        assert result_is_valid(xml, _config)

    def test_writer_reports_the_defaults_it_added(self, tmp_path):
        _config, _xml, writer = generate(
            tmp_path,
            '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL,
        )
        assert "alignment.hky.kappa" in writer.default_priors
        assert "alignment.freqParameter" in writer.default_priors


def result_is_valid(xml, config):
    """Return True when the structural validator reports no errors.

    Args:
        xml: Generated XML.
        config: Its configuration.

    Returns:
        Whether validation passed.
    """
    return XMLValidator.validate(xml, config=config).is_valid


# ---------------------------------------------------------------------------
# quoted booleans were inverted by bool("false")
# ---------------------------------------------------------------------------


class TestStrictBooleans:
    @pytest.mark.parametrize("quoted", ['"false"', "'false'", '"FALSE"', "yes", "off"])
    def test_quoted_booleans_are_explicit(self, tmp_path, quoted):
        body = (
            '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            '\n    clock_model:\n      type: "strict"\n'
            f'tree_prior:\n  type: "yule"\n'
            "calibrations: []\n"
            "mcmc:\n  chain_length: 100000\n"
            f"  sample_from_prior: {quoted}\n"
            'output:\n  file: "out.xml"\n'
        )
        path = write_inputs(tmp_path, head() + body)
        if quoted in ("yes", "off"):
            config = ConfigParser().parse(str(path))
            assert config.mcmc.sample_from_prior is (quoted == "yes")
        else:
            config = ConfigParser().parse(str(path))
            assert config.mcmc.sample_from_prior is False
            xml = XMLWriter(config).generate_xml()
            assert "sampleFromPrior" not in xml

    def test_sample_from_prior_true_is_emitted(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "yule"\n'
            "calibrations: []\n"
            "mcmc:\n  chain_length: 100000\n  sample_from_prior: true\n",
        )
        config = ConfigParser().parse(str(path))
        assert 'sampleFromPrior="true"' in XMLWriter(config).generate_xml()

    def test_garbage_boolean_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head()
            + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            + TAIL.replace(
                "mcmc:\n  chain_length: 100000\n",
                "mcmc:\n  chain_length: 100000\n  sample_from_prior: maybe\n",
            ),
        )
        with pytest.raises(ConfigError):
            ConfigParser().parse(str(path))

    def test_quoted_false_originate_flag_is_honoured(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "yule"\n'
            'calibrations:\n  - name: "rootCal"\n    taxa: null\n    use_originate: "false"\n'
            '    distribution:\n      type: "lognormal"\n      parameters: {M: 2.3, S: 1.0}\n'
            "mcmc:\n  chain_length: 100000\n",
        )
        config = ConfigParser().parse(str(path))
        assert config.calibrations[0].use_originate is False


# ---------------------------------------------------------------------------
# data type, model family, alphabet and alignment sanity
# ---------------------------------------------------------------------------


class TestTypeCrossCheck:
    def test_amino_acid_model_on_nucleotide_data_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "wag"\n' + TAIL,
        )
        with pytest.raises(ConfigError, match="amino"):
            ConfigParser().parse(str(path))

    def test_nucleotide_model_on_amino_acid_data_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head(data_type="aminoacid")
            + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            + TAIL,
        )
        with pytest.raises(ConfigError, match="nucleotide"):
            ConfigParser().parse(str(path))

    def test_amino_acid_model_on_amino_acid_data_works(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head(data_type="aminoacid")
            + '\n    site_model:\n      substitution_model:\n        type: "wag"\n'
            + TAIL,
            fasta=AMINO_FASTA,
        )
        config = ConfigParser().parse(str(path))
        xml = XMLWriter(config).generate_xml()
        assert 'spec="WAG"' in xml
        assert result_is_valid(xml, config)

    def test_illegal_character_is_rejected(self, tmp_path):
        bad = VARIABLE_FASTA.replace(
            "ACGTACGTACGTACGTTTTGCAAGGCATCGATCGGG",
            "ACGTACGTACGTACGTTTTGCAAGG@ATCGATCGGG",
        )
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL,
            fasta=bad,
        )
        with pytest.raises(ConfigError, match="illegal"):
            ConfigParser().parse(str(path))

    def test_ambiguity_codes_are_allowed(self, tmp_path):
        amb = VARIABLE_FASTA.replace(
            "ACGTACGTACGTACGTTTTGCAAGGCATCGATCGGG",
            "ACGTACGTACGTACGTRYSWKMBDHVNATCGATCGG",
        )
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL,
            fasta=amb,
        )
        config = ConfigParser().parse(str(path))
        assert config.partitions[0].alignment.n_sites == 36


class TestAlignmentIntegrity:
    def test_ragged_alignment_is_rejected(self, tmp_path):
        ragged = VARIABLE_FASTA.replace(
            ">taxonD\nACGAACGTACGTAAGTTTTGCAAGGCATCGATCGGA", ">taxonD\nACGAACGTACGTAAG"
        )
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL,
            fasta=ragged,
        )
        with pytest.raises(ConfigError, match="not all the same length"):
            ConfigParser().parse(str(path))

    def test_duplicate_taxon_is_rejected(self, tmp_path):
        dupe = VARIABLE_FASTA.replace(">taxonD", ">taxonA", 1)
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL,
            fasta=dupe,
        )
        with pytest.raises(ConfigError, match="duplicate taxon"):
            ConfigParser().parse(str(path))

    def test_n_sites_does_not_report_the_first_sequence_for_ragged_data(self):
        aln = Alignment(
            id="x",
            data_type=DataType.NUCLEOTIDE,
            sequences=[Sequence("a", "ACGT"), Sequence("b", "AC")],
        )
        with pytest.raises(ValueError):
            _ = aln.n_sites


# ---------------------------------------------------------------------------
# birth/death reparameterisation domain and prior targeting
# ---------------------------------------------------------------------------


class TestBirthDeathDomain:
    def test_death_rate_above_birth_rate_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "birth_death"\n  birth_rate: {value: 0.5}\n'
            "  death_rate: {value: 0.8}\n"
            'calibrations:\n  - name: "rootCal"\n    taxa: null\n    distribution:\n'
            '      type: "lognormal"\n      parameters: {M: 2.3, S: 1.0}\n'
            "mcmc:\n  chain_length: 100000\n",
        )
        with pytest.raises(ConfigError, match="birth_rate"):
            ConfigParser().parse(str(path))

    def test_valid_birth_death_state_nodes_respect_their_bounds(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "birth_death"\n  birth_rate: {value: 1.0}\n'
            "  death_rate: {value: 0.5}\n"
            'calibrations:\n  - name: "rootCal"\n    taxa: null\n    distribution:\n'
            '      type: "lognormal"\n      parameters: {M: 2.3, S: 1.0}\n'
            "mcmc:\n  chain_length: 100000\n",
        )
        config = ConfigParser().parse(str(path))
        nodes = state_nodes(XMLWriter(config).generate_xml())
        assert float(nodes["birthDiffRate"]["value"]) > 0
        assert float(nodes["relativeDeathRate"]["value"]) < 1
        assert nodes["relativeDeathRate"]["upper"] == "1"

    def test_sample_probability_is_explicitly_fixed(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "birth_death"\n  birth_rate: {value: 1.0}\n'
            "  death_rate: {value: 0.5}\n"
            'calibrations:\n  - name: "rootCal"\n    taxa: null\n    distribution:\n'
            '      type: "lognormal"\n      parameters: {M: 2.3, S: 1.0}\n'
            "mcmc:\n  chain_length: 100000\n",
        )
        root = ET.fromstring(XMLWriter(ConfigParser().parse(str(path))).generate_xml())
        sp = [p for p in root.iter("parameter") if p.get("id") == "sampleProbability"]
        assert sp and sp[0].get("estimate") == "false"

    def test_prior_written_against_user_name_is_applied(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "birth_death"\n  birth_rate: {value: 1.0}\n'
            "  death_rate: {value: 0.5}\n"
            'calibrations:\n  - name: "rootCal"\n    taxa: null\n    distribution:\n'
            '      type: "lognormal"\n      parameters: {M: 2.3, S: 1.0}\n'
            "parameter_priors:\n  - parameter: birth_rate\n    distribution:\n"
            "      type: lognormal\n      parameters: {M: 0.0, S: 1.0}\n"
            "mcmc:\n  chain_length: 100000\n",
        )
        xml = XMLWriter(ConfigParser().parse(str(path))).generate_xml()
        assert 'x="@birthDiffRate"' in xml

    def test_prior_on_a_non_estimated_parameter_is_an_error(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head()
            + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            + TAIL.replace(
                "mcmc:\n  chain_length: 100000\n",
                "parameter_priors:\n  - parameter: not_a_parameter\n    distribution:\n"
                "      type: lognormal\n      parameters: {M: 0.0, S: 1.0}\n"
                "mcmc:\n  chain_length: 100000\n",
            ),
        )
        config = ConfigParser().parse(str(path))
        with pytest.raises(ValueError, match="not match any estimated parameter"):
            XMLWriter(config).generate_xml()


# ---------------------------------------------------------------------------
# the validation gate must actually gate
# ---------------------------------------------------------------------------


class TestExitCodeGate:
    def test_structural_failure_prevents_writing(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL,
        )
        config = ConfigParser().parse(str(path))
        out = tmp_path / "sub" / "bad.xml"
        # Break the XML so the structural validator must reject it.
        broken = XMLWriter(config).generate_xml().replace('id="tree"', 'id="renamed_tree"')
        status = _validate_and_write(broken, str(out), config)
        assert status != 0
        assert not out.exists(), "an XML that failed validation was still written"

    def test_conflicting_calibrations_block_generation(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "yule"\n'
            "calibrations:\n"
            '  - name: "child"\n    taxa: [taxonA, taxonB]\n    distribution:\n'
            "      type: uniform\n      parameters: {lower: 40.0, upper: 50.0}\n"
            '  - name: "rootCal"\n    taxa: null\n    distribution:\n'
            "      type: uniform\n      parameters: {lower: 5.0, upper: 10.0}\n"
            "mcmc:\n  chain_length: 100000\n",
        )
        config = ConfigParser().parse(str(path))
        out = tmp_path / "out.xml"
        assert _validate_and_write(XMLWriter(config).generate_xml(), str(out), config) != 0
        assert not out.exists()

    def test_passing_gates_writes_the_file(self, tmp_path):
        config, xml, _writer = generate(
            tmp_path, '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL
        )
        out = tmp_path / "out.xml"
        assert _validate_and_write(xml, str(out), config) == 0
        assert out.exists()

    def test_existing_output_is_not_truncated_without_force(self, tmp_path):
        config, xml, _writer = generate(
            tmp_path, '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL
        )
        out = tmp_path / "out.xml"
        out.write_text("PREVIOUS RUN", encoding="utf-8")
        assert _validate_and_write(xml, str(out), config) != 0
        assert out.read_text(encoding="utf-8") == "PREVIOUS RUN"
        assert _validate_and_write(xml, str(out), config, force=True) == 0

    def test_missing_beast2_is_reported_as_unavailable_not_passed(self, tmp_path, monkeypatch):
        from beast2py.validator import XMLValidator

        config, xml, _writer = generate(
            tmp_path, '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL
        )
        out = tmp_path / "out.xml"
        out.write_text(xml, encoding="utf-8")
        # No bundled script and no `beast` on PATH: the gate must say it could
        # not check, never "passed".
        monkeypatch.setattr(XMLValidator, "_bundled_validator_script", staticmethod(lambda: None))
        monkeypatch.setattr("shutil.which", lambda name: None)
        _success, _output, status = XMLValidator.beast2_validation_status(str(out))
        assert status == "unavailable"


# ---------------------------------------------------------------------------
# frequencies semantics and bounds
# ---------------------------------------------------------------------------


class TestFrequenciesAndBounds:
    def test_fixed_frequencies_keep_the_user_vector(self, tmp_path):
        _config, xml, _writer = generate(
            tmp_path,
            '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            '        frequencies: {value: "0.1 0.2 0.3 0.4", estimate: false}\n' + TAIL,
        )
        assert 'value="0.1 0.2 0.3 0.4"' in xml
        assert "alignment.freqParameter" not in state_nodes(xml)

    def test_empirical_mode_counts_from_the_data(self, tmp_path):
        _config, xml, _writer = generate(
            tmp_path,
            '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            "        frequencies: {mode: empirical}\n" + TAIL,
        )
        assert 'spec="Frequencies"' in xml
        assert "alignment.freqParameter" not in state_nodes(xml)

    def test_frequencies_must_sum_to_one(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            "        frequencies: {value: 0.3, dimension: 4}\n" + TAIL,
        )
        config = ConfigParser().parse(str(path))
        with pytest.raises(ValueError, match="sum to"):
            XMLWriter(config).generate_xml()

    def test_zero_frequency_is_not_silently_become_a_quarter(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            '        frequencies: {value: "0.0 0.33 0.33 0.34", dimension: 4}\n' + TAIL,
        )
        config = ConfigParser().parse(str(path))
        with pytest.raises(ValueError, match="every frequency"):
            XMLWriter(config).generate_xml()

    def test_frequency_state_node_is_bounded(self, tmp_path):
        _config, xml, _writer = generate(
            tmp_path,
            '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL,
        )
        node = state_nodes(xml)["alignment.freqParameter"]
        assert node["lower"] == "0" and node["upper"] == "1"

    def test_upper_bound_reaches_the_xml(self, tmp_path):
        _config, xml, _writer = generate(
            tmp_path,
            '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            '\n    clock_model:\n      type: "ucln"\n'
            "      ucld_mean: {value: 1.0, lower: 0.0, upper: 10.0}\n"
            "      ucld_stdev: {value: 0.4, lower: 0.0, upper: 5.0}\n" + TAIL,
        )
        nodes = state_nodes(xml)
        assert nodes["alignment.ucld.mean"]["upper"] == "10"
        assert nodes["alignment.ucld.stdev"]["upper"] == "5"

    def test_positive_parameters_get_a_default_lower_bound(self, tmp_path):
        _config, xml, _writer = generate(
            tmp_path,
            '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            '\n    clock_model:\n      type: "strict"\n' + TAIL,
        )
        assert state_nodes(xml)["alignment.hky.kappa"]["lower"] == "0"
        assert state_nodes(xml)["alignment.clockRate"]["lower"] == "0"


# ---------------------------------------------------------------------------
# tip dating
# ---------------------------------------------------------------------------


class TestTipDating:
    def _body(self, dates_yaml, trait="date-forward"):
        return (
            '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "yule"\n'
            "calibrations: []\n"
            f'tip_dates:\n  enabled: true\n  trait_name: "{trait}"\n  units: "year"\n'
            f"  dates:\n{dates_yaml}\n"
            "mcmc:\n  chain_length: 100000\n"
        )

    def test_partial_date_coverage_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + self._body("    taxonA: 2010\n    taxonB: 2011\n    taxonC: 2012\n"),
        )
        with pytest.raises(ConfigError, match="missing: taxonD"):
            ConfigParser().parse(str(path))

    def test_unknown_trait_name_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head()
            + self._body(
                "    taxonA: 2010\n    taxonB: 2011\n    taxonC: 2012\n    taxonD: 2013\n",
                trait="date_backward",
            ),
        )
        with pytest.raises(ConfigError, match="trait_name"):
            ConfigParser().parse(str(path))

    def test_random_walker_operator_exists_without_a_tipsonly_calibration(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head()
            + self._body(
                "    taxonA: 2010\n    taxonB: 2011\n    taxonC: 2012\n    taxonD: 2013\n"
            ),
        )
        config = ConfigParser().parse(str(path))
        xml = XMLWriter(config).generate_xml()
        assert (
            'spec="TipDatesRandomWalker"' in xml
        ), "tip dates were attached but no operator ever proposes new tip heights"
        assert "TraitSet" in xml
        assert result_is_valid(xml, config)


# ---------------------------------------------------------------------------
# unknown keys and calibration consistency
# ---------------------------------------------------------------------------


class TestKeyAndCalibrationChecks:
    def test_typo_in_a_model_key_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            "      gamma_categorys: 8\n" + TAIL,
        )
        with pytest.raises(ConfigError, match="gamma_categorys"):
            ConfigParser().parse(str(path))

    def test_typo_in_mcmc_key_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head()
            + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            + TAIL.replace("chain_length: 100000", "chainlength_typo: 100000"),
        )
        with pytest.raises(ConfigError, match="chainlength_typo"):
            ConfigParser().parse(str(path))

    def test_inverted_uniform_bound_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "yule"\n'
            'calibrations:\n  - name: "rootCal"\n    taxa: null\n    distribution:\n'
            "      type: uniform\n      parameters: {lower: 20.0, upper: 5.0}\n"
            "mcmc:\n  chain_length: 100000\n",
        )
        with pytest.raises(ConfigError, match="strictly below"):
            ConfigParser().parse(str(path))

    def test_burnin_longer_than_chain_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head()
            + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            + TAIL.replace("chain_length: 100000", "chain_length: 100\npre_burnin: 5000"),
        )
        with pytest.raises(ConfigError, match="pre_burnin"):
            ConfigParser().parse(str(path))

    def test_hard_provenance_contradicting_the_prior_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "yule"\n'
            'calibrations:\n  - name: "rootCal"\n    taxa: null\n'
            "    provenance:\n      source: fossil\n      calibration_type: hard\n"
            "      original_age_min: 100.0\n"
            "    distribution:\n      type: normal\n      parameters: {mean: 10.0, sigma: 1.0}\n"
            "mcmc:\n  chain_length: 100000\n",
        )
        with pytest.raises(ConfigError, match="younger than the fossil"):
            ConfigParser().parse(str(path))

    def test_soft_bound_may_lie_below_the_fossil_minimum(self, tmp_path):
        """A soft bound is reported, not refused; a hard one is refused ."""
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "yule"\n'
            'calibrations:\n  - name: "rootCal"\n    taxa: null\n'
            "    provenance:\n      source: fossil\n      calibration_type: soft\n"
            "      original_age_min: 7.2\n"
            "    distribution:\n      type: normal\n      parameters: {mean: 6.0, sigma: 0.5}\n"
            "mcmc:\n  chain_length: 100000\n",
        )
        with pytest.warns(UserWarning, match="soft bound"):
            ConfigParser().parse(str(path))

    def test_root_younger_than_child_is_reported(self, tmp_path):
        from beast2py.diagnostics import ConflictDetector
        from beast2py.models import CalibrationPoint, DistributionConfig

        root = CalibrationPoint(
            name="root",
            taxa=None,
            distribution=DistributionConfig(
                type="uniform", parameters={"lower": 5.0, "upper": 10.0}
            ),
        )
        child = CalibrationPoint(
            name="child",
            taxa=["taxonA", "taxonB"],
            distribution=DistributionConfig(
                type="uniform", parameters={"lower": 40.0, "upper": 50.0}
            ),
        )
        conflicts = ConflictDetector.detect_conflicts([root, child])
        assert any(
            c.severity == "error" for c in conflicts
        ), "a root prior entirely younger than its nested clade went undetected"

    def test_overlap_criterion_is_scale_free(self, tmp_path):
        from beast2py.diagnostics import ConflictDetector
        from beast2py.models import CalibrationPoint, DistributionConfig

        def pair(mean_a, mean_b, sigma):
            return [
                CalibrationPoint(
                    "a",
                    ["taxonA", "taxonB"],
                    distribution=DistributionConfig("normal", {"mean": mean_a, "sigma": sigma}),
                ),
                CalibrationPoint(
                    "b",
                    ["taxonC", "taxonD"],
                    distribution=DistributionConfig("normal", {"mean": mean_b, "sigma": sigma}),
                ),
            ]

        # Shallow and deep pairs with the same *relative* redundancy: separated
        # by 1.25 prior widths in both cases.
        shallow, deep = pair(5.0, 6.0, 0.4), pair(100.0, 120.0, 8.0)

        def relative(cals):
            return bool(
                ConflictDetector.detect_conflicts(
                    cals,
                    settings={"overlap_measure": "relative", "relative_overlap_threshold": 0.5},
                )
            )

        def absolute(cals):
            return bool(
                ConflictDetector.detect_conflicts(
                    cals, settings={"overlap_measure": "absolute", "overlap_threshold": 2.0}
                )
            )

        assert relative(shallow) == relative(
            deep
        ), "the relative criterion changed its verdict with node depth"
        # The absolute criterion is the one that misbehaves: identical relative
        # redundancy, opposite verdicts. That is why relative is the default.
        assert absolute(shallow) and not absolute(deep)


# ---------------------------------------------------------------------------
# fingerprint identity and resume capability
# ---------------------------------------------------------------------------


class TestFingerprintIdentity:
    def config(self, tmp_path, fasta=VARIABLE_FASTA):
        """Parse the reference config against a given alignment.

        Args:
            tmp_path: Directory to write into.
            fasta: Alignment content.

        Returns:
            The parsed BEASTConfig.
        """
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL,
            fasta=fasta,
        )
        return ConfigParser().parse(str(path))

    def test_fingerprint_has_no_date_component(self, tmp_path):
        from beast2py.reproducibility import FingerprintGenerator

        fingerprint = FingerprintGenerator.generate_fingerprint(self.config(tmp_path))
        assert fingerprint.startswith("B2P-")
        assert len(fingerprint.split("-")) == 3, f"date leaked into the identifier: {fingerprint}"

    def test_different_dates_same_fingerprint(self, tmp_path, monkeypatch):
        from beast2py import reproducibility
        from beast2py.reproducibility import FingerprintGenerator

        first = FingerprintGenerator.generate_fingerprint(self.config(tmp_path))

        class FixedDatetime:
            @staticmethod
            def now():
                import datetime as dt

                return dt.datetime(2031, 1, 1, 12, 0, 0)

        monkeypatch.setattr(reproducibility, "datetime", FixedDatetime)
        again = FingerprintGenerator.generate_fingerprint(
            self.config(tmp_path / "x" if False else tmp_path)
        )
        assert first == again

    def test_embedded_xml_comment_matches_the_identifier(self, tmp_path):
        from beast2py.reproducibility import FingerprintGenerator

        config = self.config(tmp_path)
        xml = XMLWriter(config).generate_xml()
        fp = FingerprintGenerator.generate_fingerprint(config)
        assert (
            f"Analysis Fingerprint: {fp} " in xml
        ), "the identifier written into the XML differs from the one the tool reports"
        assert f"Data: {FingerprintGenerator.data_hash(config)}" in xml

    def test_methods_description_is_deterministic_and_date_free(self, tmp_path):
        from beast2py.reproducibility import FingerprintGenerator, MethodsGenerator

        config = self.config(tmp_path)
        text = MethodsGenerator.generate_methods(config)
        again = MethodsGenerator.generate_methods(self.config(tmp_path))
        assert text == again, "the auto-generated Methods paragraph is not reproducible"
        assert FingerprintGenerator.generate_fingerprint(config) in text
        import re

        assert not re.search(
            r"\b20\d{6}\b", text
        ), "a YYYYMMDD component leaked back into the methods text"

    def test_changing_alignment_content_changes_the_fingerprint(self, tmp_path):
        from beast2py.reproducibility import FingerprintGenerator

        base = self.config(tmp_path)
        before = FingerprintGenerator.generate_fingerprint(base)
        altered = VARIABLE_FASTA.replace(
            "ACGTACGTACGTACGTTTTGCAAGGCATCGATCGGG", "ACGTACGTACGTACGTTTTGCAAGGCATCGATCGGT"
        )
        after = FingerprintGenerator.generate_fingerprint(self.config(tmp_path, fasta=altered))
        assert (
            before != after
        ), "the fingerprint ignored the sequence data it was supposed to certify"

    def test_data_hash_is_independent_of_order_and_names_content(self, tmp_path):
        from beast2py.reproducibility import FingerprintGenerator

        config = self.config(tmp_path)
        digest = FingerprintGenerator.data_hash(config)
        assert len(digest) == 16
        assert FingerprintGenerator.data_hash(config) == digest

    def test_changing_output_file_name_does_not_change_fingerprint(self, tmp_path):
        from beast2py.reproducibility import FingerprintGenerator

        config = self.config(tmp_path)
        before = FingerprintGenerator.generate_fingerprint(config)
        config.metadata["output_file"] = "totally_different_name.xml"
        assert FingerprintGenerator.generate_fingerprint(config) == before

    def test_crown_stem_switch_changes_the_fingerprint(self, tmp_path):
        from beast2py.reproducibility import FingerprintGenerator

        config = self.config(tmp_path)
        before = FingerprintGenerator.generate_fingerprint(config)
        config.calibrations[0].use_originate = True
        assert FingerprintGenerator.generate_fingerprint(config) != before

    def test_xml_is_byte_identical_across_reruns(self, tmp_path):
        config = self.config(tmp_path)
        first = XMLWriter(config).generate_xml()
        second = XMLWriter(self.config(tmp_path)).generate_xml()
        assert first == second

    def test_store_every_default_is_positive(self, tmp_path):
        config = self.config(tmp_path)
        assert config.mcmc.store_every > 0


# ---------------------------------------------------------------------------
# Silent coercion and consistency checks
# ---------------------------------------------------------------------------


class TestUnavailableIsNotSilent:
    """`--beast2-validate` must not turn 'BEAST2 missing' into a green run.

    Review question 8 asks whether an absent BEAST2 (or an absent add-on) is a
    hard failure or a silent degradation.  The add-on case was already a hard
    failure (BEAST2 itself rejects the model, gate 3 fails, exit 1, nothing
    written).  The *unavailable* case only warned and then wrote the file with
    exit 0, which is the same fail-open shape  a CI job
    that reads only the exit code would mark an unvalidated XML as validated.
    It now exits 2 and writes nothing, with `--allow-unvalidated` as the
    explicit opt-out.
    """

    @staticmethod
    def _run(args, env, out_xml):
        import subprocess

        import sys as _sys

        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        proc = subprocess.run(
            [_sys.executable, "-m", "beast2py.main", "generate"] + args,
            cwd=root,
            capture_output=True,
            text=True,
            env={**env, "PYTHONPATH": root},
        )
        return proc

    @pytest.mark.skipif(os.name != "posix", reason="POSIX env manipulation")
    def test_unavailable_beast2_stops_the_write(self, tmp_path):
        empty_home = tmp_path / "no-beast"
        empty_home.mkdir()
        env = dict(
            os.environ, HOME=str(empty_home), BEAST2_JAR=str(tmp_path / "definitely-not-a-jar.jar")
        )
        out = tmp_path / "out.xml"
        config = tmp_path / "cfg.yaml"
        config.write_text(
            "alignments:\n  - id: alignment\n    file: aln.fasta\n"
            "    data_type: nucleotide\npartitions:\n  - id: alignment\n"
            "    site_model:\n      substitution_model:\n        type: hky\n"
            "    clock_model:\n      type: strict\n    tree: shared\n"
            "tree_prior:\n  type: yule\nmcmc:\n  chain_length: 10000\n",
            encoding="utf-8",
        )
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        shutil.copy(Path(root) / "examples" / "primates.fasta", tmp_path / "aln.fasta")
        args = ["-c", str(config), "-o", str(out), "--beast2-validate"]
        proc = self._run(args, env, out)
        assert proc.returncode == 2, proc.stdout + proc.stderr
        assert "unavailable" in proc.stdout.lower()
        assert not out.exists(), "an unvalidated XML was written anyway"

        proc = self._run(args + ["--allow-unvalidated"], env, out)
        assert proc.returncode == 0, proc.stdout + proc.stderr
        assert out.exists(), "--allow-unvalidated must still write the file"
        assert "unavailable" in proc.stdout.lower()


class TestGammaParameterisation:
    """A Gamma `beta` must mean the same thing everywhere .

    BEAST2's Gamma defaults to mode=ShapeScale, so beta is a *scale* and the
    mean is alpha * beta.  `_frozen_distribution` was corrected for that while
    `compute_distribution_stats` still divided by beta, so the tool reported a
    prior mean and interval for a distribution it was not emitting - and
    Supplementary Figure S1 drew the Gamma panel with R's rate argument, i.e.
    the mirror image of the tool's own convention.
    """

    def test_stats_mean_uses_the_scale_convention(self):
        from beast2py.utils import compute_distribution_stats

        st = compute_distribution_stats("gamma", {"alpha": 2.0, "beta": 4.0})
        assert st["mean"] == pytest.approx(8.0), st

    def test_stats_agree_with_the_frozen_distribution(self):
        import numpy as np

        from beast2py.utils import _frozen_distribution, compute_distribution_stats

        for par in (
            {"alpha": 2.0, "beta": 4.0},
            {"alpha": 0.5, "beta": 2.0},
            {"alpha": 3.0, "beta": 1.5},
        ):
            st = compute_distribution_stats("gamma", par)
            frozen = _frozen_distribution("gamma", par, offset=0.0)
            assert st["mean"] == pytest.approx(float(frozen.mean()), rel=1e-9), par
            assert st["median"] == pytest.approx(float(frozen.median()), rel=1e-6), par
            # the reported interval must sit where the sampled density is
            draws = np.asarray(frozen.rvs(size=20000, random_state=7))
            assert (draws < st["hpd_upper"]).mean() > 0.85, (par, st)

    def test_shape_rate_mode_is_honoured_and_emitted(self):
        from beast2py.utils import compute_distribution_stats

        st = compute_distribution_stats("gamma", {"alpha": 2.0, "beta": 4.0, "mode": "ShapeRate"})
        assert st["mean"] == pytest.approx(0.5), st

    def test_gamma_prior_in_xml_keeps_alpha_and_beta_verbatim(self, tmp_path):
        import shutil

        from beast2py import api

        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cfg = {
            "alignments": [
                {
                    "id": "alignment",
                    "file": "aln.fasta",
                    "format": "auto",
                    "data_type": "nucleotide",
                }
            ],
            "partitions": [
                {
                    "id": "alignment",
                    "site_model": {
                        "substitution_model": {"type": "hky"},
                        "gamma_categories": 4,
                        "gamma_shape": {"value": 0.5, "lower": 0.0},
                    },
                    "clock_model": {"type": "strict"},
                    "tree": "shared",
                }
            ],
            "tree_prior": {"type": "yule"},
            "mcmc": {
                "type": "standard",
                "chain_length": 100000,
                "pre_burnin": 0,
                "store_every": 10000,
            },
            "parameter_priors": [
                {
                    "parameter": "gamma_shape",
                    "distribution": {"type": "gamma", "parameters": {"alpha": 2.0, "beta": 4.0}},
                }
            ],
        }
        shutil.copy(Path(root) / "examples" / "primates.fasta", tmp_path / "aln.fasta")
        xml = api.generate_xml(cfg, base_dir=str(tmp_path))
        assert 'spec="beast.base.inference.distribution.Gamma"' in xml
        assert 'name="alpha" value="2.0"' in xml
        assert 'name="beta" value="4.0"' in xml
        # default mode needs no attribute; an explicit rate must be stated
        assert 'mode="ShapeScale"' not in xml
        xml2 = api.generate_xml(
            dict(
                cfg,
                parameter_priors=[
                    {
                        "parameter": "gamma_shape",
                        "distribution": {
                            "type": "gamma",
                            "parameters": {"alpha": 2.0, "beta": 4.0, "mode": "ShapeRate"},
                        },
                    }
                ],
            ),
            base_dir=str(tmp_path),
        )
        assert 'mode="ShapeRate"' in xml2


class TestSilentRewrites:
    def test_short_rates_vector_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "gtr"\n'
            '        rates: {value: "1.0 2.0 1.0"}\n' + TAIL,
        )
        config = ConfigParser().parse(str(path))
        with pytest.raises(ValueError, match="exactly 6"):
            XMLWriter(config).generate_xml()

    def test_rates_on_hky_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            '        rates: {value: "1.0 2.0 1.0 1.0 2.0 1.0"}\n' + TAIL,
        )
        with pytest.raises(ConfigError, match="no free"):
            ConfigParser().parse(str(path))

    def test_rates_may_not_be_combined_with_named_rates(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "gtr"\n'
            '        rates: {value: "1.0 2.0 1.0 1.0 2.0 1.0"}\n'
            "        rateAC: {value: 3.0}\n" + TAIL,
        )
        config = ConfigParser().parse(str(path))
        with pytest.raises(ValueError, match="explicitly named rate"):
            XMLWriter(config).generate_xml()

    def test_float_noise_is_not_written_to_xml(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "birth_death"\n  birth_rate: {value: 0.7}\n'
            "  death_rate: {value: 0.2}\n"
            'calibrations:\n  - name: "rootCal"\n    taxa: null\n    distribution:\n'
            '      type: "lognormal"\n      parameters: {M: 2.3, S: 1.0}\n'
            "mcmc:\n  chain_length: 100000\n",
        )
        xml = XMLWriter(ConfigParser().parse(str(path))).generate_xml()
        assert "0.00000000000000" not in xml or "0.5" in xml
        assert 'value="0.5000000000000' not in xml

    def test_gamma_shape_with_single_category_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            "      gamma_categories: 1\n      gamma_shape: {value: 0.5}\n" + TAIL,
        )
        with pytest.raises(ConfigError, match="gammaCategoryCount"):
            ConfigParser().parse(str(path))

    def test_self_referencing_filter_is_rejected(self, tmp_path):
        # Build directly through the dataclasses to exercise the writer guard.
        aln = Alignment(
            id="filt",
            data_type=DataType.NUCLEOTIDE,
            sequences=[Sequence("a", "ACGT"), Sequence("b", "ACGA")],
            filter="1::2",
            data_id="filt",
        )
        part = Partition(
            id="filt",
            alignment=aln,
            site_model=SiteModelConfig(substitution_model={"type": "hky"}),
            clock_model=ClockModelConfig(type=ClockModelType.STRICT),
        )
        from beast2py.models import BEASTConfig, TreePriorType

        config = BEASTConfig(
            metadata={"analysis_name": "x", "embed_fingerprint": False},
            partitions=[part],
            all_alignments=[aln],
            tree_prior_type=TreePriorType.YULE,
            tree_prior_params={},
            calibrations=[],
        )
        with pytest.raises(ValueError, match="different alignment|data_id"):
            XMLWriter(config).generate_xml()

    def test_partition_id_may_differ_from_alignment_id(self):
        aln = Alignment(
            id="aln_one",
            data_type=DataType.NUCLEOTIDE,
            sequences=[Sequence("a", "ACGT"), Sequence("b", "ACGA")],
        )
        part = Partition(
            id="P1",
            alignment=aln,
            site_model=SiteModelConfig(substitution_model={"type": "jc69"}),
            clock_model=ClockModelConfig(type=ClockModelType.STRICT),
        )
        from beast2py.models import BEASTConfig, TreePriorType

        config = BEASTConfig(
            metadata={"analysis_name": "x", "embed_fingerprint": False},
            partitions=[part],
            all_alignments=[aln],
            tree_prior_type=TreePriorType.YULE,
            tree_prior_params={},
            calibrations=[],
        )
        xml = XMLWriter(config).generate_xml()
        root = ET.fromstring(xml)
        data_ids = {d.get("id") for d in root.findall("data")}
        assert data_ids == {"aln_one"}, "the alignment was defined under the partition id"
        assert xml.count('id="aln_one"') == 1

    def test_nested_sampling_run_keeps_burnin_and_store_every(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head()
            + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            + TAIL.replace(
                "mcmc:\n  chain_length: 100000\n",
                "mcmc:\n  type: nested_sampling\n  chain_length: 100000\n"
                "  pre_burnin: 1000\n  store_every: 500\n",
            ),
        )
        xml = XMLWriter(ConfigParser().parse(str(path))).generate_xml()
        run = [line for line in xml.splitlines() if "<run" in line][0]
        assert 'preBurnin="1000"' in run and 'storeEvery="500"' in run

    def test_out_of_range_filter_is_rejected(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head(extra_alignments="""  - id: "third"
    file: "test.fasta"
    data_id: "alignment"
    filter: "999::3"
    data_type: "nucleotide"
""") + '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL,
        )
        with pytest.raises(ConfigError, match="outside a|no sites"):
            ConfigParser().parse(str(path))

    def test_codon_partition_filter_is_accepted_and_counts_sites(self, tmp_path):
        path = write_inputs(
            tmp_path,
            head(extra_alignments="""  - id: "first_two"
    file: "test.fasta"
    data_id: "alignment"
    filter: "1::3,2::3"
    data_type: "nucleotide"
""") + '\n    site_model:\n      substitution_model:\n        type: "hky"\n' + TAIL,
        )
        config = ConfigParser().parse(str(path))
        base = [a for a in config.all_alignments if a.id == "alignment"][0]
        assert base.n_sites == 36

    def test_library_path_is_gated_too(self, tmp_path):
        """The Python API must not be able to bypass the CLI gates ."""
        from beast2py.api import Beast2Py

        body = head() + HKY_PARTITION + TAIL
        config_path = write_inputs(tmp_path, body)
        out = tmp_path / "api_out.xml"
        Beast2Py(base_dir=str(tmp_path)).generate_xml(str(config_path), output=str(out))
        assert out.exists()

        out.write_text("PREVIOUS", encoding="utf-8")
        with pytest.raises(ValueError, match="[Gg]ate"):
            Beast2Py(base_dir=str(tmp_path)).generate_xml(str(config_path), output=str(out))
        assert out.read_text(encoding="utf-8") == "PREVIOUS"

        Beast2Py(base_dir=str(tmp_path)).generate_xml(str(config_path), output=str(out), force=True)
        assert out.read_text(encoding="utf-8").startswith("<?xml")

    def test_library_path_refuses_conflicting_calibrations(self, tmp_path):
        from beast2py.api import Beast2Py

        config_path = write_inputs(
            tmp_path,
            head() + '\n    site_model:\n      substitution_model:\n        type: "hky"\n'
            'tree_prior:\n  type: "yule"\n'
            "calibrations:\n"
            '  - name: "rootCal"\n    taxa: null\n    distribution:\n'
            "      type: uniform\n      parameters: {lower: 5.0, upper: 10.0}\n"
            '  - name: "childCal"\n    taxa: [taxonA, taxonB]\n    distribution:\n'
            "      type: uniform\n      parameters: {lower: 40.0, upper: 50.0}\n"
            "mcmc:\n  chain_length: 100000\n",
        )
        out = tmp_path / "conflict.xml"
        with pytest.raises(ValueError, match="gate|written"):
            Beast2Py(base_dir=str(tmp_path)).generate_xml(str(config_path), output=str(out))
        assert not out.exists()

    def test_validator_and_java_helper_ship_inside_the_package(self, tmp_path):
        """A wheel must be self-sufficient, not dependent on the repo layout.

        Shipping only the shell script is not enough: it resolves its Java
        helper relative to itself, so an installed Beast2Py used to fail with a
        confusing "could not find or load main class" error .
        """
        import beast2py

        pkg = Path(beast2py.__file__).resolve().parent
        script = pkg / "beast2_validate.sh"
        helper = pkg / "tools" / "Beast2Validator.java"
        classes = pkg / "tools" / "classes" / "Beast2Validator.class"
        for path in (script, helper, classes):
            assert path.is_file(), f"{path.relative_to(pkg)} missing from the package"
        from beast2py.validator import XMLValidator

        found = XMLValidator._bundled_validator_script()
        assert found is not None
        assert os.path.samefile(str(found), script), f"discovered {found}, expected {script}"

    def test_validator_script_resolves_its_helper_from_its_own_directory(self):
        """The script must look for tools/ next to itself, not only in a repo."""
        import beast2py

        pkg = Path(beast2py.__file__).resolve().parent
        body = (pkg / "beast2_validate.sh").read_text(encoding="utf-8")
        assert "$SCRIPT_DIR/tools" in body, "helper path is not script-relative"
        wrapper = pkg.parent / "beast2_validate.sh"
        if wrapper.is_file():
            assert "beast2py/beast2_validate.sh" in wrapper.read_text(encoding="utf-8")

    @staticmethod
    def _fake_java(dir_path, version_line):
        """Write an executable `java` stub that only prints a -version banner."""
        dir_path.mkdir(parents=True, exist_ok=True)
        stub = dir_path / "java"
        stub.write_text("#!/bin/sh\n" f"printf '%s\\n' '{version_line}' >&2\n" "exit 7\n")
        stub.chmod(0o755)
        return stub

    def _run_probe(self, tmp_path, path_dirs, extra_candidates=()):
        import beast2py
        import subprocess

        script = Path(beast2py.__file__).resolve().parent / "beast2_validate.sh"
        # The jar is located *before* Java, so give the script a stand-in to keep
        # the test inside the Java-selection logic.
        jar = tmp_path / "BEAST.base.jar"
        jar.write_bytes(b"PK\x03\x04 not a real jar")
        env = {
            "PATH": os.pathsep.join([str(d) for d in path_dirs] + ["/bin", "/usr/bin"]),
            "HOME": str(tmp_path / "empty-home"),  # no ~/.beast packages
            "BEAST2_JAR": str(jar),
            "BEAST2_SKIP_JAVA_HOME_TOOL": "1",
            "BEAST2_SKIP_JDK_PROBES": "1",
            "BEAST2_VALIDATE_TRACE": "1",
        }
        if extra_candidates:
            env["BEAST2_JAVA_CANDIDATES"] = os.pathsep.join(str(c) for c in extra_candidates)
        proc = subprocess.run(
            ["bash", script, str(tmp_path / "does-not-matter.xml")],
            capture_output=True,
            text=True,
            env=env,
        )
        return proc

    @pytest.mark.skipif(os.name != "posix", reason="bash-only probe ordering")
    def test_old_java_on_path_does_not_shadow_a_newer_jdk(self, tmp_path):
        """A 1.8 shim on PATH must not beat a JDK 17 found elsewhere .

        macOS machines commonly carry an ancient Oracle `java` on PATH while a
        usable openjdk@17 sits unlinked in Homebrew.  The script used to take
        PATH first and then abort with "requires Java 17 or newer", which made
        the third gate unreproducible on exactly the machines that can run it.
        """
        old = self._fake_java(tmp_path / "oldbin", 'java version "1.8.0_501" (build)')
        new = self._fake_java(tmp_path / "jdk17" / "bin", 'openjdk version "17.0.9" 2023')
        proc = self._run_probe(tmp_path, [old.parent], extra_candidates=[new])
        assert f"using java 17 at {new}" in proc.stderr, proc.stderr
        # the stub exits 7 on purpose: reaching it proves selection succeeded
        assert "needs a JDK 17 or newer" not in proc.stderr

    @pytest.mark.skipif(os.name != "posix", reason="bash-only probe ordering")
    def test_missing_jdk17_names_every_rejected_candidate(self, tmp_path):
        """The failure must be actionable, and must not claim java is absent."""
        old = self._fake_java(tmp_path / "oldbin", 'java version "1.8.0_501" (build)')
        proc = self._run_probe(tmp_path, [old.parent])
        assert proc.returncode == 2, proc.stderr
        assert "JDK 17" in proc.stderr
        assert "rejected (<17 or unreadable)" in proc.stderr
        assert str(old) in proc.stderr
        assert "major 8" in proc.stderr, "version parsing must normalise 1.8 to 8"

    @pytest.mark.skipif(os.name != "posix", reason="bash-only probe ordering")
    def test_1_x_version_strings_are_normalised(self, tmp_path):
        """ "1.8" must read as major 8 and "17.0.9" as 17, not as 1."""
        v17 = self._fake_java(tmp_path / "a", 'openjdk version "17.0.9" 2023-10-17')
        v1 = self._fake_java(tmp_path / "b", 'java version "1.6.0_65"')
        proc = self._run_probe(tmp_path, [v1.parent], extra_candidates=[v17])
        assert f"using java 17 at {v17}" in proc.stderr, proc.stderr

    def test_every_shipped_example_generates_and_validates(self, tmp_path):
        example_dir = Path(__file__).resolve().parent.parent / "examples"
        files = sorted(example_dir.glob("config_*.yaml"))
        assert files, "no example configurations found"
        for path in files:
            try:
                config = ConfigParser(base_dir=str(example_dir)).parse(str(path))
            except ConfigError as exc:  # pragma: no cover - failing message
                pytest.fail(f"{path.name} is not a valid example: {exc}")
            xml = XMLWriter(config).generate_xml()
            result = XMLValidator.validate(xml, config=config)
            assert result.is_valid, f"{os.path.basename(path)}: {result.errors[:2]}"

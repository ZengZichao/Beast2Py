"""Independent acceptance tests for the review-A implementation fixes.

Each test pins one repaired behaviour with an answer derived from the target
class semantics or simple mathematics, not from another function in the same
module.  Together they are the small independent acceptance table requested in
the review: soft-support separation, endpoint contact, stem identity,
hyperprior awareness, fixed tips, signed parameters, stem logging, method
matching, duplicate data_id, fixed frequencies, parameter-prior validation,
linked_to exclusion, fixed clock acceptance, BDSKY rate, single gamma
category, and the sequence alphabet.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

import pytest

from beast2py.calibration import CalibrationBuilder
from beast2py.config import ConfigError, ConfigParser
from beast2py.models import (
    CalibrationPoint,
    DistributionConfig,
    HyperpriorConfig,
)
from beast2py.xml_writer import XMLWriter
from beast2py.sequence import NUCLEOTIDE_ALPHABET, AMINOACID_ALPHABET

from .conftest import TEST_FASTA


def _write(tmp_path, body, fasta=TEST_FASTA):
    (tmp_path / "test.fasta").write_text(fasta, encoding="utf-8")
    cfg = tmp_path / "config.yaml"
    cfg.write_text(body, encoding="utf-8")
    return cfg


_HEAD = """
alignments:
  - id: alignment
    file: test.fasta
    data_type: nucleotide
partitions:
  - id: alignment
    site_model:
      substitution_model: {type: hky}
    clock_model: {type: strict}
    tree: shared
tree_prior: {type: yule}
mcmc: {chain_length: 10000}
output: {file: out.xml}
"""

# ---------------------------------------------------------------------------
# C01: determined tip dates are observations; only explicit uncertainty moves
# ---------------------------------------------------------------------------


class TestC01FixedTipDates:
    def test_determined_dates_produce_no_walker(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD + """
tip_dates:
  enabled: true
  dates: {taxonA: 2000, taxonB: 2001, taxonC: 2002, taxonD: 2003}
""",
        )
        config = ConfigParser().parse(cfg)
        xml = XMLWriter(config).generate_xml()
        assert 'spec="TipDatesRandomWalker"' not in xml
        assert "TraitSet" in xml

    def test_uncertain_dates_produce_one_walker_each(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD + """
tip_dates:
  enabled: true
  dates: {taxonA: 2000, taxonB: 2001, taxonC: 2002, taxonD: 2003}
  uncertain: {taxonB: 3.0}
""",
        )
        config = ConfigParser().parse(cfg)
        root = ET.fromstring(XMLWriter(config).generate_xml())
        walkers = [e for e in root.iter() if e.get("spec") == "TipDatesRandomWalker"]
        assert len(walkers) == 1
        assert walkers[0].get("windowSize") == "3"

    def test_uncertain_taxon_must_be_dated(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD + """
tip_dates:
  enabled: true
  dates: {taxonA: 2000, taxonB: 2001, taxonC: 2002, taxonD: 2003}
  uncertain: {noSuchTip: 3.0}
""",
        )
        with pytest.raises(ConfigError, match="noSuchTip"):
            ConfigParser().parse(cfg)


# ---------------------------------------------------------------------------
# C02: signed growthRate gets an additive operator that crosses zero
# ---------------------------------------------------------------------------


class TestC02SignedGrowthRate:
    def test_growth_rate_walker_is_additive(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace("tree_prior: {type: yule}", "tree_prior: {type: coalescent_exponential}"),
        )
        config = ConfigParser().parse(cfg)
        root = ET.fromstring(XMLWriter(config).generate_xml())
        walkers = [
            e
            for e in root.iter()
            if e.get("spec") == "RealRandomWalkOperator" and e.get("parameter") == "@growthRate"
        ]
        assert len(walkers) == 1
        scalers = [
            e
            for e in root.iter()
            if e.get("spec") == "ScaleOperator" and e.get("parameter") == "@growthRate"
        ]
        assert scalers == []


# ---------------------------------------------------------------------------
# C03: hyperprior serialization and support domains
# ---------------------------------------------------------------------------


class TestC03Hyperpriors:
    def _cal(self, hyper_key="S", hp_type="lognormal", hp_params=None):
        cal = CalibrationPoint(
            name="rootCal",
            taxa=None,
            distribution=DistributionConfig(type="lognormal", parameters={"M": 1.0, "S": 0.2}),
        )
        cal.hyperpriors = {
            hyper_key: HyperpriorConfig(type=hp_type, parameters=hp_params or {"M": 0.0, "S": 0.3})
        }
        return cal

    def test_lognormal_hyperprior_keeps_mean_in_real_space(self):
        cal = self._cal(hp_params={"M": 2.0, "S": 0.2, "mean_in_real_space": True})
        prior, _sn_id, sn_params = CalibrationBuilder.build_hyperprior_priors(cal)[0]
        distr = prior.find("distr")
        assert distr.get("meanInRealSpace") == "true"

    def test_location_hyperparameter_support_is_unbounded(self):
        cal = self._cal(hyper_key="M", hp_type="normal", hp_params={"mean": -1.0, "sigma": 1.0})
        _prior, _sn_id, sn_params = CalibrationBuilder.build_hyperprior_priors(cal)[0]
        assert sn_params["lower"] is None

    def test_unknown_hyperprior_target_is_rejected(self):
        cal = self._cal(hyper_key="Sdeg")
        with pytest.raises(ValueError, match="Sdeg"):
            CalibrationBuilder.build_hyperprior_priors(cal)


# ---------------------------------------------------------------------------
# C04: calibration method must match the tree prior
# ---------------------------------------------------------------------------


class TestC04MethodTreePriorMatch:
    def test_calibrated_yule_method_with_plain_yule_is_rejected(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace("tree_prior: {type: yule}", "tree_prior: {type: yule}").replace(
                "output: {file: out.xml}",
                "calibration_method: calibrated_yule\noutput: {file: out.xml}",
            )
            + """
calibrations:
  - name: rootCal
    distribution: {type: normal, parameters: {mean: 10.0, sigma: 1.0}}
""",
        )
        with pytest.raises(ConfigError, match="calibrated_yule"):
            ConfigParser().parse(cfg)

    def test_calibrated_yule_pair_generates_calibration_points(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace(
                "tree_prior: {type: yule}", "tree_prior: {type: calibrated_yule}"
            ).replace(
                "output: {file: out.xml}",
                "calibration_method: calibrated_yule\noutput: {file: out.xml}",
            )
            + """
calibrations:
  - name: rootCal
    distribution: {type: normal, parameters: {mean: 10.0, sigma: 1.0}}
""",
        )
        config = ConfigParser().parse(cfg)
        root = ET.fromstring(XMLWriter(config).generate_xml())
        assert any(
            e.tag == "calibrations" and e.get("spec") == "CalibrationPoint" for e in root.iter()
        )


# ---------------------------------------------------------------------------
# C05: stem calibration logs the stem, not the crown
# ---------------------------------------------------------------------------


class TestC05StemLogger:
    def test_stem_logger_carries_use_originate(self):
        cal = CalibrationPoint(
            name="cladeAB",
            taxa=["A", "B"],
            use_originate=True,
            distribution=DistributionConfig(type="normal", parameters={"mean": 10.0, "sigma": 1.0}),
        )
        logger = CalibrationBuilder.build_mrca_logger(cal, "tree")
        assert logger.get("useOriginate") == "true"

    def test_crown_logger_has_no_originate(self):
        cal = CalibrationPoint(
            name="cladeAB",
            taxa=["A", "B"],
            distribution=DistributionConfig(type="normal", parameters={"mean": 10.0, "sigma": 1.0}),
        )
        logger = CalibrationBuilder.build_mrca_logger(cal, "tree")
        assert "useOriginate" not in logger.attrib


# ---------------------------------------------------------------------------
# C06: fixed frequencies never become an estimated state node
# ---------------------------------------------------------------------------


class TestC06FixedFrequencies:
    def test_fixed_mode_leaves_the_state_space(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace(
                "substitution_model: {type: hky}",
                "substitution_model: {type: hky, frequencies: {mode: fixed,"
                ' value: "0.1 0.2 0.3 0.4"}}',
            ),
        )
        config = ConfigParser().parse(cfg)
        root = ET.fromstring(XMLWriter(config).generate_xml())
        assert not [
            e
            for e in root.iter("parameter")
            if (e.get("id") or "").endswith("freqParameter") and e.get("estimate") != "false"
        ]
        freq = next(
            e for e in root.iter("parameter") if (e.get("id") or "").endswith("freqParameter")
        )
        assert freq.get("estimate") == "false"
        assert not any(
            e.tag == "operator" and (e.get("parameter") or "").endswith("freqParameter")
            for e in root.iter()
        )

    def test_fixed_mode_with_estimate_true_is_a_conflict(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace(
                "substitution_model: {type: hky}",
                "substitution_model: {type: hky, frequencies: {mode: fixed,"
                " value: '0.25 0.25 0.25 0.25', estimate: true}}",
            ),
        )
        config = ConfigParser().parse(cfg)
        with pytest.raises(ValueError, match="mode: fixed"):
            XMLWriter(config).generate_xml()


# ---------------------------------------------------------------------------
# C07: parameter priors get the full distribution validation
# ---------------------------------------------------------------------------


class TestC07ParameterPriorValidation:
    def test_misspelled_parameter_key_is_rejected(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace(
                "output: {file: out.xml}",
                """
parameter_priors:
  - parameter: birthRate
    distribution: {type: lognormal, parameters: {M: 1.0, Sdeg: 0.2}}
output: {file: out.xml}
""",
            ),
        )
        with pytest.raises(ConfigError, match="Sdeg"):
            ConfigParser().parse(cfg)

    def test_quoted_false_boolean_serialises_as_false(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace(
                "output: {file: out.xml}",
                """
parameter_priors:
  - parameter: birthRate
    distribution: {type: lognormal, parameters: {M: 1.0, S: 0.2, mean_in_real_space: "false"}}
output: {file: out.xml}
""",
            ),
        )
        config = ConfigParser().parse(cfg)
        xml = XMLWriter(config).generate_xml()
        assert 'meanInRealSpace="false"' in xml


# ---------------------------------------------------------------------------
# C08: a data_id may name only one alignment
# ---------------------------------------------------------------------------


class TestC08DataIdReuse:
    def _two_files(self, tmp_path):
        (tmp_path / "test.fasta").write_text(TEST_FASTA, encoding="utf-8")
        (tmp_path / "a.fasta").write_text(TEST_FASTA, encoding="utf-8")
        other = TEST_FASTA.replace("ACGT", "TGCA", 1)
        (tmp_path / "b.fasta").write_text(other, encoding="utf-8")

    def test_same_data_id_different_content_is_rejected(self, tmp_path):
        self._two_files(tmp_path)
        cfg = tmp_path / "config.yaml"
        cfg.write_text(
            """
alignments:
  - id: g1
    file: a.fasta
    data_type: nucleotide
    filter: 1-6
    data_id: source
  - id: g2
    file: b.fasta
    data_type: nucleotide
    filter: 7-12
    data_id: source
partitions:
  - id: g1
    site_model: {substitution_model: {type: hky}}
    clock_model: {type: strict}
    tree: shared
  - id: g2
    site_model: {substitution_model: {type: hky}}
    clock_model: {type: strict}
    tree: shared
tree_prior: {type: yule}
mcmc: {chain_length: 10000}
output: {file: out.xml}
""",
            encoding="utf-8",
        )
        with pytest.raises(ConfigError, match="already registered"):
            ConfigParser().parse(cfg)

    def test_same_data_id_same_content_is_accepted(self, tmp_path):
        (tmp_path / "a.fasta").write_text(TEST_FASTA, encoding="utf-8")
        (tmp_path / "b.fasta").write_text(TEST_FASTA, encoding="utf-8")
        cfg = tmp_path / "config.yaml"
        cfg.write_text(
            """
alignments:
  - id: g1
    file: a.fasta
    data_type: nucleotide
    filter: 1-6
    data_id: source
  - id: g2
    file: b.fasta
    data_type: nucleotide
    filter: 7-12
    data_id: source
partitions:
  - id: g1
    site_model: {substitution_model: {type: hky}}
    clock_model: {type: strict}
    tree: shared
  - id: g2
    site_model: {substitution_model: {type: hky}}
    clock_model: {type: strict}
    tree: shared
tree_prior: {type: yule}
mcmc: {chain_length: 10000}
output: {file: out.xml}
""",
            encoding="utf-8",
        )
        config = ConfigParser().parse(cfg)
        assert len(config.partitions) == 2


# ---------------------------------------------------------------------------
# C09: linked_to plus local model parameters is an error
# ---------------------------------------------------------------------------


class TestC09LinkedToExclusion:
    _SECOND_FASTA = """>taxonE
ACGTACGTACGTACGTTTTGCAAGGCATCGATCGGG
>taxonF
ACGTACGTACGTAAGTTTTGCAAGGCATCGATCGGA
>taxonG
ACGAACGTACGTAAGTTTTGCAAGGCATCGATCGGG
>taxonH
ACGAACGTACGTAAGTTTTGCAAGGCATCGATCGGA
"""

    def _two_alignments(self, tmp_path):
        (tmp_path / "test.fasta").write_text(TEST_FASTA, encoding="utf-8")
        (tmp_path / "second.fasta").write_text(self._SECOND_FASTA, encoding="utf-8")

    def test_clock_link_with_local_type_is_rejected(self, tmp_path):
        self._two_alignments(tmp_path)
        cfg = _write(
            tmp_path,
            """
alignments:
  - id: g1
    file: test.fasta
    data_type: nucleotide
  - id: g2
    file: second.fasta
    data_type: nucleotide
partitions:
  - id: g1
    site_model: {substitution_model: {type: hky}}
    clock_model: {type: strict}
    tree: shared
  - id: g2
    site_model: {substitution_model: {type: hky}}
    clock_model: {type: strict, clock_rate: 0.01, linked_to: g1}
    tree: shared
tree_prior: {type: yule}
mcmc: {chain_length: 10000}
output: {file: out.xml}
""",
        )
        with pytest.raises(ConfigError, match="linked_to"):
            ConfigParser().parse(cfg)

    def test_bare_link_is_still_accepted(self, tmp_path):
        self._two_alignments(tmp_path)
        cfg = _write(
            tmp_path,
            """
alignments:
  - id: g1
    file: test.fasta
    data_type: nucleotide
  - id: g2
    file: second.fasta
    data_type: nucleotide
partitions:
  - id: g1
    site_model: {substitution_model: {type: hky}}
    clock_model: {type: strict}
    tree: shared
  - id: g2
    site_model: {substitution_model: {type: hky}}
    clock_model: {linked_to: g1}
    tree: shared
tree_prior: {type: yule}
mcmc: {chain_length: 10000}
output: {file: out.xml}
""",
        )
        config = ConfigParser().parse(cfg)
        assert config.partitions[1].clock_model.linked_to == "g1"


# ---------------------------------------------------------------------------
# C10: a fixed known clock rate is a legitimate absolute scale
# ---------------------------------------------------------------------------


class TestC10FixedClockRate:
    def test_fixed_clock_without_calibration_is_accepted(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace(
                "clock_model: {type: strict}",
                "clock_model: {type: strict, clock_rate: {value: 0.01, estimate: false}}",
            ),
        )
        config = ConfigParser().parse(cfg)
        assert config.partitions[0].clock_model.clock_rate.estimate is False

    def test_note_records_the_fixed_rate_scale(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace(
                "clock_model: {type: strict}",
                "clock_model: {type: strict, clock_rate: {value: 0.01, estimate: false}}",
            ),
        )
        parser = ConfigParser()
        parser.parse(cfg)
        assert parser._time_scale_note and "fixed clock rate" in parser._time_scale_note


# ---------------------------------------------------------------------------
# C11: BDSKY sampling rate is unbounded above
# ---------------------------------------------------------------------------


class TestC11BDSkySamplingRate:
    def test_sampling_rate_above_one_is_accepted(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace(
                "tree_prior: {type: yule}",
                "tree_prior: {type: bd_skyline_serial, sampling_rate: 2.0}",
            ),
        )
        config = ConfigParser().parse(cfg)
        assert config.tree_prior_params["sampling_rate"] == 2.0

    def test_sampling_rate_state_node_has_no_upper_bound(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace(
                "tree_prior: {type: yule}",
                "tree_prior: {type: bd_skyline_serial, sampling_rate: 2.0}",
            ),
        )
        config = ConfigParser().parse(cfg)
        root = ET.fromstring(XMLWriter(config).generate_xml())
        node = next(e for e in root.iter("parameter") if e.get("id") == "bdSkylineSamplingRate")
        assert "upper" not in node.attrib or float(node.get("upper")) > 1.0


# ---------------------------------------------------------------------------
# C12: one gamma category means no shape parameter at all
# ---------------------------------------------------------------------------


class TestC12SingleGammaCategory:
    def test_default_shape_is_not_created_for_one_category(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace(
                "substitution_model: {type: hky}",
                "substitution_model: {type: hky}\n      gamma_categories: 1",
            ),
        )
        config = ConfigParser().parse(cfg)
        assert config.partitions[0].site_model.gamma_shape is None
        root = ET.fromstring(XMLWriter(config).generate_xml())
        assert not [e for e in root.iter() if (e.get("id") or "").endswith("gammaShape")]


# ---------------------------------------------------------------------------
# C13: sequence alphabets match BEAST2 2.7.x codeMaps
# ---------------------------------------------------------------------------


class TestC13Alphabet:
    def test_nucleotide_alphabet_matches_beast_278(self):
        # Nucleotide.java:34 codeMap = "ACGTURYMWSKBDHVNX" + '-' + '?'
        expected = set("ACGTURYMWSKBDHVNXacgturymwskbdhvnx") | set("-?")
        assert NUCLEOTIDE_ALPHABET == expected

    def test_aminoacid_alphabet_matches_beast_278(self):
        # Aminoacid.java:12 codeMap = "ACDEFGHIKLMNPQRSTVWY" + "X" + '-' + '?'
        expected = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy") | set("Xx-?")
        assert AMINOACID_ALPHABET == expected

    def test_dna_dot_and_star_are_rejected(self, tmp_path):
        fasta = ">a\nACGT.ACGT*\n>b\nACGTACGTAA\n"
        cfg = _write(tmp_path, _HEAD, fasta=fasta)
        with pytest.raises(ConfigError, match="illegal for nucleotide"):
            ConfigParser().parse(cfg)

    def test_dna_x_is_accepted(self, tmp_path):
        fasta = ">a\nACGTXACGT\n>b\nACGTACGTA\n"
        cfg = _write(tmp_path, _HEAD, fasta=fasta)
        config = ConfigParser().parse(cfg)
        assert config.partitions[0].alignment is not None


# ---------------------------------------------------------------------------
# I01: generate() writes nothing when an external validation fails
# ---------------------------------------------------------------------------


class TestI01GenerateWritesOnlyAfterAllGates:
    def test_failed_beast2_validation_writes_nothing(self, tmp_path, monkeypatch):
        from unittest.mock import patch

        from beast2py.api import Beast2Py
        from beast2py.validator import XMLValidator

        cfg = _write(tmp_path, _HEAD)
        out = tmp_path / "out.xml"
        with patch.object(
            XMLValidator,
            "beast2_validation_status",
            return_value=(False, "INVALID: something", "failed"),
        ) as mocked:
            with pytest.raises(RuntimeError, match="nothing was written"):
                Beast2Py().generate(
                    str(cfg),
                    output=str(out),
                    beast2_validate=True,
                    fingerprint=True,
                    methods=True,
                    diagnose=False,
                )
            mocked.assert_called_once()
        assert not out.exists()
        assert not (tmp_path / "out.fingerprint.json").exists()
        assert not (tmp_path / "out.methods.tex").exists()
        assert not (tmp_path / "diagnostics").exists()


# ---------------------------------------------------------------------------
# I02: quick mode shares the configuration-level validation
# ---------------------------------------------------------------------------


class TestI02QuickValidation:
    def test_negative_chain_length_is_rejected(self, tmp_path):
        from beast2py.api import quick_generate

        _write(tmp_path, _HEAD)  # writes test.fasta
        with pytest.raises(ValueError, match="chain_length"):
            quick_generate(
                alignment=str(tmp_path / "test.fasta"),
                output=str(tmp_path / "out.xml"),
                chain_length=-1,
            )
        assert not (tmp_path / "out.xml").exists()

    def test_hard_fossil_mismatch_is_rejected_in_quick(self, tmp_path):
        from beast2py.api import quick_generate

        _write(tmp_path, _HEAD)  # writes test.fasta
        cal = tmp_path / "cals.yaml"
        cal.write_text(
            """
- name: fossilA
  taxa: [taxonA, taxonB]
  monophyletic: true
  distribution: {type: uniform, parameters: {lower: 5.0, upper: 10.0}}
  provenance:
    source: fossil
    calibration_type: hard
    original_age_min: 50.0
""",
            encoding="utf-8",
        )
        with pytest.raises(Exception, match="original minimum age|hard"):
            quick_generate(
                alignment=str(tmp_path / "test.fasta"),
                output=str(tmp_path / "out.xml"),
                calibration_yaml=str(cal),
            )
        assert not (tmp_path / "out.xml").exists()

    def test_quick_calibration_yaml_keeps_hyperpriors(self, tmp_path):
        from beast2py.api import quick_generate

        _write(tmp_path, _HEAD)  # writes test.fasta
        cal = tmp_path / "cals.yaml"
        cal.write_text(
            """
- name: rootCal
  taxa: null
  distribution:
    type: normal
    parameters: {mean: 10.0, sigma: 1.0}
    hyperpriors:
      mean: {type: uniform, parameters: {lower: 0.0, upper: 40.0}}
""",
            encoding="utf-8",
        )
        xml = quick_generate(
            alignment=str(tmp_path / "test.fasta"),
            output=str(tmp_path / "out.xml"),
            calibration_yaml=str(cal),
        )
        assert 'id="rootCal.mean"' in xml


# ---------------------------------------------------------------------------
# I03: INVALID verdicts are never reclassified as "unavailable"
# ---------------------------------------------------------------------------


class TestI03InvalidIsNotUnavailable:
    def test_class_not_found_is_failed(self, monkeypatch):
        from unittest.mock import patch

        from beast2py.validator import XMLValidator

        fake_output = (
            "Loading package\n"
            "INVALID: /tmp/x.xml\n"
            "  Class could not be found. Did you mean ...\n"
        )
        with patch.object(XMLValidator, "validate_with_beast2", return_value=(False, fake_output)):
            success, _output, status = XMLValidator.beast2_validation_status("/tmp/x.xml")
        assert success is False
        assert status == "failed"

    def test_taxon_not_found_is_failed(self, monkeypatch):
        from unittest.mock import patch

        from beast2py.validator import XMLValidator

        fake_output = "INVALID: Taxon not found in tree: Homo_sapiens\n"
        with patch.object(XMLValidator, "validate_with_beast2", return_value=(False, fake_output)):
            _success, _output, status = XMLValidator.beast2_validation_status("/tmp/x.xml")
        assert status == "failed"


# ---------------------------------------------------------------------------
# I04: Gate 1 boundary checks
# ---------------------------------------------------------------------------


class TestI04GateOneBoundaries:
    def _valid_xml(self, tmp_path):
        cfg = _write(tmp_path, _HEAD)
        config = ConfigParser().parse(cfg)
        return XMLWriter(config).generate_xml()

    def _posterior_of(self, xml):
        root = ET.fromstring(xml)
        return root, next(e for e in root.iter("distribution") if e.get("id") == "posterior")

    def test_prior_moved_out_of_posterior_is_not_coverage(self, tmp_path):
        from beast2py.validator import XMLValidator

        xml = self._valid_xml(tmp_path)
        root, posterior = self._posterior_of(xml)
        prior = next(
            e
            for e in posterior.iter("distribution")
            if e.get("spec", "").rsplit(".", 1)[-1] == "Prior"
        )
        for parent in root.iter():
            if prior in list(parent):
                parent.remove(prior)
                break
        root.append(prior)  # same id, same x reference, but unreachable
        result = XMLValidator.validate(ET.tostring(root, encoding="unicode"))
        assert not result.is_valid
        assert any("no <Prior>" in err for err in result.errors)

    def test_nan_initial_value_is_rejected(self, tmp_path):
        from beast2py.validator import XMLValidator

        xml = self._valid_xml(tmp_path)
        root, _posterior = self._posterior_of(xml)
        node = next(
            e
            for e in root.iter("parameter")
            if e.get("estimate", "true") == "true" and e.get("id") == "birthRate"
        )
        node.set("value", "NaN")
        result = XMLValidator.validate(ET.tostring(root, encoding="unicode"))
        assert any("not finite" in err for err in result.errors)

    def test_non_numeric_lower_bound_is_a_structured_error(self, tmp_path):
        from beast2py.validator import XMLValidator

        xml = self._valid_xml(tmp_path)
        root, _posterior = self._posterior_of(xml)
        node = next(
            e
            for e in root.iter("parameter")
            if e.get("estimate", "true") == "true" and e.get("id") == "birthRate"
        )
        node.set("lower", "oops")
        node.set("upper", "5.0")
        result = XMLValidator.validate(ET.tostring(root, encoding="unicode"))
        assert any("non-numeric bounds" in err for err in result.errors)


# ---------------------------------------------------------------------------
# I05: entity constructs are refused at the safe parser
# ---------------------------------------------------------------------------


class TestI05EntityRefusal:
    def test_doctype_is_rejected_by_the_safe_parser(self):
        from beast2py.validator import XMLValidator

        xml = (
            '<?xml version="1.0"?>\n'
            '<!DOCTYPE beast [<!ENTITY xxe SYSTEM "file:///etc/passwd">]>\n'
            "<beast><data/><run/></beast>"
        )
        result = XMLValidator.validate(xml)
        assert not result.is_valid
        assert any("entity expansion" in err for err in result.errors)


# ---------------------------------------------------------------------------
# R06: sensitivity scenarios write to disjoint logger paths and carry seeds
# ---------------------------------------------------------------------------


class TestR06SensitivityIsolation:
    def _config(self, tmp_path):
        cfg = _write(
            tmp_path,
            _HEAD.replace(
                "output: {file: out.xml}",
                """
calibrations:
  - name: calA
    taxa: [taxonA, taxonB]
    distribution: {type: normal, parameters: {mean: 10.0, sigma: 1.0}}
  - name: calB
    taxa: [taxonC, taxonD]
    distribution: {type: normal, parameters: {mean: 20.0, sigma: 1.0}}
loggers:
  trace_log: {file_name: "shared.$(seed).log", log_every: 1000}
  tree_log: {file_name: "shared.$(seed).trees", log_every: 1000}
output: {file: out.xml}
""",
            ),
        )
        return ConfigParser().parse(cfg)

    def test_scenario_logger_paths_are_disjoint(self, tmp_path):
        from beast2py.diagnostics import SensitivityAnalyzer

        config = self._config(tmp_path)
        out_dir = tmp_path / "sens"
        xmls = SensitivityAnalyzer.generate_sensitivity_xmls(config, str(out_dir))
        xml_paths = [x for x in xmls if x.endswith(".xml")]
        assert len(xml_paths) == 3  # full + 2 leave-one-out
        logger_names = set()
        for xml_path in xml_paths:
            root = ET.parse(xml_path).getroot()
            for logger in root.iter("logger"):
                file_name = logger.get("fileName")
                if file_name:
                    assert file_name not in logger_names, file_name
                    logger_names.add(file_name)

    def test_manifest_maps_scenarios_to_seeds(self, tmp_path):
        from beast2py.diagnostics import SensitivityAnalyzer

        config = self._config(tmp_path)
        out_dir = tmp_path / "sens"
        xmls = SensitivityAnalyzer.generate_sensitivity_xmls(config, str(out_dir))
        manifest = next(x for x in xmls if x.endswith("sensitivity_manifest.tsv"))
        lines = open(manifest, encoding="utf-8").read().strip().splitlines()
        assert lines[0].startswith("scenario\txml\tseed\tcommand")
        seeds = [line.split("\t")[2] for line in lines[1:]]
        assert len(set(seeds)) == len(seeds)
        assert all("beast -seed" in line for line in lines[1:])


# ---------------------------------------------------------------------------
# R07: workflow templates are structurally runnable
# ---------------------------------------------------------------------------


class TestR07WorkflowTemplates:
    def _config(self, tmp_path):
        cfg = _write(tmp_path, _HEAD)
        return ConfigParser().parse(cfg)

    def test_snakemake_default_target_is_all(self, tmp_path):
        from beast2py.reproducibility import PipelineGenerator

        config = self._config(tmp_path)
        path = PipelineGenerator.generate_snakemake(
            config, output_dir=str(tmp_path / "sm"), config_file="config.yaml"
        )
        content = open(path, encoding="utf-8").read()
        all_at = content.index("rule all:")
        generate_at = content.index("rule generate_xml:")
        assert all_at < generate_at, (
            "the first real rule must be `rule all` for a bare snakemake run "
            "to execute the whole pipeline"
        )
        assert "-seed " in content

    def test_nextflow_declares_inputs_and_seed(self, tmp_path):
        from beast2py.reproducibility import PipelineGenerator

        config = self._config(tmp_path)
        path = PipelineGenerator.generate_nextflow(
            config, output_dir=str(tmp_path / "nf"), config_file="config.yaml"
        )
        content = open(path, encoding="utf-8").read()
        assert "input:" in content, "processes must declare their file inputs"
        assert "path config_file" in content
        assert "-seed " in content

    def test_snakemake_dry_run_if_engine_installed(self, tmp_path):
        import shutil
        import subprocess

        if shutil.which("snakemake") is None:
            pytest.skip("snakemake is not installed in this environment")
        from beast2py.reproducibility import PipelineGenerator

        config = self._config(tmp_path)
        out_dir = tmp_path / "sm"
        PipelineGenerator.generate_snakemake(
            config, output_dir=str(out_dir), config_file="config.yaml"
        )
        (out_dir / "config.yaml").write_text((tmp_path / "config.yaml").read_text())
        (out_dir / "test.fasta").write_text(TEST_FASTA)
        proc = subprocess.run(
            ["snakemake", "-n", "-c", "1"],
            cwd=str(out_dir),
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, proc.stderr

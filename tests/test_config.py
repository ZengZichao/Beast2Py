"""Tests for the config module."""

import tempfile
import warnings
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


class TestNameSpellingsAreCaseInsensitive:
    """A name typed as the paper prints it must parse to the tool's own key.

    ``_norm_name`` documents case-insensitive matching for substitution models,
    clock models and distributions, and three of the four name tables honoured
    it while the clock model and the calibration distribution did not: a reader
    copying ``Normal`` or ``ChiSquare`` out of Supplementary Table S1 was
    refused, and ``type: Strict`` failed although ``type: HKY`` worked.
    """

    # (spelling printed in Supplementary Table S1, canonical key, parameters)
    TABLE_S1 = [
        ("Normal", "normal", {"mean": 10.0, "sigma": 1.0}),
        ("LogNormal", "lognormal", {"M": 0.35, "S": 0.4}),
        ("Uniform", "uniform", {"lower": 5.0, "upper": 15.0}),
        ("Exponential", "exponential", {"mean": 10.0}),
        ("Gamma", "gamma", {"alpha": 2.0, "beta": 2.0}),
        ("Beta", "beta", {"alpha": 2.0, "beta": 2.0}),
        ("OneOnX", "one_on_x", {}),
        ("Laplace", "laplace", {"mu": 10.0, "scale": 1.0}),
        ("InverseGamma", "inverse_gamma", {"alpha": 3.0, "beta": 1.0}),
        ("Poisson", "poisson", {"lambda": 5.0}),
        ("ChiSquare", "chi_square", {"dof": 4.0}),
    ]

    def setup_method(self, method):
        """Create a temp directory holding the test alignment."""
        self.tmpdir = Path(tempfile.mkdtemp())
        (self.tmpdir / "test.fasta").write_text(TEST_FASTA)

    def _parse(self, mutate):
        """Parse TEST_CONFIG_YAML after applying *mutate* to its mapping."""
        import yaml

        cfg = yaml.safe_load(TEST_CONFIG_YAML)
        mutate(cfg)
        path = self.tmpdir / "c.yaml"
        path.write_text(yaml.safe_dump(cfg))
        return ConfigParser().parse(path)

    @pytest.mark.parametrize("spelling,canonical,params", TABLE_S1)
    def test_calibration_distribution_accepts_the_printed_name(self, spelling, canonical, params):
        """Every Table S1 name resolves, including the three with underscores."""

        def mutate(cfg):
            cfg["calibrations"][0]["distribution"] = {"type": spelling, "parameters": params}

        config = self._parse(mutate)
        assert config.calibrations[0].distribution.type == canonical

    @pytest.mark.parametrize("spelling", ["Strict", "STRICT", "ucln", "UCLN", "RLC", "UCE"])
    def test_clock_model_accepts_any_case(self, spelling):
        """The clock followed a different rule from the site model."""

        def mutate(cfg):
            cfg["partitions"][0]["clock_model"]["type"] = spelling

        config = self._parse(mutate)
        assert config.partitions[0].clock_model.type.value in ("strict", "ucln", "rlc", "uce")

    def test_emitted_file_is_the_same_whichever_spelling_is_used(self):
        """Folding must not change what gets written for the same prior."""
        from beast2py.api import generate_xml

        def emit(spelling):
            import yaml

            cfg = yaml.safe_load(TEST_CONFIG_YAML)
            cfg["alignments"][0]["file"] = str(self.tmpdir / "test.fasta")
            cfg["calibrations"][0]["distribution"] = {"type": spelling, "parameters": {"dof": 4.0}}
            path = self.tmpdir / ("%s.yaml" % spelling)
            path.write_text(yaml.safe_dump(cfg))
            out = self.tmpdir / ("%s.xml" % spelling)
            generate_xml(str(path), str(out), force=True)
            return out.read_text()

        assert emit("ChiSquare") == emit("chi_square")

    def test_unknown_names_are_still_rejected_with_the_valid_list(self):
        """The gate is not widened into accepting anything."""

        def mutate(cfg):
            cfg["calibrations"][0]["distribution"] = {"type": "bananasplit", "parameters": {}}

        with pytest.raises(ConfigError) as exc:
            self._parse(mutate)
        assert "valid:" in str(exc.value)
        assert "chi_square" in str(exc.value)


class TestProvenanceValueWhitelist:
    """Provenance *values* are branched on, so a misspelling must not choose a branch.

    The key names of a ``provenance`` block are whitelisted, but the values were
    free text: ``calibration_type`` is compared against the string ``"hard"``, so
    ``hardbound``, ``Hard-bound`` or ``h`` fell into the soft branch -- turning the
    one check that stops a file from being written into a warning that a notebook,
    an ignored warning filter or a discarded stderr swallows completely.
    """

    def setup_method(self, method):
        """Create a temp directory holding the test alignment."""
        self.tmpdir = Path(tempfile.mkdtemp())
        (self.tmpdir / "test.fasta").write_text(TEST_FASTA)

    def _parse_with_provenance(self, provenance):
        """Parse TEST_CONFIG_YAML with *provenance* attached to the one calibration."""
        import yaml

        cfg = yaml.safe_load(TEST_CONFIG_YAML)
        cfg["alignments"][0]["file"] = str(self.tmpdir / "test.fasta")
        cfg["calibrations"][0]["distribution"] = {
            "type": "normal",
            "parameters": {"mean": 6.0, "sigma": 0.5},
        }
        cfg["calibrations"][0]["provenance"] = provenance
        path = self.tmpdir / "c.yaml"
        path.write_text(yaml.safe_dump(cfg))
        return ConfigParser().parse(path)

    @pytest.mark.parametrize("spelling", ["soft", "hard", "Hard", "HARD", " hard "])
    def test_calibration_type_accepts_the_documented_values(self, spelling):
        """Case and padding fold onto the two documented values."""
        config = self._parse_with_provenance(
            {"source": "fossil", "calibration_type": spelling, "reference": "DOI:10.1038/x"}
        )
        assert config.calibrations[0].provenance.calibration_type == spelling.strip().lower()

    @pytest.mark.parametrize(
        "spelling", ["hardbound", "hard-bound", "h", "hardish", "definite", ""]
    )
    def test_calibration_type_refuses_anything_else(self, spelling):
        """The value is refused, and the message lists what is allowed."""
        with pytest.raises(ConfigError) as exc:
            self._parse_with_provenance(
                {"source": "fossil", "calibration_type": spelling, "reference": "DOI:10.1038/x"}
            )
        message = str(exc.value)
        assert "calibration_type" in message
        assert "hard" in message and "soft" in message

    @pytest.mark.parametrize("spelling", ["fossil", "biogeographic", "secondary", "other"])
    def test_source_accepts_the_documented_values(self, spelling):
        """The four values Supplementary Table S6 documents all parse."""
        config = self._parse_with_provenance(
            {"source": spelling, "calibration_type": "soft", "reference": "DOI:10.1038/x"}
        )
        assert config.calibrations[0].provenance.source == spelling

    @pytest.mark.parametrize("spelling", ["Fossils", "molecular", "secondary-cal", "unclear"])
    def test_source_refuses_anything_else(self, spelling):
        """A misspelt or invented source is named, not quietly kept."""
        with pytest.raises(ConfigError) as exc:
            self._parse_with_provenance(
                {"source": spelling, "calibration_type": "soft", "reference": "DOI:10.1038/x"}
            )
        message = str(exc.value)
        assert "source" in message
        assert "biogeographic" in message

    def test_an_omitted_source_stays_legal(self):
        """``source`` is optional; only a *given* value is whitelisted."""
        config = self._parse_with_provenance(
            {"calibration_type": "soft", "reference": "DOI:10.1038/x"}
        )
        assert config.calibrations[0].provenance.source == ""

    def test_a_hard_type_that_contradicts_the_prior_still_errors(self):
        """The whitelist is not a widening: a real hard bound is still enforced."""
        with pytest.raises(ConfigError, match="younger than the fossil"):
            self._parse_with_provenance(
                {
                    "source": "fossil",
                    "calibration_type": "hard",
                    "reference": "DOI:10.1038/x",
                    "original_age_min": 7.0,
                    "original_age_max": 9.0,
                }
            )

    def test_a_misspelt_hard_type_is_refused_before_the_branch(self):
        """The fail-open path this closes, pinned end to end.

        Before the whitelist this configuration only warned: the prior's 95%
        interval lies entirely below the recorded fossil range, the author wrote
        down that the bound is hard, and the tool read the misspelling as soft.
        Refusal happens while the value is parsed, so the soft branch is no longer
        reachable by a typo, and the message still names the fossil contradiction
        it was there to prevent.
        """
        with pytest.raises(ConfigError) as exc:
            self._parse_with_provenance(
                {
                    "source": "fossil",
                    "calibration_type": "hardbound",
                    "reference": "DOI:10.1038/x",
                    "original_age_min": 7.0,
                    "original_age_max": 9.0,
                }
            )
        message = str(exc.value)
        assert "invalid calibration_type 'hardbound'" in message
        assert "(valid: hard, soft)" in message

    def test_a_genuine_soft_bound_still_only_warns(self):
        """Soft remains soft: the convention is kept, not redefined."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self._parse_with_provenance(
                {
                    "source": "fossil",
                    "calibration_type": "soft",
                    "reference": "DOI:10.1038/x",
                    "original_age_min": 7.0,
                    "original_age_max": 9.0,
                }
            )
        assert any("younger than the fossil" in str(w.message) for w in caught), [
            str(w.message) for w in caught
        ]

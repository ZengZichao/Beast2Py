"""Coverage fixtures for the 36 declared model capabilities.

The matrix in ``docs/capability_coverage.md`` maps every clock model, tree
prior, substitution model and calibration distribution to a fixture, the XML
class it must produce, and optional external Gate 3 evidence.
These tests check the fixtures exist, parse, and emit the expected class; Gate 3
records are not bundled with the repository.
"""

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from beast2py.config import ConfigParser
from beast2py.xml_writer import XMLWriter

FIXTURES = Path(__file__).parent / "fixtures" / "coverage"
DISTRIBUTIONS = {
    "normal": "beast.base.inference.distribution.Normal",
    "beta": "beast.base.inference.distribution.Beta",
    "beta_offset": "beast.base.inference.distribution.Beta",
    "poisson": "beast.base.inference.distribution.Poisson",
    "one_on_x": "beast.base.inference.distribution.OneOnX",
    "exponential": "beast.base.inference.distribution.Exponential",
    "gamma": "beast.base.inference.distribution.Gamma",
    "laplace": "beast.base.inference.distribution.LaplaceDistribution",
    "inverse_gamma": "beast.base.inference.distribution.InverseGamma",
    "chi_square": "beast.base.inference.distribution.ChiSquare",
}
AA_MATRICES = {
    "dayhoff": "Dayhoff",
    "blosum62": "Blosum62",
    "cprev": "CPREV",
    "mtrev": "MTREV",
}


def _generate(name):
    config = ConfigParser().parse(FIXTURES / f"{name}.yaml")
    return config, XMLWriter(config).generate_xml()


@pytest.mark.parametrize("name,spec", sorted(DISTRIBUTIONS.items()))
def test_distribution_fixture_emits_expected_class(name, spec):
    _, xml = _generate(f"config_dist_{name}")
    root = ET.fromstring(xml)
    specs = [e.get("spec") for e in root.iter() if e.tag == "distr"]
    assert spec in specs


@pytest.mark.parametrize("name,spec", sorted(AA_MATRICES.items()))
def test_aminoacid_matrix_fixture_emits_expected_class(name, spec):
    config, xml = _generate(f"config_aa_{name}")
    assert config.partitions[0].alignment.data_type.value == "aminoacid"
    root = ET.fromstring(xml)
    specs = [e.get("spec") for e in root.iter() if e.tag == "substModel"]
    assert spec in specs


def test_coalescent_constant_fixture_initialises_a_population_model():
    config, xml = _generate("config_coalescent_constant")
    assert config.tree_prior_type.value == "coalescent_constant"
    root = ET.fromstring(xml)
    specs = {e.get("spec") for e in root.iter()}
    assert "Coalescent" in specs
    assert "ConstantPopulation" in specs
    assert "beast.base.evolution.tree.TreeIntervals" in specs
    assert any(
        e.tag == "operator"
        and e.get("spec") == "ScaleOperator"
        and e.get("parameter") == "@popSize"
        for e in root.iter()
    )


def test_nonzero_tip_dates_operator_targets_only_dated_taxa():
    _, xml = _generate("config_tip_dates_nonzero")
    root = ET.fromstring(xml)
    walkers = [
        e for e in root.iter() if e.tag == "operator" and e.get("spec") == "TipDatesRandomWalker"
    ]
    assert len(walkers) == 1
    taxonset_id = walkers[0].get("taxonset").lstrip("@")
    taxonset = next(e for e in root.iter() if e.tag == "taxonset" and e.get("id") == taxonset_id)
    refs = {child.get("idref") or child.get("id") for child in taxonset}
    assert refs == {"taxonB", "taxonC", "taxonD"}


def test_all_zero_tip_dates_do_not_emit_a_walker():
    from tests.conftest import EXAMPLES_DIR

    config = ConfigParser(base_dir=str(EXAMPLES_DIR)).parse(EXAMPLES_DIR / "config_tip_dates.yaml")
    root = ET.fromstring(XMLWriter(config).generate_xml())
    assert not [
        e for e in root.iter() if e.tag == "operator" and e.get("spec") == "TipDatesRandomWalker"
    ]


def test_chi_square_df_is_an_integer_attribute():
    _, xml = _generate("config_dist_chi_square")
    root = ET.fromstring(xml)
    chi = next(
        e
        for e in root.iter()
        if e.tag == "distr" and e.get("spec") == "beast.base.inference.distribution.ChiSquare"
    )
    # The input is an IntegerParameter; "4.0" is rejected by BEAST2 and a
    # <parameter> child is parsed as a RealParameter.
    assert chi.get("df") == "4"
    assert not [c for c in chi if c.tag == "parameter"]


def test_beta_offset_is_legal_and_beta_without_offset_warns():
    import warnings

    with pytest.warns(UserWarning, match="Beta without offset"):
        ConfigParser().parse(FIXTURES / "config_dist_beta.yaml")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ConfigParser().parse(FIXTURES / "config_dist_beta_offset.yaml")
    assert not [w for w in caught if "Beta without offset" in str(w.message)]


def test_one_on_x_and_poisson_warn_but_parse():
    with pytest.warns(UserWarning, match="OneOnX"):
        ConfigParser().parse(FIXTURES / "config_dist_one_on_x.yaml")
    with pytest.warns(UserWarning, match="Poisson"):
        ConfigParser().parse(FIXTURES / "config_dist_poisson.yaml")


def test_strict_distributions_upgrades_the_warning(tmp_path):
    import shutil

    from beast2py.config import ConfigError

    shutil.copy(FIXTURES / "nucleotide_4.fasta", tmp_path / "nucleotide_4.fasta")
    text = (FIXTURES / "config_dist_beta.yaml").read_text()
    text = text.replace(
        'tool_version: "0.2.0"',
        'tool_version: "0.2.0"\ndiagnostics:\n  strict_distributions: true',
    )
    cfg = tmp_path / "config.yaml"
    cfg.write_text(text)
    with pytest.raises(ConfigError, match="Beta without offset"):
        ConfigParser().parse(cfg)

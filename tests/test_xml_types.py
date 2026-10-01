"""Verify Beast2Py can generate all XML types required by
BEAST.app.package.v2.7.8 and BEAST.base.package.v2.7.8.

This test does not require BEAST2 itself; it only checks that XMLWriter
produces well-formed XML for every combination of model component supported
by the BEAST 2.7.8 base and app packages.
"""

import defusedxml.ElementTree as ET
import pytest

from beast2py.models import (
    BEASTConfig,
    CalibrationMethod,
    CalibrationPoint,
    ClockModelConfig,
    ClockModelType,
    DistributionConfig,
    MCMCConfig,
    Partition,
    RealParameter,
    SiteModelConfig,
    TreePriorType,
)
from beast2py.sequence import SequenceReader
from beast2py.xml_writer import XMLWriter


@pytest.fixture(scope="module")
def aln():
    return SequenceReader.read(
        "examples/primates.fasta",
        format="auto",
        data_type="nucleotide",
        alignment_id="alignment",
    )


def _make_config(
    aln,
    subst_type="hky",
    clock_type="strict",
    tree_prior="yule",
    tree_params=None,
    cal_dist="normal",
    cal_params=None,
    gamma_cats=4,
):
    subst = {"type": subst_type}
    if subst_type in ("gtr", "sym"):
        subst["rates"] = "1.0 1.0 1.0 1.0 1.0 1.0"
    site_model = SiteModelConfig(
        substitution_model=subst,
        gamma_categories=gamma_cats,
        gamma_shape=RealParameter(value=0.5, lower=0.0) if gamma_cats > 0 else None,
    )
    clock = ClockModelConfig(type=ClockModelType(clock_type))
    part = Partition(
        id="alignment", alignment=aln, site_model=site_model, clock_model=clock, tree="shared"
    )
    params = tree_params or {}
    cal = CalibrationPoint(
        name="rootCal",
        taxa=None,
        monophyletic=False,
        distribution=DistributionConfig(
            type=cal_dist,
            parameters=cal_params or {"mean": 10.0, "sigma": 1.0},
        ),
    )
    return BEASTConfig(
        metadata={"analysis_name": "test", "beast2_version": "2.7.8", "tool_version": "0.1.0"},
        partitions=[part],
        tree_prior_type=TreePriorType(tree_prior),
        tree_prior_params=params,
        calibrations=[cal],
        calibration_method=CalibrationMethod.MRCA_PRIOR,
        mcmc=MCMCConfig(chain_length=1000000),
    )


# --- Substitution models (nucleotide) ---

NUC_MODELS = ["jc69", "hky", "tn93", "gtr", "sym", "tim", "tvm"]


@pytest.mark.parametrize("model", NUC_MODELS)
def test_nucleotide_substitution_model(aln, model):
    config = _make_config(aln, subst_type=model)
    xml = XMLWriter(config).generate_xml()
    ET.fromstring(xml)


# --- Amino acid substitution models ---

AA_MODELS = ["wag", "jtt", "dayhoff", "blosum62", "cprev", "mtrev"]


@pytest.mark.parametrize("model", AA_MODELS)
def test_aa_substitution_model(model):
    aln = SequenceReader.read(
        "examples/aminoacid.fasta",
        format="auto",
        data_type="aminoacid",
        alignment_id="alignment",
    )
    config = _make_config(aln, subst_type=model)
    xml = XMLWriter(config).generate_xml()
    ET.fromstring(xml)


# --- Clock models ---

CLOCK_MODELS = ["strict", "ucln", "uce", "rlc"]


@pytest.mark.parametrize("model", CLOCK_MODELS)
def test_clock_model(aln, model):
    config = _make_config(aln, clock_type=model)
    xml = XMLWriter(config).generate_xml()
    ET.fromstring(xml)


# --- Tree priors ---

TREE_PRIORS = [
    ("yule", {}),
    (
        "birth_death",
        {"birth_rate": {"value": 1.0, "lower": 0.0}, "death_rate": {"value": 0.5, "lower": 0.0}},
    ),
    ("calibrated_yule", {}),
    ("coalescent_constant", {}),
    ("coalescent_exponential", {"pop_size": {"value": 1.0, "lower": 0.0}, "growth_rate": 0.0}),
    ("bayesian_skyline", {}),
    ("ebsp", {}),
]


@pytest.mark.parametrize("prior,params", TREE_PRIORS)
def test_tree_prior(aln, prior, params):
    cal_method = (
        CalibrationMethod.CALIBRATED_YULE
        if prior == "calibrated_yule"
        else CalibrationMethod.MRCA_PRIOR
    )
    config = _make_config(aln, tree_prior=prior, tree_params=params)
    config.calibration_method = cal_method
    xml = XMLWriter(config).generate_xml()
    ET.fromstring(xml)


# --- Calibration distributions ---

DISTRIBUTIONS = [
    ("normal", {"mean": 10.0, "sigma": 1.0}),
    ("lognormal", {"M": 2.0, "S": 0.5}),
    ("uniform", {"lower": 5.0, "upper": 15.0}),
    ("exponential", {"mean": 10.0}),
    ("gamma", {"alpha": 2.0, "beta": 0.2}),
    ("beta", {"alpha": 2.0, "beta": 2.0}),
    ("laplace", {"mu": 10.0, "scale": 1.0}),
    ("inverse_gamma", {"alpha": 2.0, "beta": 2.0}),
    ("one_on_x", {}),
    ("poisson", {"lambda": 10.0}),
    ("chi_square", {"nu": 5.0}),
]


@pytest.mark.parametrize("dist,params", DISTRIBUTIONS)
def test_calibration_distribution(aln, dist, params):
    config = _make_config(aln, cal_dist=dist, cal_params=params)
    xml = XMLWriter(config).generate_xml()
    ET.fromstring(xml)


# --- Multi-partition ---


def test_multi_partition(aln):
    site_model_a = SiteModelConfig(
        substitution_model={"type": "hky"},
        gamma_categories=4,
        gamma_shape=RealParameter(value=0.5, lower=0.0),
    )
    site_model_b = SiteModelConfig(
        substitution_model={"type": "gtr", "rates": "1.0 1.0 1.0 1.0 1.0 1.0"},
        gamma_categories=4,
        gamma_shape=RealParameter(value=0.5, lower=0.0),
    )
    clock_a = ClockModelConfig(type=ClockModelType.UCLN)
    clock_b = ClockModelConfig(type=ClockModelType.UCLN, linked_to="alignment")
    part_a = Partition(
        id="alignment", alignment=aln, site_model=site_model_a, clock_model=clock_a, tree="shared"
    )
    part_b = Partition(
        id="alignment2", alignment=aln, site_model=site_model_b, clock_model=clock_b, tree="shared"
    )
    cal = CalibrationPoint(
        name="root",
        taxa=None,
        monophyletic=False,
        distribution=DistributionConfig(type="normal", parameters={"mean": 10.0, "sigma": 1.0}),
    )
    config = BEASTConfig(
        metadata={"analysis_name": "multi", "beast2_version": "2.7.8", "tool_version": "0.1.0"},
        partitions=[part_a, part_b],
        tree_prior_type=TreePriorType.BIRTH_DEATH,
        tree_prior_params={
            "birth_rate": {"value": 1.0, "lower": 0.0},
            "death_rate": {"value": 0.5, "lower": 0.0},
        },
        calibrations=[cal],
        calibration_method=CalibrationMethod.MRCA_PRIOR,
        mcmc=MCMCConfig(chain_length=1000000),
    )
    xml = XMLWriter(config).generate_xml()
    ET.fromstring(xml)


# --- Tip dates ---


def test_tip_dates(aln):
    from beast2py.models import TipDatesConfig

    config = _make_config(aln)
    config.tip_dates = TipDatesConfig(
        enabled=True,
        trait_name="date-backward",
        units="year",
        dates={t: str(i * 0.1) for i, t in enumerate(aln.taxa_names)},
    )
    xml = XMLWriter(config).generate_xml()
    ET.fromstring(xml)


# --- Nested sampling MCMC ---


def test_nested_sampling_mcmc(aln):
    from beast2py.models import MCMCType

    config = _make_config(aln)
    config.mcmc.type = MCMCType.NESTED_SAMPLING
    xml = XMLWriter(config).generate_xml()
    ET.fromstring(xml)

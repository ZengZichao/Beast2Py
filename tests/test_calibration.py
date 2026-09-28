"""Tests for the CalibrationBuilder module."""

from beast2py.calibration import CalibrationBuilder
from beast2py.models import (
    CalibrationPoint,
    DistributionConfig,
)


def test_build_taxonset():
    taxa = ["human", "chimp", "gorilla"]
    elem = CalibrationBuilder.build_taxonset(taxa, "testClade")
    assert elem.tag == "taxonset"
    assert elem.get("spec") == "TaxonSet"
    assert elem.get("id") == "testClade.taxonset"
    children = list(elem)
    assert len(children) == 3
    assert all(c.tag == "taxon" for c in children)
    assert children[0].get("spec") == "Taxon"
    assert children[0].get("id") == "human"


def test_build_distribution_normal():
    dist = DistributionConfig(
        type="normal",
        parameters={"mean": 6.0, "sigma": 0.5},
    )
    elem = CalibrationBuilder.build_distribution(dist, dist_id="test.distr")
    assert elem.tag == "distr"
    assert elem.get("id") == "test.distr"
    params = list(elem)
    assert len(params) == 2
    names = {p.get("name") for p in params}
    assert "mean" in names
    assert "sigma" in names


def test_build_distribution_uniform():
    dist = DistributionConfig(
        type="uniform",
        parameters={"lower": 0.0, "upper": 10.0},
    )
    elem = CalibrationBuilder.build_distribution(dist)
    assert elem.tag == "distr"
    assert elem.get("lower") == "0.0"
    assert elem.get("upper") == "10.0"


def test_build_distribution_with_offset():
    dist = DistributionConfig(
        type="normal",
        parameters={"mean": 6.0, "sigma": 0.5},
        offset=12.0,
    )
    elem = CalibrationBuilder.build_distribution(dist)
    assert elem.get("offset") == "12.0"


def test_build_mrcaprior():
    cal = CalibrationPoint(
        name="testCal",
        taxa=["human", "chimp"],
        monophyletic=True,
        distribution=DistributionConfig(
            type="normal",
            parameters={"mean": 6.0, "sigma": 0.5},
        ),
    )
    elem = CalibrationBuilder.build_mrcaprior(cal, "tree")
    assert elem.tag == "distribution"
    assert elem.get("spec") == "beast.base.evolution.tree.MRCAPrior"
    assert elem.get("tree") == "@tree"
    assert elem.get("monophyletic") == "true"


def test_build_mrcaprior_use_originate():
    cal = CalibrationPoint(
        name="testCal",
        taxa=["human", "chimp"],
        monophyletic=True,
        distribution=DistributionConfig(
            type="uniform",
            parameters={"lower": 0.0, "upper": 10.0},
        ),
        use_originate=True,
    )
    elem = CalibrationBuilder.build_mrcaprior(cal, "tree")
    assert elem.get("useOriginate") == "true"


def test_build_calibrated_yule_point():
    cal = CalibrationPoint(
        name="testCal",
        taxa=["human", "chimp"],
        monophyletic=True,
        distribution=DistributionConfig(
            type="uniform",
            parameters={"lower": 4.0, "upper": 7.0},
        ),
    )
    elem = CalibrationBuilder.build_calibrated_yule_point(cal)
    assert elem.tag == "calibrations"
    assert elem.get("spec") == "CalibrationPoint"


def test_build_mrca_logger():
    cal = CalibrationPoint(
        name="testCal",
        taxa=["human", "chimp"],
        monophyletic=True,
    )
    elem = CalibrationBuilder.build_mrca_logger(cal, "tree")
    assert elem.tag == "log"
    assert elem.get("spec") == "beast.base.evolution.tree.MRCAPrior"
    assert elem.get("tree") == "@tree"


def test_build_tip_dates_traitset():
    dates = {"human": 0, "chimp": 5}
    elem = CalibrationBuilder.build_tip_dates_traitset(dates, "date-forward", "year", "alignment")
    assert elem.tag == "trait"
    assert elem.get("traitname") == "date-forward"
    assert elem.get("units") == "year"
    assert "human = 0" in elem.get("value")
    assert "chimp = 5" in elem.get("value")

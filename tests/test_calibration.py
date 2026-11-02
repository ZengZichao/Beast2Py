"""Tests for the CalibrationBuilder module."""

import xml.etree.ElementTree as ET

import pytest

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


def _yule_point(name="testCal", taxa=("human", "chimp"), use_originate=False):
    return CalibrationPoint(
        name=name,
        taxa=list(taxa) if taxa is not None else None,
        monophyletic=True,
        distribution=DistributionConfig(
            type="uniform",
            parameters={"lower": 4.0, "upper": 7.0},
        ),
        use_originate=use_originate,
    )


def test_build_calibrated_yule_point():
    elem = CalibrationBuilder.build_calibrated_yule_point(_yule_point())
    assert elem.tag == "calibrations"
    assert elem.get("spec") == "CalibrationPoint"


@pytest.mark.parametrize(
    "use_originate,expected",
    [(False, None), (True, "true")],
    ids=["crown", "stem"],
)
def test_calibrated_yule_writes_forparent_only_for_a_stem_point(use_originate, expected):
    """A stem calibration reaches BEAST2, and only a stem one does.

    ``forParent`` is how ``CalibrationPoint`` selects the parent of the MRCA, so
    the flag has to land in the XML when the YAML asked for it: without this
    assertion a `use_originate` change would move only the fingerprint, which is
    exactly the gap the review flagged.
    """
    elem = CalibrationBuilder.build_calibrated_yule_point(_yule_point(use_originate=use_originate))
    assert elem.get("forParent") == expected


def test_calibrated_yule_refuses_a_stem_reading_of_the_root():
    with pytest.raises(ValueError, match="root has no parent"):
        CalibrationBuilder.build_calibrated_yule_point(_yule_point(taxa=None, use_originate=True))


def test_mrcaprior_refuses_a_stem_reading_of_the_root():
    with pytest.raises(ValueError, match="root has no parent"):
        CalibrationBuilder.build_mrcaprior(_yule_point(taxa=None, use_originate=True), "tree")


def test_toggling_use_originate_changes_the_emitted_xml():
    """Not just the fingerprint: the two readings must be different files."""
    crown = ET.tostring(
        CalibrationBuilder.build_calibrated_yule_point(_yule_point()), encoding="unicode"
    )
    stem = ET.tostring(
        CalibrationBuilder.build_calibrated_yule_point(_yule_point(use_originate=True)),
        encoding="unicode",
    )
    assert crown != stem
    assert 'forParent="true"' in stem


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


def test_calibrated_yule_stem_point_reaches_the_generated_xml():
    """End to end: a stem calibration in the tree prior is a ``forParent`` in XML.

    The unit test above only builds one element.  This one goes through the real
    writer on the committed ``calibrated_yule`` example, because the review's
    complaint was precisely that a YAML flag could change the fingerprint without
    changing what BEAST2 reads.
    """
    import xml.etree.ElementTree as ET

    from beast2py.config import ConfigParser
    from beast2py.xml_writer import XMLWriter

    from .conftest import EXAMPLES_DIR

    config = ConfigParser().parse(EXAMPLES_DIR / "config_calibrated_yule.yaml")
    assert config.calibrations, "the example carries at least one calibration"
    target = config.calibrations[0]
    assert target.use_originate is False, "the committed example is a crown reading"
    assert target.taxa is not None, "a root point cannot be a stem point"

    crown_xml = XMLWriter(config).generate_xml()
    target.use_originate = True
    stem_xml = XMLWriter(config).generate_xml()

    assert 'forParent="true"' in stem_xml
    assert 'forParent="true"' not in crown_xml
    points = [
        elem
        for elem in ET.fromstring(stem_xml).iter("calibrations")
        if elem.get("spec") == "CalibrationPoint"
    ]
    flagged = [elem for elem in points if elem.get("forParent") == "true"]
    assert len(flagged) == 1
    assert flagged[0].get("id").startswith(target.name)

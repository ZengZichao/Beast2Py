"""Tests for the utils module."""

from beast2py.utils import (
    make_xml_id,
    ref,
    serialize_xml,
    compute_distribution_stats,
    compute_wasserstein_distance,
    sanitize_id,
    build_namespace,
)
import xml.etree.ElementTree as ET


class TestUtils:
    """Test utility functions."""

    def test_make_xml_id(self):
        assert make_xml_id("hky.kappa") == "hky.kappa"
        assert make_xml_id("hky.kappa", "gene1") == "gene1.hky.kappa"
        assert make_xml_id("birthRate") == "birthRate"

    def test_ref(self):
        assert ref("tree") == "@tree"
        assert ref("alignment") == "@alignment"

    def test_serialize_xml(self):
        root = ET.Element("test", {"attr": "value"})
        child = ET.SubElement(root, "child")
        child.text = "content"
        xml_str = serialize_xml(root)
        assert "<test" in xml_str
        assert "<child>content</child>" in xml_str

    def test_sanitize_id(self):
        assert sanitize_id("test_name") == "test_name"
        assert sanitize_id("test name") == "test_name"
        assert sanitize_id("test-name") == "test-name"
        assert sanitize_id("test@name!") == "test_name_"

    def test_build_namespace(self):
        ns = build_namespace()
        assert "beast.base.core" in ns
        assert "beast.base.inference" in ns
        assert ":" in ns


class TestDistributionStats:
    """Test distribution statistics computation."""

    def test_normal_stats(self):
        stats = compute_distribution_stats("normal", {"mean": 5.0, "sigma": 1.0})
        assert stats["mean"] == 5.0
        assert stats["hpd_lower"] < 5.0
        assert stats["hpd_upper"] > 5.0
        assert abs(stats["hpd_lower"] - (5.0 - 1.96)) < 0.01
        assert abs(stats["hpd_upper"] - (5.0 + 1.96)) < 0.01

    def test_uniform_stats(self):
        stats = compute_distribution_stats("uniform", {"lower": 0.0, "upper": 10.0})
        assert stats["mean"] == 5.0
        assert stats["hpd_lower"] == 0.0
        assert stats["hpd_upper"] == 10.0

    def test_lognormal_stats(self):
        stats = compute_distribution_stats("lognormal", {"M": 1.0, "S": 0.5})
        assert stats["mean"] > 0
        assert stats["median"] > 0
        assert stats["hpd_lower"] < stats["hpd_upper"]

    def test_exponential_stats(self):
        stats = compute_distribution_stats("exponential", {"mean": 5.0})
        assert stats["mean"] == 5.0
        assert stats["hpd_lower"] >= 0.0

    def test_gamma_stats(self):
        # BEAST2's Gamma defaults to mode=ShapeScale, so beta is a SCALE and
        # mean = alpha * beta (Gamma.java:22-31).  This test used to assert
        # alpha / beta, i.e. the rate convention, which is what caused the
        # distribution from the one the XML emits.
        stats = compute_distribution_stats("gamma", {"alpha": 2.0, "beta": 0.5})
        assert abs(stats["mean"] - 1.0) < 0.01  # alpha * beta = 2.0 * 0.5
        rate = compute_distribution_stats("gamma", {"alpha": 2.0, "beta": 0.5, "mode": "ShapeRate"})
        assert abs(rate["mean"] - 4.0) < 0.01  # alpha / beta when asked for

    def test_stats_with_offset(self):
        stats = compute_distribution_stats("normal", {"mean": 5.0, "sigma": 1.0}, offset=2.0)
        assert stats["mean"] == 7.0


class TestWassersteinDistance:
    """Test Wasserstein distance computation."""

    def test_identical_distributions(self):
        dist = compute_wasserstein_distance(
            "normal",
            {"mean": 5.0, "sigma": 1.0},
            0.0,
            "normal",
            {"mean": 5.0, "sigma": 1.0},
            0.0,
        )
        assert dist < 0.5  # Should be very small

    def test_different_distributions(self):
        dist = compute_wasserstein_distance(
            "normal",
            {"mean": 5.0, "sigma": 1.0},
            0.0,
            "normal",
            {"mean": 15.0, "sigma": 1.0},
            0.0,
        )
        assert dist > 5.0  # Should be large

"""Tests for the utils module."""

import numpy as np

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

    def test_hpd_bounds_use_no_sampling(self):
        """The docstring claims the HPD bounds are analytic; hold it to that.

        ``compute_distribution_stats`` advertises that no sampling is involved
        and that ``seed`` cannot change the result. Both claims are checked
        here against every skewed family the docstring names, so a future
        switch back to a Monte Carlo estimator fails this test instead of only
        making the docs wrong again (CODE-005).
        """
        from scipy import stats as sp_stats

        from beast2py.utils import DEFAULT_SEED, _analytic_hpd

        cases = [
            ("gamma", {"alpha": 2.0, "beta": 0.5}, sp_stats.gamma(2.0, scale=0.5)),
            ("lognormal", {"M": 0.0, "S": 0.5}, sp_stats.lognorm(0.5, scale=1.0)),
            ("beta", {"alpha": 2.0, "beta": 5.0}, sp_stats.beta(2.0, 5.0)),
            ("inverse_gamma", {"alpha": 3.0, "beta": 2.0}, sp_stats.invgamma(3.0, scale=2.0)),
            ("chi_square", {"dof": 4.0}, sp_stats.chi2(4.0)),
        ]
        for dist_type, params, frozen in cases:
            reference = _analytic_hpd(frozen)
            baseline = compute_distribution_stats(dist_type, params)
            # Exact equality: the bounds are the quantile-minimiser output, not
            # an estimate of it.
            assert (baseline["hpd_lower"], baseline["hpd_upper"]) == reference, dist_type
            # A seed cannot move the answer, because nothing is drawn.
            for other_seed in (0, 1, DEFAULT_SEED + 1, 20261008):
                seeded = compute_distribution_stats(dist_type, params, seed=other_seed)
                assert seeded == baseline, (dist_type, other_seed)

        # Neither can ambient random state: burn the global RNG and the legacy
        # sampled estimator's own RNG, then re-check against the pre-burn result.
        np.random.seed(20261008)
        rng = np.random.default_rng(20261008)
        rng.standard_normal(10_000)
        np.random.standard_normal(10_000)
        for dist_type, params, _frozen in cases:
            after = compute_distribution_stats(dist_type, params)
            assert after["hpd_lower"] == _analytic_hpd(_frozen)[0], dist_type
            assert after["hpd_upper"] == _analytic_hpd(_frozen)[1], dist_type


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

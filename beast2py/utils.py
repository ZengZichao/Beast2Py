"""Utility functions for Beast2Py.

Provides ID generation, XML serialization helpers, distribution calculations,
and common utility functions used across modules.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from xml.dom import minidom
from typing import Any, Dict, Tuple

import numpy as np


def make_xml_id(prefix: str, partition_id: str = "") -> str:
    """Generate an XML element ID with optional partition prefix.

    Args:
        prefix: The base ID (e.g. 'hky.kappa').
        partition_id: Optional partition prefix (e.g. 'gene1').

    Returns:
        An ID string like 'gene1.hky.kappa' or 'hky.kappa'.
    """
    if partition_id:
        return f"{partition_id}.{prefix}"
    return prefix


def ref(id_val: str) -> str:
    """Format an ID reference (e.g. '@tree' or '#tree')."""
    return f"@{id_val}"


# ---------------------------------------------------------------------------
# Strict scalar coercion
# ---------------------------------------------------------------------------

_TRUE_STRINGS = {"true", "yes", "on", "1"}
_FALSE_STRINGS = {"false", "no", "off", "0"}


def as_strict_bool(value: Any, field: str) -> bool:
    """Coerce a configuration value to a bool without Python truthiness traps.

    ``bool("false")`` is ``True`` in Python, so casting a quoted ``"false"``
    through Python truthiness would invert the user's intent silently. This
    helper only accepts real booleans and an explicit set of spellings;
    everything else is a configuration error.

    Args:
        value: The raw value.
        field: Field name, used in the error message.

    Returns:
        The parsed boolean.

    Raises:
        ValueError: If the value cannot be interpreted as a boolean. Callers
            translate this into ``ConfigError``.
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if value in (0, 1):
            return bool(value)
        raise ValueError(
            f"'{field}': {value!r} is not a boolean; use true/false "
            f"(unquoted in YAML) or one of {sorted(_TRUE_STRINGS | _FALSE_STRINGS)}"
        )
    if isinstance(value, str):
        token = value.strip().lower()
        if token in _TRUE_STRINGS:
            return True
        if token in _FALSE_STRINGS:
            return False
    raise ValueError(
        f"'{field}': {value!r} is not a boolean; use true/false "
        f"(unquoted in YAML) or one of {sorted(_TRUE_STRINGS | _FALSE_STRINGS)}"
    )


def as_number(value: Any, field: str) -> float:
    """Coerce a configuration value to float, rejecting booleans and strings.

    Args:
        value: The raw value.
        field: Field name, used in the error message.

    Returns:
        The parsed float.

    Raises:
        ValueError: If the value is not a finite number.
    """
    if isinstance(value, bool) or value is None:
        raise ValueError(f"'{field}': expected a number, got {value!r}")
    try:
        num = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"'{field}': expected a number, got {value!r}") from exc
    if num != num or num in (float("inf"), float("-inf")):
        raise ValueError(f"'{field}': expected a finite number, got {value!r}")
    return num


def format_number(value: Any) -> str:
    """Render a numeric value for XML without floating-point noise.

    ``0.1 - 0.4`` must not be serialised as
    ``-0.30000000000000004``.

    Args:
        value: A number or a whitespace-separated vector of numbers.

    Returns:
        A compact decimal string.
    """
    if isinstance(value, str):
        tokens = value.split()
        if len(tokens) > 1:
            return " ".join(format_number(float(tok)) for tok in tokens)
        try:
            return format_number(float(value))
        except (TypeError, ValueError):
            return value
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    num = float(value)
    if num == int(num) and abs(num) < 1e15:
        return str(int(num))
    return f"{num:.10g}"


def serialize_xml(root: ET.Element, indent: bool = True) -> str:
    """Serialize an ElementTree Element to a formatted XML string.

    Args:
        root: The root XML Element.
        indent: Whether to pretty-print with indentation.

    Returns:
        XML string with a normalized UTF-8 XML declaration.
    """
    rough_string = ET.tostring(root, encoding="unicode")
    if indent:
        parsed = minidom.parseString(rough_string)
        pretty = parsed.toprettyxml(indent="    ", encoding=None)
        # Remove extra blank lines
        lines = [line for line in pretty.split("\n") if line.strip()]
        pretty = "\n".join(lines)
    else:
        pretty = rough_string
    # Normalize whatever declaration minidom emitted to the documented form
    body = re.sub(r"^\s*<\?xml[^?]*\?>\s*", "", pretty)
    return '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n' + body


# Deterministic default seed for Monte Carlo summaries. Reproducibility of
# the diagnostics output must not depend on ambient random state.
DEFAULT_SEED = 42

# Monte Carlo sample count used for shortest-interval HPD estimation
_HPD_N_SAMPLES = 100_000


def _analytic_hpd(dist, mass: float = 0.95) -> Tuple[float, float]:
    """Shortest interval carrying `mass` for a scipy (frozen) continuous dist.

    A sampled HPD would carry Monte Carlo discretisation error into a number
    that then *gates* a hard ``calibration_type: hard`` check: a fossil close
    to an interval edge would be classified by sampling noise rather than by
    the density.
    For a unimodal density the HPD is the width-minimising pair
    ``(ppf(q), ppf(q + mass))``, which this solves directly.

    Args:
        dist: A scipy frozen distribution exposing ``ppf``.
        mass: Probability the interval must carry.

    Returns:
        ``(lower, upper)`` endpoint pair.
    """
    from scipy.optimize import minimize_scalar

    def width(q: float) -> float:
        return float(dist.ppf(q + mass) - dist.ppf(q))

    res = minimize_scalar(
        width,
        bounds=(0.0, 1.0 - mass),
        method="bounded",
        options={"xatol": 1e-12},
    )
    q = float(res.x)
    return float(dist.ppf(q)), float(dist.ppf(q + mass))


def _shortest_interval(samples: "np.ndarray", mass: float = 0.95) -> Tuple[float, float]:
    """Return the shortest interval containing `mass` of the samples.

    For a unimodal distribution this converges to the 95% HPD interval; for
    monotone densities it converges to the one-sided HPD region.
    """
    x = np.sort(np.asarray(samples, dtype=float))
    n = len(x)
    k = int(np.floor(mass * n))
    if k >= n:
        k = n - 1
    widths = x[k:] - x[: n - k]
    i = int(np.argmin(widths))
    return float(x[i]), float(x[i + k])


def gamma_scale(parameters: Dict[str, Any]) -> Tuple[float, float, str]:
    """Return (alpha, scale, mode) for a BEAST2 Gamma, honouring ``mode``.

    BEAST2's ``beast.base.inference.distribution.Gamma`` defaults to
    ``mode=ShapeScale`` (Gamma.java:22-31), i.e. ``beta`` is a *scale* and the
    mean is ``alpha * beta``.  ``mode=ShapeRate`` inverts it.  Both
    ``_frozen_distribution`` and ``compute_distribution_stats`` go through this
    one function, so both describe the same distribution BEAST2 samples --
    treating ``beta`` as a rate where BEAST2 means a scale would silently
    shift the prior.
    """
    alpha = float(parameters.get("alpha", 1.0))
    beta = float(parameters.get("beta", 1.0))
    mode = str(parameters.get("mode", "ShapeScale")).lower()
    scale = 1.0 / beta if mode == "shaperate" else beta
    return alpha, scale, mode


def _frozen_distribution(dist_type: str, parameters: Dict[str, Any], offset: float):
    """Return a scipy frozen distribution (without offset shift) or None.

    Raises:
        ValueError: for distributions whose support cannot be normalised
            (OneOnX is improper without an upper bound).
    """
    from scipy import stats as sp_stats

    dt = dist_type.lower()

    if dt == "normal":
        return sp_stats.norm(
            loc=float(parameters.get("mean", 0.0)),
            scale=float(parameters.get("sigma", 1.0)),
        )
    if dt in ("lognormal", "log_normal"):
        M = float(parameters.get("M", 0.0))
        S = float(parameters.get("S", 1.0))
        log_mean = np.log(M) - S**2 / 2 if parameters.get("mean_in_real_space", False) else M
        return sp_stats.lognorm(s=S, scale=float(np.exp(log_mean)))
    if dt == "uniform":
        return sp_stats.uniform(
            loc=float(parameters.get("lower", 0.0)),
            scale=float(parameters.get("upper", 1.0)) - float(parameters.get("lower", 0.0)),
        )
    if dt == "exponential":
        return sp_stats.expon(scale=float(parameters.get("mean", 1.0)))
    if dt == "gamma":
        # See gamma_scale(): beta is a SCALE unless mode=ShapeRate
        alpha, scale, _ = gamma_scale(parameters)
        return sp_stats.gamma(a=alpha, scale=scale)
    if dt == "beta":
        return sp_stats.beta(
            float(parameters.get("alpha", 1.0)),
            float(parameters.get("beta", 1.0)),
        )
    if dt in ("laplace", "laplacedistribution"):
        return sp_stats.laplace(
            loc=float(parameters.get("mu", 0.0)),
            scale=float(parameters.get("scale", 1.0)),
        )
    if dt in ("inverse_gamma", "inversegamma"):
        return sp_stats.invgamma(
            a=float(parameters.get("alpha", 1.0)),
            scale=float(parameters.get("beta", 1.0)),
        )
    if dt == "chi_square":
        return sp_stats.chi2(df=float(parameters.get("dof", 1.0)))
    if dt == "poisson":
        return sp_stats.poisson(mu=float(parameters.get("lambda", 1.0)))
    if dt in ("one_on_x", "oneonx", "one_div_x"):
        raise ValueError(
            "OneOnX is an improper density (f(x) = 1/x) without an upper "
            "bound, so its mean and HPD interval are undefined. Specify "
            "explicit bounds via a uniform distribution instead."
        )
    raise ValueError(f"Unsupported distribution type: {dist_type}")


def compute_distribution_stats(
    dist_type: str,
    parameters: Dict[str, Any],
    offset: float = 0.0,
    seed: int = DEFAULT_SEED,
) -> Dict[str, Any]:
    """Compute summary statistics (mean, median, 95% HPD bounds) for a distribution.

    Closed-form 95% HPD intervals are used where available (normal, uniform,
    exponential, Laplace, one-sided gamma); for the remaining skewed
    distributions the HPD interval is estimated from the shortest interval of
    a large deterministic Monte Carlo sample. Poisson (discrete) uses the
    central 95% probability interval.

    Args:
        dist_type: Distribution type name.
        parameters: Distribution parameters.
        offset: Distribution offset.
        seed: Seed for deterministic Monte Carlo summaries.

    Returns:
        Dict with 'mean', 'median', 'hpd_lower', 'hpd_upper' and
        'interval_type' keys.

    Raises:
        ValueError: If the distribution type is unsupported or its summary
            statistics are undefined (OneOnX).
    """
    from scipy import stats as sp_stats

    dt = dist_type.lower()
    offset_f = float(offset)

    if dt == "normal":
        mean = float(parameters.get("mean", 0.0)) + offset_f
        sigma = float(parameters.get("sigma", 1.0))
        return {
            "mean": mean,
            "median": mean,
            "hpd_lower": mean - 1.959964 * sigma,
            "hpd_upper": mean + 1.959964 * sigma,
            "interval_type": "hpd",
        }

    if dt == "uniform":
        lower = float(parameters.get("lower", 0.0)) + offset_f
        upper = float(parameters.get("upper", 1.0)) + offset_f
        mean = (lower + upper) / 2
        return {
            "mean": mean,
            "median": mean,
            "hpd_lower": lower,
            "hpd_upper": upper,
            "interval_type": "hpd",
        }

    if dt == "exponential":
        mean_val = float(parameters.get("mean", 1.0))
        # Monotone decreasing density: the 95% HPD region is one-sided
        hpd_lower = 0.0 + offset_f
        hpd_upper = sp_stats.expon.ppf(0.95, scale=mean_val) + offset_f
        return {
            "mean": mean_val + offset_f,
            "median": sp_stats.expon.ppf(0.5, scale=mean_val) + offset_f,
            "hpd_lower": hpd_lower,
            "hpd_upper": hpd_upper,
            "interval_type": "hpd",
        }

    if dt == "gamma":
        # Same parameterisation as the emitted XML: mean = alpha * scale, so a
        # ShapeScale Gamma reports the mean BEAST2 will actually sample
        # (; the rate convention is only used for mode=ShapeRate).
        alpha, scale, _mode = gamma_scale(parameters)
        median = sp_stats.gamma.ppf(0.5, alpha, scale=scale) + offset_f
        mean_val = alpha * scale + offset_f
        if alpha <= 1.0:
            # Monotone decreasing density: one-sided HPD
            return {
                "mean": mean_val,
                "median": median,
                "hpd_lower": 0.0 + offset_f,
                "hpd_upper": sp_stats.gamma.ppf(0.95, alpha, scale=scale) + offset_f,
                "interval_type": "hpd",
            }
        lo, hi = _analytic_hpd(sp_stats.gamma(alpha, scale=scale))
        return {
            "mean": mean_val,
            "median": median,
            "hpd_lower": lo + offset_f,
            "hpd_upper": hi + offset_f,
            "interval_type": "hpd",
        }

    if dt in ("laplace", "laplacedistribution"):
        mu = float(parameters.get("mu", 0.0)) + offset_f
        scale = float(parameters.get("scale", 1.0))
        # Exact symmetric 95% HPD: P(|X - mu| <= d) = 1 - exp(-d/b) = 0.95
        d = scale * np.log(20.0)
        return {
            "mean": mu,
            "median": mu,
            "hpd_lower": mu - d,
            "hpd_upper": mu + d,
            "interval_type": "hpd",
        }

    if dt == "poisson":
        lam = float(parameters.get("lambda", 1.0))
        mean_val = lam + offset_f
        lower = float(sp_stats.poisson.ppf(0.025, lam)) + offset_f
        upper = float(sp_stats.poisson.ppf(0.975, lam)) + offset_f
        return {
            "mean": mean_val,
            "median": float(sp_stats.poisson.ppf(0.5, lam)) + offset_f,
            "hpd_lower": lower,
            "hpd_upper": upper,
            "interval_type": "central",
        }

    if dt in ("lognormal", "log_normal"):
        M = float(parameters.get("M", 0.0))
        S = float(parameters.get("S", 1.0))
        mean_in_real_space = parameters.get("mean_in_real_space", False)
        log_mean = np.log(M) - S**2 / 2 if mean_in_real_space else M
        real_mean = np.exp(log_mean + S**2 / 2)
        lo, hi = _analytic_hpd(sp_stats.lognorm(S, scale=np.exp(log_mean)))
        return {
            "mean": float(real_mean) + offset_f,
            "median": float(np.exp(log_mean)) + offset_f,
            "hpd_lower": lo + offset_f,
            "hpd_upper": hi + offset_f,
            "interval_type": "hpd",
        }

    if dt == "beta":
        alpha = float(parameters.get("alpha", 1.0))
        beta_p = float(parameters.get("beta", 1.0))
        lo, hi = _analytic_hpd(sp_stats.beta(alpha, beta_p))
        return {
            "mean": alpha / (alpha + beta_p) + offset_f,
            "median": float(sp_stats.beta.ppf(0.5, alpha, beta_p)) + offset_f,
            "hpd_lower": lo + offset_f,
            "hpd_upper": hi + offset_f,
            "interval_type": "hpd",
        }

    if dt in ("inverse_gamma", "inversegamma"):
        alpha = float(parameters.get("alpha", 1.0))
        beta = float(parameters.get("beta", 1.0))
        mean_val = beta / (alpha - 1) + offset_f if alpha > 1 else float("inf")
        lo, hi = _analytic_hpd(sp_stats.invgamma(alpha, scale=beta))
        return {
            "mean": mean_val,
            "median": float(sp_stats.invgamma.ppf(0.5, alpha, scale=beta)) + offset_f,
            "hpd_lower": lo + offset_f,
            "hpd_upper": hi + offset_f,
            "interval_type": "hpd",
        }

    if dt == "chi_square":
        dof = float(parameters.get("dof", 1.0))
        lo, hi = _analytic_hpd(sp_stats.chi2(dof))
        return {
            "mean": dof + offset_f,
            "median": float(sp_stats.chi2.ppf(0.5, dof)) + offset_f,
            "hpd_lower": lo + offset_f,
            "hpd_upper": hi + offset_f,
            "interval_type": "hpd",
        }

    # Remaining known types raise a descriptive error (one_on_x and unknown)
    _frozen_distribution(dist_type, parameters, offset_f)
    raise ValueError(f"Unsupported distribution type: {dist_type}")


def compute_wasserstein_distance(
    dist1_type: str,
    dist1_params: Dict[str, Any],
    dist1_offset: float,
    dist2_type: str,
    dist2_params: Dict[str, Any],
    dist2_offset: float,
    n_samples: int = 10000,
    seed: int = DEFAULT_SEED,
) -> float:
    """Compute the Wasserstein (Earth Mover's) distance between two distributions.

    Uses scipy.stats.wasserstein_distance on deterministic sampled data.

    Args:
        dist1_type: First distribution type.
        dist1_params: First distribution parameters.
        dist1_offset: First distribution offset.
        dist2_type: Second distribution type.
        dist2_params: Second distribution parameters.
        dist2_offset: Second distribution offset.
        n_samples: Number of samples for Monte Carlo approximation.
        seed: Seed for deterministic sampling.

    Returns:
        The Wasserstein distance.
    """
    from scipy import stats as sp_stats

    samples1 = _sample_distribution(dist1_type, dist1_params, dist1_offset, n_samples, seed=seed)
    samples2 = _sample_distribution(dist2_type, dist2_params, dist2_offset, n_samples, seed=seed)

    return float(sp_stats.wasserstein_distance(samples1, samples2))


def _sample_distribution(
    dist_type: str,
    parameters: Dict[str, Any],
    offset: float,
    n: int,
    seed: int = DEFAULT_SEED,
) -> "np.ndarray":
    """Sample n values from a distribution (deterministic given seed).

    Raises:
        ValueError: If the distribution cannot be sampled (OneOnX).
    """
    dist = _frozen_distribution(dist_type, parameters, offset=0.0)
    rng = np.random.default_rng(seed)
    return dist.rvs(size=n, random_state=rng) + offset


def sanitize_id(name: str) -> str:
    """Sanitize a name for use as an XML ID.

    Replaces spaces and special characters with underscores.

    Args:
        name: The input name.

    Returns:
        A sanitized ID string.
    """
    return re.sub(r"[^a-zA-Z0-9_.\-]", "_", name)


def build_namespace() -> str:
    """Build the default BEAST2 v2.7.x namespace string.

    Returns:
        Colon-separated namespace string.
    """
    return ":".join(
        [
            "beast.base.evolution.alignment",
            "beast.base.evolution.speciation",
            "beast.pkgmgmt",
            "beast.base.core",
            "beast.base.inference",
            "beast.base.evolution.tree.coalescent",
            "beast.base.evolution.tree",
            "beast.base.evolution.likelihood",
            "beast.base.inference.distribution",
            "beast.base.evolution.operator",
            "beast.base.inference.operator",
            "beast.base.evolution.sitemodel",
            "beast.base.evolution.substitutionmodel",
            "beast.base.evolution.branchratemodel",
            "beast.base.inference.parameter",
            "beast.base.math",
            "beast.base.util",
        ]
    )

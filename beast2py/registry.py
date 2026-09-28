"""Model registry for Beast2Py.

Implements a Registry Pattern that maps configuration model names to
their XML spec values and parameter specifications. This enables
extensible model support without modifying core code.
"""

from __future__ import annotations

from typing import Any, Dict, Optional


class ModelRegistry:
    """Registry mapping model names to BEAST2 XML spec values.

    Supports registration of substitution models, clock models, tree priors,
    and calibration distributions.
    """

    # ---- Substitution Models ----
    # Maps config name -> (XML spec, parameter specs)
    # parameter specs: list of (param_name, xml_attribute, dimension, default_value)
    _substitution_models: Dict[str, Dict[str, Any]] = {
        "hky": {
            "spec": "HKY",
            "params": {"kappa": {"attr": "kappa", "dim": 1, "default": 2.0, "lower": 0.0}},
            "frequencies": True,
            "is_aminoacid": False,
        },
        "gtr": {
            "spec": "GTR",
            "params": {
                "rateAC": {"attr": "rateAC", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateAG": {"attr": "rateAG", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateAT": {"attr": "rateAT", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateCG": {"attr": "rateCG", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateCT": {"attr": "rateCT", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateGT": {"attr": "rateGT", "dim": 1, "default": 1.0, "lower": 0.0},
            },
            "frequencies": True,
            "is_aminoacid": False,
        },
        "jc69": {
            "spec": "JukesCantor",
            "params": {},
            "frequencies": False,
            "is_aminoacid": False,
        },
        "tn93": {
            "spec": "TN93",
            "params": {
                "kappa1": {"attr": "kappa1", "dim": 1, "default": 1.0, "lower": 0.0},
                "kappa2": {"attr": "kappa2", "dim": 1, "default": 1.0, "lower": 0.0},
            },
            "frequencies": True,
            "is_aminoacid": False,
        },
        "sym": {
            # SYM: GTR with equal base frequencies. All six rates are free
            # (BEAST normalises them); frequencies are forced uniform via
            # <frequencies estimate="false">, see model.py.
            "spec": "SYM",
            "params": {
                "rateAC": {"attr": "rateAC", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateAG": {"attr": "rateAG", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateAT": {"attr": "rateAT", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateCG": {"attr": "rateCG", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateCT": {"attr": "rateCT", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateGT": {"attr": "rateGT", "dim": 1, "default": 1.0, "lower": 0.0},
            },
            "frequencies": "uniform",
            "is_aminoacid": False,
        },
        "tim": {
            # TIM2: rateAG, rateCT and two transversion rates (BEAST2 TIM inputs)
            "spec": "TIM",
            "params": {
                "rateAG": {"attr": "rateAG", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateCT": {"attr": "rateCT", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateTransversions1": {
                    "attr": "rateTransversions1",
                    "dim": 1,
                    "default": 1.0,
                    "lower": 0.0,
                },
                "rateTransversions2": {
                    "attr": "rateTransversions2",
                    "dim": 1,
                    "default": 1.0,
                    "lower": 0.0,
                },
            },
            "frequencies": True,
            "is_aminoacid": False,
        },
        "tvm": {
            # TVM: four transversion rates plus a shared transition rate (BEAST2 TVM inputs)
            "spec": "TVM",
            "params": {
                "rateAC": {"attr": "rateAC", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateAT": {"attr": "rateAT", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateCG": {"attr": "rateCG", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateGT": {"attr": "rateGT", "dim": 1, "default": 1.0, "lower": 0.0},
                "rateTransitions": {
                    "attr": "rateTransitions",
                    "dim": 1,
                    "default": 1.0,
                    "lower": 0.0,
                },
            },
            "frequencies": True,
            "is_aminoacid": False,
        },
        # Amino acid models: EmpiricalSubstitutionModel instances with built-in
        # empirical equilibrium frequencies, so no <frequencies> element is
        # generated (their frequencies input is optional in BEAST2).
        "wag": {
            "spec": "WAG",
            "params": {},
            "frequencies": False,
            "is_aminoacid": True,
        },
        "jtt": {
            "spec": "JTT",
            "params": {},
            "frequencies": False,
            "is_aminoacid": True,
        },
        "dayhoff": {
            "spec": "Dayhoff",
            "params": {},
            "frequencies": False,
            "is_aminoacid": True,
        },
        "blosum62": {
            "spec": "Blosum62",
            "params": {},
            "frequencies": False,
            "is_aminoacid": True,
        },
        "cprev": {
            "spec": "CPREV",
            "params": {},
            "frequencies": False,
            "is_aminoacid": True,
        },
        "mtrev": {
            "spec": "MTREV",
            "params": {},
            "frequencies": False,
            "is_aminoacid": True,
        },
    }

    # ---- Calibration Distributions ----
    # Maps config name -> (XML spec, parameter mapping)
    _distributions: Dict[str, Dict[str, Any]] = {
        "normal": {
            "spec": "beast.base.inference.distribution.Normal",
            "params": {"mean": "mean", "sigma": "sigma"},
            "use_param_elements": True,
        },
        "lognormal": {
            "spec": "beast.base.inference.distribution.LogNormalDistributionModel",
            "params": {"M": "M", "S": "S"},
            "use_param_elements": True,
            "bool_attrs": {"mean_in_real_space": "meanInRealSpace"},
        },
        "uniform": {
            "spec": "beast.base.inference.distribution.Uniform",
            "params": {"lower": "lower", "upper": "upper"},
            "use_param_elements": False,
        },
        "exponential": {
            "spec": "beast.base.inference.distribution.Exponential",
            "params": {"mean": "mean"},
            "use_param_elements": True,
        },
        "gamma": {
            "spec": "beast.base.inference.distribution.Gamma",
            "params": {"alpha": "alpha", "beta": "beta"},
            "use_param_elements": True,
        },
        "beta": {
            "spec": "beast.base.inference.distribution.Beta",
            "params": {"alpha": "alpha", "beta": "beta"},
            "use_param_elements": True,
        },
        "one_on_x": {
            "spec": "beast.base.inference.distribution.OneOnX",
            "params": {},
            "use_param_elements": False,
        },
        "laplace": {
            "spec": "beast.base.inference.distribution.LaplaceDistribution",
            "params": {"mu": "mu", "scale": "scale"},
            "use_param_elements": True,
        },
        "inverse_gamma": {
            "spec": "beast.base.inference.distribution.InverseGamma",
            "params": {"alpha": "alpha", "beta": "beta"},
            "use_param_elements": True,
        },
        "poisson": {
            "spec": "beast.base.inference.distribution.Poisson",
            "params": {"lambda": "lambda"},
            "use_param_elements": True,
        },
        "chi_square": {
            # BEAST2 ChiSquare input name is "df" (not "dof")
            "spec": "beast.base.inference.distribution.ChiSquare",
            "params": {"dof": "df"},
            "use_param_elements": True,
        },
    }

    @classmethod
    def get_substitution_model(cls, name: str) -> Optional[Dict[str, Any]]:
        """Get substitution model specification by name.

        Args:
            name: Model name (case-insensitive).

        Returns:
            Model specification dict, or None if not found.
        """
        return cls._substitution_models.get(name.lower())

    @classmethod
    def get_distribution_spec(cls, name: str) -> Optional[Dict[str, Any]]:
        """Get calibration distribution specification by name.

        Args:
            name: Distribution type name (case-insensitive).

        Returns:
            Distribution specification dict, or None if not found.
        """
        # Try direct lookup
        result = cls._distributions.get(name.lower())
        if result:
            return result
        # Try with underscores converted
        result = cls._distributions.get(name.lower().replace(" ", "_"))
        return result

    @classmethod
    def register_substitution_model(cls, name: str, spec: Dict[str, Any]) -> None:
        """Register a custom substitution model.

        Args:
            name: Model name.
            spec: Model specification dict.
        """
        cls._substitution_models[name.lower()] = spec

    @classmethod
    def register_distribution(cls, name: str, spec: Dict[str, Any]) -> None:
        """Register a custom distribution.

        Args:
            name: Distribution name.
            spec: Distribution specification dict.
        """
        cls._distributions[name.lower()] = spec

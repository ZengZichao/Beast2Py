"""Model assembly module for Beast2Py.

Builds XML elements for substitution models, site models, clock models,
tree priors, operators, and state space using xml.etree.ElementTree.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional, Tuple

from .models import (
    BEASTConfig,
    CalibrationMethod,
    ClockModelConfig,
    ClockModelType,
    DataType,
    InitTreeType,
    RealParameter,
    SiteModelConfig,
    TreePriorType,
)
from .partition import PartitionManager
from .registry import ModelRegistry
from .utils import format_number, make_xml_id, ref, sanitize_id

# --- Parameter roles and their intrinsic support ---------------------------
#
# bounds came from whatever the YAML happened to contain, so a
# user who forgot ``lower`` on a frequency vector let an additive
# DeltaExchangeOperator propose negative frequencies (likelihood = NaN), and
# ``bdSkylineSamplingRate`` had no upper bound although a sampling proportion
# must lie in [0, 1]. The support below is a property of the *statistical
# role*, not of the file, so it is applied unless the user overrides it.
PARAM_SUPPORT = {
    "rate": (0.0, None),  # strictly positive scale/ratio
    "probability": (0.0, 1.0),
    "frequency": (0.0, 1.0),
    "integer": (1.0, None),  # counts / group sizes
    "boolean": (None, None),
    "signed": (None, None),  # may be negative (exponential growth rate)
    "tree": (None, None),
}


def _state_node(
    param_id: str,
    raw: Any,
    *,
    role: str = "rate",
    default_value: float = 1.0,
    dimension: int = 1,
) -> Tuple[str, Dict[str, Any]]:
    """Build one state-node record, honouring the user's own settings.

    ``estimate`` was read nowhere in this module and every state
    node was hard-coded to ``"estimate": True``, so ``kappa: {estimate: false}``
    still entered the state space, got an operator, and appeared in the trace.
    ``upper`` was dropped on the clock and gamma-shape nodes.

    Args:
        param_id: XML id of the state node.
        raw: The user's value: a number, a string vector, or a mapping with
            ``value``/``lower``/``upper``/``dimension``/``estimate``.
        role: Key into :data:`PARAM_SUPPORT` supplying default bounds.
        default_value: Value used when the user gives none.
        dimension: Default dimension (overridden by a vector value or mapping).

    Returns:
        ``(param_id, params)`` as consumed by ``build_state``. ``params``
        carries ``estimate`` so callers can define the parameter inline when it
        is fixed, and ``role`` so the writer can attach a default prior.
    """
    lower_default, upper_default = PARAM_SUPPORT.get(role, (None, None))
    value = default_value
    lower = lower_default
    upper = upper_default
    dim = dimension
    estimate = True

    if isinstance(raw, RealParameter):
        source: Dict[str, Any] = {
            "value": raw.value,
            "lower": raw.lower,
            "upper": raw.upper,
            "dimension": raw.dimension,
            "estimate": raw.estimate,
        }
    elif isinstance(raw, dict):
        source = raw
    elif raw is None:
        source = {}
    else:
        source = {"value": raw}

    if source.get("value") is not None:
        value = source["value"]
    if source.get("lower") is not None:
        lower = source["lower"]
    if source.get("upper") is not None:
        upper = source["upper"]
    if source.get("dimension") is not None:
        dim = int(source["dimension"])
    if "estimate" in source and source["estimate"] is not None:
        estimate = _as_bool(source["estimate"], param_id)

    # A space-separated vector implies its own dimension.
    if isinstance(value, str) and len(value.split()) > 1:
        tokens = value.split()
        if source.get("dimension") is None and dim == 1:
            dim = len(tokens)
        elif dim != len(tokens):
            raise ValueError(
                f"{param_id}: value lists {len(tokens)} entries but dimension is {dim}"
            )

    return (
        param_id,
        {
            "value": value,
            "dimension": dim,
            "lower": lower,
            "upper": upper,
            "estimate": estimate,
            "role": role,
        },
    )


def _as_bool(value: Any, where: str) -> bool:
    """Strict boolean coercion with a :class:`ValueError` for bad input.

    Args:
        value: Raw value.
        where: Location used in the error message.

    Returns:
        The boolean.

    Raises:
        ValueError: If the value is not an unambiguous boolean.
    """
    from .utils import as_strict_bool

    try:
        return as_strict_bool(value, where)
    except ValueError as exc:
        raise ValueError(str(exc)) from exc


def _normalised_frequencies(value: Any, dimension: int, where: str) -> str:
    """Validate and normalise an equilibrium-frequency vector.

    A vector that does not sum to 1 is renormalised or rejected rather than
    written out verbatim, and an explicit ``0.0`` entry stays ``0.0`` instead
    of falling back to a default through Python falsiness.

    Args:
        value: Scalar, list, or whitespace-separated string of frequencies.
        dimension: Expected number of frequencies.
        where: Location used in error messages.

    Returns:
        A space-separated vector of ``dimension`` frequencies summing to 1.

    Raises:
        ValueError: On a wrong length, a non-positive entry, or a total too far
            from 1 to be a normalisation rounding artefact.
    """
    if isinstance(value, (list, tuple)):
        tokens = [float(v) for v in value]
    else:
        tokens = [float(tok) for tok in str(value).split()]
    if len(tokens) == 1 and dimension != 1:
        tokens = tokens * dimension
    if len(tokens) != dimension:
        raise ValueError(
            f"{where}: expected {dimension} frequencies, got {len(tokens)} ({value!r})"
        )
    if any(tok <= 0 for tok in tokens):
        raise ValueError(f"{where}: every frequency must be > 0, got {tokens}")
    total = sum(tokens)
    if abs(total - 1.0) > 0.02:
        raise ValueError(
            f"{where}: frequencies sum to {total:.6g}, not 1.0. Supply a proper "
            f"probability vector (a rounding tolerance of 0.02 is accepted and "
            f"the vector is then renormalised)."
        )
    return " ".join(format_number(tok / total) for tok in tokens)


def _rate_value(raw: Any, where: str) -> float:
    """Read a numeric value out of a tree-prior parameter setting.

    Args:
        raw: Number, or mapping carrying ``value``.
        where: Location used in the error message.

    Returns:
        The float value.

    Raises:
        ValueError: If the value is not numeric.
    """
    if isinstance(raw, dict):
        raw = raw.get("value", 1.0)
    if raw is None or isinstance(raw, bool):
        raise ValueError(f"tree_prior.{where}: expected a number, got {raw!r}")
    return float(raw)


def _rate_estimate(raw: Any, where: str) -> bool:
    """Read the ``estimate`` flag of a tree-prior parameter setting.

    Args:
        raw: Number, or mapping carrying ``estimate``.
        where: Location used in the error message.

    Returns:
        True when the derived state node should be estimated.
    """
    if isinstance(raw, dict) and raw.get("estimate") is not None:
        return _as_bool(raw["estimate"], f"tree_prior.{where}.estimate")
    return True


def _estimate_calibration_age(dist) -> float:
    """Rough point age of a calibration distribution (for initial-tree sizing).

    Returns 0.0 when no reasonable point estimate exists. Only used to pick a
    starting root height, never for the model itself.
    """
    if dist is None:
        return 0.0
    p = dist.parameters
    t = str(dist.type).lower()
    try:
        if t == "normal":
            age = float(p.get("mean", 0.0))
        elif t in ("lognormal", "log_normal"):
            age = float(p.get("M", 0.0))
        elif t == "uniform":
            age = float(p.get("upper", 0.0))
        elif t == "exponential":
            age = 2.0 * float(p.get("mean", 0.0))
        elif t == "gamma":
            age = float(p.get("alpha", 1.0)) * float(p.get("beta", 1.0))
        else:
            age = 0.0
    except (TypeError, ValueError):
        return 0.0
    return age + float(dist.offset or 0.0)


class ModelBuilder:
    """Assemble XML elements for all model components."""

    @staticmethod
    def build_substitution_model(
        partition_id: str,
        site_model: SiteModelConfig,
        data_id: Optional[str] = None,
        data_type: Optional[DataType] = None,
    ) -> Tuple[ET.Element, List[Tuple[str, Dict[str, Any]]]]:
        """Build substitution model XML element.

        Args:
            partition_id: Partition identifier.
            site_model: Site model configuration.
            data_id: Alignment ID (used for uniform/fitted frequencies).
            data_type: Alignment data type (nucleotide or aminoacid).

        Returns:
            Tuple of (XML element, list of (stateNode ID, params) tuples).
        """
        subst_config = dict(site_model.substitution_model)  # copy: never mutate the config
        model_type = str(subst_config.get("type", "hky")).lower()
        model_spec = ModelRegistry.get_substitution_model(model_type)

        if model_spec is None:
            raise ValueError(f"Unknown substitution model: {model_type}")

        # --- Data-type guard  ---
        # The previous guard tested `freqs_mode is not True` on models whose
        # registry entry stores frequencies as the literal True, so it could
        # never fire and WAG on nucleotide data was emitted unchecked.
        is_aa_model = bool(model_spec.get("is_aminoacid"))
        if data_type is not None:
            if is_aa_model and data_type != DataType.AMINOACID:
                raise ValueError(
                    f"Substitution model '{model_type}' is an amino-acid matrix and is "
                    f"undefined on {data_type.value} data"
                )
            if not is_aa_model and data_type == DataType.AMINOACID:
                raise ValueError(
                    f"Substitution model '{model_type}' is a nucleotide model and cannot be "
                    f"applied to amino-acid data"
                )

        # Legacy convenience: a "rates" vector is expanded into the individual
        # rate parameters for models whose rates are free (GTR, SYM). Models
        # with constrained rate structure (TIM, TVM) must use their own named
        # parameters, since a 6-free-rate vector does not implement their
        # rate equalities.
        if "rates" in subst_config:
            rates_raw = subst_config.pop("rates")
            if model_type not in ("gtr", "sym"):
                raise ValueError(
                    f"Substitution model '{model_type}' does not take a 'rates' vector; "
                    f"rateAC..rateGT do not exist as free parameters for it. Use its own "
                    f"named rate parameters (see Supplementary Table S2)."
                )
            if isinstance(rates_raw, dict):
                rates_str = str(rates_raw.get("value", "1.0 1.0 1.0 1.0 1.0 1.0"))
                rates_lower = rates_raw.get("lower", 0.0)
                rates_upper = rates_raw.get("upper")
                rates_estimate = rates_raw.get("estimate", True)
            else:
                rates_str = str(rates_raw)
                rates_lower = 0.0
                rates_upper = None
                rates_estimate = True
            rate_values = rates_str.split()
            rate_names = ["rateAC", "rateAG", "rateAT", "rateCG", "rateCT", "rateGT"]
            # a short vector used to be silently padded with 1.0
            # and a long one truncated.
            if len(rate_values) != 6:
                raise ValueError(
                    f"'rates' for a {model_type.upper()} model must list exactly 6 values "
                    f"(rateAC rateAG rateAT rateCG rateCT rateGT); got {len(rate_values)}: "
                    f"{rates_str!r}"
                )
            clash = [name for name in rate_names if name in subst_config]
            if clash:
                raise ValueError(
                    f"'rates' cannot be combined with explicitly named rate parameters "
                    f"{', '.join(clash)}: the expansion would silently overwrite one of them"
                )
            for i, rate_name in enumerate(rate_names):
                subst_config[rate_name] = {
                    "value": float(rate_values[i]),
                    "lower": rates_lower,
                    "upper": rates_upper,
                    "estimate": rates_estimate,
                }

        spec = model_spec["spec"]
        state_nodes: List[Tuple[str, Dict[str, Any]]] = []

        # Create substitution model element
        subst_elem = ET.Element(
            "substModel", {"spec": spec, "id": make_xml_id(f"{model_type}", partition_id)}
        )

        # Add model-specific parameters; every parameter is wired to the
        # substitution model through the XML attribute named in the registry
        # (kappa, kappa1, rateAC, rateTransversions1, ...). A *fixed* parameter
        # is instead nested as an inline <parameter estimate="false"> so it
        # never enters the state space.
        for param_name, param_spec in model_spec.get("params", {}).items():
            attr_name = param_spec["attr"]
            dim = param_spec.get("dim", 1)
            default_val = param_spec.get("default", 1.0)
            lower = param_spec.get("lower")

            user_val = subst_config.get(param_name)
            merged: Dict[str, Any] = {"value": default_val, "dimension": dim}
            if lower is not None:
                merged["lower"] = lower
            if isinstance(user_val, dict):
                merged.update(user_val)
            elif user_val is not None:
                merged["value"] = user_val

            param_id = make_xml_id(f"{model_type}.{param_name}", partition_id)
            node_id, node = _state_node(param_id, merged, role="rate", dimension=dim)
            if node["estimate"]:
                subst_elem.set(attr_name, ref(param_id))
                state_nodes.append((node_id, node))
            else:
                inline = ET.SubElement(
                    subst_elem,
                    "parameter",
                    {
                        "id": param_id,
                        "name": attr_name,
                        "value": format_number(node["value"]),
                        "estimate": "false",
                        "dimension": str(node["dimension"]),
                    },
                )
                if node.get("lower") is not None:
                    inline.set("lower", format_number(node["lower"]))
                if node.get("upper") is not None:
                    inline.set("upper", format_number(node["upper"]))

        # --- Equilibrium frequencies ---
        # Registry values: True = estimated vector, "uniform" = forced uniform
        # (SYM), False = the matrix carries its own empirical frequencies.
        freqs_mode = model_spec.get("frequencies", False)
        if freqs_mode and data_type == DataType.AMINOACID and freqs_mode is not True:
            raise ValueError(f"Substitution model '{model_type}' does not support amino acid data")

        if freqs_mode and data_type != DataType.AMINOACID:
            freq_raw = subst_config.get("frequencies", {})
            freq_opts: Dict[str, Any] = dict(freq_raw) if isinstance(freq_raw, dict) else {}
            if not isinstance(freq_raw, dict) and freq_raw is not None:
                freq_opts = {"value": freq_raw}
            freq_mode = str(freq_opts.pop("mode", "estimated")).lower()
            if freq_mode not in ("estimated", "uniform", "empirical", "fixed"):
                raise ValueError(
                    f"frequencies.mode must be estimated, uniform, empirical or fixed; "
                    f"got {freq_mode!r}"
                )
            freq_estimate = _as_bool(freq_opts.get("estimate", True), f"{partition_id}.frequencies")
            if freq_mode == "estimated":
                freq_mode = "fixed" if not freq_estimate else "estimated"

            if freq_mode in ("uniform", "empirical"):
                # <frequencies spec="Frequencies" data="@aln" estimate=...>:
                # estimate=false is *uniform over characters* and estimate=true
                # is empirical counting (Frequencies.java).
                if not data_id:
                    raise ValueError(
                        f"data_id is required for {freq_mode} frequencies on '{model_type}'"
                    )
                ET.SubElement(
                    subst_elem,
                    "frequencies",
                    {
                        "id": make_xml_id("freqs", partition_id),
                        "spec": "Frequencies",
                        "data": ref(data_id),
                        "estimate": "true" if freq_mode == "empirical" else "false",
                    },
                )
            else:
                # An explicit user vector, either estimated or written out as
                # given. estimate=false used to be routed to
                # <frequencies estimate="false">, which BEAST2 documents as
                # "uniform over characters" — the user's numbers were discarded.
                freq_id = make_xml_id("freqParameter", partition_id)
                freqs_child = ET.SubElement(
                    subst_elem,
                    "frequencies",
                    {
                        "id": make_xml_id("freqs", partition_id),
                        "spec": "Frequencies",
                    },
                )
                # NB: Frequencies accepts *either* data= (counting or uniform)
                # *or* an explicit frequencies vector — supplying both is a
                # validation error in BEAST2.
                default_dim = 20 if data_type == DataType.AMINOACID else 4
                freq_node_id, freq_node = _state_node(
                    freq_id,
                    {**freq_opts, "dimension": freq_opts.get("dimension", default_dim)},
                    role="frequency",
                    default_value=1.0 / default_dim,
                    dimension=default_dim,
                )
                freq_value = _normalised_frequencies(
                    freq_node["value"], int(freq_node["dimension"]), f"{partition_id}.frequencies"
                )
                freq_node["value"] = freq_value
                if freq_node["estimate"]:
                    freqs_child.set("frequencies", ref(freq_id))
                    state_nodes.append((freq_node_id, freq_node))
                else:
                    ET.SubElement(
                        freqs_child,
                        "parameter",
                        {
                            "id": freq_id,
                            "name": "frequencies",
                            "value": freq_value,
                            "dimension": str(freq_node["dimension"]),
                            "lower": format_number(freq_node["lower"]),
                            "upper": format_number(freq_node["upper"]),
                            "estimate": "false",
                        },
                    )

        return subst_elem, state_nodes

    @staticmethod
    def build_site_model(
        partition_id: str,
        site_model: SiteModelConfig,
        data_id: Optional[str] = None,
        data_type: Optional[DataType] = None,
    ) -> Tuple[ET.Element, List[Tuple[str, Dict[str, Any]]]]:
        """Build site model XML element.

        Args:
            partition_id: Partition identifier.
            site_model: Site model configuration.
            data_id: Alignment ID (passed through to the substitution model).
            data_type: Alignment data type.

        Returns:
            Tuple of (XML element, list of (stateNode ID, params) tuples).
        """
        site_model_id = make_xml_id("siteModel", partition_id)
        state_nodes: List[Tuple[str, Dict[str, Any]]] = []

        sm_elem = ET.Element("siteModel", {"id": site_model_id, "spec": "SiteModel"})

        # Gamma categories
        gamma_cats = site_model.gamma_categories
        sm_elem.set("gammaCategoryCount", str(gamma_cats))

        if gamma_cats > 0 and site_model.gamma_shape:
            gs = site_model.gamma_shape
            shape_id = make_xml_id("gammaShape", partition_id)
            _, node = _state_node(shape_id, gs, role="rate", default_value=0.5)
            if gs.estimate:
                ET.SubElement(
                    sm_elem,
                    "parameter",
                    {"idref": shape_id, "name": "shape"},
                )
                state_nodes.append((shape_id, node))
            else:
                # Fixed parameter: define inline so it stays out of the state
                attrs = {
                    "id": shape_id,
                    "name": "shape",
                    "value": format_number(node["value"]),
                    "estimate": "false",
                }
                if node.get("lower") is not None:
                    attrs["lower"] = format_number(node["lower"])
                if node.get("upper") is not None:
                    attrs["upper"] = format_number(node["upper"])
                ET.SubElement(sm_elem, "parameter", attrs)

        # Proportion invariant
        if site_model.proportion_invariant:
            pi = site_model.proportion_invariant
            pi_id = make_xml_id("proportionInvariant", partition_id)
            _, node = _state_node(pi_id, pi, role="probability", default_value=0.0)
            if pi.estimate:
                ET.SubElement(
                    sm_elem,
                    "parameter",
                    {"idref": pi_id, "name": "proportionInvariant"},
                )
                state_nodes.append((pi_id, node))
            else:
                attrs = {
                    "id": pi_id,
                    "name": "proportionInvariant",
                    "value": format_number(node["value"]),
                    "estimate": "false",
                }
                if node.get("lower") is not None:
                    attrs["lower"] = format_number(node["lower"])
                if node.get("upper") is not None:
                    attrs["upper"] = format_number(node["upper"])
                ET.SubElement(sm_elem, "parameter", attrs)

        # Add substitution model
        subst_elem, subst_state_nodes = ModelBuilder.build_substitution_model(
            partition_id, site_model, data_id=data_id, data_type=data_type
        )
        sm_elem.append(subst_elem)
        state_nodes.extend(subst_state_nodes)

        return sm_elem, state_nodes

    @staticmethod
    def build_clock_model(
        partition_id: str,
        clock: ClockModelConfig,
        tree_id: str,
        n_taxa: int,
    ) -> Tuple[ET.Element, List[Tuple[str, Dict[str, Any]]]]:
        """Build clock model XML element.

        Args:
            partition_id: Partition identifier.
            clock: Clock model configuration.
            tree_id: Tree ID.
            n_taxa: Number of taxa (for rateCategories dimension).

        Returns:
            Tuple of (XML element, list of (stateNode ID, params) tuples).
        """
        clock_id = make_xml_id("clockModel", partition_id)
        state_nodes: List[Tuple[str, Dict[str, Any]]] = []

        if clock.type == ClockModelType.STRICT:
            clock_elem = ET.Element(
                "branchRateModel",
                {"id": clock_id, "spec": "beast.base.evolution.branchratemodel.StrictClockModel"},
            )

            # Clock rate parameter
            rate_id = make_xml_id("clockRate", partition_id)
            cr = clock.clock_rate or RealParameter(value=1.0)
            _, node = _state_node(rate_id, cr, role="rate", default_value=1.0)
            if node["estimate"]:
                # Estimated: reference state node via idref
                ET.SubElement(
                    clock_elem,
                    "parameter",
                    {
                        "idref": rate_id,
                        "name": "clock.rate",
                    },
                )
                state_nodes.append((rate_id, node))
            else:
                # Fixed: define inline with id
                ET.SubElement(
                    clock_elem,
                    "parameter",
                    {
                        "id": rate_id,
                        "name": "clock.rate",
                        "value": format_number(node["value"]),
                        "lower": format_number(node["lower"]),
                        "estimate": "false",
                    },
                )

        elif clock.type in (ClockModelType.UCLN, ClockModelType.UCE):
            clock_elem = ET.Element(
                "branchRateModel",
                {
                    "id": clock_id,
                    "spec": "beast.base.evolution.branchratemodel.UCRelaxedClockModel",
                    "normalize": "true",
                },
            )

            # Clock rate (ucld.mean)
            mean_id = make_xml_id("ucld.mean", partition_id)
            um = clock.ucld_mean or RealParameter(value=1.0)
            _, mean_node = _state_node(mean_id, um, role="rate", default_value=1.0)
            if mean_node["estimate"]:
                ET.SubElement(
                    clock_elem,
                    "parameter",
                    {
                        "idref": mean_id,
                        "name": "clock.rate",
                    },
                )
                state_nodes.append((mean_id, mean_node))
            else:
                ET.SubElement(
                    clock_elem,
                    "parameter",
                    {
                        "id": mean_id,
                        "name": "clock.rate",
                        "value": format_number(mean_node["value"]),
                        "lower": format_number(mean_node["lower"]),
                        "estimate": "false",
                    },
                )

            # Distribution for relaxed clock
            distr_id = make_xml_id("ucld.distr", partition_id)
            if clock.type == ClockModelType.UCLN:
                distr_elem = ET.SubElement(
                    clock_elem,
                    "distr",
                    {
                        "id": distr_id,
                        "spec": "beast.base.inference.distribution.LogNormalDistributionModel",
                    },
                )
                ET.SubElement(distr_elem, "parameter", {"name": "M", "value": "1.0"})
                # ucld.stdev
                #
                # "probability" clamped this to (0, 1], which rejected perfectly
                # ordinary UCLD settings: S is a log-scale standard deviation of
                # among-branch rates, not a proportion, and S > 1 is common for
                # strongly rate-varying phylogenies.
                stdev_id = make_xml_id("ucld.stdev", partition_id)
                us = clock.ucld_stdev or RealParameter(value=0.333)
                _, stdev_node = _state_node(stdev_id, us, role="rate", default_value=0.333)
                if stdev_node["estimate"]:
                    ET.SubElement(
                        distr_elem,
                        "parameter",
                        {
                            "idref": stdev_id,
                            "name": "S",
                        },
                    )
                    state_nodes.append((stdev_id, stdev_node))
                else:
                    ET.SubElement(
                        distr_elem,
                        "parameter",
                        {
                            "id": stdev_id,
                            "name": "S",
                            "value": format_number(stdev_node["value"]),
                            "lower": format_number(stdev_node["lower"]),
                            "upper": format_number(stdev_node["upper"]),
                            "estimate": "false",
                        },
                    )
            else:  # UCE
                distr_elem = ET.SubElement(
                    clock_elem,
                    "distr",
                    {"id": distr_id, "spec": "beast.base.inference.distribution.Exponential"},
                )
                ET.SubElement(distr_elem, "parameter", {"name": "mean", "value": "1.0"})

            # Tree input (required by UCRelaxedClockModel)
            ET.SubElement(clock_elem, "tree", {"idref": tree_id})

            # Rate categories: one per branch of the rooted tree (2n-2).
            # BEAST2 auto-expands the dimension, but emitting the correct
            # dimension avoids a runtime warning and initialisation ambiguity.
            branch_count = max(1, 2 * n_taxa - 2)
            cats_id = make_xml_id("rateCategories", partition_id)
            ET.SubElement(
                clock_elem,
                "parameter",
                {
                    "idref": cats_id,
                    "name": "rateCategories",
                },
            )
            state_nodes.append(
                (
                    cats_id,
                    {
                        "value": 1,
                        "dimension": branch_count,
                        "estimate": True,
                        "lower": 0.0,
                        "spec": "parameter.IntegerParameter",
                        "role": "integer",
                    },
                )
            )

        elif clock.type == ClockModelType.RLC:
            branch_count = max(1, 2 * n_taxa - 2)
            clock_elem = ET.Element(
                "branchRateModel",
                {
                    "id": clock_id,
                    "spec": "beast.base.evolution.branchratemodel.RandomLocalClockModel",
                    "ratesAreMultipliers": "false",
                    "tree": ref(tree_id),
                },
            )

            # Rates: defined once in the state, referenced here
            rates_id = make_xml_id("rates", partition_id)
            ET.SubElement(
                clock_elem,
                "parameter",
                {"idref": rates_id, "name": "rates"},
            )
            state_nodes.append(
                (
                    rates_id,
                    {
                        "value": 1,
                        "dimension": branch_count,
                        "lower": 0.0,
                        "estimate": True,
                        "role": "rate",
                    },
                )
            )

            # Indicators: boolean state parameter, referenced here
            indicators_id = make_xml_id("indicators", partition_id)
            ET.SubElement(
                clock_elem,
                "parameter",
                {"idref": indicators_id, "name": "indicators"},
            )
            state_nodes.append(
                (
                    indicators_id,
                    {
                        "value": "false",
                        "dimension": branch_count,
                        "estimate": True,
                        "spec": "parameter.BooleanParameter",
                        "role": "boolean",
                    },
                )
            )

        else:
            raise ValueError(f"Unknown clock model type: {clock.type}")

        return clock_elem, state_nodes

    @staticmethod
    def build_tree_prior(
        config: BEASTConfig,
        tree_id: str,
        partition_mgr: PartitionManager,
    ) -> Tuple[ET.Element, List[Tuple[str, Dict[str, Any]]], List[str]]:
        """Build tree prior XML element.

        Args:
            config: BEAST configuration.
            tree_id: Tree ID.
            partition_mgr: Partition manager.

        Returns:
            Tuple of (tree prior XML element, list of (stateNode ID, params) tuples,
                       list of prior distribution IDs to include in prior compound).
        """
        tp_type = config.tree_prior_type
        params = config.tree_prior_params
        state_nodes: List[Tuple[str, Dict[str, Any]]] = []
        prior_ids: List[str] = []

        if tp_type == TreePriorType.YULE:
            prior_id = "yule"
            birth_rate_id = "birthRate"

            yule_elem = ET.Element(
                "distribution",
                {
                    "id": prior_id,
                    "spec": "beast.base.evolution.speciation.YuleModel",
                    "tree": ref(tree_id),
                },
            )

            ET.SubElement(
                yule_elem,
                "parameter",
                {
                    "name": "birthDiffRate",
                    "idref": birth_rate_id,
                },
            )

            _, br_node = _state_node(
                birth_rate_id, params.get("birth_rate", 1.0), role="rate", default_value=1.0
            )
            if not br_node["estimate"]:
                raise ValueError(
                    "tree_prior.birth_rate: the birth (or birth-difference) rate is the tree "
                    "prior's own parameter; fixing it leaves the prior without a value to "
                    "evaluate. Estimate it, or choose a different tree prior."
                )
            state_nodes.append((birth_rate_id, br_node))
            prior_ids.append(prior_id)

            return yule_elem, state_nodes, prior_ids

        elif tp_type == TreePriorType.BIRTH_DEATH:
            prior_id = "birthDeath"
            bdr_id = "birthDiffRate"
            rdr_id = "relativeDeathRate"
            sp_id = "sampleProbability"

            bd_elem = ET.Element(
                "distribution",
                {
                    "id": prior_id,
                    "spec": "beast.base.evolution.speciation.BirthDeathGernhard08Model",
                    "tree": ref(tree_id),
                },
            )

            # Get raw birth/death rates from config and compute derived params.
            #
            # the model is parameterised as birth - death and
            # death / birth, both of which live on a bounded domain. A death
            # rate >= birth rate used to be converted anyway, producing
            # birthDiffRate = -0.3 with lower="0.0" and relativeDeathRate = 1.6
            # with upper="1.0" — an initial state outside both of its own
            # bounds, which BEAST2 refuses at start-up while the generator
            # exited 0.
            br_raw = params.get("birth_rate", {"value": 1.0, "lower": 0.0})
            dr_raw = params.get("death_rate", {"value": 0.5, "lower": 0.0})
            birth_val = _rate_value(br_raw, "birth_rate")
            death_val = _rate_value(dr_raw, "death_rate")
            birth_estimate = _rate_estimate(br_raw, "birth_rate")
            death_estimate = _rate_estimate(dr_raw, "death_rate")
            if birth_estimate != death_estimate:
                raise ValueError(
                    "tree_prior: birth_rate and death_rate must be estimated or fixed together. "
                    "The BirthDeathGernhard08Model state nodes are the derived quantities "
                    "birth - death and death / birth, so fixing only one of the two inputs "
                    "leaves an estimated state node whose value depends on a constant that "
                    "carries no prior."
                )
            if birth_val <= 0:
                raise ValueError(f"tree_prior.birth_rate must be > 0, got {birth_val}")
            if death_val >= birth_val:
                raise ValueError(
                    f"tree_prior: death_rate ({death_val}) must be strictly below birth_rate "
                    f"({birth_val}) — see the BirthDeathGernhard08Model parameterisation; use "
                    f"tree_prior.type: bd_skyline_serial (BDSKY) for a critical or "
                    f"sub-critical process."
                )
            bdr_value = birth_val - death_val
            rdr_value = death_val / birth_val

            _, bdr_node = _state_node(bdr_id, bdr_value, role="rate", default_value=bdr_value)
            bdr_node["estimate"] = birth_estimate
            bdr_node["value"] = round(bdr_value, 10)
            _, rdr_node = _state_node(
                rdr_id, rdr_value, role="probability", default_value=rdr_value
            )
            rdr_node["estimate"] = birth_estimate
            rdr_node["value"] = round(rdr_value, 10)

            if birth_estimate:
                # Estimated: reference state nodes from the prior element.
                ET.SubElement(bd_elem, "parameter", {"name": "birthDiffRate", "idref": bdr_id})
                ET.SubElement(bd_elem, "parameter", {"name": "relativeDeathRate", "idref": rdr_id})
                state_nodes.append((bdr_id, bdr_node))
                state_nodes.append((rdr_id, rdr_node))
            else:
                # Both fixed: define the derived parameters inline so they never
                # become state nodes at all.
                for name, node, node_id in (
                    ("birthDiffRate", bdr_node, bdr_id),
                    ("relativeDeathRate", rdr_node, rdr_id),
                ):
                    attrs = {
                        "id": node_id,
                        "name": name,
                        "value": format_number(node["value"]),
                        "estimate": "false",
                    }
                    if node.get("lower") is not None:
                        attrs["lower"] = format_number(node["lower"])
                    if node.get("upper") is not None:
                        attrs["upper"] = format_number(node["upper"])
                    ET.SubElement(bd_elem, "parameter", attrs)

            # sampleProbability: BEAST2's Parameter.estimate defaults to true,
            # so this inline parameter used to be an estimated quantity that was
            # absent from <state> — an orphan the operator-coverage invariant
            # could not see. It is a data property (taxon sampling
            # fraction), so it is now explicitly fixed.
            sp_raw = params.get("sample_probability", params.get("sampling_rate", 1.0))
            sp_value = _rate_value(sp_raw, "sample_probability")
            if not 0.0 < sp_value <= 1.0:
                raise ValueError(
                    f"tree_prior.sample_probability must lie in (0, 1], got {sp_value}"
                )
            ET.SubElement(
                bd_elem,
                "parameter",
                {
                    "name": "sampleProbability",
                    "id": sp_id,
                    "value": format_number(sp_value),
                    "lower": "0.0",
                    "upper": "1.0",
                    "estimate": "false",
                },
            )

            prior_ids.append(prior_id)
            return bd_elem, state_nodes, prior_ids

        elif tp_type == TreePriorType.COALESCENT_CONSTANT:
            prior_id = "coalescent"
            pop_size_id = "popSize"

            coal_elem = ET.Element(
                "distribution",
                {"id": prior_id, "spec": "Coalescent"},
            )

            # Tree intervals
            ET.SubElement(
                coal_elem,
                "treeIntervals",
                {
                    "id": "TreeIntervals",
                    "spec": "beast.base.evolution.tree.TreeIntervals",
                    "tree": ref(tree_id),
                },
            )

            # Population model
            pop_model = ET.SubElement(
                coal_elem,
                "populationModel",
                {"spec": "ConstantPopulation", "id": "ConstantPopulation"},
            )

            ET.SubElement(
                pop_model,
                "parameter",
                {"name": "popSize", "idref": pop_size_id},
            )

            _, ps_node = _state_node(
                pop_size_id, params.get("pop_size", 1.0), role="rate", default_value=1.0
            )
            if not ps_node["estimate"]:
                raise ValueError(
                    "tree_prior.pop_size: the coalescent population size is the tree prior's "
                    "own parameter and cannot be fixed while the prior is evaluated."
                )
            state_nodes.append((pop_size_id, ps_node))
            prior_ids.append(prior_id)

            return coal_elem, state_nodes, prior_ids

        elif tp_type == TreePriorType.COALESCENT_EXPONENTIAL:
            prior_id = "coalescent"
            pop_size_id = "popSize"
            growth_rate_id = "growthRate"

            coal_elem = ET.Element(
                "distribution",
                {"id": prior_id, "spec": "Coalescent"},
            )

            ET.SubElement(
                coal_elem,
                "treeIntervals",
                {
                    "id": "TreeIntervals",
                    "spec": "beast.base.evolution.tree.TreeIntervals",
                    "tree": ref(tree_id),
                },
            )

            pop_model = ET.SubElement(
                coal_elem,
                "populationModel",
                {"spec": "ExponentialGrowth", "id": "ExponentialGrowth"},
            )

            ET.SubElement(
                pop_model,
                "parameter",
                {"name": "popSize", "idref": pop_size_id},
            )

            _, ps_node = _state_node(
                pop_size_id, params.get("pop_size", 1.0), role="rate", default_value=1.0
            )
            if not ps_node["estimate"]:
                raise ValueError("tree_prior.pop_size cannot be fixed for an estimated prior")
            state_nodes.append((pop_size_id, ps_node))

            ET.SubElement(
                pop_model,
                "parameter",
                {"name": "growthRate", "idref": growth_rate_id},
            )

            # growthRate is signed: an unbounded flat prior on it is improper,
            # so it carries role "signed" and gets a default Normal prior.
            _, gr_node = _state_node(
                growth_rate_id, params.get("growth_rate", 0.0), role="signed", default_value=0.0
            )
            state_nodes.append((growth_rate_id, gr_node))

            prior_ids.append(prior_id)
            return coal_elem, state_nodes, prior_ids

        elif tp_type == TreePriorType.CALIBRATED_YULE:
            prior_id = "calibratedYule"
            birth_rate_id = "birthRate"

            cyule_elem = ET.Element(
                "distribution",
                {
                    "spec": "beast.base.evolution.speciation.CalibratedYuleModel",
                    "type": "full",
                    "id": prior_id,
                },
            )

            ET.SubElement(
                cyule_elem,
                "parameter",
                {"name": "birthRate", "idref": birth_rate_id},
            )
            ET.SubElement(cyule_elem, "tree", {"idref": tree_id})

            _, br_node = _state_node(
                birth_rate_id, params.get("birth_rate", 1.0), role="rate", default_value=1.0
            )
            if not br_node["estimate"]:
                raise ValueError(
                    "tree_prior.birth_rate: CalibratedYuleModel evaluates the birth rate as its "
                    "own parameter and cannot be fixed."
                )
            state_nodes.append((birth_rate_id, br_node))
            prior_ids.append(prior_id)

            # Add calibration points as children of CalibratedYuleModel
            from .calibration import CalibrationBuilder

            first_aln_id = partition_mgr.get_first_alignment_id()
            for cal in config.calibrations:
                cal_elem = CalibrationBuilder.build_calibrated_yule_point(
                    cal, all_taxa_ref=first_aln_id
                )
                cyule_elem.append(cal_elem)

            return cyule_elem, state_nodes, prior_ids

        elif tp_type == TreePriorType.BAYESIAN_SKYLINE:
            prior_id = "bayesianSkyline"
            bsp_id = "bPopSizes"
            bsg_id = "bGroupSizes"

            bsp_elem = ET.Element(
                "distribution",
                {"id": prior_id, "spec": "BayesianSkyline"},
            )

            # Population sizes and group sizes; input names are "popSizes" and
            # "groupSizes" (child tag must match the input name).
            n_taxa = partition_mgr.get_n_taxa()
            n_groups = min(5, max(1, n_taxa - 1))

            ps_raw = params.get("pop_sizes", {})
            ps_dim = (
                int(ps_raw.get("dimension", n_groups)) if isinstance(ps_raw, dict) else n_groups
            )
            gs_raw = params.get("group_sizes", {})
            gs_dim = int(gs_raw.get("dimension", ps_dim)) if isinstance(gs_raw, dict) else ps_dim

            ET.SubElement(
                bsp_elem,
                "parameter",
                {"idref": bsp_id, "name": "popSizes"},
            )
            _, ps_node = _state_node(bsp_id, ps_raw, role="rate", dimension=ps_dim)
            state_nodes.append((bsp_id, ps_node))

            ET.SubElement(
                bsp_elem,
                "parameter",
                {"idref": bsg_id, "name": "groupSizes"},
            )
            # Group sizes are integers; BEAST2 auto-expands all-1 initial
            # values to cover the coalescent events.
            _, gs_node = _state_node(
                bsg_id, gs_raw, role="integer", default_value=1, dimension=gs_dim
            )
            gs_node["spec"] = "parameter.IntegerParameter"
            state_nodes.append((bsg_id, gs_node))

            # Tree intervals
            ET.SubElement(
                bsp_elem,
                "treeIntervals",
                {
                    "id": "TreeIntervals",
                    "spec": "beast.base.evolution.tree.TreeIntervals",
                    "tree": ref(tree_id),
                },
            )

            prior_ids.append(prior_id)
            return bsp_elem, state_nodes, prior_ids

        elif tp_type == TreePriorType.EBSP:
            # The Extended Bayesian Skyline Plot is a PopulationFunction, so it
            # must be wrapped in a Coalescent distribution. In BEAST 2.7 the
            # CompoundPopulationFunction is part of the base package.
            prior_id = "ebsp"
            pop_sizes_id = "populationSizes"
            pop_indicators_id = "populationIndicators"

            ebsp_elem = ET.Element(
                "distribution",
                {"id": prior_id, "spec": "Coalescent"},
            )

            n_taxa = partition_mgr.get_n_taxa()
            n_pop = max(2, n_taxa - 1)

            # Tree intervals: defined here, referenced from the population function
            ET.SubElement(
                ebsp_elem,
                "treeIntervals",
                {
                    "id": "TreeIntervals",
                    "spec": "beast.base.evolution.tree.TreeIntervals",
                    "tree": ref(tree_id),
                },
            )

            cpf_elem = ET.SubElement(
                ebsp_elem,
                "populationModel",
                {
                    "id": "ebsp.popFunction",
                    "spec": "beast.base.evolution.tree.coalescent.CompoundPopulationFunction",
                    "type": "linear",
                    "useIntervalsMiddle": "false",
                },
            )

            ET.SubElement(
                cpf_elem,
                "parameter",
                {"idref": pop_sizes_id, "name": "populationSizes"},
            )
            _, pop_node = _state_node(
                pop_sizes_id, params.get("pop_sizes", 1.0), role="rate", dimension=n_pop
            )
            pop_node["dimension"] = n_pop
            state_nodes.append((pop_sizes_id, pop_node))

            ET.SubElement(
                cpf_elem,
                "parameter",
                {"idref": pop_indicators_id, "name": "populationIndicators"},
            )
            state_nodes.append(
                (
                    pop_indicators_id,
                    {
                        "value": "true",
                        "dimension": n_pop - 1,
                        "estimate": True,
                        "spec": "parameter.BooleanParameter",
                        "role": "boolean",
                    },
                )
            )

            ET.SubElement(cpf_elem, "itree", {"idref": "TreeIntervals"})

            prior_ids.append(prior_id)
            return ebsp_elem, state_nodes, prior_ids

        elif tp_type == TreePriorType.BD_SKYLINE_SERIAL:
            # Serial-sampling birth-death skyline prior. Requires the bdsky
            # add-on package (github.com/BEAST2-Dev/bdsky); conditioned on the
            # root so no origin parameter is needed. All rates are relative
            # to the first interval (one skyline interval by default).
            prior_id = "bdSkylineSerial"
            birth_rate_id = "bdSkylineBirthRate"
            death_rate_id = "bdSkylineDeathRate"
            sampling_rate_id = "bdSkylineSamplingRate"
            rho_id = "bdSkylineRho"

            bds_elem = ET.Element(
                "distribution",
                {
                    "id": prior_id,
                    "spec": "bdsky.evolution.speciation.BirthDeathSkylineModel",
                    "tree": ref(tree_id),
                    "contemp": "false",
                    "conditionOnRoot": "true",
                },
            )

            _, br_node = _state_node(
                birth_rate_id, params.get("birth_rate", 1.0), role="rate", default_value=1.0
            )
            ET.SubElement(bds_elem, "parameter", {"name": "birthRate", "idref": birth_rate_id})
            state_nodes.append((birth_rate_id, br_node))

            _, dr_node = _state_node(
                death_rate_id, params.get("death_rate", 0.5), role="rate", default_value=0.5
            )
            ET.SubElement(bds_elem, "parameter", {"name": "deathRate", "idref": death_rate_id})
            state_nodes.append((death_rate_id, dr_node))

            # The sampling proportion psi must lie in [0, 1]; found
            # it emitted with a lower bound only.
            _, sr_node = _state_node(
                sampling_rate_id,
                params.get("sampling_rate", 0.1),
                role="probability",
                default_value=0.1,
            )
            ET.SubElement(
                bds_elem, "parameter", {"name": "samplingRate", "idref": sampling_rate_id}
            )
            state_nodes.append((sampling_rate_id, sr_node))

            # Serial sampling: no contemporary (rho) sampling
            rho_raw = params.get("rho", 0.0)
            _, rho_node = _state_node(rho_id, rho_raw, role="probability", default_value=0.0)
            ET.SubElement(
                bds_elem,
                "parameter",
                {
                    "name": "rho",
                    "id": rho_id,
                    "value": format_number(rho_node["value"]),
                    "lower": format_number(rho_node["lower"]),
                    "upper": format_number(rho_node["upper"]),
                    "estimate": "false",
                },
            )

            prior_ids.append(prior_id)
            return bds_elem, state_nodes, prior_ids

        else:
            raise ValueError(f"Tree prior type {tp_type} not yet implemented")

    @staticmethod
    def build_operators(
        config: BEASTConfig,
        partition_mgr: PartitionManager,
        state_node_ids: List[str],
    ) -> List[ET.Element]:
        """Build operator XML elements based on model configuration.

        Args:
            config: BEAST configuration.
            partition_mgr: Partition manager.
            state_node_ids: List of stateNode IDs that need operators.

        Returns:
            List of operator XML elements.
        """
        operators: List[ET.Element] = []
        tree_id = partition_mgr.get_tree_id()
        weights = config.operator_weights

        def get_weight(name: str, default: float) -> float:
            return float(weights.get(name, default))

        # --- Tree operators (shared) ---
        # ScaleOperator for tree
        op = ET.Element(
            "operator",
            {
                "id": "treeScaler",
                "spec": "ScaleOperator",
                "scaleFactor": "0.75",
                "weight": str(get_weight("tree_scaler", 3.0)),
                "tree": ref(tree_id),
            },
        )
        operators.append(op)

        # Uniform operator (fully qualified to avoid clash with distribution.Uniform)
        op = ET.Element(
            "operator",
            {
                "spec": "beast.base.evolution.operator.Uniform",
                "weight": str(get_weight("uniform", 30.0)),
                "tree": ref(tree_id),
            },
        )
        operators.append(op)

        # SubtreeSlide
        op = ET.Element(
            "operator",
            {
                "spec": "SubtreeSlide",
                "weight": str(get_weight("subtree_slide", 15.0)),
                "gaussian": "true",
                "size": "1.7",
                "tree": ref(tree_id),
            },
        )
        operators.append(op)

        # Exchange (narrow)
        op = ET.Element(
            "operator",
            {
                "id": "narrow",
                "spec": "Exchange",
                "isNarrow": "true",
                "weight": str(get_weight("exchange_narrow", 15.0)),
                "tree": ref(tree_id),
            },
        )
        operators.append(op)

        # Exchange (wide)
        op = ET.Element(
            "operator",
            {
                "id": "wide",
                "spec": "Exchange",
                "isNarrow": "false",
                "weight": str(get_weight("exchange_wide", 3.0)),
                "tree": ref(tree_id),
            },
        )
        operators.append(op)

        # WilsonBalding
        op = ET.Element(
            "operator",
            {
                "spec": "WilsonBalding",
                "weight": str(get_weight("wilson_balding", 1.0)),
                "tree": ref(tree_id),
            },
        )
        operators.append(op)

        # UpDownOperator scales the tree down and the clock rate up. Pairing
        # them is essential when the clock rate is estimated; without an "up"
        # item the operator only disturbs the tree (BEAST2 logs a warning).
        up_rate_id: Optional[str] = None
        for pid in partition_mgr.partition_ids:
            p = partition_mgr.get_partition(pid)
            if p.clock_model.linked_to:
                continue
            if p.clock_model.type == ClockModelType.STRICT:
                cr = p.clock_model.clock_rate
                if cr is None or cr.estimate:
                    up_rate_id = make_xml_id("clockRate", pid)
                    break
            elif p.clock_model.type in (ClockModelType.UCLN, ClockModelType.UCE):
                um = p.clock_model.ucld_mean
                if um is None or um.estimate:
                    up_rate_id = make_xml_id("ucld.mean", pid)
                    break

        updown_attrs = {
            "spec": "UpDownOperator",
            "scaleFactor": "0.75",
            "weight": str(get_weight("up_down", 3.0)),
            "down": ref(tree_id),
        }
        if up_rate_id and up_rate_id in state_node_ids:
            updown_attrs["up"] = ref(up_rate_id)
        elif up_rate_id:
            # The only clock rate the scan found is fixed (estimate: false), so
            # it is not a state node. Drop the id rather than reference a
            # parameter that does not exist in the state.
            up_rate_id = None
        op = ET.Element("operator", updown_attrs)
        operators.append(op)

        # --- Partition-specific operators ---
        for pid in partition_mgr.partition_ids:
            p = partition_mgr.get_partition(pid)

            # Substitution model parameter operators
            sm_type = p.site_model.substitution_model.get("type", "hky").lower()
            model_spec = ModelRegistry.get_substitution_model(sm_type)

            if model_spec:
                for param_name, param_spec in model_spec.get("params", {}).items():
                    param_id = make_xml_id(f"{sm_type}.{param_name}", pid)
                    if param_id in state_node_ids:
                        op = ET.Element(
                            "operator",
                            {
                                "id": f"{pid}.{param_name}Scaler",
                                "spec": "ScaleOperator",
                                "scaleFactor": "0.75",
                                "weight": str(get_weight(f"{pid}_{param_name}_scaler", 0.1)),
                                "parameter": ref(param_id),
                            },
                        )
                        operators.append(op)

                # Frequencies operator: emitted whenever an estimated frequency
                # vector exists in the state, regardless of which registry mode
                # created it (an explicit user vector on a "uniform" model such
                # as SYM is still an estimated parameter that must be proposed).
                freq_id = make_xml_id("freqParameter", pid)
                if freq_id in state_node_ids:
                    op = ET.Element(
                        "operator",
                        {
                            "id": f"{pid}.frequenciesDelta",
                            "spec": "DeltaExchangeOperator",
                            "delta": "0.01",
                            "weight": str(get_weight(f"{pid}_frequencies_delta", 0.1)),
                            "parameter": ref(freq_id),
                        },
                    )
                    operators.append(op)

            # Gamma shape operator
            if p.site_model.gamma_categories > 0 and p.site_model.gamma_shape:
                shape_id = make_xml_id("gammaShape", pid)
                if shape_id in state_node_ids:
                    op = ET.Element(
                        "operator",
                        {
                            "id": f"{pid}.gammaShapeScaler",
                            "spec": "ScaleOperator",
                            "scaleFactor": "0.75",
                            "weight": str(get_weight(f"{pid}_gamma_shape_scaler", 0.1)),
                            "parameter": ref(shape_id),
                        },
                    )
                    operators.append(op)

            # Proportion invariant operator: an estimated state node without
            # an operator is never sampled (BEAST2 only prints a warning and
            # the parameter stays at its initial value).
            pi = p.site_model.proportion_invariant
            if pi is not None and pi.estimate:
                pi_id = make_xml_id("proportionInvariant", pid)
                if pi_id in state_node_ids:
                    op = ET.Element(
                        "operator",
                        {
                            "id": f"{pid}.proportionInvariantScaler",
                            "spec": "ScaleOperator",
                            "scaleFactor": "0.75",
                            "weight": str(get_weight(f"{pid}_prop_invariant_scaler", 0.1)),
                            "parameter": ref(pi_id),
                        },
                    )
                    operators.append(op)

            # Clock model operators
            clock = p.clock_model
            if not clock.linked_to:
                if clock.type in (ClockModelType.UCLN, ClockModelType.UCE):
                    # ucld.stdev scaler
                    stdev_id = make_xml_id("ucld.stdev", pid)
                    if stdev_id in state_node_ids:
                        op = ET.Element(
                            "operator",
                            {
                                "id": f"{pid}.stdevScaler",
                                "spec": "ScaleOperator",
                                "scaleFactor": "0.5",
                                "weight": str(get_weight(f"{pid}_stdev_scaler", 1.0)),
                                "parameter": ref(stdev_id),
                            },
                        )
                        operators.append(op)

                    # rateCategories uniform operator
                    cats_id = make_xml_id("rateCategories", pid)
                    if cats_id in state_node_ids:
                        op = ET.Element(
                            "operator",
                            {
                                "id": f"{pid}.categoriesRandomWalk",
                                "spec": "UniformOperator",
                                "weight": str(get_weight(f"{pid}_categories", 1.0)),
                                "parameter": ref(cats_id),
                            },
                        )
                        operators.append(op)

                    # ucld.mean scaler
                    mean_id = make_xml_id("ucld.mean", pid)
                    if mean_id in state_node_ids:
                        op = ET.Element(
                            "operator",
                            {
                                "id": f"{pid}.meanScaler",
                                "spec": "ScaleOperator",
                                "scaleFactor": "0.75",
                                "weight": str(get_weight(f"{pid}_mean_scaler", 1.0)),
                                "parameter": ref(mean_id),
                            },
                        )
                        operators.append(op)

                elif clock.type == ClockModelType.RLC:
                    # indicators bit flip
                    ind_id = make_xml_id("indicators", pid)
                    if ind_id in state_node_ids:
                        op = ET.Element(
                            "operator",
                            {
                                "id": f"{pid}.indicatorFlip",
                                "spec": "BitFlipOperator",
                                "weight": str(get_weight(f"{pid}_indicator_flip", 1.0)),
                                "parameter": ref(ind_id),
                            },
                        )
                        operators.append(op)

                    # rates scaler
                    rates_id = make_xml_id("rates", pid)
                    if rates_id in state_node_ids:
                        op = ET.Element(
                            "operator",
                            {
                                "id": f"{pid}.rateScaler",
                                "spec": "ScaleOperator",
                                "scaleFactor": "0.5",
                                "weight": str(get_weight(f"{pid}_rate_scaler", 1.0)),
                                "parameter": ref(rates_id),
                            },
                        )
                        operators.append(op)

                elif clock.type == ClockModelType.STRICT:
                    # clock rate scaler
                    rate_id = make_xml_id("clockRate", pid)
                    if rate_id in state_node_ids:
                        op = ET.Element(
                            "operator",
                            {
                                "id": f"{pid}.clockRateScaler",
                                "spec": "ScaleOperator",
                                "scaleFactor": "0.75",
                                "weight": str(get_weight(f"{pid}_clock_rate_scaler", 1.0)),
                                "parameter": ref(rate_id),
                            },
                        )
                        operators.append(op)

        # --- Tree prior parameter operators ---
        tp_type = config.tree_prior_type
        if tp_type in (TreePriorType.YULE, TreePriorType.CALIBRATED_YULE):
            if "birthRate" in state_node_ids:
                op = ET.Element(
                    "operator",
                    {
                        "spec": "ScaleOperator",
                        "weight": str(get_weight("birth_rate_scaler", 1.0)),
                        "scaleFactor": "0.25",
                    },
                )
                ET.SubElement(op, "parameter", {"idref": "birthRate"})
                operators.append(op)

        elif tp_type == TreePriorType.BIRTH_DEATH:
            for pid_name, op_weight_name in [
                ("birthDiffRate", "birth_diff_rate_scaler"),
                ("relativeDeathRate", "relative_death_rate_scaler"),
            ]:
                if pid_name in state_node_ids:
                    op = ET.Element(
                        "operator",
                        {
                            "spec": "ScaleOperator",
                            "weight": str(get_weight(op_weight_name, 1.0)),
                            "scaleFactor": "0.5",
                        },
                    )
                    ET.SubElement(op, "parameter", {"idref": pid_name})
                    operators.append(op)

        elif tp_type in (TreePriorType.COALESCENT_CONSTANT, TreePriorType.COALESCENT_EXPONENTIAL):
            if "popSize" in state_node_ids:
                op = ET.Element(
                    "operator",
                    {
                        "id": "popSizeScaler",
                        "spec": "ScaleOperator",
                        "scaleFactor": "0.75",
                        "weight": str(get_weight("pop_size_scaler", 3.0)),
                        "parameter": ref("popSize"),
                    },
                )
                operators.append(op)

            if tp_type == TreePriorType.COALESCENT_EXPONENTIAL and "growthRate" in state_node_ids:
                op = ET.Element(
                    "operator",
                    {
                        "id": "growthRateScaler",
                        "spec": "ScaleOperator",
                        "scaleFactor": "0.75",
                        "weight": str(get_weight("growth_rate_scaler", 1.0)),
                        "parameter": ref("growthRate"),
                    },
                )
                operators.append(op)

        elif tp_type == TreePriorType.BAYESIAN_SKYLINE:
            if "bPopSizes" in state_node_ids:
                op = ET.Element(
                    "operator",
                    {
                        "id": "bPopSizesScaler",
                        "spec": "ScaleOperator",
                        "scaleFactor": "0.75",
                        "weight": str(get_weight("b_pop_sizes_scaler", 1.0)),
                        "parameter": ref("bPopSizes"),
                    },
                )
                operators.append(op)
            if "bGroupSizes" in state_node_ids:
                # IntegerParameter requires the intparameter input
                op = ET.Element(
                    "operator",
                    {
                        "id": "bGroupSizesDelta",
                        "spec": "DeltaExchangeOperator",
                        "delta": "1.0",
                        "integer": "true",
                        "weight": str(get_weight("b_group_sizes_delta", 1.0)),
                    },
                )
                op.set("intparameter", ref("bGroupSizes"))
                operators.append(op)

        elif tp_type == TreePriorType.EBSP:
            if "populationSizes" in state_node_ids:
                op = ET.Element(
                    "operator",
                    {
                        "id": "populationSizesScaler",
                        "spec": "ScaleOperator",
                        "scaleFactor": "0.75",
                        "weight": str(get_weight("pop_sizes_scaler", 1.0)),
                        "parameter": ref("populationSizes"),
                    },
                )
                operators.append(op)
            if "populationIndicators" in state_node_ids:
                op = ET.Element(
                    "operator",
                    {
                        "id": "populationIndicatorsBitFlip",
                        "spec": "BitFlipOperator",
                        "weight": str(get_weight("pop_indicators_bitflip", 1.0)),
                        "parameter": ref("populationIndicators"),
                    },
                )
                operators.append(op)

        elif tp_type == TreePriorType.BD_SKYLINE_SERIAL:
            for pid_name, op_weight_name in [
                ("bdSkylineBirthRate", "bd_birth_rate_scaler"),
                ("bdSkylineDeathRate", "bd_death_rate_scaler"),
                ("bdSkylineSamplingRate", "bd_sampling_rate_scaler"),
            ]:
                if pid_name in state_node_ids:
                    op = ET.Element(
                        "operator",
                        {
                            "id": f"{pid_name}Scaler",
                            "spec": "ScaleOperator",
                            "weight": str(get_weight(op_weight_name, 1.0)),
                            "scaleFactor": "0.5",
                            "parameter": ref(pid_name),
                        },
                    )
                    operators.append(op)

        # --- Calibration hyperprior operators ---
        # build_hyperprior_priors registers each hyperprior parameter as an
        # estimated state node; without an operator it would never be sampled.
        for cal in config.calibrations:
            if not cal.distribution or not cal.hyperpriors:
                continue
            for key in cal.hyperpriors:
                if key not in cal.distribution.parameters:
                    continue
                param_id = sanitize_id(f"{cal.name}.{key}")
                if param_id in state_node_ids:
                    op = ET.Element(
                        "operator",
                        {
                            "id": f"{param_id}.scaler",
                            "spec": "ScaleOperator",
                            "scaleFactor": "0.5",
                            "weight": str(get_weight(f"{cal.name}_{key}_scaler", 1.0)),
                            "parameter": ref(param_id),
                        },
                    )
                    operators.append(op)

        # --- Tip dates operator ---
        # TipDatesRandomWalker used to be emitted only for a
        # calibration that declared `tipsonly: true`. That is *not* how tip
        # dating normally works — in BEAUti the dates are a trait on the tree,
        # not a calibration — so the common configuration produced a TraitSet
        # that no operator ever moved: the tip heights were never proposed, the
        # dates contributed nothing, and the run returned an undated tree while
        # the user read a dated one. The operator is now emitted whenever
        # tip_dates is enabled, over the taxa that actually carry a date.
        if config.tip_dates and config.tip_dates.enabled:
            dated_taxa = sorted(
                taxon for taxon, value in config.tip_dates.dates.items() if float(value) != 0.0
            )
            if dated_taxa:
                # The TaxonSet itself is emitted at the <beast> root by
                # XMLWriter, sharing the taxon-definition bookkeeping with the
                # calibration taxon sets: a TaxonSet may carry either `alignment`
                # or nested <taxon> children, never both, and a nested
                # <taxonset> inside the operator would collide with the
                # taxonset= attribute (both are BEAST2 parser errors).
                op = ET.Element(
                    "operator",
                    {
                        "id": "tipDatesRandomWalker",
                        "spec": "TipDatesRandomWalker",
                        "windowSize": "1",
                        "weight": str(get_weight("tip_dates", 1.0)),
                        "taxonset": ref("tipdates.taxonset"),
                        "tree": ref(tree_id),
                    },
                )
                operators.append(op)

        return operators

    @staticmethod
    def build_state(
        config: BEASTConfig,
        partition_mgr: PartitionManager,
        state_node_ids_with_params: List[Tuple[str, Dict[str, Any]]],
    ) -> ET.Element:
        """Build the state space XML element.

        Args:
            config: BEAST configuration.
            partition_mgr: Partition manager.
            state_node_ids_with_params: List of (id, params) tuples for state nodes.

        Returns:
            State XML element.
        """
        state_elem = ET.Element(
            "state", {"id": "state", "storeEvery": str(config.mcmc.store_every)}
        )

        # Tree state node
        tree_id = partition_mgr.get_tree_id()
        tree_state = ET.SubElement(
            state_elem,
            "tree",
            {"id": tree_id, "name": "stateNode", "estimate": "true"},
        )

        # BEAST2 resolves MRCAPrior/TipDatesRandomWalker taxon names via
        # tree.getTaxaNames(), which requires the tree's taxonset input.
        # Without it, any MRCA calibration crashes at startup with
        # "Cannot find taxon X in data" (BEAUti emits the same shape).
        ET.SubElement(
            tree_state,
            "taxonset",
            {
                "id": f"{tree_id}.taxonset",
                "spec": "TaxonSet",
                "alignment": ref(partition_mgr.get_first_alignment_id()),
            },
        )

        # Tip dates: the TraitSet must be attached to the tree itself (as the
        # tree's "trait" input) for tip ages to take effect.
        if config.tip_dates and config.tip_dates.enabled:
            from .calibration import CalibrationBuilder

            trait_elem = CalibrationBuilder.build_tip_dates_traitset(
                config.tip_dates.dates,
                config.tip_dates.trait_name,
                config.tip_dates.units,
                partition_mgr.get_first_alignment_id(),
            )
            tree_state.append(trait_elem)

        # Add all parameter state nodes
        for node_id, params in state_node_ids_with_params:
            estimate = params.get("estimate", True)
            if isinstance(estimate, str):
                estimate = _as_bool(estimate, node_id)
            param_attrs = {
                "id": node_id,
                "name": "stateNode",
                "value": format_number(params.get("value", 1.0)),
                "dimension": str(params.get("dimension", 1)),
                "estimate": "true" if estimate else "false",
            }
            # Support custom spec (e.g., IntegerParameter)
            if params.get("spec"):
                param_attrs["spec"] = params["spec"]
            param_elem = ET.SubElement(
                state_elem,
                "parameter",
                param_attrs,
            )
            if params.get("lower") is not None:
                param_elem.set("lower", format_number(params["lower"]))
            if params.get("upper") is not None:
                param_elem.set("upper", format_number(params["upper"]))

        return state_elem

    @staticmethod
    def build_init(
        config: BEASTConfig,
        partition_mgr: PartitionManager,
        tree_id: str,
    ) -> ET.Element:
        """Build the init (initial tree) XML element.

        Args:
            config: BEAST configuration.
            partition_mgr: Partition manager.
            tree_id: Tree ID.

        Returns:
            Init XML element.
        """
        init_type = config.initialization.tree_type
        first_aln_id = partition_mgr.get_first_alignment_id()

        if init_type == InitTreeType.RANDOM:
            # Check if CalibratedYule method is used
            if config.calibration_method == CalibrationMethod.CALIBRATED_YULE:
                # Use CalibratedYuleInitialTree for CalibratedYule method
                init_elem = ET.Element(
                    "init",
                    {
                        "spec": "beast.base.evolution.speciation.CalibratedYuleInitialTree",
                        "id": "calibratedYuleTree",
                        "initial": ref(tree_id),
                    },
                )
                # TaxonSet referencing the alignment
                ET.SubElement(
                    init_elem, "taxonset", {"spec": "TaxonSet", "alignment": ref(first_aln_id)}
                )
                # Reference calibration points
                for cal in config.calibrations:
                    ET.SubElement(
                        init_elem, "calibrations", {"idref": f"{sanitize_id(cal.name)}.cal"}
                    )
                return init_elem

            init_elem = ET.Element(
                "init",
                {
                    "spec": "beast.base.evolution.tree.coalescent.RandomTree",
                    "id": "randomTree",
                    "initial": ref(tree_id),
                    "taxa": ref(first_aln_id),
                },
            )

            pop_model = ET.SubElement(
                init_elem,
                "populationModel",
                {"spec": "ConstantPopulation"},
            )
            ET.SubElement(
                pop_model,
                "parameter",
                {"name": "popSize", "value": "1"},
            )

            # Constrain the initial tree by all MRCAPrior calibrations:
            # monophyletic clades, and root calibrations (taxa=None) in
            # particular. RandomTree resamples the initial tree until every
            # constraint holds; without them a random starting tree violates
            # hard-bounded calibrations and the chain starts at an infinite
            # posterior, which BEAST2 refuses to run.
            for cal in config.calibrations:
                if cal.taxa is None or cal.monophyletic:
                    ET.SubElement(init_elem, "constraint", {"idref": sanitize_id(cal.name)})

            # Scale the initial tree so deep calibrations fit inside it: a
            # clade calibrated at, say, 6 Ma cannot sit below a coalescent-
            # scale root of ~2 Ma (invalid state, -inf likelihood). BEAST2's
            # RandomTree rescales the sampled tree to this root height when
            # the constraints allow it.
            if config.calibrations and config.calibration_method == CalibrationMethod.MRCA_PRIOR:
                max_age = max(
                    (_estimate_calibration_age(cal.distribution) for cal in config.calibrations),
                    default=0.0,
                )
                if max_age > 0:
                    init_elem.set("rootHeight", str(round(max_age * 1.5, 3)))

            return init_elem

        elif init_type == InitTreeType.UPGMA:
            init_elem = ET.Element(
                "init",
                {
                    "spec": "beast.base.evolution.tree.ClusterTree",
                    "id": "upgmaTree",
                    "initial": ref(tree_id),
                    "clusterType": "upgma",
                    "taxa": ref(first_aln_id),
                },
            )
            return init_elem

        elif init_type == InitTreeType.NEWICK:
            newick_file = config.initialization.newick_file
            if newick_file:
                with open(newick_file, "r", encoding="utf-8") as f:
                    newick_str = f.read().strip()
                init_elem = ET.Element(
                    "init",
                    {
                        "spec": "beast.base.evolution.tree.TreeParser",
                        "id": "newickTree",
                        "initial": ref(tree_id),
                        "taxa": ref(first_aln_id),
                        "newick": newick_str,
                    },
                )
                return init_elem
            else:
                raise ValueError("newick_file must be specified for newick initialization")

        else:
            raise ValueError(f"Unknown init tree type: {init_type}")

"""Calibration assembly module for Beast2Py.

Builds XML elements for MRCAPrior, TaxonSet, calibration distributions,
and CalibratedYuleModel calibration points.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from typing import List, Optional

from .models import CalibrationPoint, DistributionConfig
from .registry import ModelRegistry
from .utils import ref, sanitize_id


class CalibrationBuilder:
    """Build calibration-related XML elements."""

    # Scale-typed distribution parameters are unbounded above in density as
    # they approach 0 (e.g. Normal sigma at t == mean), so a hyperprior state
    # node with lower=0 combined with a scale operator drives the posterior
    # to +inf. A small positive floor keeps the model proper.
    SCALE_PARAM_KEYS = frozenset({"sigma", "S", "scale"})
    SCALE_PARAM_FLOOR = 1e-3

    @staticmethod
    def build_taxonset(
        taxa: List[str],
        cal_name: str,
        defined_taxa: Optional[set] = None,
    ) -> ET.Element:
        """Build a TaxonSet XML element.

        Uses spec='Taxon' with id for first occurrence of each taxon,
        and idref for subsequent references to avoid duplicate IDs.

        Args:
            taxa: List of taxon names.
            cal_name: Calibration name (used for TaxonSet ID).
            defined_taxa: Set of already-defined taxon names (updated in place).

        Returns:
            TaxonSet XML element.
        """
        if defined_taxa is None:
            defined_taxa = set()
        ts_id = sanitize_id(f"{cal_name}.taxonset")
        ts_elem = ET.Element("taxonset", {"spec": "TaxonSet", "id": ts_id})
        for taxon in taxa:
            if taxon in defined_taxa:
                ET.SubElement(ts_elem, "taxon", {"spec": "Taxon", "idref": taxon})
            else:
                ET.SubElement(ts_elem, "taxon", {"spec": "Taxon", "id": taxon})
                defined_taxa.add(taxon)
        return ts_elem

    @staticmethod
    def build_distribution(
        dist: DistributionConfig,
        dist_id: Optional[str] = None,
        hyperpriors: Optional[dict] = None,
        state_nodes: Optional[List[tuple]] = None,
        owner_name: str = "",
    ) -> ET.Element:
        """Build a calibration distribution XML element.

        Args:
            dist: Distribution configuration.
            dist_id: Optional ID for the distribution element.
            hyperpriors: Optional dict mapping parameter name to
                HyperpriorConfig; listed parameters become estimated state
                nodes referenced from the distribution.
            state_nodes: Optional output list collecting
                (stateNode ID, params) tuples for hyperprior-estimated
                parameters.
            owner_name: Name prefix (e.g. calibration name) for state IDs.

        Returns:
            Distribution XML element.
        """
        dist_spec = ModelRegistry.get_distribution_spec(dist.type)
        if dist_spec is None:
            raise ValueError(f"Unknown distribution type: {dist.type}")

        spec = dist_spec["spec"]
        attrs = {"spec": spec}
        if dist_id:
            attrs["id"] = dist_id

        # Add offset if present
        if dist.offset and dist.offset != 0.0:
            attrs["offset"] = str(dist.offset)

        # BEAST2's Gamma selects its parameterisation with a `mode` *attribute*
        # (ShapeScale | ShapeRate | ShapeMean | OneParameter; Gamma.java:22-31),
        # not with a child parameter.  Beast2Py's summary statistics follow the
        # requested mode, so the XML has to state it as well - otherwise a
        # `mode: ShapeRate` prior is sampled as a scale prior while the reported
        # mean and interval use the rate convention.
        if dist.type.lower() in ("gamma", "inverse_gamma", "inversegamma"):
            mode = (dist.parameters or {}).get("mode")
            if mode is not None:
                allowed = (
                    {"shapescale", "shaperate", "shapemean", "oneparameter"}
                    if dist.type.lower() == "gamma"
                    else {"shapescale", "shaperate"}
                )
                if str(mode).lower().replace("_", "") not in allowed:
                    raise ValueError(
                        f"Gamma mode {mode!r} is not one of "
                        + ", ".join(sorted(allowed))
                        + " (BEAST2 Gamma.java:22-31)."
                    )
                attrs["mode"] = str(mode)

        dist_elem = ET.Element("distr", attrs)

        # Add boolean attributes (e.g., meanInRealSpace for LogNormal)
        bool_attrs = dist_spec.get("bool_attrs", {})
        for config_key, xml_attr in bool_attrs.items():
            if config_key in dist.parameters:
                val = dist.parameters[config_key]
                dist_elem.set(xml_attr, str(bool(val)).lower())

        # Add parameters
        use_param_elements = dist_spec.get("use_param_elements", False)
        param_mapping = dist_spec.get("params", {})

        for config_key, xml_attr in param_mapping.items():
            if config_key not in dist.parameters:
                continue
            val = dist.parameters[config_key]

            hp = (hyperpriors or {}).get(config_key)
            if hp is not None:
                # Hyperprior on this parameter: make it an estimated state
                # node and reference it from the distribution.
                param_id = (
                    sanitize_id(f"{owner_name}.{config_key}")
                    if owner_name
                    else sanitize_id(config_key)
                )
                if use_param_elements:
                    ET.SubElement(
                        dist_elem,
                        "parameter",
                        {"name": xml_attr, "idref": param_id},
                    )
                else:
                    dist_elem.set(xml_attr, ref(param_id))
                if state_nodes is not None:
                    state_nodes.append(
                        (
                            param_id,
                            {
                                "value": float(val) if not isinstance(val, bool) else str(val),
                                "lower": (
                                    CalibrationBuilder.SCALE_PARAM_FLOOR
                                    if config_key in CalibrationBuilder.SCALE_PARAM_KEYS
                                    else 0.0
                                ),
                                "estimate": True,
                            },
                        )
                    )
                continue

            if use_param_elements:
                ET.SubElement(
                    dist_elem,
                    "parameter",
                    {"name": xml_attr, "value": str(val)},
                )
            else:
                dist_elem.set(xml_attr, str(val))

        return dist_elem

    @staticmethod
    def build_hyperprior_priors(cal: CalibrationPoint) -> List[tuple]:
        """Build Prior distribution elements for a calibration's hyperpriors.

        Each hyperprior turns one calibration distribution parameter into an
        estimated state node (referenced from the calibration distribution by
        build_distribution) and adds a Prior on it to the posterior's prior
        compound distribution.

        Args:
            cal: Calibration point with distribution and hyperpriors.

        Returns:
            List of (Prior XML element, stateNode ID, stateNode params)
            tuples. Empty if the calibration has no hyperpriors.
        """
        if not cal.distribution or not cal.hyperpriors:
            return []

        results: List[tuple] = []
        for config_key, hp in cal.hyperpriors.items():
            if config_key not in cal.distribution.parameters:
                continue
            hp_spec = ModelRegistry.get_distribution_spec(hp.type)
            if hp_spec is None:
                raise ValueError(f"Unknown hyperprior distribution type: {hp.type}")

            param_id = sanitize_id(f"{cal.name}.{config_key}")
            val = cal.distribution.parameters[config_key]

            prior_elem = ET.Element(
                "distribution",
                {
                    "id": f"{param_id}.prior",
                    "spec": "beast.base.inference.distribution.Prior",
                    "x": ref(param_id),
                },
            )
            distr_attrs = {"spec": hp_spec["spec"]}
            if hp.offset:
                distr_attrs["offset"] = str(hp.offset)
            distr_elem = ET.SubElement(prior_elem, "distr", distr_attrs)

            use_param_elements = hp_spec.get("use_param_elements", False)
            for hp_key, xml_attr in hp_spec.get("params", {}).items():
                if hp_key in hp.parameters:
                    if use_param_elements:
                        ET.SubElement(
                            distr_elem,
                            "parameter",
                            {"name": xml_attr, "value": str(hp.parameters[hp_key])},
                        )
                    else:
                        distr_elem.set(xml_attr, str(hp.parameters[hp_key]))

            results.append(
                (
                    prior_elem,
                    param_id,
                    {
                        "value": float(val) if not isinstance(val, bool) else str(val),
                        "lower": (
                            CalibrationBuilder.SCALE_PARAM_FLOOR
                            if config_key in CalibrationBuilder.SCALE_PARAM_KEYS
                            else 0.0
                        ),
                        "estimate": True,
                    },
                )
            )

        return results

    @staticmethod
    def build_mrcaprior(
        cal: CalibrationPoint,
        tree_id: str,
        alignment_id: str = None,
    ) -> ET.Element:
        """Build an MRCAPrior XML element.

        Args:
            cal: Calibration point.
            tree_id: Tree ID.

        Returns:
            MRCAPrior XML element.
        """
        attrs = {
            "id": sanitize_id(cal.name),
            "spec": "beast.base.evolution.tree.MRCAPrior",
            "tree": ref(tree_id),
            "monophyletic": str(cal.monophyletic).lower(),
        }

        if cal.use_originate and cal.taxa is not None:
            attrs["useOriginate"] = "true"

        if cal.tipsonly:
            attrs["tipsonly"] = "true"

        mrca_elem = ET.Element("distribution", attrs)

        # Add TaxonSet reference
        if cal.taxa is not None:
            taxonset_id = sanitize_id(f"{cal.name}.taxonset")
            ET.SubElement(mrca_elem, "taxonset", {"idref": taxonset_id})
        else:
            # All taxa (root calibration) - define inline with alignment reference
            ts_attrs = {"spec": "TaxonSet", "id": sanitize_id(f"{cal.name}.taxonset")}
            if alignment_id:
                ts_attrs["alignment"] = ref(alignment_id)
            taxonset_elem = ET.Element("taxonset", ts_attrs)
            mrca_elem.append(taxonset_elem)

        # Add distribution
        if cal.distribution:
            distr_elem = CalibrationBuilder.build_distribution(
                cal.distribution,
                dist_id=f"{sanitize_id(cal.name)}.distr",
                hyperpriors=cal.hyperpriors,
                owner_name=sanitize_id(cal.name),
            )
            mrca_elem.append(distr_elem)

        return mrca_elem

    @staticmethod
    def build_calibrated_yule_point(
        cal: CalibrationPoint,
        all_taxa_ref: str = None,
    ) -> ET.Element:
        """Build a CalibratedYuleModel CalibrationPoint XML element.

        Args:
            cal: Calibration point.

        Returns:
            CalibrationPoint XML element.
        """
        cal_elem = ET.Element(
            "calibrations",
            {
                "spec": "CalibrationPoint",
                "id": f"{sanitize_id(cal.name)}.cal",
            },
        )

        # TaxonSet
        if cal.taxa is not None:
            ts_elem = ET.Element("taxonset", {"spec": "TaxonSet"})
            for taxon in cal.taxa:
                ET.SubElement(ts_elem, "taxon", {"idref": taxon})
            cal_elem.append(ts_elem)
        elif all_taxa_ref:
            # Root calibration - reference all taxa via alignment
            ET.SubElement(
                cal_elem, "taxonset", {"spec": "TaxonSet", "alignment": ref(all_taxa_ref)}
            )

        # Distribution
        if cal.distribution:
            distr_elem = CalibrationBuilder.build_distribution(
                cal.distribution,
                dist_id=f"{sanitize_id(cal.name)}.distr",
                hyperpriors=cal.hyperpriors,
                owner_name=sanitize_id(cal.name),
            )
            cal_elem.append(distr_elem)

        return cal_elem

    @staticmethod
    def build_mrca_logger(
        cal: CalibrationPoint,
        tree_id: str,
    ) -> ET.Element:
        """Build an MRCA time logger element.

        Args:
            cal: Calibration point.
            tree_id: Tree ID.

        Returns:
            Logger XML element for MRCA time.
        """
        taxonset_id = sanitize_id(f"{cal.name}.taxonset")
        log_elem = ET.Element(
            "log",
            {
                "spec": "beast.base.evolution.tree.MRCAPrior",
                "tree": ref(tree_id),
                "id": f"{sanitize_id(cal.name)}.age",
                "taxonset": ref(taxonset_id),
            },
        )
        return log_elem

    @staticmethod
    def build_tip_dates_traitset(
        taxon_dates: dict,
        trait_name: str,
        units: str,
        alignment_id: str,
    ) -> ET.Element:
        """Build a TraitSet XML element for tip dates.

        Args:
            taxon_dates: Dict mapping taxon name to date value.
            trait_name: Trait name (date-forward or date-backward).
            units: Time units (year, month, day).
            alignment_id: Alignment ID for taxa reference.

        Returns:
            TraitSet XML element.
        """
        # Build trait value string
        parts = [f"{taxon} = {date}" for taxon, date in taxon_dates.items()]
        value = ",\n\t\t\t".join(parts)

        trait_elem = ET.Element(
            "trait",
            {
                "spec": "beast.base.evolution.tree.TraitSet",
                "traitname": trait_name,
                "units": units,
                "value": value,
            },
        )

        ET.SubElement(
            trait_elem,
            "taxa",
            {
                "spec": "TaxonSet",
                "alignment": ref(alignment_id),
            },
        )

        return trait_elem

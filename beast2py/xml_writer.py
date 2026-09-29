"""XML generation module for Beast2Py.

Assembles complete BEAST2 XML using xml.etree.ElementTree, combining
data, partitions, models, calibrations, operators, loggers, and MCMC
configuration into a single valid XML document.
"""

from __future__ import annotations

import math as _math
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional, Tuple

from .calibration import CalibrationBuilder
from .models import (
    BEASTConfig,
    CalibrationMethod,
    MCMCType,
)
from .model import ModelBuilder
from .partition import PartitionManager
from .utils import build_namespace, ref, serialize_xml


class XMLWriter:
    """Generate complete BEAST2 XML using xml.etree.ElementTree."""

    def __init__(self, config: BEASTConfig):
        """Initialize the XML writer.

        Args:
            config: BEAST configuration object.
        """
        self.config = config
        self.partition_mgr = PartitionManager(config.partitions)
        self.tree_id = self.partition_mgr.get_tree_id()
        self._state_node_params: List[Tuple[str, Dict[str, Any]]] = []
        self._clock_model_ids: Dict[str, str] = {}
        # State nodes that already carry a density in the prior compound.
        self._prior_covered: set = set()
        # State nodes that received an automatically supplied default prior.
        self._default_prior_ids: List[str] = []

    @property
    def default_priors(self) -> List[str]:
        """Ids of parameters whose prior was supplied by the generator.

        Returns:
            List of state node ids, empty until ``generate_xml`` has run.
        """
        return list(self._default_prior_ids)

    @property
    def add_on_dependencies(self) -> List[str]:
        """Add-on packages this XML needs at run time.

        ``bdsky.*`` and ``nestedsampling.*`` spec classes are not part of the
        BEAST2 core, and their exact inputs could not be verified against the
        core sources. Reporting them lets the CLI tell the user that a missing
        add-on is a hard start-up failure, not a silent downgrade.

        Returns:
            Sorted list of package names inferred from the emitted specs.
        """
        root = getattr(self, "_root", None)
        if root is None:
            return []
        packages = set()
        for elem in root.iter():
            spec = elem.get("spec", "")
            if spec.startswith("bdsky."):
                packages.add("bdsky")
            elif spec.startswith("nestedsampling."):
                packages.add("nested-sampling")
        return sorted(packages)

    def generate_xml(self) -> str:
        """Generate the complete BEAST2 XML string.

        Returns:
            A formatted XML string.
        """
        # Create root element
        root = ET.Element(
            "beast",
            {
                "version": "2.0",
                "namespace": build_namespace(),
            },
        )

        # 1. Data elements (alignments)
        # Every declared alignment is emitted exactly once, in dependency
        # order (a filter's source must precede the filter). The old loop
        # compared alignment ids against *partition* ids and re-emitted data
        # per partition, so a shared alignment was defined twice and a filter
        # could reference itself.
        for data_elem in self._build_data_elements():
            root.append(data_elem)

        # 2. TaxonSets for calibrations
        defined_taxa: set = set()
        for cal in self.config.calibrations:
            if cal.taxa is not None:
                root.append(CalibrationBuilder.build_taxonset(cal.taxa, cal.name, defined_taxa))

        # Tip-dating needs the set of dated tips as a TaxonSet for the
        # TipDatesRandomWalker operator. It is emitted here, sharing the taxon
        # definition bookkeeping so no taxon is declared twice.
        if self.config.tip_dates and self.config.tip_dates.enabled:
            dated_taxa = sorted(
                taxon for taxon, value in self.config.tip_dates.dates.items() if float(value) != 0.0
            )
            if dated_taxa:
                root.append(CalibrationBuilder.build_taxonset(dated_taxa, "tipdates", defined_taxa))

        # 3. Posterior distribution
        root.append(self._build_posterior())

        # 4. MCMC run
        root.append(self._build_mcmc())

        # 5. Embed fingerprint comment
        if self.config.metadata.get("embed_fingerprint", True):
            from .reproducibility import FingerprintGenerator

            root.insert(0, ET.Comment(FingerprintGenerator.generate_comment(self.config)))

        # Add XML declaration (normalized UTF-8 form, see utils.serialize_xml)
        xml_str = serialize_xml(root, indent=True)
        self._root = root

        # Fail closed on identifier collisions before the file is handed to a
        # user: sanitize_id is many-to-one (two CJK taxon or calibration names
        # can collapse to the same ASCII id), and the standalone validator is
        # the only thing that used to notice — a caller going straight to
        # XMLWriter got a silently ambiguous model.
        seen_ids: Dict[str, int] = {}
        for elem in root.iter():
            elem_id = elem.get("id")
            if elem_id:
                seen_ids[elem_id] = seen_ids.get(elem_id, 0) + 1
        dupes = sorted(k for k, v in seen_ids.items() if v > 1)
        if dupes:
            raise ValueError(
                f"Duplicate XML ids generated: {', '.join(dupes)}. Names must differ after "
                f" sanitisation; rename the conflicting taxa or calibration points."
            )

        return xml_str

    def _build_data_elements(self) -> List[ET.Element]:
        """Build every ``<data>`` element exactly once, in dependency order.

        Returns:
            Ordered list of data elements; a filtered alignment always follows
            the alignment it filters.

        Raises:
            ValueError: If a filtered alignment points at itself.
        """
        ordered: List[Any] = []
        for part in self.config.partitions:
            ordered.append(part.alignment)
        for aln in self.config.all_alignments:
            if aln.id not in {a.id for a in ordered}:
                ordered.append(aln)

        unique: List[Any] = []
        seen: set = set()
        for aln in ordered:
            if aln.id in seen:
                continue
            seen.add(aln.id)
            unique.append(aln)

        emitted: List[ET.Element] = []
        # Plain alignments first, so a FilteredAlignment's source is always
        # defined before it is referenced.
        for aln in [a for a in unique if not a.filter] + [a for a in unique if a.filter]:
            if aln.filter:
                source_id = aln.data_id
                if not source_id:
                    raise ValueError(
                        f"Alignment '{aln.id}' has a filter but no data_id: it would "
                        f"reference itself and lose all sequence data"
                    )
                if source_id == aln.id:
                    raise ValueError(
                        f"Alignment '{aln.id}': data_id must name a different alignment"
                    )
                if source_id not in seen:
                    raise ValueError(
                        f"Alignment '{aln.id}': data_id references unknown alignment "
                        f"'{source_id}'"
                    )
                data_elem = ET.Element(
                    "data",
                    {
                        "data": ref(source_id),
                        "dataType": aln.data_type.value,
                        "filter": aln.filter,
                        "id": aln.id,
                        "spec": "beast.base.evolution.alignment.FilteredAlignment",
                    },
                )
            else:
                data_elem = ET.Element(
                    "data",
                    {
                        "id": aln.id,
                        "dataType": aln.data_type.value,
                    },
                )
                for seq in aln.sequences:
                    seq_elem = ET.SubElement(data_elem, "sequence", {"taxon": seq.taxon})
                    seq_elem.text = seq.sequence
            emitted.append(data_elem)

        return emitted

    def _build_posterior(self) -> ET.Element:
        """Build the posterior CompoundDistribution element.

        Contains prior (tree prior + calibrations + parameter priirs) and
        likelihood (TreeLikelihood for each partition).
        """
        posterior = ET.Element(
            "distribution",
            {
                "id": "posterior",
                "spec": "CompoundDistribution",
            },
        )

        # Build prior
        prior = ET.SubElement(
            posterior,
            "distribution",
            {
                "id": "prior",
                "spec": "CompoundDistribution",
            },
        )

        # Tree prior
        tree_prior_elem, tp_state_nodes, _tp_prior_ids = ModelBuilder.build_tree_prior(
            self.config, self.tree_id, self.partition_mgr
        )
        prior.append(tree_prior_elem)
        for sn_id, sn_params in tp_state_nodes:
            self._state_node_params.append((sn_id, sn_params))

        # MRCA calibration priors
        if self.config.calibration_method == CalibrationMethod.MRCA_PRIOR:
            first_aln_id = self.partition_mgr.get_first_alignment_id()
            for cal in self.config.calibrations:
                mrca_elem = CalibrationBuilder.build_mrcaprior(
                    cal, self.tree_id, alignment_id=first_aln_id
                )
                prior.append(mrca_elem)
                # Hyperpriors: estimated calibration distribution parameters
                for prior_elem, sn_id, sn_params in CalibrationBuilder.build_hyperprior_priors(cal):
                    prior.append(prior_elem)
                    self._state_node_params.append((sn_id, sn_params))
                    # The hyperprior density *is* the prior for that node.
                    self._prior_covered.add(sn_id)
        elif self.config.calibration_method == CalibrationMethod.CALIBRATED_YULE:
            # Calibration points are children of the CalibratedYuleModel; only
            # their hyperprior Priors go into the prior compound here.
            for cal in self.config.calibrations:
                for prior_elem, sn_id, sn_params in CalibrationBuilder.build_hyperprior_priors(cal):
                    prior.append(prior_elem)
                    self._state_node_params.append((sn_id, sn_params))
                    self._prior_covered.add(sn_id)

        # Parameter priors are resolved after the model has been assembled, so
        # that a user-facing name can be matched against the real XML id
        # (see _add_user_priors).
        user_priors = list(self.config.parameter_priors)

        # Build likelihood
        likelihood = ET.SubElement(
            posterior,
            "distribution",
            {
                "id": "likelihood",
                "spec": "CompoundDistribution",
            },
        )

        for p in self.config.partitions:
            tl_elem = self._build_tree_likelihood(p)
            likelihood.append(tl_elem)

        self._add_user_priors(prior, user_priors)
        self._add_default_priors(prior)

        return posterior

    # Aliases from the names a user naturally writes in ``parameter_priors`` to
    # the XML ids the generator actually creates. renaming the
    # birth-death state nodes to birthDiffRate / relativeDeathRate made a prior
    # written against ``birth_rate`` match nothing, silently leaving that
    # parameter without the prior the user believed they had set.
    PRIOR_ALIASES = {
        "birth_rate": ("birthRate", "birthDiffRate"),
        "death_rate": ("relativeDeathRate",),
        "sample_probability": ("sampleProbability",),
        "sampling_rate": ("sampleProbability", "bdSkylineSamplingRate"),
        "pop_size": ("popSize",),
        "pop_sizes": ("bPopSizes", "populationSizes"),
        "growth_rate": ("growthRate",),
        "clock_rate": ("clockRate",),
        "ucld_mean": ("ucld.mean",),
        "ucld_stdev": ("ucld.stdev",),
        "gamma_shape": ("gammaShape",),
        "proportion_invariant": ("proportionInvariant",),
        "kappa": ("kappa",),
    }

    def _resolve_prior_target(self, name: str, state_ids: List[str]) -> Optional[str]:
        """Map a user-supplied parameter name onto a real state node id.

        Args:
            name: The name written in the YAML prior block.
            state_ids: Ids of the state nodes actually generated.

        Returns:
            The matching state node id, or None when nothing matches.
        """
        candidates = [name, *self.PRIOR_ALIASES.get(name, ())]
        for pid in self.partition_mgr.partition_ids:
            candidates += [f"{pid}.{c}" for c in list(candidates)]
        for cand in candidates:
            if cand in state_ids:
                return cand
        # Suffix match, e.g. "kappa" -> "gene1.hky.kappa"
        for sid in state_ids:
            if sid.endswith("." + name):
                return sid
            for alias in self.PRIOR_ALIASES.get(name, ()):
                if sid.endswith("." + alias) or sid == alias:
                    return sid
        return None

    def _add_user_priors(self, prior: ET.Element, user_priors: List[Dict[str, Any]]) -> None:
        """Append user-declared parameter priors, resolving their targets.

        Args:
            prior: The prior CompoundDistribution element.
            user_priors: Raw ``parameter_priors`` entries.

        Raises:
            ValueError: If a prior names a parameter that is not estimated.
        """
        state_ids = [sn[0] for sn in self._state_node_params]
        for pp in user_priors:
            param_name = pp.get("parameter", "")
            target = self._resolve_prior_target(param_name, state_ids)
            if target is None:
                raise ValueError(
                    f"parameter_priors: '{param_name}' does not match any estimated parameter "
                    f"of this analysis. Estimated state nodes are: "
                    f"{', '.join(sorted(state_ids)) if state_ids else '(none)'}. A prior whose "
                    f"target does not exist is silently never evaluated."
                )
            resolved = dict(pp)
            resolved["parameter"] = target
            elem = self._build_parameter_prior(resolved)
            if elem is not None:
                prior.append(elem)
                self._prior_covered.add(target)

    def _default_prior_element(self, node_id: str, node: Dict[str, Any]) -> Optional[ET.Element]:
        """Build a weakly informative default prior for one estimated parameter.

        Continuous parameters on an unbounded or half-bounded support are not
        proper without a density: a flat prior on (0, inf) puts infinite mass
        everywhere, so the "posterior" is not a probability distribution and the
        analysis cannot be reproduced from the XML because no prior was written
        down. The defaults below follow the conventions BEAUti
        ships: a diffuse lognormal on positive rates, a truncated normal
        (mean 0.4, sigma 0.3334) on the UCLN stdev, a uniform on the bounded
        probability/frequency parameters, and a broad normal on signed rates.

        Args:
            node_id: The state node id.
            node: Its parameters, including ``role``, ``lower`` and ``upper``.

        Returns:
            A Prior distribution element, or None when no prior is required
            (bounded discrete parameters are uniform by construction).
        """
        from .models import DistributionConfig

        role = node.get("role", "rate")
        if role in ("integer", "boolean", "tree"):
            return None

        value = node.get("value", 1.0)
        prior_id = f"{node_id}.prior"

        def distr(config: DistributionConfig, dist_id: str) -> ET.Element:
            return CalibrationBuilder.build_distribution(config, dist_id=dist_id)

        if role in ("probability", "frequency"):
            lower = node.get("lower") if node.get("lower") is not None else 0.0
            upper = node.get("upper") if node.get("upper") is not None else 1.0
            config = DistributionConfig(
                type="uniform", parameters={"lower": float(lower), "upper": float(upper)}
            )
        elif role == "signed":
            config = DistributionConfig(type="normal", parameters={"mean": 0.0, "sigma": 10.0})
        else:  # "rate" and anything else that must be positive
            try:
                start = float(value)
            except (TypeError, ValueError):
                start = 1.0
            log_mean = _math.log(start) if start > 0 else 0.0
            if node_id.endswith("ucld.stdev"):
                config = DistributionConfig(
                    type="normal", parameters={"mean": 0.4, "sigma": 0.3334}
                )
            else:
                # Diffuse lognormal centred on the initial value: proper on
                # (0, inf) and a ~500x spread either way, so it informs the
                # support without informing the answer.
                config = DistributionConfig(
                    type="lognormal", parameters={"M": round(log_mean, 6), "S": 2.0}
                )

        prior_elem = ET.Element(
            "distribution",
            {
                "id": prior_id,
                "spec": "beast.base.inference.distribution.Prior",
                "x": ref(node_id),
            },
        )
        prior_elem.append(distr(config, f"{node_id}.distr"))
        return prior_elem

    def _add_default_priors(self, prior: ET.Element) -> None:
        """Attach default priors to every estimated parameter that lacks one.

        Args:
            prior: The prior CompoundDistribution element.
        """
        self._default_prior_ids = []
        for node_id, node in self._state_node_params:
            estimate = node.get("estimate", True)
            if isinstance(estimate, str):
                estimate = estimate.lower() not in ("false", "no", "0")
            if not estimate:
                continue
            if node_id in self._prior_covered:
                continue
            elem = self._default_prior_element(node_id, node)
            if elem is not None:
                prior.append(elem)
                self._default_prior_ids.append(node_id)
                self._prior_covered.add(node_id)

    def _build_tree_likelihood(self, partition) -> ET.Element:
        """Build a TreeLikelihood element for a partition."""
        pid = partition.id
        data_id = partition.alignment.id
        data_type = partition.alignment.data_type

        # Site model: built for unlinked partitions, referenced when linked
        site_model_ref = self.partition_mgr.get_site_ref(pid)
        if partition.site_model.linked_to:
            site_model_elem = None
        else:
            site_model_elem, sm_state_nodes = ModelBuilder.build_site_model(
                pid, partition.site_model, data_id=data_id, data_type=data_type
            )
            for sn_id, sn_params in sm_state_nodes:
                self._state_node_params.append((sn_id, sn_params))

        # Build or reference clock model
        clock = partition.clock_model
        clock_elem = None
        if clock.linked_to:
            clock_ref = self.partition_mgr.get_clock_ref(pid)
        else:
            clock_elem, clock_state_nodes = ModelBuilder.build_clock_model(
                pid, clock, self.tree_id, partition.alignment.n_taxa
            )
            clock_ref = f"{pid}.clockModel"
            self._clock_model_ids[pid] = clock_ref
            for sn_id, sn_params in clock_state_nodes:
                self._state_node_params.append((sn_id, sn_params))

        # Build TreeLikelihood
        tl_id = f"{pid}.treeLikelihood"
        tl_elem = ET.Element(
            "distribution",
            {
                "id": tl_id,
                "spec": "TreeLikelihood",
                "data": ref(data_id),
                "tree": ref(self.tree_id),
            },
        )

        if site_model_elem is not None:
            tl_elem.append(site_model_elem)
        else:
            tl_elem.append(ET.Element("siteModel", {"idref": site_model_ref}))

        if clock.linked_to:
            ET.SubElement(
                tl_elem, "branchRateModel", {"idref": self.partition_mgr.get_clock_ref(pid)}
            )
        elif clock_elem is not None:
            tl_elem.append(clock_elem)

        return tl_elem

    def _build_parameter_prior(self, pp: Dict[str, Any]) -> Optional[ET.Element]:
        """Build a parameter prior distribution element.

        Args:
            pp: Parameter prior config dict with 'parameter' and 'distribution' keys.

        Returns:
            Prior distribution XML element, or None.
        """
        param_name = pp.get("parameter", "")
        dist_config = pp.get("distribution", {})
        if not param_name or not dist_config:
            return None

        prior_id = f"{param_name}.prior"
        prior_elem = ET.Element(
            "distribution",
            {
                "id": prior_id,
                "spec": "beast.base.inference.distribution.Prior",
                "x": ref(param_name),
            },
        )

        # Build the distribution
        from .models import DistributionConfig

        dist = DistributionConfig(
            type=dist_config.get("type", "uniform"),
            parameters=dist_config.get("parameters", {}),
            offset=float(dist_config.get("offset", 0.0)),
        )
        distr_elem = CalibrationBuilder.build_distribution(dist, dist_id=f"{param_name}.distr")
        prior_elem.append(distr_elem)

        return prior_elem

    def _build_mcmc(self) -> ET.Element:
        """Build the MCMC run element.

        Contains state, distribution reference, operators, loggers, and init.
        """
        mcmc = self.config.mcmc

        if mcmc.mcmc_type == MCMCType.NESTED_SAMPLING:
            # Nested sampling is provided by the NS add-on package
            # (github.com/BEAST2-Dev/nested-sampling, class nestedsampling.gss.NS
            # since v1.2.0); it is not part of the BEAST2 core packages. The
            # add-on's own settings are written here, and the shared ones
            # (preBurnin, storeEvery) are no longer dropped for this run type.
            # Note that MCMC has no "seed" *input* — the RNG seed
            # is a command-line argument (`beast -seed N`) — so mcmc.seed is
            # deliberately not serialised here; doing so makes the XML
            # unparseable.
            run_attrs = {
                "id": "mcmc",
                "spec": "nestedsampling.gss.NS",
                "chainLength": str(mcmc.chain_length),
                "preBurnin": str(mcmc.pre_burnin),
                "storeEvery": str(mcmc.store_every),
                "particleCount": str(mcmc.particle_count),
                "subChainLength": str(mcmc.sub_chain_length),
            }
        else:
            run_attrs = {
                "id": "mcmc",
                "spec": "MCMC",
                "chainLength": str(mcmc.chain_length),
                "preBurnin": str(mcmc.pre_burnin),
                "storeEvery": str(mcmc.store_every),
            }

        run_elem = ET.Element("run", run_attrs)

        # Sample from prior
        if mcmc.sample_from_prior:
            run_elem.set("sampleFromPrior", "true")

        # Init element
        init_elem = ModelBuilder.build_init(self.config, self.partition_mgr, self.tree_id)
        run_elem.append(init_elem)

        # State element
        state_elem = ModelBuilder.build_state(
            self.config,
            self.partition_mgr,
            self._state_node_params,
        )
        run_elem.append(state_elem)

        # Distribution reference
        ET.SubElement(run_elem, "distribution", {"idref": "posterior"})

        # Operators
        all_state_ids = [sn[0] for sn in self._state_node_params]
        operators = ModelBuilder.build_operators(self.config, self.partition_mgr, all_state_ids)
        for op in operators:
            run_elem.append(op)

        # Invariant: every estimated state node must be targeted by at least
        # one operator. BEAST2 only prints a warning for orphans and pins the
        # parameter to its initial value, producing silently wrong results.
        targeted: set = set()
        for op in operators:
            for attr in ("parameter", "intparameter", "tree", "taxonset"):
                v = op.get(attr)
                if v and v.startswith("@"):
                    targeted.add(v[1:])
            for child in op:
                idref = child.get("idref")
                if idref:
                    targeted.add(idref)

        def _estimated(node: Dict[str, Any]) -> bool:
            estimate = node.get("estimate", True)
            if isinstance(estimate, str):
                return estimate.lower() not in ("false", "no", "0")
            return bool(estimate)

        orphans = [
            sn_id
            for sn_id, sn_params in self._state_node_params
            if sn_id not in targeted and sn_id != self.tree_id and _estimated(sn_params)
        ]
        if orphans:
            raise ValueError(f"Estimated state nodes without operators: {orphans}")
        # The mirror invariant: an operator must not target a parameter the user
        # asked to keep fixed.
        fixed = {sn_id for sn_id, sn in self._state_node_params if not _estimated(sn)}
        if fixed & targeted:
            raise ValueError(
                f"Operators target fixed (estimate: false) parameters: "
                f"{', '.join(sorted(fixed & targeted))}"
            )

        # Loggers
        # Trace log
        trace_log = self._build_trace_log()
        run_elem.append(trace_log)

        # Tree log
        tree_log = self._build_tree_log()
        run_elem.append(tree_log)

        # Screen log
        screen_log = self._build_screen_log()
        run_elem.append(screen_log)

        return run_elem

    def _build_trace_log(self) -> ET.Element:
        """Build the trace log element."""
        trace_cfg = self.config.loggers.trace_log
        file_name = trace_cfg.get("file_name", "output.$(seed).log")
        log_every = str(trace_cfg.get("log_every", 1000))

        logger = ET.Element(
            "logger",
            {
                "id": "tracelog",
                "fileName": file_name,
                "logEvery": log_every,
                "mode": "autodetect",
                "model": ref("posterior"),
            },
        )

        # Standard log entries
        ET.SubElement(logger, "log", {"idref": "posterior"})
        ET.SubElement(logger, "log", {"idref": "prior"})
        ET.SubElement(logger, "log", {"idref": "likelihood"})

        # Tree height
        ET.SubElement(
            logger,
            "log",
            {
                "id": "treeHeight",
                "spec": "beast.base.evolution.tree.TreeHeightLogger",
                "tree": ref(self.tree_id),
            },
        )

        # MRCA times (skip for CalibratedYule method - no MRCAPrior elements)
        if self.config.calibration_method != CalibrationMethod.CALIBRATED_YULE:
            for cal in self.config.calibrations:
                if not cal.tipsonly:
                    mrca_log = CalibrationBuilder.build_mrca_logger(cal, self.tree_id)
                    logger.append(mrca_log)

        # All estimated parameters: log exactly the state nodes that were
        # collected while assembling the model. This covers substitution
        # model rates, frequencies, gamma shape, clock parameters, and tree
        # prior parameters, and automatically skips linked partitions and
        # fixed (non-estimated) parameters.
        seen: set = set()
        for sn_id, _sn_params in self._state_node_params:
            if sn_id in seen:
                continue
            seen.add(sn_id)
            ET.SubElement(logger, "parameter", {"idref": sn_id, "name": "log"})

        return logger

    def _build_tree_log(self) -> ET.Element:
        """Build the tree log element."""
        tree_cfg = self.config.loggers.tree_log
        file_name = tree_cfg.get("file_name", "output.$(seed).trees")
        log_every = str(tree_cfg.get("log_every", 1000))

        logger = ET.Element(
            "logger",
            {
                "id": "treelog",
                "fileName": file_name,
                "logEvery": log_every,
                "mode": "tree",
            },
        )

        ET.SubElement(logger, "log", {"idref": self.tree_id})

        return logger

    def _build_screen_log(self) -> ET.Element:
        """Build the screen log element."""
        screen_cfg = self.config.loggers.screen_log
        log_every = str(screen_cfg.get("log_every", 10000))

        logger = ET.Element(
            "logger",
            {
                "id": "screenlog",
                "logEvery": log_every,
                "mode": "autodetect",
                "model": ref("posterior"),
            },
        )

        ET.SubElement(logger, "log", {"idref": "posterior"})
        ET.SubElement(logger, "log", {"idref": "prior"})
        ET.SubElement(logger, "log", {"idref": "likelihood"})

        return logger

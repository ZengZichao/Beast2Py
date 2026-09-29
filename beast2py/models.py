"""Data model definitions for Beast2Py.

This module defines all dataclasses used throughout the tool, representing
alignments, partitions, calibrations, models, and configuration objects.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Union


# Metadata keys that document an analysis without describing it. They are kept
# out of the fingerprint so that crediting a co-author or updating a prose
# summary cannot make one analysis look like two.
NON_SCIENCE_METADATA = frozenset(
    {"author", "date", "output_file", "analysis_description", "description", "notes"}
)


class DataType(Enum):
    """Sequence data type."""

    NUCLEOTIDE = "nucleotide"
    AMINOACID = "aminoacid"


class TreePriorType(Enum):
    """Tree prior model types."""

    YULE = "yule"
    CALIBRATED_YULE = "calibrated_yule"
    BIRTH_DEATH = "birth_death"
    COALESCENT_CONSTANT = "coalescent_constant"
    COALESCENT_EXPONENTIAL = "coalescent_exponential"
    BAYESIAN_SKYLINE = "bayesian_skyline"
    EBSP = "ebsp"
    BD_SKYLINE_SERIAL = "bd_skyline_serial"


class ClockModelType(Enum):
    """Clock model types."""

    STRICT = "strict"
    UCLN = "ucln"  # UC Relaxed LogNormal
    UCE = "uce"  # UC Relaxed Exponential
    RLC = "rlc"  # Random Local Clock


class MCMCType(Enum):
    """MCMC run types."""

    STANDARD = "standard"
    NESTED_SAMPLING = "nested_sampling"


class CalibrationMethod(Enum):
    """Calibration method selection."""

    MRCA_PRIOR = "mrca_prior"
    CALIBRATED_YULE = "calibrated_yule"


class InitTreeType(Enum):
    """Initial tree generation strategy."""

    RANDOM = "random"
    UPGMA = "upgma"
    NEWICK = "newick"


@dataclass
class Sequence:
    """A single sequence entry."""

    taxon: str
    sequence: str


@dataclass
class Alignment:
    """A sequence alignment.

    Attributes:
        id: Alignment identifier.
        data_type: Nucleotide or amino acid.
        sequences: List of Sequence objects.
        filter: Optional FilteredAlignment filter expression (e.g. "3::3").
        data_id: Source alignment ID for a filtered alignment.
        source_file: Path the alignment was read from, if any.
        observed_characters: Non-gap characters actually present in the data,
            filled in by ``sequence.validate_alignment`` so the generator can
            cross-check the declared data type against the model family.
    """

    id: str
    data_type: DataType
    sequences: List[Sequence]
    filter: Optional[str] = None
    data_id: Optional[str] = None
    source_file: Optional[str] = None
    observed_characters: Optional[str] = None

    @property
    def taxa_names(self) -> List[str]:
        """Return list of taxon names."""
        return [seq.taxon for seq in self.sequences]

    @property
    def n_taxa(self) -> int:
        """Return number of taxa."""
        return len(self.sequences)

    @property
    def site_lengths(self) -> List[int]:
        """Return the length of every sequence, in file order."""
        return [len(seq.sequence) for seq in self.sequences]

    @property
    def n_sites(self) -> int:
        """Return number of sites.

        Raises:
            ValueError: If the alignment is ragged. Reporting the first
                sequence's length for ragged data made the printed and
                fingerprinted site count disagree with the file.
        """
        lengths = set(self.site_lengths)
        if not lengths:
            return 0
        if len(lengths) > 1:
            raise ValueError(
                f"Alignment '{self.id}' is ragged: sequence lengths are "
                f"{sorted(lengths)}. Use site_lengths after validating the alignment."
            )
        return lengths.pop()

    def content_digest(self) -> str:
        """Return a SHA-256 digest prefix of the alignment's *content*.

        Taxon names and sequence characters both feed the digest, so two
        datasets of the same shape but different sites produce different
        fingerprints ( /.

        Returns:
            First 16 hex characters of the digest.
        """
        import hashlib

        hasher = hashlib.sha256()
        for seq in sorted(self.sequences, key=lambda s: s.taxon):
            hasher.update(seq.taxon.encode("utf-8"))
            hasher.update(b"\0")
            hasher.update(seq.sequence.encode("utf-8"))
            hasher.update(b"\0")
        return hasher.hexdigest()[:16]


@dataclass
class RealParameter:
    """A real-valued parameter with optional bounds.

    Attributes:
        value: Parameter value (float or space-separated string for multi-dim).
        lower: Lower bound.
        upper: Upper bound.
        dimension: Parameter dimension.
        estimate: Whether to estimate this parameter.
    """

    value: Union[float, str] = 1.0
    lower: Optional[float] = None
    upper: Optional[float] = None
    dimension: int = 1
    estimate: bool = True


@dataclass
class DistributionConfig:
    """Configuration for a calibration distribution.

    Attributes:
        type: Distribution type name (normal, lognormal, uniform, etc.).
        parameters: Distribution-specific parameters.
        offset: Offset for the distribution (default 0.0).
    """

    type: str
    parameters: Dict[str, Union[float, bool, str]] = field(default_factory=dict)
    offset: float = 0.0


@dataclass
class HyperpriorConfig:
    """Hyperprior configuration for a distribution parameter.

    Attributes:
        type: Distribution type for the hyperprior.
        parameters: Parameters for the hyperprior distribution.
        offset: Offset for the hyperprior distribution.
    """

    type: str
    parameters: Dict[str, Union[float, bool, str]] = field(default_factory=dict)
    offset: float = 0.0


@dataclass
class Provenance:
    """Calibration provenance metadata (core innovation).

    Attributes:
        source: Source type (fossil, biogeographic, secondary, other).
        reference: DOI or URL reference.
        calibration_type: Soft or hard boundary.
        original_age_min: Original minimum age.
        original_age_max: Original maximum age.
        notes: Free-text notes.
    """

    source: str = ""
    reference: str = ""
    calibration_type: str = "soft"
    original_age_min: Optional[float] = None
    original_age_max: Optional[float] = None
    notes: str = ""


@dataclass
class CalibrationPoint:
    """A node calibration point.

    Attributes:
        name: Calibration point name (used as ID).
        taxa: List of taxon names, or None for all taxa (root).
        monophyletic: Whether to enforce monophyly.
        distribution: Calibration distribution configuration.
        use_originate: Whether to use the originate (stem) node.
        tipsonly: Whether this is a tip-date calibration.
        provenance: Provenance metadata.
        hyperpriors: Hyperprior configurations keyed by parameter name.
    """

    name: str
    taxa: Optional[List[str]]
    monophyletic: bool = True
    distribution: Optional[DistributionConfig] = None
    use_originate: bool = False
    tipsonly: bool = False
    provenance: Optional[Provenance] = None
    hyperpriors: Dict[str, HyperpriorConfig] = field(default_factory=dict)


@dataclass
class SiteModelConfig:
    """Site model configuration.

    Attributes:
        substitution_model: Dict with 'type' and model-specific parameters.
        gamma_categories: Number of gamma rate categories (0 = no gamma).
        gamma_shape: Gamma shape parameter.
        proportion_invariant: Proportion of invariant sites.
        linked_to: Partition ID whose site model is shared, or None.
    """

    substitution_model: Dict[str, Any]
    gamma_categories: int = 0
    gamma_shape: Optional[RealParameter] = None
    proportion_invariant: Optional[RealParameter] = None
    linked_to: Optional[str] = None


@dataclass
class ClockModelConfig:
    """Clock model configuration.

    Attributes:
        type: Clock model type.
        clock_rate: Clock rate parameter.
        ucld_mean: UCLD mean parameter (for relaxed clocks).
        ucld_stdev: UCLD stdev parameter (for relaxed clocks).
        linked_to: Partition ID to share clock model with, or None.
    """

    type: ClockModelType = ClockModelType.STRICT
    clock_rate: Optional[RealParameter] = None
    ucld_mean: Optional[RealParameter] = None
    ucld_stdev: Optional[RealParameter] = None
    linked_to: Optional[str] = None


@dataclass
class Partition:
    """A partition combining alignment, site model, and clock model.

    Attributes:
        id: Partition identifier.
        alignment: Alignment object.
        site_model: Site model configuration.
        clock_model: Clock model configuration.
        tree: Tree sharing setting. Only ``"shared"`` is supported today:
            Beast2Py assembles one linked topology for all partitions, and the
            parser rejects any other value instead of silently ignoring it
            . Per-partition unlinked trees are on the roadmap.
    """

    id: str
    alignment: Alignment
    site_model: SiteModelConfig
    clock_model: ClockModelConfig
    tree: str = "shared"


@dataclass
class MCMCConfig:
    """MCMC run configuration.

    Attributes:
        chain_length: Number of MCMC generations.
        pre_burnin: Pre-burnin length.
        store_every: State storage interval **in logged samples** (MCMC.java
            tests ``(sampleNr + 1) % storeEvery == 0``), not in generations; 0
            or -1 disables checkpointing. The default scales with the chain
            length so a long run gets ~10 resume points instead of one at the
            very end.
        sample_from_prior: Whether to sample from prior only.
        mcmc_type: Standard or nested sampling.
        particle_count: NS particle count.
        sub_chain_length: NS sub-chain length.
        seed: RNG seed written to the XML, so log file names built from the
            ``$(seed)`` macro and the run itself are reproducible.
    """

    chain_length: int = 10000000
    pre_burnin: int = 0
    store_every: int = 10000
    sample_from_prior: bool = False
    mcmc_type: MCMCType = MCMCType.STANDARD
    particle_count: int = 1
    sub_chain_length: int = 10000
    seed: int = 1


@dataclass
class LoggerConfig:
    """Logger configuration.

    Attributes:
        trace_log: Trace log file name and log interval.
        tree_log: Tree log file name and log interval.
        screen_log: Screen log interval.
    """

    trace_log: Dict[str, Any] = field(
        default_factory=lambda: {"file_name": "output.$(seed).log", "log_every": 1000}
    )
    tree_log: Dict[str, Any] = field(
        default_factory=lambda: {"file_name": "output.$(seed).trees", "log_every": 1000}
    )
    screen_log: Dict[str, Any] = field(default_factory=lambda: {"log_every": 10000})


@dataclass
class TipDatesConfig:
    """Tip dates (TraitSet) configuration.

    Attributes:
        enabled: Whether tip dates are enabled.
        trait_name: Trait name (date-forward or date-backward).
        units: Time units (year, month, day).
        dates: Dict mapping taxon name to date value.
        file: Optional CSV file with tip dates.
    """

    enabled: bool = False
    trait_name: str = "date-forward"
    units: str = "year"
    dates: Dict[str, float] = field(default_factory=dict)
    file: Optional[str] = None


@dataclass
class InitializationConfig:
    """Initial tree configuration.

    Attributes:
        tree_type: Initialization strategy (random, upgma, newick).
        newick_file: Path to Newick tree file (for tree_type=newick).
    """

    tree_type: InitTreeType = InitTreeType.RANDOM
    newick_file: Optional[str] = None
    # SHA-256 of the starting-tree file's bytes. Only the *type* of initial tree
    # used to reach the fingerprint, so two analyses seeded with completely
    # different user topologies shared one identity.
    newick_digest: Optional[str] = None


@dataclass
class BEASTConfig:
    """Complete BEAST2 analysis configuration.

    This is the top-level configuration object that aggregates all sub-configurations.

    Attributes:
        metadata: Analysis metadata dict.
        partitions: List of Partition objects.
        tree_prior_type: Tree prior type.
        tree_prior_params: Tree prior parameters.
        calibrations: List of calibration points.
        calibration_method: Calibration method (mrca_prior or calibrated_yule).
        parameter_priors: List of parameter prior configurations.
        mcmc: MCMC configuration.
        loggers: Logger configuration.
        tip_dates: Tip dates configuration.
        initialization: Initialization configuration.
        operator_weights: Operator weight overrides.
    """

    metadata: Dict[str, Any]
    partitions: List[Partition]
    tree_prior_type: TreePriorType
    tree_prior_params: Dict[str, Any]
    calibrations: List[CalibrationPoint]
    all_alignments: List[Alignment] = field(default_factory=list)
    calibration_method: CalibrationMethod = CalibrationMethod.MRCA_PRIOR
    parameter_priors: List[Dict[str, Any]] = field(default_factory=list)
    mcmc: MCMCConfig = field(default_factory=MCMCConfig)
    loggers: LoggerConfig = field(default_factory=LoggerConfig)
    tip_dates: Optional[TipDatesConfig] = None
    initialization: InitializationConfig = field(default_factory=InitializationConfig)
    operator_weights: Dict[str, float] = field(default_factory=dict)
    diagnostics: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize the *scientific* content of the config for fingerprinting.

        Everything that changes the analysis must change the fingerprint, and
        nothing that does not may. Accordingly the alignment content digest,
        per-partition tree reference, calibration stem/crown and tipsonly
        flags, hyperpriors, explicit parameter priors and tip dates are all
        included, while the output file name is excluded.

        Returns:
            A JSON-compatible dict.
        """
        import json

        def _default(obj: Any) -> Any:
            if isinstance(obj, Enum):
                return obj.value
            if hasattr(obj, "__dict__"):
                return {k: v for k, v in obj.__dict__.items()}
            return str(obj)

        # Only the fields that describe the *analysis* identify it. The
        # documentary ones (who ran it, when, prose about it) used to change the
        # fingerprint, contradicting this module's own "deliberately contains no
        # calendar date" claim and splitting one analysis across several ids.
        metadata = {
            k: v for k, v in self.metadata.items() if k not in NON_SCIENCE_METADATA
        }

        def _param_repr(param: Any) -> Any:
            if param is None:
                return None
            if hasattr(param, "__dict__"):
                return dict(param.__dict__)
            return param

        return json.loads(
            json.dumps(
                {
                    "metadata": metadata,
                    "partitions": [
                        {
                            "id": p.id,
                            "alignment_id": p.alignment.id,
                            "data_type": p.alignment.data_type.value,
                            "n_taxa": p.alignment.n_taxa,
                            "n_sites": p.alignment.n_sites,
                            "alignment_digest": p.alignment.content_digest(),
                            "filter": p.alignment.filter,
                            "tree": p.tree,
                            "site_model": {
                                "substitution_model": p.site_model.substitution_model,
                                "gamma_categories": p.site_model.gamma_categories,
                                "gamma_shape": _param_repr(p.site_model.gamma_shape),
                                "proportion_invariant": _param_repr(
                                    p.site_model.proportion_invariant
                                ),
                                "linked_to": p.site_model.linked_to,
                            },
                            "clock_model": {
                                "type": p.clock_model.type.value,
                                "clock_rate": _param_repr(p.clock_model.clock_rate),
                                "ucld_mean": _param_repr(p.clock_model.ucld_mean),
                                "ucld_stdev": _param_repr(p.clock_model.ucld_stdev),
                                "linked_to": p.clock_model.linked_to,
                            },
                        }
                        for p in self.partitions
                    ],
                    "tree_prior_type": self.tree_prior_type.value,
                    "tree_prior_params": self.tree_prior_params,
                    "calibrations": [
                        {
                            "name": c.name,
                            "taxa": c.taxa,
                            "monophyletic": c.monophyletic,
                            "use_originate": c.use_originate,
                            "tipsonly": c.tipsonly,
                            "distribution": c.distribution.__dict__ if c.distribution else None,
                            "hyperpriors": {
                                key: hp.__dict__ for key, hp in sorted(c.hyperpriors.items())
                            },
                            "provenance": c.provenance.__dict__ if c.provenance else None,
                        }
                        for c in self.calibrations
                    ],
                    "calibration_method": self.calibration_method.value,
                    "parameter_priors": self.parameter_priors,
                    "mcmc": self.mcmc.__dict__,
                    "loggers": {
                        "trace_log": self.loggers.trace_log,
                        "tree_log": self.loggers.tree_log,
                        "screen_log": self.loggers.screen_log,
                    },
                    "tip_dates": self.tip_dates.__dict__ if self.tip_dates else None,
                    "init_tree_type": self.initialization.tree_type.value,
                    "init_tree_newick_digest": self.initialization.newick_digest,
                    "operator_weights": self.operator_weights,
                },
                default=_default,
                sort_keys=True,
            )
        )

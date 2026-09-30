"""Configuration file parsing and validation for Beast2Py.

Parses YAML configuration files, validates parameters, fills defaults,
and constructs BEASTConfig objects.

Validation is deliberately *fail-closed*: a key that is not understood, a
value that is not a boolean, a ragged alignment or a malformed calibration
raises :class:`ConfigError` instead of silently falling back to a default.
Silent fallback was the single largest source of wrong-but-valid XML in the

"""

from __future__ import annotations

import hashlib
import math
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import yaml

from .models import (
    BEASTConfig,
    CalibrationMethod,
    CalibrationPoint,
    ClockModelConfig,
    ClockModelType,
    DistributionConfig,
    HyperpriorConfig,
    InitializationConfig,
    InitTreeType,
    LoggerConfig,
    MCMCConfig,
    MCMCType,
    Partition,
    Provenance,
    RealParameter,
    SiteModelConfig,
    TipDatesConfig,
    TreePriorType,
)
from .registry import ModelRegistry
from .utils import as_number, as_strict_bool


class ConfigError(Exception):
    """Configuration parsing or validation error."""

    pass


# Sections whose keys are enumerated exhaustively. Any key not listed here is a
# configuration error (: a typo such as ``gamma_categorys`` used to
# be silently replaced by its default).
KNOWN_KEYS: Dict[str, Set[str]] = {
    "": {
        "metadata",
        "alignments",
        "partitions",
        "tree_prior",
        "calibrations",
        "calibration_method",
        "parameter_priors",
        "priors",
        "mcmc",
        "loggers",
        "tip_dates",
        "initialization",
        "operator_weights",
        "diagnostics",
        "output",
    },
    "metadata": {
        "analysis_name",
        "analysis_description",
        "author",
        "date",
        "beast2_version",
        "tool_version",
        "output_file",
        "embed_fingerprint",
    },
    "alignment": {
        "id",
        "file",
        "format",
        "data_type",
        "filter",
        "data_id",
    },
    "partition": {"id", "site_model", "clock_model", "tree"},
    "site_model": {
        "substitution_model",
        "gamma_categories",
        "gamma_shape",
        "proportion_invariant",
        "linked_to",
    },
    "clock_model": {
        "type",
        "clock_rate",
        "ucld_mean",
        "ucld_stdev",
        "linked_to",
    },
    "parameter": {"value", "lower", "upper", "dimension", "estimate", "mode"},
    "tree_prior": {
        "type",
        "birth_rate",
        "death_rate",
        "sampling_rate",
        "sample_probability",
        "pop_size",
        "growth_rate",
        "pop_sizes",
        "group_sizes",
        "rho",
    },
    "calibration": {
        "name",
        "taxa",
        "monophyletic",
        "distribution",
        "use_originate",
        "tipsonly",
        "provenance",
    },
    "distribution": {"type", "parameters", "offset", "hyperpriors"},
    "provenance": {
        "source",
        "reference",
        "calibration_type",
        "original_age_min",
        "original_age_max",
        "notes",
    },
    "hyperprior": {"type", "parameters", "offset"},
    "mcmc": {
        "chain_length",
        "pre_burnin",
        "store_every",
        "sample_from_prior",
        "type",
        "particle_count",
        "sub_chain_length",
        "seed",
    },
    "loggers": {"trace_log", "tree_log", "screen_log"},
    "logger": {"file_name", "log_every"},
    "tip_dates": {"enabled", "trait_name", "units", "dates", "file"},
    "initialization": {"tree_type", "newick_file"},
    "diagnostics": {
        "overlap_measure",
        "overlap_threshold",
        "relative_overlap_threshold",
    },
    "parameter_prior": {"parameter", "distribution"},
    "output": {"file", "embed_fingerprint"},
}

VALID_TRAIT_NAMES = {"date-forward", "date-backward"}
VALID_DATE_UNITS = {"year", "month", "day", "coordinate"}


class _RegistryNames:
    """A live read-only view over a :class:`ModelRegistry` name table.

    A plain ``frozenset(registry._distributions)`` would snapshot the table when
    this module is first imported, so a later ``register_*`` call still would not
    be visible to the config validator -- which is the bug being fixed. Membership
    and iteration resolve against the registry on every use instead
    .

    Args:
        table: Attribute name on ``ModelRegistry`` holding the dict to mirror.
    """

    __slots__ = ("_table",)

    def __init__(self, table: str) -> None:
        self._table = table

    def _names(self):
        return getattr(ModelRegistry, self._table).keys()

    def __contains__(self, item) -> bool:
        return item in self._names()

    def __iter__(self):
        return iter(self._names())

    def __len__(self) -> int:
        return len(self._names())

    def __repr__(self) -> str:
        return f"_RegistryNames({sorted(self._names())!r})"


class ConfigParser:
    """Parse and validate YAML configuration files."""

    # Supported model type sets.
    #
    # Substitution models and calibration distributions are *derived from
    # ModelRegistry* rather than repeated here. They used to be a second,
    # hard-coded copy of the same names, which made the advertised extension
    # points dead: `ModelRegistry.register_substitution_model("mygtr", ...)`
    # taught the registry but every config still failed with "invalid
    # substitution model", while the registration leaked into later analyses in
    # the same process.
    VALID_SUBSTITUTION_MODELS = _RegistryNames("_substitution_models")
    VALID_CLOCK_MODELS = {"strict", "ucln", "uce", "rlc"}
    VALID_TREE_PRIORS = {
        "yule",
        "calibrated_yule",
        "birth_death",
        "coalescent_constant",
        "coalescent_exponential",
        "bayesian_skyline",
        "ebsp",
        "bd_skyline_serial",
    }
    VALID_DISTRIBUTIONS = _RegistryNames("_distributions")
    VALID_DATA_TYPES = {"nucleotide", "aminoacid"}

    # Parameters a distribution cannot be sampled without.  BEAST2 fills in its
    # own defaults for anything it does not receive, so a missing required key
    # silently changes which prior is fitted.
    DISTRIBUTION_REQUIRED = {
        "normal": {"mean", "sigma"},
        "lognormal": {"M", "S"},
        "uniform": {"lower", "upper"},
        "exponential": {"mean"},
        "gamma": {"alpha", "beta"},
        "beta": {"alpha", "beta"},
        "laplace": {"mu", "scale"},
        "inverse_gamma": {"alpha", "beta"},
        "poisson": {"lambda"},
        "chi_square": {"dof"},
        "one_on_x": set(),
    }
    # Keys accepted on top of the registry's parameter names: real BEAST2
    # boolean attributes and the min/max spelling of a Uniform's bounds.
    DISTRIBUTION_EXTRA_KEYS = {
        "lognormal": {"mean_in_real_space"},
        "uniform": {"min", "max"},
    }

    @staticmethod
    def _norm_name(value: Any) -> Any:
        """Fold a model/prior name to lower case with surrounding space removed.

        Substitution models, clock models, distributions and tree priors are
        all matched case-insensitively, so a name resolves to the same key
        however it is capitalised. Authors should not have to guess which rule
        applies.

        Args:
            value: The configured name; non-strings pass through untouched.

        Returns:
            The normalised name, or ``value`` when it is not a string.
        """
        return value.strip().lower() if isinstance(value, str) else value

    @staticmethod
    def _canonical_name(value: Any, names, field: str, where: str) -> str:
        """Fold a user-typed model/distribution name onto the tool's own key.

        Every named choice in a configuration file is meant to be matched
        case-insensitively, but three call sites compared the raw string against
        a lower-case name table: ``type: Strict`` was refused for the clock while
        ``type: HKY`` was accepted for the site model, and the calibration
        distributions a reader copies straight out of Supplementary Table S1
        (``Normal``, ``LogNormal``, ``ChiSquare``) all failed.  Separators are
        ignored so the printed forms ``OneOnX`` / ``InverseGamma`` / ``ChiSquare``
        resolve to ``one_on_x`` / ``inverse_gamma`` / ``chi_square``.

        Args:
            value: The configured name.
            names: Iterable of canonical names to resolve against.
            field: Name used in the error message.
            where: Location used in the error message.

        Returns:
            The canonical key, so downstream lookups keep seeing lower-case names.

        Raises:
            ConfigError: If no canonical name matches, listing the valid ones.
        """
        if not isinstance(value, str):
            raise ConfigError(
                f"{where}: {field} must be a string, got {value!r}"
            )
        folded = "".join(ch for ch in value if ch not in " _-").lower()
        for name in names:
            if "".join(ch for ch in name if ch not in " _-").lower() == folded:
                return name
        raise ConfigError(
            f"{where}: invalid {field} '{value}' "
            f"(valid: {', '.join(sorted(names))})"
        )

    @staticmethod
    def _parse_enum(enum_cls, value: Any, field: str) -> Any:
        """Parse an enum value, raising ConfigError with the list of valid
        options instead of a bare ValueError when the value is invalid.

        Matching is case-insensitive so that every named choice in a config
        file follows the same rule ."""
        try:
            return enum_cls(value)
        except ValueError:
            pass
        if isinstance(value, str):
            wanted = value.strip().lower()
            for member in enum_cls:
                if isinstance(member.value, str) and member.value.lower() == wanted:
                    return member
            valid = ", ".join(m.value for m in enum_cls)
            raise ConfigError(f"Invalid {field}: '{value}' (valid: {valid})")

    def __init__(self, base_dir: Optional[str] = None):
        """Initialize the config parser.

        Args:
            base_dir: Base directory for resolving relative file paths.
        """
        self.base_dir = Path(base_dir) if base_dir else Path.cwd()

    # -- key checking -------------------------------------------------------

    @staticmethod
    def _check_keys(raw: Dict[str, Any], section: str, where: str) -> None:
        """Reject unknown keys in a configuration section.

        Args:
            raw: The mapping to check.
            section: Key into :data:`KNOWN_KEYS` naming the section.
            where: Human-readable location used in the error message.

        Raises:
            ConfigError: If any key is not part of the section's schema.
        """
        known = KNOWN_KEYS.get(section)
        if known is None or not isinstance(raw, dict):
            return
        unknown = [k for k in raw if k not in known]
        if unknown:
            raise ConfigError(
                f"{where}: unknown key(s) {', '.join(sorted(map(str, unknown)))}. "
                f"Allowed keys for this section: {', '.join(sorted(known))}. "
                f"Beast2Py refuses unknown keys rather than silently applying a "
                f"default, because a typo in a model or MCMC setting changes the "
                f"analysis without changing the XML's validity."
            )

    @staticmethod
    def _as_bool(raw: Dict[str, Any], key: str, default: bool, where: str) -> bool:
        """Read a boolean from a mapping with strict coercion.

        Args:
            raw: The mapping to read.
            key: Key name.
            default: Value used when the key is absent.
            where: Location used in error messages.

        Returns:
            The boolean value.

        Raises:
            ConfigError: If the value is present but not a boolean.
        """
        if key not in raw or raw[key] is None:
            return default
        try:
            return as_strict_bool(raw[key], f"{where}.{key}")
        except ValueError as exc:
            raise ConfigError(str(exc)) from exc

    @staticmethod
    def _as_num(raw: Dict[str, Any], key: str, where: str, default: Optional[float] = None):
        """Read a number from a mapping, rejecting booleans and stray strings.

        Args:
            raw: The mapping to read.
            key: Key name.
            where: Location used in error messages.
            default: Returned when the key is absent; ``None`` means "absent".

        Returns:
            The float value, or ``default``.

        Raises:
            ConfigError: If the value is not a finite number.
        """
        if key not in raw or raw[key] is None:
            return default
        try:
            value = as_number(raw[key], f"{where}.{key}")
        except ValueError as exc:
            raise ConfigError(str(exc)) from exc
        # The docstring always promised to reject non-finite numbers but never
        # did, so `kappa: .inf` escaped as a raw ValueError from a later int()
        # cast and reached the user as a Python traceback.
        if not math.isfinite(value):
            raise ConfigError(
                f"'{where}.{key}' must be a finite number, got {raw[key]!r}"
            )
        return value

    def parse(self, config_path: str) -> BEASTConfig:
        """Parse a YAML configuration file and return a BEASTConfig object.

        Args:
            config_path: Path to the YAML configuration file.

        Returns:
            A BEASTConfig object.

        Raises:
            ConfigError: If the configuration is invalid.
        """
        config_path = Path(config_path)
        if not config_path.is_absolute():
            config_path = self.base_dir / config_path

        if not config_path.exists():
            raise ConfigError(f"Configuration file not found: {config_path}")

        try:
            with open(config_path, "r", encoding="utf-8") as f:
                raw = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ConfigError(f"Invalid YAML in {config_path}: {e}")

        if not isinstance(raw, dict):
            raise ConfigError("Configuration root must be a mapping")

        # Set base dir to config file's directory for resolving relative paths
        self.base_dir = config_path.parent

        return self._build_config(raw)

    def parse_dict(self, raw: Dict[str, Any]) -> BEASTConfig:
        """Parse a configuration dictionary directly.

        Args:
            raw: Configuration dictionary.

        Returns:
            A BEASTConfig object.
        """
        return self._build_config(raw)

    def _build_config(self, raw: Dict[str, Any]) -> BEASTConfig:
        """Build a BEASTConfig from a raw dictionary."""
        errors: List[str] = []
        self._check_keys(raw, "", "Configuration")

        # --- Metadata ---
        metadata = dict(raw.get("metadata", {}) or {})
        self._check_keys(metadata, "metadata", "metadata")
        if not metadata.get("analysis_name"):
            metadata["analysis_name"] = "unnamed_analysis"
        metadata.setdefault("beast2_version", "2.7.8")
        metadata.setdefault("tool_version", "0.1.0")

        # --- Alignments and Partitions ---
        alignments_raw = raw.get("alignments", [])
        partitions_raw = raw.get("partitions", [])

        if not alignments_raw:
            errors.append("At least one alignment must be specified")
        if not partitions_raw:
            errors.append("At least one partition must be specified")

        if errors:
            raise ConfigError("; ".join(errors))

        # Read sequence files
        from .sequence import AlignmentError, SequenceReader, validate_alignment

        alignments: Dict[str, Any] = {}
        # Ids that were auto-registered as the unfiltered parent of a filtered
        # alignment.  A later *explicit* declaration of one of these is the
        # author naming that source, not a duplicate (0).
        implicit_alignment_ids: Set[str] = set()
        for index, aln_raw in enumerate(alignments_raw):
            where = f"alignments[{index}]"
            if not isinstance(aln_raw, dict) or "id" not in aln_raw:
                raise ConfigError(f"{where}: each alignment must be a mapping with an 'id'")
            self._check_keys(aln_raw, "alignment", where)
            aln_id = aln_raw["id"]
            aln_file = aln_raw.get("file", "")
            aln_format = aln_raw.get("format", "auto")
            aln_data_type = aln_raw.get("data_type", "nucleotide")
            aln_filter = aln_raw.get("filter")

            if aln_data_type not in self.VALID_DATA_TYPES:
                raise ConfigError(
                    f"Alignment '{aln_id}': invalid data_type '{aln_data_type}' "
                    f"(valid: {', '.join(sorted(self.VALID_DATA_TYPES))})"
                )

            # Resolve file path
            if aln_file:
                aln_path = Path(aln_file)
                if not aln_path.is_absolute():
                    aln_path = self.base_dir / aln_path
                if not aln_path.exists():
                    raise ConfigError(f"Alignment '{aln_id}': file not found: {aln_path}")
                base_id = aln_raw.get("data_id") or f"{aln_id}.unfiltered"
                try:
                    base_alignment = SequenceReader.read(
                        str(aln_path),
                        format=aln_format,
                        data_type=aln_data_type,
                        filter=None,
                        # Without a filter the file *is* the alignment, so it
                        # keeps the declared id.
                        alignment_id=base_id if aln_filter else aln_id,
                    )
                except AlignmentError as exc:
                    raise ConfigError(f"Alignment '{aln_id}': {exc}") from exc
                base_alignment.source_file = str(aln_path)
                # Compareability checks run on the real data, before any filter.
                try:
                    validate_alignment(base_alignment)
                except AlignmentError as exc:
                    raise ConfigError(str(exc)) from exc

                if aln_filter:
                    # An alignment declaring both a file and a filter is a
                    # FilteredAlignment over that file. The unfiltered source is
                    # registered under its own id so the filter never points at
                    # itself and the sequence data is not dropped.
                    if base_id == aln_id:
                        raise ConfigError(
                            f"Alignment '{aln_id}': data_id must differ from the alignment's "
                            f"own id when a filter is declared"
                        )
                    from .models import Alignment as AlignmentModel

                    if base_id not in alignments:
                        alignments[base_id] = base_alignment
                        implicit_alignment_ids.add(base_id)
                    alignment = AlignmentModel(
                        id=aln_id,
                        data_type=base_alignment.data_type,
                        sequences=base_alignment.sequences,
                        filter=aln_filter,
                        observed_characters=base_alignment.observed_characters,
                    )
                    alignment.data_id = base_id
                    self._validate_filter_sites(alignment, aln_filter, aln_id)
                else:
                    alignment = base_alignment
            else:
                # Alignment without file: a FilteredAlignment over a base one
                if aln_filter:
                    base_id = aln_raw.get("data_id", "")
                    if not base_id:
                        raise ConfigError(
                            f"Alignment '{aln_id}': 'filter' requires 'data_id' naming the "
                            f"source alignment; without it the filter would reference itself"
                        )
                    if base_id == aln_id:
                        raise ConfigError(
                            f"Alignment '{aln_id}': 'data_id' must name a different alignment"
                        )
                    if base_id not in alignments:
                        raise ConfigError(
                            f"Alignment '{aln_id}': 'data_id' references unknown alignment "
                            f"'{base_id}' (declare the source alignment before its filter)"
                        )
                    from .models import Alignment as AlignmentModel

                    base_aln = alignments[base_id]
                    alignment = AlignmentModel(
                        id=aln_id,
                        data_type=base_aln.data_type,
                        sequences=base_aln.sequences,
                        filter=aln_filter,
                        observed_characters=base_aln.observed_characters,
                    )
                    self._validate_filter_sites(alignment, aln_filter, aln_id)
                else:
                    raise ConfigError(f"Alignment '{aln_id}': no file or filter specified")

            explicit_data_id = aln_raw.get("data_id")
            if explicit_data_id is not None:
                alignment.data_id = explicit_data_id
            # Otherwise keep the data_id the branch above already resolved: a
            # file+filter alignment owns an implicit unfiltered parent, and
            # clobbering it with None deferred the failure to generation, where
            # it surfaced as an unrelated "no matching alignment found".

            if aln_id in alignments and aln_id not in implicit_alignment_ids:
                raise ConfigError(f"Duplicate alignment id: {aln_id}")
            # An implicit parent may be re-declared explicitly further down the
            # list; declaration order must not decide whether the config is
            # valid (0).
            implicit_alignment_ids.discard(aln_id)
            alignments[aln_id] = alignment

        # --- Partitions ---
        partitions: List[Partition] = []
        partition_ids: Set[str] = set()
        for index, p_raw in enumerate(partitions_raw):
            where = f"partitions[{index}]"
            if not isinstance(p_raw, dict) or "id" not in p_raw:
                raise ConfigError(f"{where}: each partition must be a mapping with an 'id'")
            self._check_keys(p_raw, "partition", where)
            p_id = p_raw["id"]
            if p_id in partition_ids:
                raise ConfigError(f"Duplicate partition ID: {p_id}")
            partition_ids.add(p_id)

            if p_id not in alignments:
                raise ConfigError(f"Partition '{p_id}': no matching alignment found")

            alignment = alignments[p_id]

            # A second partition reusing the same sites would count them in the
            # likelihood twice, which raises the posterior to a power instead of
            # sharpening it. Enforced after the loop by
            # _validate_partition_site_coverage, which needs every partition
            # present before it can compare coverage pairwise.
            site_model = self._parse_site_model(
                p_raw.get("site_model", {}) or {}, errors, p_id, alignment
            )
            clock_model = self._parse_clock_model(p_raw.get("clock_model", {}) or {}, errors, p_id)

            tree_ref = p_raw.get("tree", "shared")
            self._validate_tree_ref(tree_ref, p_id, partition_ids)

            partitions.append(
                Partition(
                    id=p_id,
                    alignment=alignment,
                    site_model=site_model,
                    clock_model=clock_model,
                    tree=tree_ref,
                )
            )

        self._validate_partition_site_coverage(partitions, alignments)

        # Validate linked_to references
        for p in partitions:
            if p.clock_model.linked_to:
                if p.clock_model.linked_to not in partition_ids:
                    errors.append(
                        f"Partition '{p.id}': clock_model.linked_to references "
                        f"non-existent partition '{p.clock_model.linked_to}'"
                    )
                elif p.clock_model.linked_to == p.id:
                    errors.append(
                        f"Partition '{p.id}': clock_model.linked_to must "
                        f"reference another partition"
                    )
            if p.site_model.linked_to:
                if p.site_model.linked_to not in partition_ids:
                    errors.append(
                        f"Partition '{p.id}': site_model.linked_to references "
                        f"non-existent partition '{p.site_model.linked_to}'"
                    )
                elif p.site_model.linked_to == p.id:
                    errors.append(
                        f"Partition '{p.id}': site_model.linked_to must reference another partition"
                    )

        self._validate_linked_model_consistency(partitions, partitions_raw, errors)

        if errors:
            raise ConfigError("; ".join(errors))
        # --- Tree Prior ---
        tp_raw = raw.get("tree_prior", {}) or {}
        self._check_keys(tp_raw, "tree_prior", "tree_prior")
        tp_type_str = ConfigParser._norm_name(tp_raw.get("type", "yule"))
        if tp_type_str not in self.VALID_TREE_PRIORS:
            raise ConfigError(
                f"Invalid tree_prior type: {tp_type_str} "
                f"(valid: {', '.join(sorted(self.VALID_TREE_PRIORS))})"
            )

        tree_prior_type = TreePriorType(tp_type_str)
        tree_prior_params = {k: v for k, v in tp_raw.items() if k != "type"}
        self._validate_tree_prior_params(tree_prior_type, tree_prior_params)

        # --- Calibrations ---
        cals_raw = raw.get("calibrations", []) or []
        calibrations = []
        all_taxa: Set[str] = set()
        for aln in alignments.values():
            all_taxa.update(aln.taxa_names)

        for cal_index, cal_raw in enumerate(cals_raw):
            parsed = self._parse_calibration(
                cal_raw, errors, all_taxa, f"calibrations[{cal_index}]"
            )
            # Calibration names become XML ids, so a repeat only surfaced much
            # later as "Duplicate XML ids generated: X, X.age, X.distr, ..."
            # with advice about renaming taxa -- sending the author to the wrong
            # place for what is a repeated `name:`.
            clash = [c.name for c in calibrations if c.name == parsed.name]
            if clash:
                raise ConfigError(
                    f"Duplicate calibration name: '{parsed.name}' is declared "
                    f"{len(clash) + 1} times (calibrations[{cal_index}]). Give every "
                    f"calibration a distinct 'name:'; the taxa may stay the same."
                )
            calibrations.append(parsed)

        # --- Calibration method ---
        cal_method_str = raw.get("calibration_method", "mrca_prior")
        calibration_method = ConfigParser._parse_enum(
            CalibrationMethod, cal_method_str, "calibration_method"
        )

        # --- Parameter priors ---
        parameter_priors = raw.get("parameter_priors", raw.get("priors", [])) or []
        for pp_index, pp in enumerate(parameter_priors):
            if not isinstance(pp, dict):
                raise ConfigError(f"parameter_priors[{pp_index}]: must be a mapping")
            self._check_keys(pp, "parameter_prior", f"parameter_priors[{pp_index}]")
            self._check_keys(
                pp.get("distribution", {}) or {},
                "distribution",
                f"parameter_priors[{pp_index}].distribution",
            )
            if "parameter" not in pp:
                raise ConfigError(f"parameter_priors[{pp_index}]: missing 'parameter'")

        # --- Identifiability: an absolute time axis needs either a clock rate
        # that is estimated, or at least one absolute calibration. With a
        # fixed clock rate and no calibration the timestamps are unidentified.
        has_calibration = bool(calibrations)
        fixed_clock = all(
            (p.clock_model.clock_rate is not None and not p.clock_model.clock_rate.estimate)
            or p.clock_model.linked_to
            for p in partitions
        ) and any(p.clock_model.clock_rate is not None for p in partitions)
        if fixed_clock and not has_calibration:
            raise ConfigError(
                "The clock rate is fixed and no calibration is supplied, so absolute "
                "time is not identifiable: node heights can only be estimated in "
                "expected-substitutions-per-site units. Either estimate the clock rate "
                "or add at least one absolute (e.g. fossil) calibration."
            )

        # --- MCMC ---
        mcmc_raw = raw.get("mcmc", {}) or {}
        self._check_keys(mcmc_raw, "mcmc", "mcmc")
        mcmc = self._parse_mcmc(mcmc_raw)

        # --- Loggers ---
        loggers_raw = raw.get("loggers", {}) or {}
        self._check_keys(loggers_raw, "loggers", "loggers")
        for logger_name in ("trace_log", "tree_log", "screen_log"):
            if logger_name in loggers_raw:
                self._check_keys(loggers_raw[logger_name], "logger", f"loggers.{logger_name}")
        loggers = LoggerConfig(
            trace_log=dict(
                loggers_raw.get("trace_log", {"file_name": "output.$(seed).log", "log_every": 1000})
            ),
            tree_log=dict(
                loggers_raw.get(
                    "tree_log", {"file_name": "output.$(seed).trees", "log_every": 1000}
                )
            ),
            screen_log=dict(loggers_raw.get("screen_log", {"log_every": 10000})),
        )

        # --- Tip dates ---
        tip_dates = None
        td_raw = raw.get("tip_dates", {}) or {}
        self._check_keys(td_raw, "tip_dates", "tip_dates")
        if self._as_bool(td_raw, "enabled", False, "tip_dates"):
            tip_dates = self._parse_tip_dates(td_raw, partitions)

        # --- Initialization ---
        init_raw = raw.get("initialization", {}) or {}
        self._check_keys(init_raw, "initialization", "initialization")
        init_type_str = init_raw.get("tree_type", "random")
        parsed_init_type = ConfigParser._parse_enum(
            InitTreeType, init_type_str, "initialization.tree_type"
        )
        newick_path = init_raw.get("newick_file")
        newick_digest = None
        if newick_path:
            # Resolve like every other file in a config, and fail here rather
            # than at XML assembly with a bare FileNotFoundError traceback.
            candidate = Path(newick_path)
            if not candidate.is_absolute():
                candidate = self.base_dir / candidate
            if not candidate.exists():
                raise ConfigError(
                    f"initialization.newick_file: file not found: {candidate}"
                )
            newick_bytes = candidate.read_bytes()
            newick_path = str(candidate)
            # The starting tree is written into the XML, so it is part of the
            # analysis; only its type used to reach the fingerprint, letting two
            # different topologies share one identity.
            newick_digest = hashlib.sha256(newick_bytes).hexdigest()
        initialization = InitializationConfig(
            tree_type=parsed_init_type,
            newick_file=newick_path,
            newick_digest=newick_digest,
        )

        # --- Operator weights ---
        operator_weights = raw.get("operator_weights", {}) or {}

        # --- Diagnostics (overlap criteria exposed to the user
        #: the absolute 2.0 Ma threshold was only editable in the
        # source) ---
        diag_raw = raw.get("diagnostics", {}) or {}
        self._check_keys(diag_raw, "diagnostics", "diagnostics")
        diagnostics = dict(diag_raw)
        measure = str(diagnostics.get("overlap_measure", "relative")).lower()
        if measure not in ("relative", "absolute", "both"):
            raise ConfigError(
                f"diagnostics.overlap_measure must be relative, absolute or both, got "
                f"{diagnostics.get('overlap_measure')!r}"
            )
        diagnostics["overlap_measure"] = measure
        for key in ("overlap_threshold", "relative_overlap_threshold"):
            if key in diagnostics:
                diagnostics[key] = self._as_num(diagnostics, key, "diagnostics")

        # --- Output ---
        output_raw = raw.get("output", {}) or {}
        self._check_keys(output_raw, "output", "output")
        metadata["output_file"] = output_raw.get("file", "output.xml")
        metadata["embed_fingerprint"] = self._as_bool(
            output_raw, "embed_fingerprint", True, "output"
        )

        if errors:
            raise ConfigError("; ".join(errors))

        config = BEASTConfig(
            metadata=metadata,
            partitions=partitions,
            all_alignments=list(alignments.values()),
            tree_prior_type=tree_prior_type,
            tree_prior_params=tree_prior_params,
            calibrations=calibrations,
            calibration_method=calibration_method,
            parameter_priors=parameter_priors,
            mcmc=mcmc,
            loggers=loggers,
            tip_dates=tip_dates,
            initialization=initialization,
            operator_weights=operator_weights,
            diagnostics=diagnostics,
        )
        self._validate_calibration_provenance(config)
        return config

    # -- section helpers ----------------------------------------------------

    def _validate_tree_ref(self, tree_ref: Any, partition_id: str, known_ids: Set[str]) -> None:
        """Reject tree settings the generator cannot honour.

        ``tree: separate`` (and any partition-specific tree) was
        parsed, stored and then never read, so multi-locus analyses wanting
        *unlinked* topologies silently produced a fully linked single-tree
        analysis. Rather than keep that fail-open middle ground, the setting is
        now refused at parse time until per-partition trees are implemented.

        Args:
            tree_ref: Value of the partition's ``tree`` key.
            partition_id: Partition being validated.
            known_ids: Partition ids seen so far.

        Raises:
            ConfigError: If the requested tree is not the shared tree.
        """
        if tree_ref == "shared":
            return
        raise ConfigError(
            f"Partition '{partition_id}': tree: {tree_ref!r} requests a tree other than the "
            f"shared one, but Beast2Py assembles every partition on a single linked topology "
            f"(per-partition / unlinked trees are not implemented). Silently falling back to "
            f"the shared tree would turn an intended *BEAST-style multi-locus analysis into "
            f"a fully linked one, so this is an error. Remove the 'tree' key to analyse the "
            f"shared topology, or run one analysis per independent partition."
        )

    def _resolve_effective_sites(self, alignment, alignments):
        """Map an alignment to (source key, root-relative 1-based site set).

        Walks the ``data_id`` chain to the file-backed root, then replays every
        filter along the way in the coordinate space it actually applies to, so
        chained filters are compared against the sites they really keep.

        Args:
            alignment: The alignment a partition refers to.
            alignments: All parsed alignments, keyed by id.

        Returns:
            Tuple of a source key identifying the underlying sequence data and
            the set of *root-relative* site positions it contributes.
        """
        chain = []
        cur = alignment
        seen = set()
        while True:
            chain.append(cur)
            parent_id = getattr(cur, "data_id", None)
            if not parent_id or parent_id in seen or parent_id not in alignments:
                break
            seen.add(cur.id)
            cur = alignments[parent_id]
        chain.reverse()
        root = chain[0]

        sites = list(range(1, root.n_sites + 1))
        for stage in chain[1:]:
            expr = getattr(stage, "filter", None)
            if not expr:
                continue
            picked = ConfigParser.resolve_filter_sites(expr, len(sites), stage.id)
            sites = sorted(sites[i - 1] for i in picked if i - 1 < len(sites))
        source_key = getattr(root, "source_file", None) or root.id
        return source_key, set(sites)

    def _validate_partition_site_coverage(self, partitions, alignments) -> None:
        """Reject partitions that feed the same site into the likelihood twice.

        a comment here claimed this was already rejected, but no
        code implemented it, so copying one 898-site alignment into two
        partitions produced a working XML whose posterior is the true posterior
        squared -- over-confident, not sharper.

        Args:
            partitions: All parsed partitions.
            alignments: All parsed alignments.

        Raises:
            ConfigError: If two partitions cover a common site.
        """
        effective = {}
        for p in partitions:
            aln = alignments.get(p.id)
            if aln is None:
                continue
            effective[p.id] = self._resolve_effective_sites(aln, alignments)

        ids = sorted(effective)
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                a, b = ids[i], ids[j]
                (src_a, sites_a), (src_b, sites_b) = effective[a], effective[b]
                if src_a != src_b:
                    continue
                shared = sites_a & sites_b
                if shared:
                    sample = ", ".join(str(s) for s in sorted(shared)[:5])
                    raise ConfigError(
                        f"Partitions '{a}' and '{b}' both draw on the same sequence "
                        f"data ({src_a}) and share {len(shared)} site(s) (e.g. {sample}). "
                        f"Every site may enter the likelihood exactly once: duplicated "
                        f"or overlapping partitions raise the posterior to a power rather "
                        f"than sharpening it. Carve the source into disjoint subsets with "
                        f"'filter' (e.g. codon positions '1::3,2::3' and '3::3')."
                    )

    def _validate_linked_model_consistency(self, partitions, partitions_raw, errors) -> None:
        """Reject a ``linked_to`` that also declares a *different* model.

        When a partition links its site or clock model to another partition the
        generator skips building its own model entirely, so anything declared
        alongside the link was discarded without a word: ``gtr`` plus
        ``linked_to: g1`` quietly produced a second copy of g1's HKY, and a
        relaxed clock linked to a strict one quietly became strict. The analysis
        that ran was not the one written.

        Only fields the author actually wrote are compared: the dataclasses fill
        in defaults (gamma_categories 0, clock type strict), and comparing those
        against the target would reject the documented, correct "link and
        declare nothing else" style used by the shipped examples.

        Args:
            partitions: All parsed partitions.
            partitions_raw: The raw ``partitions`` list, for declared-key checks.
            errors: Accumulator list, extended in place.
        """
        by_id = {p.id: p for p in partitions}
        raw_by_id = {
            r["id"]: r for r in partitions_raw if isinstance(r, dict) and "id" in r
        }

        for p in partitions:
            raw = raw_by_id.get(p.id, {})
            raw_sm = raw.get("site_model", {}) or {}
            raw_cm = raw.get("clock_model", {}) or {}

            target_id = p.site_model.linked_to
            if target_id and target_id in by_id:
                target = by_id[target_id].site_model
                declared_sm = raw_sm.get("substitution_model")
                if isinstance(declared_sm, dict) and "type" in declared_sm:
                    mine = ConfigParser._norm_name(declared_sm.get("type"))
                    theirs = ConfigParser._norm_name(
                        (target.substitution_model or {}).get("type")
                    )
                    if mine and theirs and mine != theirs:
                        errors.append(
                            f"Partition '{p.id}': site_model.linked_to '{target_id}' "
                            f"replaces this partition's whole site model, but the two "
                            f"declare different substitution models ({mine} vs "
                            f"{theirs}). Linking silently discards {mine}; either drop "
                            f"'substitution_model' here or stop linking."
                        )
                if "gamma_categories" in raw_sm:
                    if p.site_model.gamma_categories != target.gamma_categories:
                        errors.append(
                            f"Partition '{p.id}': site_model.linked_to '{target_id}' "
                            f"ignores this partition's gamma_categories="
                            f"{p.site_model.gamma_categories}; the linked target uses "
                            f"{target.gamma_categories}. Make them agree or stop linking."
                        )
                if "proportion_invariant" in raw_sm:
                    if bool(p.site_model.proportion_invariant) != bool(
                        target.proportion_invariant
                    ):
                        missing = (
                            "does not " if not target.proportion_invariant else ""
                        )
                        errors.append(
                            f"Partition '{p.id}': site_model.linked_to '{target_id}' "
                            f"ignores this partition's proportion_invariant setting; "
                            f"the linked target {missing}use one. Make them agree or "
                            f"stop linking."
                        )

            clock_target_id = p.clock_model.linked_to
            if clock_target_id and clock_target_id in by_id:
                if "type" in raw_cm:
                    mine = getattr(p.clock_model.type, "value", p.clock_model.type)
                    target_clock = by_id[clock_target_id].clock_model
                    theirs = getattr(target_clock.type, "value", target_clock.type)
                    if mine and theirs and mine != theirs:
                        errors.append(
                            f"Partition '{p.id}': clock_model.linked_to "
                            f"'{clock_target_id}' reuses that partition's clock, but "
                            f"this partition declares type '{mine}' while "
                            f"'{clock_target_id}' uses '{theirs}'. Linking discards "
                            f"'{mine}' and runs '{theirs}' instead."
                        )

    @staticmethod
    def parse_filter_terms(expr: str, n_sites: int, aln_id: str) -> List[range]:
        """Parse a BEAST2 filter expression into 1-based inclusive ranges.

        Grammar (as implemented by FilteredAlignment.parseFilterSpec): a
        comma-separated list of singletons ``N``, ranges ``from-to`` with an
        optional ``\\step`` suffix, or iterators ``from:to:step`` where an empty
        ``to`` means "last site".

        Args:
            expr: The BEAST2 filter expression.
            n_sites: Site count the expression is validated against.
            aln_id: Alignment id used in error messages.

        Returns:
            List of ``range`` objects (upper bound exclusive) over 1-based
            site positions.

        Raises:
            ConfigError: If a term is malformed or addresses absent sites.
        """
        spans: List[range] = []
        for raw_token in expr.split(","):
            token = raw_token.strip()
            if not token:
                continue
            step = 1
            try:
                if "\\" in token:
                    head, tail = token.split("\\", 1)
                    step = int(tail)
                else:
                    head = token
                if "-" in head:
                    parts = head.split("-")
                    start = int(parts[0]) if parts[0] else 1
                    end = int(parts[1]) if len(parts) > 1 and parts[1] else n_sites
                elif head.count(":") >= 2:
                    parts = head.split(":")
                    start = int(parts[0]) if parts[0] else 1
                    end = int(parts[1]) if parts[1] else n_sites
                    step = int(parts[2]) if len(parts) > 2 and parts[2] else 1
                elif head.isdigit():
                    start = end = int(head)
                else:
                    raise ValueError
            except ValueError:
                raise ConfigError(
                    f"Alignment '{aln_id}': unparsable filter term '{raw_token}' in '{expr}'. "
                    f"Expected a site number, a 'from-to[\\step]' range, or a "
                    f"'from:to:step' iterator (empty parts take their defaults)."
                )
            if step < 1:
                raise ConfigError(f"Alignment '{aln_id}': filter step must be >= 1: '{expr}'")
            if start < 1 or end > n_sites or start > end:
                raise ConfigError(
                    f"Alignment '{aln_id}': filter '{expr}' addresses sites {start}-{end} "
                    f"outside a {n_sites}-site alignment"
                )
            spans.append(range(start, end + 1, step))
        return spans

    @classmethod
    def resolve_filter_sites(cls, expr: str, n_sites: int, aln_id: str) -> Set[int]:
        """Materialise the set of 1-based site positions a filter selects.

        Used both to validate an expression and, crucially, to compare the
        coverage of two partitions over a shared source alignment.
        """
        sites: Set[int] = set()
        for span in cls.parse_filter_terms(expr, n_sites, aln_id):
            sites.update(span)
        return sites

    def _validate_filter_sites(self, alignment, expr: str, aln_id: str) -> None:
        """Check a FilteredAlignment expression selects at least one site.

        Args:
            alignment: The filtered alignment (sequences come from its source).
            expr: The BEAST2 filter expression.
            aln_id: Alignment id for error messages.

        Raises:
            ConfigError: If the expression is malformed or selects no sites.
        """
        if not self.resolve_filter_sites(expr, alignment.n_sites, aln_id):
            raise ConfigError(
                f"Alignment '{aln_id}': filter '{expr}' selects no sites from a "
                f"{alignment.n_sites}-site alignment"
            )

    def _validate_tree_prior_params(self, tp_type: TreePriorType, params: Dict[str, Any]) -> None:
        """Range-check tree-prior parameters against their definitions.

        a birth-death prior with death_rate > birth_rate produced a
        ``birthDiffRate`` state node below its own lower bound and a
        ``relativeDeathRate`` above its upper bound, which BEAST2 rejects at
        run time while the generator exited 0.

        Args:
            tp_type: The tree prior type.
            params: Raw (unconverted) tree prior parameters.

        Raises:
            ConfigError: If a parameter is outside its defined domain.
        """

        def value_of(key: str, default: float) -> float:
            raw = params.get(key, default)
            if isinstance(raw, dict):
                return ConfigParser._num_or_default(
                    raw, "value", f"tree_prior.{key}", default
                )
            try:
                return as_number(raw, f"tree_prior.{key}")
            except ValueError as exc:
                raise ConfigError(str(exc)) from exc

        if tp_type == TreePriorType.BIRTH_DEATH:
            birth = value_of("birth_rate", 1.0)
            death = value_of("death_rate", 0.5)
            if birth <= 0:
                raise ConfigError(
                    f"tree_prior.birth_rate: {birth} must be > 0 for a birth-death prior"
                )
            if death < 0:
                raise ConfigError(
                    f"tree_prior.death_rate: {death} must be >= 0 for a birth-death prior"
                )
            if death >= birth:
                raise ConfigError(
                    f"tree_prior: death_rate ({death}) must be smaller than birth_rate "
                    f"({birth}) for the BirthDeathGernhard08Model, whose parameters are "
                    f"birth - death (>= 0) and death / birth (< 1). A super-critical "
                    f"process (death >= birth) needs the BDSKY add-on "
                    f"(tree_prior.type: bd_skyline_serial), which parameterises "
                    f"birth, death and sampling rates separately."
                )
            sampling = value_of("sample_probability", value_of("sampling_rate", 1.0))
            if not 0.0 < sampling <= 1.0:
                raise ConfigError(
                    f"tree_prior.sample_probability: {sampling} must lie in (0, 1] "
                    f"(it is a sampling probability)"
                )
        if tp_type == TreePriorType.BD_SKYLINE_SERIAL:
            # (lower bound, upper bound, is the lower bound exclusive)
            for key, lo, hi, excl in (
                # A zero speciation rate makes the diversification process
                # degenerate, so birth_rate excludes 0. The others legitimately
                # take 0: death_rate 0 is pure birth, and BDSKY can parameterise
                # sampling through rho *or* sampling_rate, leaving the other at 0.
                ("birth_rate", 0.0, None, True),
                ("death_rate", 0.0, None, False),
                ("sampling_rate", 0.0, 1.0, False),
                ("rho", 0.0, 1.0, False),
            ):
                if key not in params:
                    continue
                val = value_of(key, 0.1)
                below = val <= lo if excl else val < lo
                if below or (hi is not None and val > hi):
                    bound = "inf" if hi is None else hi
                    interval = f"({lo}, {bound}]" if excl else f"[{lo}, {bound}]"
                    raise ConfigError(f"tree_prior.{key}: {val} must lie in {interval}")
        if tp_type == TreePriorType.COALESCENT_CONSTANT:
            if value_of("pop_size", 1.0) <= 0:
                raise ConfigError("tree_prior.pop_size must be > 0")
        if tp_type == TreePriorType.COALESCENT_EXPONENTIAL:
            if value_of("pop_size", 1.0) <= 0:
                raise ConfigError("tree_prior.pop_size must be > 0")
        if tp_type in (TreePriorType.YULE, TreePriorType.CALIBRATED_YULE):
            if value_of("birth_rate", 1.0) <= 0:
                raise ConfigError("tree_prior.birth_rate must be > 0")

    def _parse_tip_dates(
        self, td_raw: Dict[str, Any], partitions: List[Partition]
    ) -> TipDatesConfig:
        """Parse and validate a ``tip_dates`` block.

        Missing taxa are rejected rather than silently treated as
        contemporaneous samples, and an unrecognised ``trait_name`` is an
        error instead of producing valid-but-not-used XML.

        Args:
            td_raw: The raw tip_dates mapping.
            partitions: Parsed partitions, used to enumerate tip names.

        Returns:
            A validated TipDatesConfig.

        Raises:
            ConfigError: On an unknown trait name, unknown units, unknown or
                missing tip names, or non-numeric dates.
        """
        trait_name = td_raw.get("trait_name", "date-forward")
        if trait_name not in VALID_TRAIT_NAMES:
            raise ConfigError(
                f"tip_dates.trait_name: {trait_name!r} is not a name BEAST2 reads as a "
                f"date trait (valid: {', '.join(sorted(VALID_TRAIT_NAMES))}). Any other "
                f"trait is carried through the analysis without ever dating the tips."
            )
        units = td_raw.get("units", "year")
        if units not in VALID_DATE_UNITS:
            raise ConfigError(
                f"tip_dates.units: {units!r} invalid (valid: {', '.join(sorted(VALID_DATE_UNITS))})"
            )

        dates: Dict[str, float] = {}
        raw_dates = td_raw.get("dates", {}) or {}
        if not isinstance(raw_dates, dict):
            raise ConfigError("tip_dates.dates must be a mapping of taxon name to date")
        for taxon, value in raw_dates.items():
            dates[str(taxon)] = self._as_num({"d": value}, "d", f"tip_dates.dates.{taxon}")

        csv_file = td_raw.get("file")
        if csv_file:
            path = Path(csv_file)
            if not path.is_absolute():
                path = self.base_dir / path
            if not path.exists():
                raise ConfigError(f"tip_dates.file not found: {path}")
            dates.update(self._read_tip_date_csv(path))

        if not dates:
            raise ConfigError("tip_dates is enabled but no dates were supplied")

        all_taxa: Set[str] = set()
        for part in partitions:
            all_taxa.update(part.alignment.taxa_names)
        unknown = sorted(set(dates) - all_taxa)
        if unknown:
            raise ConfigError(
                f"tip_dates.dates names taxa absent from the alignments: {', '.join(unknown)}"
            )
        missing = sorted(all_taxa - set(dates))
        if missing:
            raise ConfigError(
                f"tip_dates.dates covers {len(dates)} of {len(all_taxa)} taxa; these are "
                f"missing: {', '.join(missing)}. BEAST2 treats an undated tip as sampled at "
                f"time zero, which silently turns those lineages into contemporaneous "
                f"samples. Supply a date for every tip, or disable tip_dates."
            )
        return TipDatesConfig(
            enabled=True, trait_name=trait_name, units=units, dates=dates, file=csv_file
        )

    @staticmethod
    def _read_tip_date_csv(path: Path) -> Dict[str, float]:
        """Read a two-column taxon,date CSV of the BEAUti TipDateTree shape.

        Args:
            path: CSV path; a header line is skipped when non-numeric.

        Returns:
            Mapping of taxon name to numeric date.

        Raises:
            ConfigError: If a data row is unparsable.
        """
        import csv

        dates: Dict[str, float] = {}
        with open(path, newline="", encoding="utf-8") as handle:
            for row_index, row in enumerate(csv.reader(handle)):
                if not row:
                    continue
                if len(row) < 2:
                    raise ConfigError(f"{path}: row {row_index + 1} needs at least 2 columns")
                try:
                    dates[row[0].strip()] = float(row[1])
                except ValueError:
                    if row_index == 0:
                        continue  # header
                    raise ConfigError(
                        f"{path}: row {row_index + 1}: {row[1]!r} is not a numeric date"
                    )
        return dates

    def _validate_calibration_provenance(self, config: BEASTConfig) -> None:
        """Cross-check recorded fossil ages against the prior that is used.

        ``provenance.original_age_min/max`` were parsed and stored
        but never compared with the distribution's support, so a calibration
        recorded as ">= 5.5 Ma" could carry a prior with ~98% of its mass below
        the fossil's age without any complaint.

        Args:
            config: The assembled configuration.

        Raises:
            ConfigError: If the prior's 95% interval lies wholly outside the
                recorded fossil range.
        """
        from .utils import compute_distribution_stats

        for cal in config.calibrations:
            prov = cal.provenance
            if prov is None or cal.distribution is None:
                continue
            if prov.original_age_min is None and prov.original_age_max is None:
                continue
            # A *soft* bound is allowed to put most of its mass outside the
            # fossil range: that is precisely the BEAST2 convention (e.g. a
            # lognormal whose 2.5% quantile sits at the hard minimum). Only a
            # calibration declared "hard" must actually contain the evidence it
            # cites, so soft cases are reported as warnings, not errors.
            hard = str(prov.calibration_type).lower() == "hard"
            try:
                stats = compute_distribution_stats(
                    cal.distribution.type, cal.distribution.parameters, cal.distribution.offset
                )
            except ValueError:
                continue  # improper distributions (OneOnX) are already excluded
            lo, hi = stats["hpd_lower"], stats["hpd_upper"]
            problem = None
            if prov.original_age_min is not None and hi < float(prov.original_age_min):
                problem = (
                    f"provenance records an original minimum age of "
                    f"{prov.original_age_min} but the 95% interval of the "
                    f"{cal.distribution.type} prior is [{lo:.3g}, {hi:.3g}], i.e. entirely "
                    f"younger than the fossil"
                )
            elif prov.original_age_max is not None and lo > float(prov.original_age_max):
                problem = (
                    f"provenance records an original maximum age of "
                    f"{prov.original_age_max} but the 95% interval of the "
                    f"{cal.distribution.type} prior is [{lo:.3g}, {hi:.3g}], i.e. entirely "
                    f"older than the fossil"
                )
            if problem:
                message = f"Calibration '{cal.name}': {problem}."
                if hard:
                    raise ConfigError(
                        message + " A hard-bound calibration must contain the age it claims to "
                        "enforce; either change calibration_type to soft (with a "
                        "soft-bound distribution) or move the prior."
                    )
                warnings.warn(
                    message + " Acceptable for a soft bound (the fossil minimum is imposed "
                    "softly), but the recorded range and the prior should be "
                    "explained together in the methods text.",
                    stacklevel=2,
                )

    def _parse_site_model(
        self,
        sm_raw: Dict[str, Any],
        errors: List[str],
        partition_id: str,
        alignment: Any = None,
    ) -> SiteModelConfig:
        """Parse site model configuration.

        Args:
            sm_raw: Raw ``site_model`` mapping.
            errors: Accumulator for recoverable messages.
            partition_id: Partition being parsed.
            alignment: Parsed alignment, used for the data-type / model-family
                cross-check.

        Returns:
            A SiteModelConfig.

        Raises:
            ConfigError: On unknown keys, an unknown model, or an inconsistent
                data type / alphabet / model family triple.
        """
        self._check_keys(sm_raw, "site_model", f"Partition '{partition_id}': site_model")
        subst_raw = sm_raw.get("substitution_model", {}) or {}
        if not isinstance(subst_raw, dict):
            raise ConfigError(f"Partition '{partition_id}': substitution_model must be a mapping")

        # Validate substitution model type
        sm_type = self._canonical_name(
            subst_raw.get("type", "hky"), self.VALID_SUBSTITUTION_MODELS,
            "substitution model", f"Partition '{partition_id}'")

        # Validate substitution-model keys against the registry's parameter list
        from .registry import ModelRegistry

        model_spec = ModelRegistry.get_substitution_model(sm_type) or {}
        allowed = {"type", "rates", "frequencies"} | set(model_spec.get("params", {}))
        unknown = sorted(set(subst_raw) - allowed)
        if unknown:
            raise ConfigError(
                f"Partition '{partition_id}': substitution_model '{sm_type}' has unknown "
                f"key(s) {', '.join(unknown)} (allowed: {', '.join(sorted(allowed))})"
            )
        for key, val in subst_raw.items():
            if isinstance(val, dict):
                where_param = f"Partition '{partition_id}': substitution_model.{key}"
                self._check_keys(val, "parameter", where_param)
                # Only the key names were checked, so `kappa: {value: .inf}`
                # survived parsing and blew up during XML assembly as a raw
                # "cannot convert float infinity to integer" traceback.
                if "value" in val and val["value"] is not None and not isinstance(
                    val["value"], str
                ):
                    try:
                        num = as_number(val["value"], f"{where_param}.value")
                    except ValueError as exc:
                        raise ConfigError(str(exc)) from exc
                    self._finite(where_param, num)
                    self._check_bound(
                        where_param, num, val.get("lower"), val.get("upper")
                    )
        # A rates vector is only meaningful where all six exchangeabilities are
        # free; catching it here turns "the generator raised later, or worse
        # ignored it" into a parse-time message naming the model.
        if "rates" in subst_raw and sm_type not in ("gtr", "sym"):
            raise ConfigError(
                f"Partition '{partition_id}': substitution_model '{sm_type}' has no free "
                f"rateAC..rateGT vector; only gtr and sym accept 'rates'. Use the model's "
                f"own named rate parameters (kappa for hky, kappa1/kappa2 for tn93, "
                f"rateAG/rateCT/rateTransversions1-2 for tim, rateAC/AT/CG/GT/"
                f"rateTransitions for tvm)."
            )

        # --- Three-way cross-check: declared data type x model family x the
        # alphabet actually observed in the alignment ---
        if model_spec.get("is_aminoacid") and (
            alignment is None or alignment.data_type.value != "aminoacid"
        ):
            raise ConfigError(
                f"Partition '{partition_id}': model '{sm_type}' is an amino-acid "
                f"substitution model (a 20 x 20 exchangeability matrix) and is undefined on "
                f"{'' if alignment is None else alignment.data_type.value} data. Set the "
                f"alignment's data_type to aminoacid or choose a nucleotide model."
            )
        if not model_spec.get("is_aminoacid") and alignment is not None:
            if alignment.data_type.value == "aminoacid":
                raise ConfigError(
                    f"Partition '{partition_id}': model '{sm_type}' is a nucleotide model "
                    f"(its parameters such as kappa are transition/transversion ratios, which "
                    f"have no meaning for 20-state amino-acid data) but the alignment is "
                    f"declared aminoacid. Use WAG/JTT/Dayhoff/Blosum62/CPREV/MTREV."
                )
        if alignment is not None and alignment.data_type.value == "aminoacid":
            observed = set(alignment.observed_characters or "")
            nucleotide_only = {"A", "C", "G", "T", "U"}
            if observed and observed <= nucleotide_only:
                raise ConfigError(
                    f"Partition '{partition_id}': the alignment is declared aminoacid but only "
                    f"contains the characters {''.join(sorted(observed))}, i.e. a nucleotide "
                    f"alphabet. Either data_type is wrong or the wrong file was referenced."
                )

        # Parse gamma categories
        gamma_cats = self._as_num(sm_raw, "gamma_categories", f"'{partition_id}'.site_model", 0)
        if gamma_cats != int(gamma_cats) or int(gamma_cats) < 0 or int(gamma_cats) > 100:
            raise ConfigError(
                f"Partition '{partition_id}': gamma_categories must be a non-negative integer, "
                f"got {sm_raw.get('gamma_categories')!r}"
            )
        gamma_cats = int(gamma_cats)
        if "gamma_shape" in sm_raw and gamma_cats <= 1:
            raise ConfigError(
                f"Partition '{partition_id}': gamma_shape was supplied but "
                f"gamma_categories is {gamma_cats}. SiteModel.java documents the shape input "
                f"as 'Ignored if gammaCategoryCount 1 or less', so an estimated shape would "
                f"be a pure prior random walk that contributes nothing to the likelihood "
                f"while still appearing as a column in the trace."
            )

        # Parse gamma shape
        gamma_shape = None
        if gamma_cats > 0:
            gs_raw = sm_raw.get("gamma_shape")
            if gs_raw:
                gamma_shape = self._parse_real_parameter(gs_raw, f"'{partition_id}'.gamma_shape")
            else:
                gamma_shape = RealParameter(value=0.5, lower=0.0, upper=None)

        # Parse proportion invariant
        prop_inv = None
        pi_raw = sm_raw.get("proportion_invariant")
        if pi_raw:
            prop_inv = self._parse_real_parameter(pi_raw, f"'{partition_id}'.proportion_invariant")

        return SiteModelConfig(
            substitution_model=subst_raw,
            gamma_categories=gamma_cats,
            gamma_shape=gamma_shape,
            proportion_invariant=prop_inv,
            linked_to=sm_raw.get("linked_to"),
        )

    def _parse_clock_model(
        self,
        cm_raw: Dict[str, Any],
        errors: List[str],
        partition_id: str,
    ) -> ClockModelConfig:
        """Parse clock model configuration."""
        self._check_keys(cm_raw, "clock_model", f"Partition '{partition_id}': clock_model")
        cm_type_str = self._canonical_name(
            cm_raw.get("type", "strict"), self.VALID_CLOCK_MODELS,
            "clock model", f"Partition '{partition_id}'")

        prefix = f"Partition '{partition_id}': clock_model"
        clock_rate = None
        if "clock_rate" in cm_raw:
            clock_rate = self._parse_real_parameter(cm_raw["clock_rate"], f"{prefix}.clock_rate")

        ucld_mean = None
        ucld_stdev = None
        if cm_type_str in ("ucln", "uce"):
            if "ucld_mean" in cm_raw:
                ucld_mean = self._parse_real_parameter(cm_raw["ucld_mean"], f"{prefix}.ucld_mean")
            if "ucld_stdev" in cm_raw:
                ucld_stdev = self._parse_real_parameter(
                    cm_raw["ucld_stdev"], f"{prefix}.ucld_stdev"
                )
        else:
            for stray in ("ucld_mean", "ucld_stdev"):
                if stray in cm_raw:
                    raise ConfigError(
                        f"{prefix}.{stray} only applies to the ucln/uce clocks; the "
                        f"{cm_type_str} clock has no such parameter"
                    )

        linked_to = cm_raw.get("linked_to")

        return ClockModelConfig(
            type=ClockModelType(cm_type_str),
            clock_rate=clock_rate,
            ucld_mean=ucld_mean,
            ucld_stdev=ucld_stdev,
            linked_to=linked_to,
        )

    def _parse_calibration(
        self,
        cal_raw: Dict[str, Any],
        errors: List[str],
        all_taxa: set,
        where: str,
        check_taxa: bool = True,
    ) -> Optional[CalibrationPoint]:
        """Parse a calibration point.

        Args:
            cal_raw: Raw calibration mapping.
            errors: Accumulator for recoverable messages.
            all_taxa: Taxon names present across the alignments.
            where: Location used in error messages.
            check_taxa: Whether to verify taxon membership. Standalone
                calibration files are parsed before any alignment is known, so
                they pass ``False`` and get the check later.

        Returns:
            A validated CalibrationPoint.

        Raises:
            ConfigError: On any structural, boolean, taxon or distribution
                problem, including an inverted ``uniform`` bound.
        """
        if not isinstance(cal_raw, dict):
            raise ConfigError(f"{where}: must be a mapping")
        self._check_keys(cal_raw, "calibration", where)
        name = cal_raw.get("name", "")
        if not name:
            raise ConfigError(f"{where}: missing 'name'")

        taxa = cal_raw.get("taxa")
        if taxa is not None and check_taxa:
            unknown = [t for t in taxa if t not in all_taxa]
            if unknown:
                raise ConfigError(
                    f"Calibration '{name}': taxon(s) {', '.join(map(str, unknown))} not found "
                    f"in any alignment"
                )

        monophyletic = self._as_bool(cal_raw, "monophyletic", True, where)
        use_originate = self._as_bool(cal_raw, "use_originate", False, where)
        tipsonly = self._as_bool(cal_raw, "tipsonly", False, where)

        # Parse distribution
        dist = None
        dist_raw = cal_raw.get("distribution")
        hp_raw: Dict[str, Any] = {}
        if dist_raw:
            if not isinstance(dist_raw, dict):
                raise ConfigError(f"Calibration '{name}': distribution must be a mapping")
            self._check_keys(dist_raw, "distribution", f"Calibration '{name}': distribution")
            dist_type = self._canonical_name(
                dist_raw.get("type", ""), self.VALID_DISTRIBUTIONS,
                "distribution type", f"Calibration '{name}'")
            parameters = dict(dist_raw.get("parameters", {}) or {})
            offset = self._as_num(dist_raw, "offset", f"Calibration '{name}'", 0.0) or 0.0
            dist = DistributionConfig(type=dist_type, parameters=parameters, offset=float(offset))
            self._validate_distribution(dist, name)
            hp_raw = dist_raw.get("hyperpriors", {}) or {}

        # Parse provenance
        provenance = None
        prov_raw = cal_raw.get("provenance")
        if prov_raw:
            if not isinstance(prov_raw, dict):
                raise ConfigError(f"Calibration '{name}': provenance must be a mapping")
            self._check_keys(prov_raw, "provenance", f"Calibration '{name}': provenance")
            provenance = Provenance(
                source=prov_raw.get("source", ""),
                reference=prov_raw.get("reference", ""),
                calibration_type=prov_raw.get("calibration_type", "soft"),
                original_age_min=self._as_num(
                    prov_raw, "original_age_min", f"Calibration '{name}': provenance"
                ),
                original_age_max=self._as_num(
                    prov_raw, "original_age_max", f"Calibration '{name}': provenance"
                ),
                notes=str(prov_raw.get("notes", "") or ""),
            )
            if provenance.original_age_min is not None and provenance.original_age_max is not None:
                if provenance.original_age_min > provenance.original_age_max:
                    raise ConfigError(
                        f"Calibration '{name}': provenance.original_age_min "
                        f"({provenance.original_age_min}) exceeds original_age_max "
                        f"({provenance.original_age_max})"
                    )

        # Parse hyperpriors
        hyperpriors: Dict[str, HyperpriorConfig] = {}
        for param_name, hp_config in hp_raw.items():
            where_hp = f"Calibration '{name}': hyperpriors.{param_name}"
            if not isinstance(hp_config, dict):
                raise ConfigError(f"{where_hp}: must be a mapping")
            self._check_keys(hp_config, "hyperprior", where_hp)
            hp_type = self._canonical_name(
                hp_config.get("type", "uniform"), self.VALID_DISTRIBUTIONS,
                "distribution type", where_hp)
            if param_name not in (dist.parameters if dist else {}):
                raise ConfigError(
                    f"{where_hp}: '{param_name}' is not a parameter of the "
                    f"{dist.type if dist else 'absent'} distribution"
                )
            hp_dist = DistributionConfig(
                type=hp_type,
                parameters=dict(hp_config.get("parameters", {}) or {}),
                offset=float(self._as_num(hp_config, "offset", where_hp, 0.0) or 0.0),
            )
            self._validate_distribution(hp_dist, f"{name} hyperprior {param_name}")
            hyperpriors[param_name] = HyperpriorConfig(
                type=hp_type, parameters=hp_dist.parameters, offset=hp_dist.offset
            )

        return CalibrationPoint(
            name=name,
            taxa=taxa,
            monophyletic=monophyletic,
            distribution=dist,
            use_originate=use_originate,
            tipsonly=tipsonly,
            provenance=provenance,
            hyperpriors=hyperpriors,
        )

    def _validate_distribution(self, dist: DistributionConfig, where: str) -> None:
        """Range-check a parametric distribution.

        ``uniform {lower: 20, upper: 5}`` produced a Uniform whose
        density is undefined (zero) everywhere, with exit code 0.

        Args:
            dist: The parsed distribution configuration.
            where: Location used in error messages.

        Raises:
            ConfigError: If a parameter is outside its admissible domain.
        """
        params = dist.parameters

        # ---- parameter key closure  ---------------------------
        # Unknown and missing keys used to pass straight through: `Sdeg: 0.4`
        # emitted a LogNormal carrying only M, BEAST2 then sampled its own
        # default S=1.0, and every gate -- BEAST2 included -- stayed green on
        # a different prior than the one the user wrote.
        from .registry import ModelRegistry

        registry_spec = ModelRegistry.get_distribution_spec(dist.type) or {}
        canonical = set(registry_spec.get("params", {})) | set(
            registry_spec.get("bool_attrs", {})
        )
        allowed = canonical | self.DISTRIBUTION_EXTRA_KEYS.get(dist.type, set())
        unknown = set(params) - allowed
        if unknown:
            raise ConfigError(
                f"Calibration '{where}': {dist.type} accepts "
                f"{sorted(canonical) or 'no parameters'}, but got unrecognised "
                f"key(s) {sorted(unknown)}. An unknown key is dropped without a "
                f"word and BEAST2 substitutes its own default, so the prior that "
                f"gets sampled is not the prior you wrote."
            )
        missing = self.DISTRIBUTION_REQUIRED.get(dist.type, set()) - set(params)
        if missing:
            raise ConfigError(
                f"Calibration '{where}': {dist.type} is missing required "
                f"parameter(s) {sorted(missing)}; BEAST2 would silently fall back "
                f"to its own defaults and sample a different prior."
            )

        spec = {
            "normal": {"mean": None, "sigma": (0.0, None)},
            "lognormal": {"S": (0.0, None)},  # M is a log-scale mean: any real value
            "uniform": {"lower": None, "upper": None},
            "exponential": {"mean": (0.0, None)},
            "gamma": {"alpha": (0.0, None), "beta": (0.0, None)},
            "beta": {"alpha": (0.0, None), "beta": (0.0, None)},
            "laplace": {"mu": None, "scale": (0.0, None)},
            "inverse_gamma": {"alpha": (0.0, None), "beta": (0.0, None)},
            "poisson": {"lambda": (0.0, None)},
            "chi_square": {"dof": (0.0, None)},
            "one_on_x": {},
        }[dist.type]

        for key, range_spec in spec.items():
            if key not in params:
                continue
            try:
                val = as_number(params[key], f"Calibration '{where}': {key}")
            except ValueError as exc:
                raise ConfigError(str(exc)) from exc
            if range_spec is None:
                continue
            lo, hi = range_spec
            if lo is not None and val <= lo:
                raise ConfigError(
                    f"Calibration '{where}': {dist.type}.{key} = {val} must be > {lo}"
                )
            if hi is not None and val >= hi:
                raise ConfigError(
                    f"Calibration '{where}': {dist.type}.{key} = {val} must be < {hi}"
                )
        if dist.type == "uniform":
            lower = params.get("lower", params.get("min"))
            upper = params.get("upper", params.get("max"))
            if lower is None or upper is None:
                raise ConfigError(
                    f"Calibration '{where}': uniform requires both 'lower' and 'upper'"
                )
            lo = as_number(lower, f"'{where}'.uniform.lower")
            hi = as_number(upper, f"'{where}'.uniform.upper")
            if lo + dist.offset >= hi + dist.offset:
                raise ConfigError(
                    f"Calibration '{where}': uniform lower ({lo}) must be strictly below "
                    f"upper ({hi}); as given, the density is zero on the whole support."
                )
        if dist.type == "lognormal":
            for flag in ("mean_in_real_space",):
                if flag in params:
                    try:
                        params[flag] = as_strict_bool(params[flag], f"'{where}'.{flag}")
                    except ValueError as exc:
                        raise ConfigError(str(exc)) from exc
            # M is the *log*-scale mean (median = exp(M)) unless
            # meanInRealSpace is set, in which case BEAST2 rescales it and M must
            # be a positive real-space mean.
            if params.get("mean_in_real_space") and "M" in params:
                if as_number(params["M"], f"'{where}'.lognormal.M") <= 0:
                    raise ConfigError(
                        f"Calibration '{where}': lognormal M must be > 0 when "
                        f"mean_in_real_space is true (M is then the real-space mean). "
                        f"With mean_in_real_space false, M is the log-scale mean and any "
                        f"real value is valid."
                    )

    @staticmethod
    def _num_or_default(raw: Dict[str, Any], key: str, where: str, default: float) -> float:
        """Read a number, substituting *default* only when the key is absent.

        ``x or default`` also swallows 0 and None, which silently overrode an
        author's explicit value and made the range guards that follow it
        unreachable (1).

        Args:
            raw: The mapping to read.
            key: The key to read.
            where: Location used in error messages.
            default: Value used when ``key`` is not present at all.

        Returns:
            The parsed number, or ``default``.

        Raises:
            ConfigError: If the value is not a finite number.
        """
        if key not in raw or raw[key] is None:
            return default
        value = as_number(raw[key], f"'{where}.{key}'")
        if not math.isfinite(value):
            raise ConfigError(
                f"'{where}.{key}' must be a finite number, got {raw[key]!r}"
            )
        return value

    def _int_or_default(self, raw: Dict[str, Any], key: str, where: str,
                        default: int) -> int:
        """Integer form of :meth:`_num_or_default`."""
        return int(self._num_or_default(raw, key, where, default))

    def _parse_mcmc(self, mcmc_raw: Dict[str, Any]) -> MCMCConfig:
        """Parse MCMC configuration.

        Args:
            mcmc_raw: Raw ``mcmc`` mapping.

        Returns:
            A validated MCMCConfig.

        Raises:
            ConfigError: On non-positive lengths, an inverted burn-in, or a
                non-boolean ``sample_from_prior``.
        """
        mcmc_type_str = mcmc_raw.get("type", "standard")
        # `int(x or DEFAULT)` maps an explicit 0/None back to the default, which
        # made every "<= 0" guard below unreachable: `chain_length: 0` produced a
        # silently 10^7-state chain. Default only on *absence*.
        chain_length = self._int_or_default(mcmc_raw, "chain_length", "mcmc", 10000000)
        pre_burnin = self._int_or_default(mcmc_raw, "pre_burnin", "mcmc", 0)
        if "store_every" in mcmc_raw:
            store_every = self._int_or_default(mcmc_raw, "store_every", "mcmc", 0)
        else:
            # MCMC.java counts this in *logged samples*, not generations
            # (`(sampleNr + 1) % storeEvery == 0`), so a single fixed number
            # would checkpoint a 10^7-generation chain once, at the very end.
            # Scale it with the chain so any run gets ~10 resume points.
            store_every = max(1, chain_length // 10000)
        particle_count = self._int_or_default(mcmc_raw, "particle_count", "mcmc", 1)
        sub_chain_length = self._int_or_default(mcmc_raw, "sub_chain_length", "mcmc", 10000)
        if chain_length <= 0:
            raise ConfigError(f"mcmc.chain_length must be > 0, got {chain_length}")
        if pre_burnin < 0:
            raise ConfigError(f"mcmc.pre_burnin must be >= 0, got {pre_burnin}")
        if pre_burnin >= chain_length:
            raise ConfigError(
                f"mcmc.pre_burnin ({pre_burnin}) must be smaller than chain_length "
                f"({chain_length}); otherwise the chain never samples."
            )
        if store_every == 0:
            raise ConfigError(
                "mcmc.store_every = 0 is not a valid interval. Use a positive number of "
                "samples to checkpoint every N logged samples, or -1 to disable "
                "checkpointing; 0 is neither, and would be written into the XML as "
                "storeEvery=\"0\"."
            )
        if store_every < -1:
            raise ConfigError(
                f"mcmc.store_every = {store_every} is not a valid interval; use a positive "
                f"number of samples (10000 by default) or -1 to disable checkpointing."
            )
        if store_every == -1:
            # Legal but reported: the validator warns that the chain cannot be
            # resumed, which is what complained about.
            pass
        if particle_count < 1:
            raise ConfigError(f"mcmc.particle_count must be >= 1, got {particle_count}")
        if sub_chain_length < 1:
            raise ConfigError(f"mcmc.sub_chain_length must be >= 1, got {sub_chain_length}")
        seed_raw = mcmc_raw.get("seed", 1)
        seed = self._int_or_default({"seed": seed_raw}, "seed", "mcmc", 1)
        if seed < 1:
            raise ConfigError(f"mcmc.seed must be >= 1, got {seed}")
        return MCMCConfig(
            chain_length=chain_length,
            pre_burnin=pre_burnin,
            store_every=store_every,
            sample_from_prior=self._as_bool(mcmc_raw, "sample_from_prior", False, "mcmc"),
            mcmc_type=ConfigParser._parse_enum(MCMCType, mcmc_type_str, "mcmc.type"),
            particle_count=particle_count,
            sub_chain_length=sub_chain_length,
            seed=seed,
        )

    def _parse_real_parameter(self, raw: Any, where: str = "parameter") -> RealParameter:
        """Parse a RealParameter from raw config (dict or scalar).

        Args:
            raw: A mapping with value/lower/upper/dimension/estimate, or a bare
                number.
            where: Location used in error messages.

        Returns:
            A RealParameter whose ``estimate`` flag and both bounds are honoured.

        Raises:
            ConfigError: On unknown keys, non-boolean ``estimate``, non-numeric
                bounds, or a value outside ``[lower, upper]``.
        """
        self._check_keys(raw, "parameter", where)
        if isinstance(raw, dict):
            value = raw.get("value", 1.0)
            lower = self._as_num(raw, "lower", where)
            upper = self._as_num(raw, "upper", where)
            dimension = int(self._as_num(raw, "dimension", where, 1) or 1)
            estimate = self._as_bool(raw, "estimate", True, where)
            if dimension < 1:
                raise ConfigError(f"{where}.dimension must be >= 1, got {dimension}")
            if isinstance(value, str):
                tokens = value.split()
                if len(tokens) > 1:
                    try:
                        numeric = [
                            self._finite(where, as_number(tok, f"{where}.value"))
                            for tok in tokens
                        ]
                    except ValueError as exc:
                        raise ConfigError(str(exc)) from exc
                    if dimension == 1:
                        dimension = len(numeric)
                    if len(numeric) != dimension:
                        raise ConfigError(
                            f"{where}: value lists {len(numeric)} entries but dimension is "
                            f"{dimension}"
                        )
                    for num in numeric:
                        self._check_bound(where, num, lower, upper)
                else:
                    num = self._finite(where, as_number(tokens[0], f"{where}.value"))
                    self._check_bound(where, num, lower, upper)
                    value = num
            elif value is not None:
                num = self._finite(where, as_number(value, f"{where}.value"))
                self._check_bound(where, num, lower, upper)
                value = num
            return RealParameter(
                value=value, lower=lower, upper=upper, dimension=dimension, estimate=estimate
            )
        return RealParameter(
            value=self._finite(where, as_number(raw, where)), estimate=True
        )

    @staticmethod
    def _finite(where: str, value: float) -> float:
        """Reject NaN/inf before it reaches XML attribute or int() conversion.

        Args:
            where: Location used in error messages.
            value: Candidate number.

        Returns:
            ``value`` unchanged.

        Raises:
            ConfigError: If the number is not finite.
        """
        if not math.isfinite(value):
            raise ConfigError(f"{where}: must be a finite number, got {value!r}")
        return value

    @staticmethod
    def _check_bound(where: str, value: float, lower: Optional[float], upper: Optional[float]):
        """Assert an initial value respects the parameter's bounds.

        Args:
            where: Location used in error messages.
            value: Candidate initial value.
            lower: Lower bound or None.
            upper: Upper bound or None.

        Raises:
            ConfigError: If the bounds are inverted or exclude the value.
        """
        if lower is not None and upper is not None and lower >= upper:
            raise ConfigError(f"{where}: lower ({lower}) must be strictly below upper ({upper})")
        if lower is not None and value < lower:
            raise ConfigError(
                f"{where}: value {value} is below its lower bound {lower}; BEAST2 rejects "
                f"an initial state outside the parameter's domain"
            )
        if upper is not None and value > upper:
            raise ConfigError(
                f"{where}: value {value} is above its upper bound {upper}; BEAST2 rejects "
                f"an initial state outside the parameter's domain"
            )

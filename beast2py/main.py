"""Command-line interface for Beast2Py.

Provides subcommands for generating BEAST2 XML, running calibration diagnostics,
validating XML, generating methods descriptions, and creating pipeline files.

Usage:
    beast2py generate --config config.yaml --output output.xml
    beast2py diagnose --config config.yaml --report report.html
    beast2py quick --alignment input.fasta --output output.xml
    beast2py validate --xml output.xml
    beast2py methods --config config.yaml --output methods.tex
    beast2py pipeline --config config.yaml --type snakemake
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import List, Optional

from . import __version__
from .config import ConfigError, ConfigParser
from .models import (
    BEASTConfig,
    CalibrationMethod,
    CalibrationPoint,
    ClockModelConfig,
    ClockModelType,
    InitTreeType,
    MCMCConfig,
    Partition,
    RealParameter,
    SiteModelConfig,
    TreePriorType,
)
from .sequence import SequenceReader

# ============================================================================
# Color / formatting helpers
# ============================================================================


class Colors:
    """ANSI color codes for terminal output."""

    RED = "\033[91m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    BLUE = "\033[94m"
    BOLD = "\033[1m"
    RESET = "\033[0m"


def _info(msg: str) -> None:
    """Print an info message."""
    print(f"{Colors.BLUE}[INFO]{Colors.RESET} {msg}")


def _success(msg: str) -> None:
    """Print a success message."""
    print(f"{Colors.GREEN}[OK]{Colors.RESET} {msg}")


def _warn(msg: str) -> None:
    """Print a warning message."""
    print(f"{Colors.YELLOW}[WARN]{Colors.RESET} {msg}", file=sys.stderr)


def _error(msg: str) -> None:
    """Print an error message."""
    print(f"{Colors.RED}[ERROR]{Colors.RESET} {msg}", file=sys.stderr)


def _verbose(msg: str, verbose: bool) -> None:
    """Print a verbose message if verbose mode is enabled."""
    if verbose:
        print(f"  {Colors.BOLD}[verbose]{Colors.RESET} {msg}")


# ============================================================================
# Subcommand: generate
# ============================================================================


def _validate_and_write(
    xml_str: str,
    output_path: str,
    config: BEASTConfig,
    beast2_validate: bool = False,
    beast2_optional: bool = False,
    beast2_path: str = "beast",
    force: bool = False,
    verbose: bool = False,
) -> int:
    """Run the three validation gates, then write the artifact.

    The gates are, in order:

    1. Beast2Py's own structural checks (ids, idrefs, duplicate ids, sequence
       sanity, prior coverage).
    2. Calibration conflict detection, so a model that contradicts itself is
       not silently emitted.
    3. The optional BEAST2 check, which parses *and* initialises the model
       (``XMLParser`` calls ``initAndValidate`` on every object).

    Nothing is written until every gate that ran has passed, and a failed gate
    returns non-zero. Previously ``generate`` printed
    "BEAST2 validation reported issues" and still returned 0 after writing the
    file, so a pipeline using exit codes recorded an invalid model as verified
    . The three outcomes are reported with distinct wording so
    that "[OK] structural checks passed" is never mistaken for "BEAST2 can run
    this analysis".

    Args:
        xml_str: The generated XML.
        output_path: Destination XML path.
        config: Parsed configuration, used for the conflict gate.
        beast2_validate: Whether to invoke the external BEAST2 check.
        beast2_optional: Downgrade an *unavailable* BEAST2 check from a hard
            stop to a warning. Only meaningful with ``beast2_validate``: asking
            for BEAST2 validation and then writing the file because BEAST2
            happens to be missing is a fail-open behaviour, so by default this
            exits 2 and writes nothing.
        beast2_path: BEAST2 executable or validation script.
        force: Overwrite an existing output file.
        verbose: Print per-gate detail.

    Returns:
        ``0`` when everything that ran passed, ``1`` otherwise.
    """
    out_path = Path(output_path)

    if out_path.exists() and not force:
        _error(
            f"Refusing to overwrite existing file: {out_path}. Pass --force to replace "
            f"it. Silently truncating an output while leaving its sibling .fingerprint.json "
            f"or .methods.tex behind mixes two analyses in one name."
        )
        return 1

    # --- Gate 1: structural checks ---
    from .validator import XMLValidator

    result = XMLValidator.validate(xml_str, config=config)
    if result.errors:
        _error("Gate 1/3 structural checks failed:")
        for err in result.errors:
            _error(f"  {err}")
        return 1
    for warning in result.warnings:
        _warn(f"  {warning}")
    _success("Gate 1/3 passed: Beast2Py structural checks (ids, idrefs, priors, operators).")
    _verbose(f"Gate 1: {len(result.errors)} errors, {len(result.warnings)} warnings", verbose)

    # --- Gate 2: calibration conflicts ---
    from .diagnostics import ConflictDetector

    conflicts = ConflictDetector.detect_conflicts(
        config.calibrations, settings=getattr(config, "diagnostics", None)
    )
    hard = [c for c in conflicts if c.severity == "error"]
    soft = [c for c in conflicts if c.severity != "error"]
    for conflict in soft:
        _warn(
            f"  [{conflict.conflict_type}] {conflict.calibration_a} <-> "
            f"{conflict.calibration_b}: {conflict.description}"
        )
    if hard:
        _error(f"Gate 2/3 failed: {len(hard)} calibration conflict(s):")
        for conflict in hard:
            _error(
                f"  [{conflict.conflict_type}] {conflict.calibration_a} <-> "
                f"{conflict.calibration_b}: {conflict.description}"
            )
            if conflict.recommendation:
                _error(f"    -> {conflict.recommendation}")
        return 1
    _success("Gate 2/3 passed: no contradictory calibration conflicts.")

    add_ons = [
        name
        for prefix, name in (("bdsky.", "bdsky"), ("nestedsampling.", "nested-sampling"))
        if prefix in xml_str
    ]
    if add_ons:
        _warn(
            "This XML requires BEAST2 add-on package(s): "
            + ", ".join(add_ons)
            + ". Their classes are not part of the core packages, so without them "
            "installed BEAST2 fails at start-up (it does not silently downgrade). "
            "Install them with `beast -pm` / the add-on manager before running."
        )

    # --- Gate 3: BEAST2 parse + initAndValidate (optional) ---
    if beast2_validate:
        import tempfile

        _info("Gate 3/3: validating with BEAST2 (parse + initialise)...")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            "w", suffix=".xml", dir=str(out_path.parent), delete=False, encoding="utf-8"
        ) as handle:
            handle.write(xml_str)
            probe_path = handle.name
        try:
            success, output, status = XMLValidator.beast2_validation_status(
                probe_path, beast2_path=beast2_path
            )
        finally:
            Path(probe_path).unlink(missing_ok=True)
        if status == "unavailable":
            message = (
                "Gate 3/3 unavailable: no BEAST2 installation was found, so the model was "
                "checked structurally only. Install BEAST2 v2.7.x (or pass --beast2-path, "
                "or set BEAST2_JAR) before treating this XML as BEAST2-verified; pass "
                "--allow-unvalidated to write the file anyway."
            )
            if not beast2_optional:
                _error(message)
                print(output.strip())
                return 2
            _warn(message)
            print(output.strip())
        elif not success:
            _error("Gate 3/3 failed: BEAST2 rejected the model.")
            print(output.strip())
            return 1
        else:
            _success("Gate 3/3 passed: BEAST2 parsed and initialised the model.")

    # --- All gates that ran have passed: write ---
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(xml_str, encoding="utf-8")
    _success(f"XML written to: {out_path}")
    return 0


def cmd_generate(args: argparse.Namespace) -> int:
    """Generate BEAST2 XML from a YAML configuration file.

    This is the primary command for creating BEAST2 XML files. It reads a
    YAML configuration, assembles the complete XML, optionally runs
    diagnostics, generates a fingerprint, and produces a methods description.
    """
    verbose = args.verbose
    config_path = args.config
    output_path = args.output

    _info(f"Reading configuration: {config_path}")
    try:
        parser = ConfigParser()
        config = parser.parse(config_path)
    except ConfigError as e:
        _error(f"Configuration error: {e}")
        return 1
    except FileNotFoundError as e:
        _error(f"File not found: {e}")
        return 1

    _verbose(f"Analysis: {config.metadata.get('analysis_name', 'N/A')}", verbose)
    _verbose(f"Partitions: {len(config.partitions)}", verbose)
    _verbose(f"Calibrations: {len(config.calibrations)}", verbose)
    _verbose(f"Tree prior: {config.tree_prior_type.value}", verbose)

    if config.calibrations and config.initialization.tree_type == InitTreeType.UPGMA:
        _warn(
            "UPGMA starting trees have genetic-distance-scale heights, which are "
            "typically far below MRCA calibration ages; BEAST2 may fail to "
            "initialise. Consider initialization.tree_type: random."
        )

    # Generate XML
    _info("Generating BEAST2 XML...")
    from .xml_writer import XMLWriter

    try:
        writer = XMLWriter(config)
        xml_str = writer.generate_xml()
    except ValueError as e:
        # e.g. cyclic/invalid partition links detected while assembling models
        _error(f"Model assembly error: {e}")
        return 1

    # Validate, then write only if every gate that ran has passed
    status = _validate_and_write(
        xml_str,
        output_path,
        config,
        beast2_validate=args.beast2_validate,
        beast2_path=args.beast2_path,
        beast2_optional=getattr(args, "allow_unvalidated", False),
        force=getattr(args, "force", False),
        verbose=verbose,
    )
    if status != 0:
        return status

    # Optional: Flag reproducibility siblings that are about to go stale: an overwritten
    # XML with an untouched .fingerprint.json / .methods.tex on disk mixes two
    # analyses under one name .
    if not (args.fingerprint and args.methods):
        for suffix in (".fingerprint.json", ".methods.tex"):
            stale = Path(output_path).with_suffix(suffix)
            if stale.exists():
                _warn(
                    f"{stale} belongs to an earlier run of this output and is not being "
                    f"refreshed. Pass --fingerprint --methods to regenerate it."
                )
    if args.fingerprint:
        _info("Generating analysis fingerprint...")
        from .reproducibility import FingerprintGenerator

        fp_path = str(Path(output_path).with_suffix(".fingerprint.json"))
        # Digest the file as it now sits on disk, so the sidecar describes the
        # artefact BEAST2 will be handed rather than the intent behind it.
        written_xml = Path(output_path).read_text(encoding="utf-8")
        FingerprintGenerator.write_fingerprint_file(config, fp_path, written_xml)
        fingerprint = FingerprintGenerator.generate_fingerprint(config)
        _success(f"Fingerprint: {fingerprint}")
        _success(f"Fingerprint file: {fp_path}")

    # Optional: Methods description
    if args.methods:
        _info("Generating methods description...")
        from .reproducibility import MethodsGenerator

        methods_text = MethodsGenerator.generate_methods(config)
        methods_path = str(Path(output_path).with_suffix(".methods.tex"))
        Path(methods_path).write_text(methods_text, encoding="utf-8")
        _success(f"Methods description: {methods_path}")

    # Optional: Diagnostics
    if args.diagnose:
        _info("Running calibration diagnostics...")
        _run_diagnostics(
            config,
            output_dir=str(Path(output_path).parent / "diagnostics"),
            run_sensitivity=True,
            run_visualization=True,
            verbose=verbose,
        )

    _success("Done!")
    return 0


# Public alias: the Python API (beast2py.api) routes its write path through this
# same gate helper, so the CLI and the library cannot disagree about what
# "validated" means ( covered only the CLI before).
validate_and_write = _validate_and_write


# ============================================================================
# Subcommand: diagnose
# ============================================================================


def cmd_diagnose(args: argparse.Namespace) -> int:
    """Run calibration diagnostics on a configuration.

    Detects calibration conflicts, generates sensitivity analysis XMLs,
    creates distribution visualizations, and produces an HTML report.
    This is the core innovation command of Beast2Py.
    """
    verbose = args.verbose
    config_path = args.config

    _info(f"Reading configuration: {config_path}")
    try:
        parser = ConfigParser()
        config = parser.parse(config_path)
    except ConfigError as e:
        _error(f"Configuration error: {e}")
        return 1
    except FileNotFoundError as e:
        _error(f"File not found: {e}")
        return 1

    report_path = args.report
    output_dir = args.sensitivity_xml_dir or str(Path(report_path).parent / "diagnostics")

    _run_diagnostics(
        config,
        output_dir=output_dir,
        run_sensitivity=args.sensitivity,
        run_visualization=args.visualize,
        report_path=report_path,
        verbose=verbose,
    )

    _success("Diagnostics complete!")
    return 0


def _run_diagnostics(
    config: BEASTConfig,
    output_dir: str,
    run_sensitivity: bool,
    run_visualization: bool,
    report_path: Optional[str] = None,
    verbose: bool = False,
) -> None:
    """Run the diagnostic workflow and print results."""
    from .diagnostics import DiagnosticEngine, ReportGenerator

    _info("Running diagnostic engine...")
    report = DiagnosticEngine.run(
        config,
        output_dir=output_dir,
        run_sensitivity=run_sensitivity,
        run_visualization=run_visualization,
    )

    # Print conflicts
    n_errors = sum(1 for c in report.conflicts if c.severity == "error")
    n_warnings = sum(1 for c in report.conflicts if c.severity == "warning")
    unchecked = [c for c in report.conflicts if c.conflict_type == "unchecked"]
    actionable = [c for c in report.conflicts if c.conflict_type != "unchecked"]
    if actionable:
        _warn(f"Found {n_errors} error(s) and {n_warnings} warning(s):")
    for c in actionable:
        icon = "✗" if c.severity == "error" else "⚠"
        print(f"  {icon} [{c.conflict_type}] {c.calibration_a} ↔ {c.calibration_b}")
        print(f"    {c.description}")
        if c.recommendation:
            print(f"    → {c.recommendation}")
    if not actionable and not unchecked:
        _success("No calibration conflicts detected.")
    if unchecked:
        # "Nothing found" and "nothing checked" are different statements, and
        # only the first one deserves a green tick .
        _warn(
            f"{len(unchecked)} calibration pair(s) could NOT be checked "
            f"({', '.join(c.calibration_a + '/' + c.calibration_b for c in unchecked)}): "
            f"a distribution with undefined quantiles has no comparable interval."
        )
        for c in unchecked:
            print(f"  ○ [unchecked] {c.calibration_a} ↔ {c.calibration_b}")
            print(f"    {c.description}")

    # Print sensitivity info
    if report.sensitivity_xmls:
        _info(f"Generated {len(report.sensitivity_xmls)} sensitivity analysis XML files:")
        for x in report.sensitivity_xmls:
            print(f"  - {x}")

    # Print visualization info
    if report.visualization_file:
        _success(f"Visualization: {report.visualization_file}")
    if report.visualization_error:
        _warn(report.visualization_error)

    # Generate HTML report
    if report_path:
        _info(f"Generating HTML report: {report_path}")
        ReportGenerator.generate_report(
            config=config,
            conflicts=report.conflicts,
            summaries=report.summaries,
            sensitivity_xmls=report.sensitivity_xmls,
            visualization_file=report.visualization_file,
            output_file=report_path,
        )
        _success(f"Report written to: {report_path}")

    _verbose(f"Summaries: {len(report.summaries)} calibration points", verbose)


# ============================================================================
# Subcommand: quick
# ============================================================================


def cmd_quick(args: argparse.Namespace) -> int:
    """Quick mode: generate XML from command-line arguments without a full YAML config.

    Supports single partition + single calibration point scenarios.
    For complex analyses, use the `generate` command with a YAML config file.
    """
    verbose = args.verbose

    # Read alignment
    _info(f"Reading alignment: {args.alignment}")
    try:
        alignment = SequenceReader.read(
            args.alignment,
            format="auto",
            data_type=getattr(args, "data_type", "auto"),
            alignment_id="alignment",
        )
    except FileNotFoundError as e:
        _error(str(e))
        return 1
    except ValueError as e:
        # AlignmentError subclasses ValueError, so header/alphabet problems from
        # the reader land here with a message that already names the record.
        _error(f"Failed to parse alignment: {e}")
        return 1

    _verbose(f"Taxa: {alignment.n_taxa}, Sites: {alignment.n_sites}", verbose)

    # Parse substitution model
    subst_type = args.subst_model.lower()
    if subst_type not in ConfigParser.VALID_SUBSTITUTION_MODELS:
        _error(
            f"Invalid substitution model: '{subst_type}' "
            f"(valid: {', '.join(sorted(ConfigParser.VALID_SUBSTITUTION_MODELS))})"
        )
        return 1
    subst_model = {"type": subst_type}

    # Parse gamma categories
    gamma_cats = args.gamma_categories or 0
    gamma_shape = None
    if gamma_cats > 0:
        gamma_shape = RealParameter(value=0.5, lower=0.0)

    site_model = SiteModelConfig(
        substitution_model=subst_model,
        gamma_categories=gamma_cats,
        gamma_shape=gamma_shape,
    )

    # Parse clock model
    clock_type_str = args.clock_model.lower()
    if clock_type_str not in ConfigParser.VALID_CLOCK_MODELS:
        _error(
            f"Invalid clock model: '{clock_type_str}' "
            f"(valid: {', '.join(sorted(ConfigParser.VALID_CLOCK_MODELS))})"
        )
        return 1
    clock_model = ClockModelConfig(type=ClockModelType(clock_type_str))

    # Build partition
    partition = Partition(
        id="alignment",
        alignment=alignment,
        site_model=site_model,
        clock_model=clock_model,
        tree="shared",
    )

    # Parse tree prior
    tp_type_str = args.tree_prior.lower()
    if tp_type_str not in ConfigParser.VALID_TREE_PRIORS:
        _error(
            f"Invalid tree prior: '{tp_type_str}' "
            f"(valid: {', '.join(sorted(ConfigParser.VALID_TREE_PRIORS))})"
        )
        return 1
    tree_prior_type = TreePriorType(tp_type_str)

    # Parse calibrations from YAML file
    calibrations: List[CalibrationPoint] = []
    if args.calibration_yaml:
        _info(f"Reading calibrations: {args.calibration_yaml}")
        calibrations = _parse_calibration_yaml(
            args.calibration_yaml, known_taxa=alignment.taxa_names
        )
        _verbose(f"Calibration points: {len(calibrations)}", verbose)

    # Parse MCMC
    mcmc = MCMCConfig(
        chain_length=args.chain_length,
        pre_burnin=args.pre_burnin,
    )

    # Build config
    metadata = {
        "analysis_name": args.name or "quick_analysis",
        "beast2_version": "2.7.8",
        "tool_version": __version__,
        "output_file": args.output,
    }

    config = BEASTConfig(
        metadata=metadata,
        partitions=[partition],
        tree_prior_type=tree_prior_type,
        tree_prior_params={},
        calibrations=calibrations,
        calibration_method=CalibrationMethod.MRCA_PRIOR,
        mcmc=mcmc,
    )

    # Generate XML
    _info("Generating BEAST2 XML...")
    from .xml_writer import XMLWriter

    writer = XMLWriter(config)
    xml_str = writer.generate_xml()

    status = _validate_and_write(
        xml_str,
        args.output,
        config,
        force=getattr(args, "force", False),
        verbose=verbose,
    )
    if status != 0:
        return status

    _success("Done!")
    return 0


def _parse_calibration_yaml(yaml_path: str, known_taxa=None) -> List[CalibrationPoint]:
    """Parse a calibration-only YAML file.

    Expected format:
        - name: humanChimpMRCA
          taxa: ["human", "chimp"]
          monophyletic: true
          distribution:
            type: normal
            parameters: {mean: 6.0, sigma: 0.5}

    The parsing is delegated to :meth:`ConfigParser._parse_calibration` so the
    ``--calibration-yaml`` shortcut gets exactly the same checks as a full YAML
    configuration. Its own lenient copy used to accept ``monophyletic:
    "false"`` (inverted by Python truthiness), an inverted ``uniform`` bound
    and an unknown distribution type, all silently.

    Args:
        yaml_path: Path to the YAML list.
        known_taxa: Taxon names of the analysis, or None to skip membership
            checking.

    Returns:
        List of validated CalibrationPoint objects.

    Raises:
        ConfigError: If any entry is invalid.
        ValueError: If the file is not a list.
    """
    import yaml

    with open(yaml_path, "r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    if not isinstance(raw, list):
        raise ValueError("Calibration YAML must be a list of calibration points")

    parser = ConfigParser()
    errors: List[str] = []
    calibrations: List[CalibrationPoint] = []
    for index, cal_raw in enumerate(raw):
        cal = parser._parse_calibration(
            cal_raw,
            errors,
            set(known_taxa or ()),
            f"calibrations[{index}]",
            check_taxa=known_taxa is not None,
        )
        if cal is not None:
            calibrations.append(cal)
    if errors:
        raise ConfigError("; ".join(errors))
    return calibrations


# ============================================================================
# Subcommand: validate
# ============================================================================


def cmd_validate(args: argparse.Namespace) -> int:
    """Validate a BEAST2 XML file, reporting each level separately.

    Three levels are distinguished on purpose: Beast2Py's structural checks,
    and (optionally) BEAST2 parsing plus model initialisation. "BEAST2 could
    read this file" is not the same statement as "this model is what the user
    specified", and the old wording blurred them.
    """
    xml_path = args.xml
    _info(f"Validating XML: {xml_path}")

    from .validator import XMLValidator

    # Structural validation
    result = XMLValidator.validate_file(xml_path)

    if result.errors:
        _error(f"Level 1/2 (structural checks) failed with {len(result.errors)} error(s):")
        for err in result.errors:
            _error(f"  {err}")
    else:
        _success("Level 1/2 passed: structural checks (ids, idrefs, priors, bounds).")

    if result.warnings:
        _warn(f"{len(result.warnings)} warning(s):")
        for w in result.warnings:
            _warn(f"  {w}")

    # BEAST2 validation
    if args.beast2:
        _info("Level 2/2: validating with BEAST2 (parse + model initialisation)...")
        beast2_path = args.beast2_path or "beast"
        success, output, status = XMLValidator.beast2_validation_status(xml_path, beast2_path)
        if status == "unavailable":
            _warn(
                "Level 2/2 skipped: BEAST2 is not installed or its validation script was "
                "not found, so this XML has only been checked structurally. Do not treat "
                "the structural pass as BEAST2 verification."
            )
            print(output.strip())
        elif success:
            _success("Level 2/2 passed: BEAST2 parsed and initialised the model.")
        else:
            _error("Level 2/2 failed: BEAST2 rejected the model:")
            print(output)
            return 1

    return 0 if result.is_valid else 1


# ============================================================================
# Subcommand: methods
# ============================================================================


def cmd_methods(args: argparse.Namespace) -> int:
    """Generate a LaTeX-formatted methods description from a configuration.

    The methods description includes substitution model, clock model, tree prior,
    calibration points, MCMC settings, and the analysis fingerprint.
    """
    config_path = args.config
    _info(f"Reading configuration: {config_path}")

    try:
        parser = ConfigParser()
        config = parser.parse(config_path)
    except ConfigError as e:
        _error(f"Configuration error: {e}")
        return 1
    except FileNotFoundError as e:
        _error(f"File not found: {e}")
        return 1

    _info("Generating methods description...")
    from .reproducibility import MethodsGenerator

    methods_text = MethodsGenerator.generate_methods(config)

    output_path = args.output
    if output_path:
        output_dir = Path(output_path).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        Path(output_path).write_text(methods_text, encoding="utf-8")
        _success(f"Methods description written to: {output_path}")
    else:
        # Print to stdout
        print()
        print(methods_text)
        print()

    return 0


# ============================================================================
# Subcommand: pipeline
# ============================================================================


def cmd_pipeline(args: argparse.Namespace) -> int:
    """Generate Snakemake or Nextflow pipeline files for reproducible analysis.

    The pipeline includes rules for XML generation, BEAST2 execution,
    diagnostics, and methods description generation.
    """
    config_path = args.config
    _info(f"Reading configuration: {config_path}")

    try:
        parser = ConfigParser()
        config = parser.parse(config_path)
    except ConfigError as e:
        _error(f"Configuration error: {e}")
        return 1
    except FileNotFoundError as e:
        _error(f"File not found: {e}")
        return 1

    output_dir = args.output_dir or "."
    os.makedirs(output_dir, exist_ok=True)

    from .reproducibility import PipelineGenerator

    pipeline_type = args.type.lower()

    if pipeline_type == "snakemake":
        _info("Generating Snakemake pipeline...")
        snakefile = PipelineGenerator.generate_snakemake(
            config,
            output_dir=output_dir,
            config_file=config_path,
        )
        _success(f"Snakefile written to: {snakefile}")

    elif pipeline_type == "nextflow":
        _info("Generating Nextflow pipeline...")
        nf_file = PipelineGenerator.generate_nextflow(
            config,
            output_dir=output_dir,
            config_file=config_path,
        )
        _success(f"Nextflow file written to: {nf_file}")

    elif pipeline_type == "both":
        _info("Generating Snakemake pipeline...")
        snakefile = PipelineGenerator.generate_snakemake(
            config,
            output_dir=output_dir,
            config_file=config_path,
        )
        _success(f"Snakefile written to: {snakefile}")

        _info("Generating Nextflow pipeline...")
        nf_file = PipelineGenerator.generate_nextflow(
            config,
            output_dir=output_dir,
            config_file=config_path,
        )
        _success(f"Nextflow file written to: {nf_file}")

    else:
        _error(f"Unknown pipeline type: {pipeline_type}. Use snakemake, nextflow, or both.")
        return 1

    _success("Done!")
    return 0


# ============================================================================
# Subcommand: fingerprint
# ============================================================================


def cmd_fingerprint(args: argparse.Namespace) -> int:
    """Generate an analysis fingerprint for reproducibility tracking.

    The fingerprint is a unique identifier based on the configuration hash,
    tool version, and generation timestamp.
    """
    config_path = args.config
    _info(f"Reading configuration: {config_path}")

    try:
        parser = ConfigParser()
        config = parser.parse(config_path)
    except ConfigError as e:
        _error(f"Configuration error: {e}")
        return 1
    except FileNotFoundError as e:
        _error(f"File not found: {e}")
        return 1

    from .reproducibility import FingerprintGenerator

    fingerprint = FingerprintGenerator.generate_fingerprint(config)
    fp_dict = FingerprintGenerator.generate_fingerprint_dict(config)

    _success(f"Fingerprint: {fingerprint}")
    print()

    if args.verbose:
        for k, v in fp_dict.items():
            print(f"  {k}: {v}")
        print()

    if args.output:
        FingerprintGenerator.write_fingerprint_file(config, args.output)
        _success(f"Fingerprint file written to: {args.output}")

    return 0


# ============================================================================
# Subcommand: list-models
# ============================================================================


def cmd_list_models(args: argparse.Namespace) -> int:
    """List all supported models (substitution, clock, tree prior, distributions)."""
    from .config import ConfigParser

    print(f"\n{Colors.BOLD}Beast2Py v{__version__} — Supported Models{Colors.RESET}\n")

    print(f"{Colors.BOLD}Substitution Models:{Colors.RESET}")
    for m in sorted(ConfigParser.VALID_SUBSTITUTION_MODELS):
        print(f"  - {m}")

    print(f"\n{Colors.BOLD}Clock Models:{Colors.RESET}")
    for m in sorted(ConfigParser.VALID_CLOCK_MODELS):
        print(f"  - {m}")

    print(f"\n{Colors.BOLD}Tree Priors:{Colors.RESET}")
    for m in sorted(ConfigParser.VALID_TREE_PRIORS):
        print(f"  - {m}")

    print(f"\n{Colors.BOLD}Calibration Distributions:{Colors.RESET}")
    for m in sorted(ConfigParser.VALID_DISTRIBUTIONS):
        print(f"  - {m}")

    print(f"\n{Colors.BOLD}Data Types:{Colors.RESET}")
    for m in sorted(ConfigParser.VALID_DATA_TYPES):
        print(f"  - {m}")

    print()
    return 0


# ============================================================================
# Argument parser construction
# ============================================================================


def build_parser() -> argparse.ArgumentParser:
    """Build the main argument parser with all subcommands.

    Returns:
        Configured ArgumentParser.
    """
    parser = argparse.ArgumentParser(
        prog="beast2py",
        description=(
            "Beast2Py: A Python framework for reproducible "
            "divergence time estimation with automated calibration prior "
            "specification, validation, and diagnostics."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  beast2py generate -c config.yaml -o output.xml\n"
            "  beast2py diagnose -c config.yaml --report report.html\n"
            "  beast2py quick -a input.fasta -o output.xml \\\n"
            "      --tree-prior yule --subst-model hky --clock-model strict \\\n"
            "      --calibration-yaml cals.yaml --chain-length 10000000\n"
            "  beast2py validate --xml output.xml\n"
            "  beast2py methods -c config.yaml -o methods.tex\n"
            "  beast2py pipeline -c config.yaml --type snakemake\n"
        ),
    )

    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"Beast2Py v{__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        title="Available commands",
        metavar="<command>",
    )

    # --- generate ---
    p_gen = subparsers.add_parser(
        "generate",
        help="Generate BEAST2 XML from a YAML configuration file",
        description="Generate BEAST2 XML from a YAML configuration file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_gen.add_argument("-c", "--config", required=True, help="YAML configuration file path")
    p_gen.add_argument("-o", "--output", required=True, help="Output XML file path")
    p_gen.add_argument(
        "--diagnose",
        action="store_true",
        help="Also run calibration diagnostics",
    )
    p_gen.add_argument(
        "--fingerprint",
        action="store_true",
        help="Generate analysis fingerprint file",
    )
    p_gen.add_argument(
        "--methods",
        action="store_true",
        help="Generate LaTeX methods description",
    )
    p_gen.add_argument(
        "--beast2-validate",
        action="store_true",
        dest="beast2_validate",
        help="Validate XML with BEAST2 -validate",
    )
    p_gen.add_argument(
        "--beast2-path",
        default="beast",
        dest="beast2_path",
        help="Path to BEAST2 executable (default: beast)",
    )
    p_gen.add_argument(
        "--allow-unvalidated",
        action="store_true",
        dest="allow_unvalidated",
        help="With --beast2-validate, write the XML even when BEAST2 is not "
        "installed. By default that case exits 2, so a CI job cannot turn "
        "'BEAST2 missing' into a green run.",
    )
    p_gen.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the output file if it already exists",
    )
    p_gen.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p_gen.set_defaults(func=cmd_generate)

    # --- diagnose ---
    p_diag = subparsers.add_parser(
        "diagnose",
        help="Run calibration diagnostics (conflict detection, sensitivity analysis)",
        description=(
            "Run calibration diagnostics: conflict detection, prior sensitivity "
            "analysis, distribution visualization, and HTML report generation. "
            "This is the core innovation of Beast2Py."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_diag.add_argument("-c", "--config", required=True, help="YAML configuration file path")
    p_diag.add_argument("--report", required=True, help="Output HTML report path")
    p_diag.add_argument(
        "--sensitivity",
        action="store_true",
        help="Generate prior sensitivity analysis XMLs",
    )
    p_diag.add_argument(
        "--visualize",
        action="store_true",
        help="Generate calibration distribution density plots",
    )
    p_diag.add_argument(
        "--sensitivity-xml-dir",
        default=None,
        dest="sensitivity_xml_dir",
        help="Output directory for sensitivity XMLs (default: <report_dir>/diagnostics)",
    )
    p_diag.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p_diag.set_defaults(func=cmd_diagnose)

    # --- quick ---
    p_quick = subparsers.add_parser(
        "quick",
        help="Quick mode: generate XML from command-line arguments (single partition)",
        description=(
            "Quick mode for single-partition + single-tree-prior analyses. "
            "Calibration points are specified via a separate YAML file. "
            "For complex configurations, use 'generate' with a full YAML config."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_quick.add_argument(
        "-a", "--alignment", required=True, help="Sequence file path (FASTA/NEXUS)"
    )
    p_quick.add_argument("-o", "--output", required=True, help="Output XML file path")
    p_quick.add_argument(
        "--data-type",
        default="auto",
        dest="data_type",
        choices=["auto", "nucleotide", "aminoacid"],
        help=(
            "Sequence type. 'auto' infers it from the residues present, which "
            "is what makes an amino-acid alignment work here at all (default: auto)"
        ),
    )
    p_quick.add_argument(
        "--tree-prior",
        default="yule",
        dest="tree_prior",
        help="Tree prior type (default: yule)",
    )
    p_quick.add_argument(
        "--subst-model",
        default="hky",
        dest="subst_model",
        help="Substitution model (default: hky)",
    )
    p_quick.add_argument(
        "--clock-model",
        default="strict",
        dest="clock_model",
        help="Clock model (default: strict)",
    )
    p_quick.add_argument(
        "--gamma-categories",
        type=int,
        default=0,
        dest="gamma_categories",
        help="Number of gamma rate categories (default: 0 = no gamma)",
    )
    p_quick.add_argument(
        "--calibration-yaml",
        default=None,
        dest="calibration_yaml",
        help="YAML file with calibration points",
    )
    p_quick.add_argument(
        "--chain-length",
        type=int,
        default=10000000,
        dest="chain_length",
        help="MCMC chain length (default: 10000000)",
    )
    p_quick.add_argument(
        "--pre-burnin",
        type=int,
        default=0,
        dest="pre_burnin",
        help="MCMC pre-burnin (default: 0)",
    )
    p_quick.add_argument("--name", default=None, help="Analysis name")
    p_quick.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the output file if it already exists",
    )
    p_quick.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p_quick.set_defaults(func=cmd_quick)

    # --- validate ---
    p_val = subparsers.add_parser(
        "validate",
        help="Validate a BEAST2 XML file",
        description="Validate a BEAST2 XML file for structural correctness and ID consistency.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_val.add_argument("--xml", required=True, help="XML file to validate")
    p_val.add_argument(
        "--beast2",
        action="store_true",
        help="Also validate with BEAST2 -validate",
    )
    p_val.add_argument(
        "--beast2-path",
        default="beast",
        dest="beast2_path",
        help="Path to BEAST2 executable (default: beast)",
    )
    p_val.set_defaults(func=cmd_validate)

    # --- methods ---
    p_methods = subparsers.add_parser(
        "methods",
        help="Generate a LaTeX methods description",
        description="Generate a LaTeX-formatted methods description from a YAML configuration.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_methods.add_argument("-c", "--config", required=True, help="YAML configuration file path")
    p_methods.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output file path (default: print to stdout)",
    )
    p_methods.set_defaults(func=cmd_methods)

    # --- pipeline ---
    p_pipe = subparsers.add_parser(
        "pipeline",
        help="Generate Snakemake or Nextflow pipeline files",
        description="Generate reproducible pipeline files (Snakemake or Nextflow).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_pipe.add_argument("-c", "--config", required=True, help="YAML configuration file path")
    p_pipe.add_argument(
        "--type",
        default="snakemake",
        choices=["snakemake", "nextflow", "both"],
        help="Pipeline type (default: snakemake)",
    )
    p_pipe.add_argument(
        "--output-dir",
        default=".",
        dest="output_dir",
        help="Output directory (default: current directory)",
    )
    p_pipe.set_defaults(func=cmd_pipeline)

    # --- fingerprint ---
    p_fp = subparsers.add_parser(
        "fingerprint",
        help="Generate an analysis fingerprint",
        description="Generate a unique analysis fingerprint for reproducibility tracking.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_fp.add_argument("-c", "--config", required=True, help="YAML configuration file path")
    p_fp.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output JSON file path (default: print to stdout only)",
    )
    p_fp.add_argument("-v", "--verbose", action="store_true", help="Show detailed fingerprint info")
    p_fp.set_defaults(func=cmd_fingerprint)

    # --- list-models ---
    p_lm = subparsers.add_parser(
        "list-models",
        help="List all supported models",
        description=(
            "List all supported substitution models, clock models, "
            "tree priors, and distributions."
        ),
    )
    p_lm.set_defaults(func=cmd_list_models)

    return parser


# ============================================================================
# Main entry point
# ============================================================================


def main(argv: Optional[List[str]] = None) -> int:
    """Main CLI entry point.

    Args:
        argv: Command-line arguments (default: sys.argv[1:]).

    Returns:
        Exit code (0 for success, non-zero for error).
    """
    parser = build_parser()

    if argv is None:
        argv = sys.argv[1:]

    # No command provided
    if not argv:
        parser.print_help()
        return 0

    args = parser.parse_args(argv)

    # Dispatch to subcommand handler
    if hasattr(args, "func"):
        try:
            return args.func(args)
        except KeyboardInterrupt:
            _error("Interrupted by user.")
            return 130
        except Exception as e:
            _error(f"Unexpected error: {e}")
            import traceback

            traceback.print_exc()
            return 1
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    sys.exit(main())

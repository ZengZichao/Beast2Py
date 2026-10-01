"""Beast2Py API: Programmatic interface for BEAST2 XML generation, calibration
diagnostics, validation, and reproducibility.

This module provides a clean, high-level API that mirrors all CLI functionality.
Users can import Beast2Py as a Python library and call these functions directly
without invoking the command-line interface.

Usage example::

    from beast2py.api import Beast2Py

    # Create an API instance
    b2p = Beast2Py()

    # Generate XML from a YAML config file
    xml_str = b2p.generate_xml("config.yaml")

    # Run calibration diagnostics
    report = b2p.diagnose("config.yaml", output_dir="diagnostics/")

    # Validate generated XML
    result = b2p.validate_xml(xml_str)

    # Generate LaTeX methods description
    methods = b2p.generate_methods("config.yaml")

Alternatively, use the module-level convenience functions::

    from beast2py.api import generate_xml, diagnose, validate_xml

    xml_str = generate_xml("config.yaml")
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from .config import ConfigError, ConfigParser
from .models import (
    BEASTConfig,
    CalibrationMethod,
    CalibrationPoint,
    ClockModelConfig,
    ClockModelType,
    DistributionConfig,
    MCMCConfig,
    Partition,
    Provenance,
    RealParameter,
    SiteModelConfig,
    TreePriorType,
)
from .sequence import SequenceReader

# ============================================================================
# Type aliases
# ============================================================================

ConfigInput = Union[str, Dict[str, Any], BEASTConfig]
"""Type for configuration input: file path, dict, or BEASTConfig object."""


# ============================================================================
# Helper functions
# ============================================================================


def _resolve_config(
    config: ConfigInput,
    base_dir: Optional[str] = None,
) -> BEASTConfig:
    """Resolve a configuration input into a BEASTConfig object.

    Args:
        config: Configuration input (file path, dict, or BEASTConfig).
        base_dir: Base directory for resolving relative file paths.

    Returns:
        A BEASTConfig object.

    Raises:
        ConfigError: If the configuration is invalid.
        FileNotFoundError: If a config file is not found.
    """
    if isinstance(config, BEASTConfig):
        return config

    if isinstance(config, str):
        parser = ConfigParser(base_dir=base_dir)
        return parser.parse(config)

    if isinstance(config, dict):
        parser = ConfigParser(base_dir=base_dir)
        return parser.parse_dict(config)

    raise TypeError(
        f"config must be a file path (str), dict, or BEASTConfig, got {type(config).__name__}"
    )


# ============================================================================
# Core API class
# ============================================================================


class Beast2Py:
    """High-level API for Beast2Py functionality.

    This class provides a programmatic interface to all Beast2Py features,
    mirroring the CLI commands. Each method accepts a configuration input
    (file path, dict, or BEASTConfig object) and returns structured results.

    Attributes:
        version: Beast2Py version string.
    """

    def __init__(self, base_dir: Optional[str] = None):
        """Initialize the Beast2Py API.

        Args:
            base_dir: Base directory for resolving relative file paths
                in configuration. Defaults to current working directory.
        """
        self.base_dir = base_dir or str(Path.cwd())
        from . import __version__

        self.version = __version__

    # ------------------------------------------------------------------------
    # XML Generation
    # ------------------------------------------------------------------------

    def generate_xml(
        self,
        config: ConfigInput,
        output: Optional[str] = None,
        force: bool = False,
    ) -> str:
        """Generate BEAST2 XML from a configuration.

        Args:
            config: Configuration input (file path, dict, or BEASTConfig).
            output: Optional output file path. When given, the XML goes through
                the same validation gates as the CLI and is written only if
                they all pass.
            force: Overwrite an existing output file (default: refuse).

        Returns:
            The generated XML string.

        Raises:
            ConfigError: If the configuration is invalid.
            FileNotFoundError: If referenced files are not found.
            ValueError: If an output path was given and a validation gate
                failed, or the file exists and ``force`` was not set. Nothing
                is written in that case: writing first and gating afterwards
                would leave a file on disk that never passed the checks.

        Example::

            b2p = Beast2Py()
            xml = b2p.generate_xml("config.yaml", output="output.xml")
        """
        beast_config = _resolve_config(config, self.base_dir)

        from .xml_writer import XMLWriter

        writer = XMLWriter(beast_config)
        xml_str = writer.generate_xml()

        if output:
            from .main import validate_and_write

            status = validate_and_write(xml_str, output, beast_config, force=force)
            if status != 0:
                raise ValueError(
                    f"Validation gates failed for {output}; nothing was written. "
                    f"See the messages above, or call XMLValidator.validate() directly "
                    f"for the structured result."
                )

        return xml_str

    # ------------------------------------------------------------------------
    # Quick Mode (single partition from command-line args)
    # ------------------------------------------------------------------------

    def quick_generate(
        self,
        alignment: str,
        output: str,
        tree_prior: str = "yule",
        subst_model: str = "hky",
        clock_model: str = "strict",
        gamma_categories: int = 0,
        calibration_yaml: Optional[str] = None,
        chain_length: int = 10_000_000,
        pre_burnin: int = 0,
        name: Optional[str] = None,
        force: bool = False,
    ) -> str:
        """Quick generation mode for single-partition analyses.

        Generates a BEAST2 XML from simple command-line-style parameters
        without requiring a full YAML configuration file. This is useful
        for rapid prototyping and simple analyses.

        Args:
            alignment: Path to the sequence file (FASTA/NEXUS).
            output: Output XML file path.
            tree_prior: Tree prior type (default: "yule").
            subst_model: Substitution model type (default: "hky").
            clock_model: Clock model type (default: "strict").
            gamma_categories: Number of gamma rate categories (default: 0).
            calibration_yaml: Optional path to a calibration YAML file.
            chain_length: MCMC chain length (default: 10,000,000).
            pre_burnin: MCMC pre-burnin length (default: 0).
            name: Analysis name (default: derived from alignment filename).
            force: Overwrite an existing output file.

        Returns:
            The generated XML string.

        Raises:
            FileNotFoundError: If the alignment file is not found.
            ValueError: If parameters are invalid.

        Example::

            b2p = Beast2Py()
            xml = b2p.quick_generate(
                alignment="primates.fasta",
                output="output.xml",
                tree_prior="yule",
                subst_model="hky",
                clock_model="strict",
            )
        """
        # Read alignment
        aln = SequenceReader.read(
            alignment,
            format="auto",
            data_type="nucleotide",
            alignment_id="alignment",
        )

        # Build substitution model
        subst_type = subst_model.lower()
        if subst_type not in ConfigParser.VALID_SUBSTITUTION_MODELS:
            raise ConfigError(
                f"Invalid substitution model: '{subst_type}' "
                f"(valid: {', '.join(sorted(ConfigParser.VALID_SUBSTITUTION_MODELS))})"
            )
        subst_model_dict = {"type": subst_type}

        gamma_shape = None
        if gamma_categories > 0:
            gamma_shape = RealParameter(value=0.5, lower=0.0)

        site_model = SiteModelConfig(
            substitution_model=subst_model_dict,
            gamma_categories=gamma_categories,
            gamma_shape=gamma_shape,
        )

        clock_type_str = clock_model.lower()
        if clock_type_str not in ConfigParser.VALID_CLOCK_MODELS:
            raise ConfigError(
                f"Invalid clock model: '{clock_type_str}' "
                f"(valid: {', '.join(sorted(ConfigParser.VALID_CLOCK_MODELS))})"
            )
        clock_model_cfg = ClockModelConfig(type=ClockModelType(clock_type_str))

        partition = Partition(
            id="alignment",
            alignment=aln,
            site_model=site_model,
            clock_model=clock_model_cfg,
            tree="shared",
        )

        tp_type_str = tree_prior.lower()
        if tp_type_str not in ConfigParser.VALID_TREE_PRIORS:
            raise ConfigError(
                f"Invalid tree prior: '{tp_type_str}' "
                f"(valid: {', '.join(sorted(ConfigParser.VALID_TREE_PRIORS))})"
            )
        tp_type = TreePriorType(tp_type_str)

        # Parse calibrations from YAML if provided
        calibrations: List[CalibrationPoint] = []
        if calibration_yaml:
            calibrations = _parse_calibration_yaml_file(calibration_yaml)

        mcmc = MCMCConfig(
            chain_length=chain_length,
            pre_burnin=pre_burnin,
        )

        metadata = {
            "analysis_name": name or "quick_analysis",
            "beast2_version": "2.7.8",
            "tool_version": self.version,
            "output_file": output,
        }

        beast_config = BEASTConfig(
            metadata=metadata,
            partitions=[partition],
            tree_prior_type=tp_type,
            tree_prior_params={},
            calibrations=calibrations,
            calibration_method=CalibrationMethod.MRCA_PRIOR,
            mcmc=mcmc,
        )

        from .xml_writer import XMLWriter

        writer = XMLWriter(beast_config)
        xml_str = writer.generate_xml()

        # Same gate as the CLI and generate_xml(): validate first, write only on
        # success.
        from .main import validate_and_write

        status = validate_and_write(xml_str, output, beast_config, force=force)
        if status != 0:
            raise ValueError(f"Validation gates failed for {output}; nothing was written.")

        return xml_str

    # ------------------------------------------------------------------------
    # Calibration Diagnostics
    # ------------------------------------------------------------------------

    def diagnose(
        self,
        config: ConfigInput,
        output_dir: str = "diagnostics",
        run_sensitivity: bool = True,
        run_visualization: bool = True,
        report_path: Optional[str] = None,
    ) -> Any:
        """Run calibration diagnostics on a configuration.

        Performs calibration conflict detection, sensitivity analysis,
        and distribution visualization.

        Args:
            config: Configuration input (file path, dict, or BEASTConfig).
            output_dir: Output directory for diagnostic files.
            run_sensitivity: Whether to generate sensitivity analysis XMLs.
            run_visualization: Whether to generate distribution plots.
            report_path: Optional path for the HTML report. If None,
                no HTML report is generated.

        Returns:
            A DiagnosticReport object containing conflicts, summaries,
            sensitivity XML paths, and visualization file path.

        Raises:
            ConfigError: If the configuration is invalid.

        Example::

            b2p = Beast2Py()
            report = b2p.diagnose(
                "config.yaml",
                output_dir="diagnostics/",
                report_path="report.html",
            )
            for c in report.conflicts:
                print(f"[{c.severity}] {c.description}")
        """
        beast_config = _resolve_config(config, self.base_dir)

        from .diagnostics import DiagnosticEngine, ReportGenerator

        report = DiagnosticEngine.run(
            beast_config,
            output_dir=output_dir,
            run_sensitivity=run_sensitivity,
            run_visualization=run_visualization,
        )

        if report_path:
            ReportGenerator.generate_report(
                config=beast_config,
                conflicts=report.conflicts,
                summaries=report.summaries,
                sensitivity_xmls=report.sensitivity_xmls,
                visualization_file=report.visualization_file,
                output_file=report_path,
            )
            report.report_path = report_path

        return report

    def detect_conflicts(
        self,
        config: ConfigInput,
    ) -> List[Any]:
        """Detect calibration conflicts in a configuration.

        Args:
            config: Configuration input (file path, dict, or BEASTConfig).

        Returns:
            List of Conflict objects, each describing a detected conflict
            with type, severity, description, and recommendation.

        Example::

            b2p = Beast2Py()
            conflicts = b2p.detect_conflicts("config.yaml")
            for c in conflicts:
                print(f"[{c.severity}] {c.conflict_type}: {c.description}")
        """
        beast_config = _resolve_config(config, self.base_dir)

        from .diagnostics import ConflictDetector

        return ConflictDetector.detect_conflicts(
            beast_config.calibrations, settings=getattr(beast_config, "diagnostics", None)
        )

    def generate_sensitivity_xmls(
        self,
        config: ConfigInput,
        output_dir: str = "sensitivity",
    ) -> List[str]:
        """Generate prior-only sampling XMLs for sensitivity analysis.

        Generates N+1 XML files (N = number of calibration points):
        one with all calibrations and one for each leave-one-out scenario.

        Args:
            config: Configuration input (file path, dict, or BEASTConfig).
            output_dir: Output directory for XML files.

        Returns:
            List of generated XML file paths.

        Example::

            b2p = Beast2Py()
            xmls = b2p.generate_sensitivity_xmls("config.yaml")
            print(f"Generated {len(xmls)} sensitivity XML files")
        """
        beast_config = _resolve_config(config, self.base_dir)

        from .diagnostics import SensitivityAnalyzer

        return SensitivityAnalyzer.generate_sensitivity_xmls(beast_config, output_dir)

    def plot_calibrations(
        self,
        config: ConfigInput,
        output_file: str,
    ) -> str:
        """Generate density plots of calibration distributions.

        Args:
            config: Configuration input (file path, dict, or BEASTConfig).
            output_file: Output image file path (PNG/PDF/SVG).

        Returns:
            Path to the generated image file.

        Example::

            b2p = Beast2Py()
            b2p.plot_calibrations("config.yaml", "calibrations.png")
        """
        beast_config = _resolve_config(config, self.base_dir)

        from .diagnostics import VisualizationGenerator

        return VisualizationGenerator.plot_calibration_distributions(
            beast_config.calibrations, output_file
        )

    # ------------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------------

    def validate_xml(
        self,
        xml: Union[str, Path],
        beast2_validate: bool = False,
        beast2_path: str = "beast",
    ) -> Any:
        """Validate a BEAST2 XML file or string.

        Performs structural validation including XML well-formedness,
        required elements, id/idref consistency, and duplicate ID detection.
        Optionally validates with BEAST2's native parser.

        Args:
            xml: XML string or file path to validate.
            beast2_validate: Whether to also validate with BEAST2 -validate.
            beast2_path: Path to BEAST2 executable (default: "beast").

        Returns:
            A ValidationResult object with errors and warnings lists.

        Example::

            b2p = Beast2Py()
            result = b2p.validate_xml("output.xml")
            if result.is_valid:
                print("XML is valid")
            else:
                for err in result.errors:
                    print(f"Error: {err}")
        """
        from .validator import XMLValidator

        # Normalize Path inputs so both file paths and XML strings work
        if isinstance(xml, Path):
            xml = str(xml)

        # Determine if input is a file path or XML string
        if isinstance(xml, str) and not xml.strip().startswith("<"):
            # It's a file path
            result = XMLValidator.validate_file(xml)
            if beast2_validate:
                success, output, status = XMLValidator.beast2_validation_status(xml, beast2_path)
                if status == "unavailable":
                    result.add_warning(
                        "BEAST2 check skipped (not installed or its validation script was "
                        "not found): this XML has only been checked structurally"
                    )
                elif not success:
                    result.add_error(f"BEAST2 validation failed: {output}")
        else:
            # It's an XML string
            result = XMLValidator.validate(xml)
            if beast2_validate:
                # Write to temp file for BEAST2 validation
                import tempfile

                with tempfile.NamedTemporaryFile(mode="w", suffix=".xml", delete=False) as tmp:
                    tmp.write(xml)
                    tmp_path = tmp.name
                try:
                    success, output, status = XMLValidator.beast2_validation_status(
                        tmp_path, beast2_path
                    )
                    if status == "unavailable":
                        result.add_warning(
                            "BEAST2 check skipped (not installed): structural checks only"
                        )
                    elif not success:
                        result.add_error(f"BEAST2 validation failed: {output}")
                finally:
                    os.unlink(tmp_path)

        return result

    # ------------------------------------------------------------------------
    # Reproducibility
    # ------------------------------------------------------------------------

    def generate_fingerprint(
        self,
        config: ConfigInput,
        output: Optional[str] = None,
    ) -> str:
        """Generate an analysis fingerprint for reproducibility tracking.

        The fingerprint is a deterministic identifier based on the SHA-256
        hash of the scientific configuration and the tool version; it carries
        no date component, so reruns over the same data reproduce it exactly.

        Args:
            config: Configuration input (file path, dict, or BEASTConfig).
            output: Optional output JSON file path. If provided, detailed
                fingerprint information is written to this file.

        Returns:
            The fingerprint string (e.g., "B2P-a3f7b2c9d1e8-0.1.1").

        Example::

            b2p = Beast2Py()
            fp = b2p.generate_fingerprint("config.yaml", output="fp.json")
            print(f"Analysis fingerprint: {fp}")
        """
        beast_config = _resolve_config(config, self.base_dir)

        from .reproducibility import FingerprintGenerator

        fingerprint = FingerprintGenerator.generate_fingerprint(beast_config)

        if output:
            FingerprintGenerator.write_fingerprint_file(beast_config, output)

        return fingerprint

    def generate_methods(
        self,
        config: ConfigInput,
        output: Optional[str] = None,
    ) -> str:
        """Generate a LaTeX-formatted methods description.

        The methods description includes substitution model, clock model,
        tree prior, calibration points (with distribution parameters and
        provenance), MCMC settings, and the analysis fingerprint.

        Args:
            config: Configuration input (file path, dict, or BEASTConfig).
            output: Optional output file path. If None, the methods text
                is returned but not written to disk.

        Returns:
            LaTeX-formatted methods description string.

        Example::

            b2p = Beast2Py()
            methods = b2p.generate_methods("config.yaml", output="methods.tex")
        """
        beast_config = _resolve_config(config, self.base_dir)

        from .reproducibility import MethodsGenerator

        methods_text = MethodsGenerator.generate_methods(beast_config)

        if output:
            output_dir = Path(output).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            Path(output).write_text(methods_text, encoding="utf-8")

        return methods_text

    def generate_pipeline(
        self,
        config: ConfigInput,
        pipeline_type: str = "snakemake",
        output_dir: str = ".",
        config_file: Optional[str] = None,
    ) -> Union[str, Tuple[str, str]]:
        """Generate Snakemake or Nextflow pipeline files.

        Args:
            config: Configuration input (file path, dict, or BEASTConfig).
            pipeline_type: Pipeline type ("snakemake", "nextflow", or "both").
            output_dir: Output directory for pipeline files.
            config_file: Path to the YAML config file (for pipeline references).
                Required if config is provided as a dict or BEASTConfig.

        Returns:
            Path to the generated pipeline file, or a tuple of
            (snakemake_path, nextflow_path) if pipeline_type="both".

        Raises:
            ValueError: If pipeline_type is invalid.

        Example::

            b2p = Beast2Py()
            path = b2p.generate_pipeline(
                "config.yaml",
                pipeline_type="snakemake",
                output_dir="pipeline/",
            )
        """
        beast_config = _resolve_config(config, self.base_dir)

        if isinstance(config, str):
            config_file = config_file or config

        from .reproducibility import PipelineGenerator

        os.makedirs(output_dir, exist_ok=True)

        if pipeline_type == "snakemake":
            return PipelineGenerator.generate_snakemake(
                beast_config,
                output_dir=output_dir,
                config_file=config_file or "config.yaml",
            )
        elif pipeline_type == "nextflow":
            return PipelineGenerator.generate_nextflow(
                beast_config,
                output_dir=output_dir,
                config_file=config_file or "config.yaml",
            )
        elif pipeline_type == "both":
            sm_path = PipelineGenerator.generate_snakemake(
                beast_config,
                output_dir=output_dir,
                config_file=config_file or "config.yaml",
            )
            nf_path = PipelineGenerator.generate_nextflow(
                beast_config,
                output_dir=output_dir,
                config_file=config_file or "config.yaml",
            )
            return (sm_path, nf_path)
        else:
            raise ValueError(
                f"Invalid pipeline_type: {pipeline_type}. "
                "Use 'snakemake', 'nextflow', or 'both'."
            )

    # ------------------------------------------------------------------------
    # Configuration and Models
    # ------------------------------------------------------------------------

    def parse_config(
        self,
        config: Union[str, Dict[str, Any]],
    ) -> BEASTConfig:
        """Parse a YAML configuration file or dictionary.

        Args:
            config: Path to YAML config file, or a configuration dictionary.

        Returns:
            A BEASTConfig object.

        Raises:
            ConfigError: If the configuration is invalid.

        Example::

            b2p = Beast2Py()
            cfg = b2p.parse_config("config.yaml")
            print(f"Partitions: {len(cfg.partitions)}")
        """
        return _resolve_config(config, self.base_dir)

    def list_models(self) -> Dict[str, List[str]]:
        """List all supported models.

        Returns:
            Dictionary with keys 'substitution_models', 'clock_models',
            'tree_priors', 'calibration_distributions', and 'data_types',
            each containing a sorted list of supported model names.

        Example::

            b2p = Beast2Py()
            models = b2p.list_models()
            print(models["substitution_models"])
        """
        return {
            "substitution_models": sorted(ConfigParser.VALID_SUBSTITUTION_MODELS),
            "clock_models": sorted(ConfigParser.VALID_CLOCK_MODELS),
            "tree_priors": sorted(ConfigParser.VALID_TREE_PRIORS),
            "calibration_distributions": sorted(ConfigParser.VALID_DISTRIBUTIONS),
            "data_types": sorted(ConfigParser.VALID_DATA_TYPES),
        }

    # ------------------------------------------------------------------------
    # Sequence Reading
    # ------------------------------------------------------------------------

    def read_sequence(
        self,
        file_path: str,
        format: str = "auto",
        data_type: str = "nucleotide",
        alignment_id: str = "alignment",
    ) -> Any:
        """Read a sequence alignment file.

        Args:
            file_path: Path to the sequence file (FASTA/NEXUS).
            format: File format ("fasta", "nexus", or "auto").
            data_type: "nucleotide" or "aminoacid".
            alignment_id: ID for the alignment.

        Returns:
            An Alignment object.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the format is unsupported or parsing fails.

        Example::

            b2p = Beast2Py()
            aln = b2p.read_sequence("primates.fasta")
            print(f"Taxa: {aln.n_taxa}, Sites: {aln.n_sites}")
        """
        return SequenceReader.read(
            file_path,
            format=format,
            data_type=data_type,
            alignment_id=alignment_id,
        )

    # ------------------------------------------------------------------------
    # All-in-one: Generate XML with optional diagnostics and reproducibility
    # ------------------------------------------------------------------------

    def generate(
        self,
        config: ConfigInput,
        output: str,
        diagnose: bool = False,
        fingerprint: bool = False,
        methods: bool = False,
        beast2_validate: bool = False,
        beast2_path: str = "beast",
        force: bool = False,
    ) -> Dict[str, Any]:
        """Generate BEAST2 XML with optional diagnostics and reproducibility artifacts.

        This is the API equivalent of the CLI ``generate`` command, combining
        XML generation, structural validation, optional BEAST2 validation,
        diagnostics, fingerprint, and methods description in one call.

        Args:
            config: Configuration input (file path, dict, or BEASTConfig).
            output: Output XML file path.
            diagnose: Whether to also run calibration diagnostics.
            fingerprint: Whether to generate an analysis fingerprint file.
            methods: Whether to generate a LaTeX methods description.
            beast2_validate: Whether to validate with BEAST2 -validate.
            beast2_path: Path to BEAST2 executable.

        Returns:
            Dictionary with keys:
                - 'xml_path': Path to generated XML.
                - 'validation': ValidationResult object.
                - 'fingerprint': Fingerprint string (if fingerprint=True).
                - 'fingerprint_path': Fingerprint file path (if fingerprint=True).
                - 'methods_path': Methods file path (if methods=True).
                - 'diagnostic_report': DiagnosticReport (if diagnose=True).

        Example::

            b2p = Beast2Py()
            result = b2p.generate(
                "config.yaml",
                output="output.xml",
                diagnose=True,
                fingerprint=True,
                methods=True,
            )
            print(f"XML: {result['xml_path']}")
            print(f"Valid: {result['validation'].is_valid}")
        """
        results: Dict[str, Any] = {}

        # Generate XML — generate_xml runs the same gates as the CLI and raises
        # without writing if any of them fails.
        xml_str = self.generate_xml(config, output=output, force=force)
        results["xml_path"] = output

        # Validate
        from .validator import XMLValidator

        result = XMLValidator.validate(xml_str)
        results["validation"] = result

        if beast2_validate:
            success, output_msg, status = XMLValidator.beast2_validation_status(
                output, beast2_path=beast2_path
            )
            results["beast2_validation"] = {
                "success": success,
                "status": status,
                "output": output_msg,
            }

        # Fingerprint
        if fingerprint:
            fp_path = str(Path(output).with_suffix(".fingerprint.json"))
            fp = self.generate_fingerprint(config, output=fp_path)
            results["fingerprint"] = fp
            results["fingerprint_path"] = fp_path

        # Methods
        if methods:
            methods_path = str(Path(output).with_suffix(".methods.tex"))
            self.generate_methods(config, output=methods_path)
            results["methods_path"] = methods_path

        # Diagnostics
        if diagnose:
            diag_dir = str(Path(output).parent / "diagnostics")
            report = self.diagnose(
                config,
                output_dir=diag_dir,
                run_sensitivity=True,
                run_visualization=True,
            )
            results["diagnostic_report"] = report

        return results


# ============================================================================
# Module-level convenience functions
# ============================================================================


def generate_xml(
    config: ConfigInput,
    output: Optional[str] = None,
    base_dir: Optional[str] = None,
    force: bool = False,
) -> str:
    """Generate BEAST2 XML from a configuration.

    Args:
        config: Configuration input (file path, dict, or BEASTConfig).
        output: Optional output file path.
        base_dir: Base directory for resolving relative paths.
        force: Overwrite an existing output file.

    Returns:
        The generated XML string.

    Example::

        from beast2py.api import generate_xml
        xml = generate_xml("config.yaml", output="output.xml")
    """
    api = Beast2Py(base_dir=base_dir)
    return api.generate_xml(config, output=output, force=force)


def quick_generate(
    alignment: str,
    output: str,
    tree_prior: str = "yule",
    subst_model: str = "hky",
    clock_model: str = "strict",
    gamma_categories: int = 0,
    calibration_yaml: Optional[str] = None,
    chain_length: int = 10_000_000,
    pre_burnin: int = 0,
    name: Optional[str] = None,
    force: bool = False,
) -> str:
    """Quick generation mode for single-partition analyses.

    Args:
        alignment: Path to the sequence file (FASTA/NEXUS).
        output: Output XML file path.
        tree_prior: Tree prior type (default: "yule").
        subst_model: Substitution model type (default: "hky").
        clock_model: Clock model type (default: "strict").
        gamma_categories: Number of gamma rate categories (default: 0).
        calibration_yaml: Optional path to a calibration YAML file.
        chain_length: MCMC chain length (default: 10,000,000).
        pre_burnin: MCMC pre-burnin length (default: 0).
        name: Analysis name.
        force: Overwrite an existing output file.

    Returns:
        The generated XML string.

    Raises:
        ValueError: If a validation gate failed, or the output exists and
            ``force`` was not set. Nothing is written in that case.
    """
    api = Beast2Py()
    return api.quick_generate(
        alignment=alignment,
        output=output,
        tree_prior=tree_prior,
        subst_model=subst_model,
        clock_model=clock_model,
        gamma_categories=gamma_categories,
        calibration_yaml=calibration_yaml,
        chain_length=chain_length,
        pre_burnin=pre_burnin,
        name=name,
        force=force,
    )


def diagnose(
    config: ConfigInput,
    output_dir: str = "diagnostics",
    run_sensitivity: bool = True,
    run_visualization: bool = True,
    report_path: Optional[str] = None,
    base_dir: Optional[str] = None,
) -> Any:
    """Run calibration diagnostics on a configuration.

    Args:
        config: Configuration input (file path, dict, or BEASTConfig).
        output_dir: Output directory for diagnostic files.
        run_sensitivity: Whether to generate sensitivity analysis XMLs.
        run_visualization: Whether to generate distribution plots.
        report_path: Optional path for the HTML report.
        base_dir: Base directory for resolving relative paths.

    Returns:
        A DiagnosticReport object.
    """
    api = Beast2Py(base_dir=base_dir)
    return api.diagnose(
        config,
        output_dir=output_dir,
        run_sensitivity=run_sensitivity,
        run_visualization=run_visualization,
        report_path=report_path,
    )


def validate_xml(
    xml: Union[str, Path],
    beast2_validate: bool = False,
    beast2_path: str = "beast",
) -> Any:
    """Validate a BEAST2 XML file or string.

    Args:
        xml: XML string or file path.
        beast2_validate: Whether to also validate with BEAST2 -validate.
        beast2_path: Path to BEAST2 executable.

    Returns:
        A ValidationResult object.
    """
    api = Beast2Py()
    return api.validate_xml(
        xml,
        beast2_validate=beast2_validate,
        beast2_path=beast2_path,
    )


def generate_fingerprint(
    config: ConfigInput,
    output: Optional[str] = None,
    base_dir: Optional[str] = None,
) -> str:
    """Generate an analysis fingerprint.

    Args:
        config: Configuration input (file path, dict, or BEASTConfig).
        output: Optional output JSON file path.
        base_dir: Base directory for resolving relative paths.

    Returns:
        The fingerprint string.
    """
    api = Beast2Py(base_dir=base_dir)
    return api.generate_fingerprint(config, output=output)


def generate_methods(
    config: ConfigInput,
    output: Optional[str] = None,
    base_dir: Optional[str] = None,
) -> str:
    """Generate a LaTeX methods description.

    Args:
        config: Configuration input (file path, dict, or BEASTConfig).
        output: Optional output file path.
        base_dir: Base directory for resolving relative paths.

    Returns:
        LaTeX-formatted methods description string.
    """
    api = Beast2Py(base_dir=base_dir)
    return api.generate_methods(config, output=output)


def generate_pipeline(
    config: ConfigInput,
    pipeline_type: str = "snakemake",
    output_dir: str = ".",
    config_file: Optional[str] = None,
    base_dir: Optional[str] = None,
) -> Union[str, Tuple[str, str]]:
    """Generate Snakemake or Nextflow pipeline files.

    Args:
        config: Configuration input (file path, dict, or BEASTConfig).
        pipeline_type: "snakemake", "nextflow", or "both".
        output_dir: Output directory for pipeline files.
        config_file: Path to the YAML config file (for pipeline references).
        base_dir: Base directory for resolving relative paths.

    Returns:
        Path to the generated pipeline file(s).
    """
    api = Beast2Py(base_dir=base_dir)
    return api.generate_pipeline(
        config,
        pipeline_type=pipeline_type,
        output_dir=output_dir,
        config_file=config_file,
    )


def list_models() -> Dict[str, List[str]]:
    """List all supported models.

    Returns:
        Dictionary of model categories and their supported types.
    """
    api = Beast2Py()
    return api.list_models()


def parse_config(
    config: Union[str, Dict[str, Any]],
    base_dir: Optional[str] = None,
) -> BEASTConfig:
    """Parse a YAML configuration file or dictionary.

    Args:
        config: Path to YAML config file, or a configuration dictionary.
        base_dir: Base directory for resolving relative paths.

    Returns:
        A BEASTConfig object.
    """
    api = Beast2Py(base_dir=base_dir)
    return api.parse_config(config)


def read_sequence(
    file_path: str,
    format: str = "auto",
    data_type: str = "nucleotide",
    alignment_id: str = "alignment",
) -> Any:
    """Read a sequence alignment file.

    Args:
        file_path: Path to the sequence file (FASTA/NEXUS).
        format: File format ("fasta", "nexus", or "auto").
        data_type: "nucleotide" or "aminoacid".
        alignment_id: ID for the alignment.

    Returns:
        An Alignment object.
    """
    api = Beast2Py()
    return api.read_sequence(
        file_path,
        format=format,
        data_type=data_type,
        alignment_id=alignment_id,
    )


# ============================================================================
# Internal helpers
# ============================================================================


def _parse_calibration_yaml_file(yaml_path: str) -> List[CalibrationPoint]:
    """Parse a calibration-only YAML file.

    Args:
        yaml_path: Path to the calibration YAML file.

    Returns:
        List of CalibrationPoint objects.
    """
    import yaml

    with open(yaml_path, "r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    if not isinstance(raw, list):
        raise ValueError("Calibration YAML must be a list of calibration points")

    calibrations: List[CalibrationPoint] = []
    for cal_raw in raw:
        name = cal_raw.get("name", "")
        taxa = cal_raw.get("taxa")
        where = f"calibrations[{len(calibrations)}]"
        # Strict coercion: this shortcut path used to read these flags raw, so
        # `monophyletic: "false"` inverted the user's intent through Python
        # truthiness even though the YAML-config path rejected it.
        monophyletic = ConfigParser._as_bool(cal_raw, "monophyletic", True, where)
        use_originate = ConfigParser._as_bool(cal_raw, "use_originate", False, where)
        tipsonly = ConfigParser._as_bool(cal_raw, "tipsonly", False, where)

        dist = None
        dist_raw = cal_raw.get("distribution")
        if dist_raw:
            dist = DistributionConfig(
                type=dist_raw.get("type", "normal"),
                parameters=dist_raw.get("parameters", {}),
                offset=float(dist_raw.get("offset", 0.0)),
            )
            # Same range and support checks as the full configuration path
            # (: an inverted `uniform` was accepted silently).
            ConfigParser()._validate_distribution(dist, name)

        provenance = None
        prov_raw = cal_raw.get("provenance")
        if prov_raw:
            provenance = Provenance(
                source=prov_raw.get("source", ""),
                reference=prov_raw.get("reference", ""),
                calibration_type=prov_raw.get("calibration_type", "soft"),
                original_age_min=prov_raw.get("original_age_min"),
                original_age_max=prov_raw.get("original_age_max"),
                notes=prov_raw.get("notes", ""),
            )

        calibrations.append(
            CalibrationPoint(
                name=name,
                taxa=taxa,
                monophyletic=monophyletic,
                distribution=dist,
                use_originate=use_originate,
                tipsonly=tipsonly,
                provenance=provenance,
            )
        )

    return calibrations

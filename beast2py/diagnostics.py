"""Calibration diagnostics engine for Beast2Py.

Core innovation module containing:
- ConflictDetector: Detects temporal, distribution, and monophyly conflicts between calibrations
- SensitivityAnalyzer: Generates prior-only sampling XMLs for sensitivity analysis
- VisualizationGenerator: Creates density plots of calibration distributions
- ReportGenerator: Generates HTML diagnostic reports
"""

from __future__ import annotations

import copy
import html
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from . import __version__
from .models import BEASTConfig, CalibrationPoint, DistributionConfig
from .utils import compute_distribution_stats, compute_wasserstein_distance

# ============================================================================
# Data Classes for Diagnostic Results
# ============================================================================


@dataclass
class Conflict:
    """A detected calibration conflict."""

    conflict_type: str  # temporal, distribution_overlap, monophyly, range
    severity: str  # error, warning, info
    calibration_a: str  # Name of first calibration
    calibration_b: str  # Name of second calibration
    description: str
    recommendation: str = ""


@dataclass
class CalibrationSummary:
    """Summary of a calibration point for reporting."""

    name: str
    taxa: Optional[List[str]]
    taxa_count: int
    distribution_type: str
    distribution_params: Dict[str, Any]
    offset: float
    stats: Dict[str, float]  # mean, median, hpd_lower, hpd_upper
    provenance_source: str = ""
    provenance_reference: str = ""
    calibration_type: str = ""


@dataclass
class DiagnosticReport:
    """Complete diagnostic report."""

    config: BEASTConfig
    conflicts: List[Conflict]
    summaries: List[CalibrationSummary]
    sensitivity_xmls: List[str] = field(default_factory=list)
    visualization_file: Optional[str] = None
    report_path: Optional[str] = None
    visualization_error: Optional[str] = None


# ============================================================================
# Conflict Detector
# ============================================================================


class ConflictDetector:
    """Detect conflicts between calibration points."""

    # Threshold for distribution overlap warning (Wasserstein distance, in the
    # analysis' time units). Applied only when ``overlap_measure`` asks for
    # it; the relative criterion below is the default because an absolute
    # cut-off is scale-dependent.
    OVERLAP_DISTANCE_THRESHOLD = 2.0
    # Dimensionless threshold: W divided by the mean width of the two priors'
    # 95% intervals. Below ~0.15 the two priors are nearly interchangeable.
    OVERLAP_RATIO_THRESHOLD = 0.15

    @staticmethod
    def detect_conflicts(
        calibrations: List[CalibrationPoint],
        overlap_threshold: Optional[float] = None,
        relative_threshold: Optional[float] = None,
        settings: Optional[Dict[str, Any]] = None,
    ) -> List[Conflict]:
        """Detect all types of conflicts between calibration points.

        Args:
            calibrations: List of calibration points.
            overlap_threshold: Absolute Wasserstein threshold in time units.
                ``None`` means "use the class default only if
                ``settings['overlap_measure']`` asks for the absolute
                criterion"; pass a number to force it.
            relative_threshold: Dimensionless threshold, or None for the class
                default.
            settings: The configuration's ``diagnostics`` block, e.g.
                ``{"overlap_measure": "relative",
                   "overlap_threshold": 2.0,
                   "relative_overlap_threshold": 0.15}``.

        Returns:
            List of detected conflicts.
        """
        conflicts: List[Conflict] = []

        if len(calibrations) < 2:
            return conflicts

        settings = settings or {}
        measure = str(settings.get("overlap_measure", "relative")).lower()
        if measure not in ("relative", "absolute", "both"):
            raise ValueError(
                "diagnostics.overlap_measure must be relative, absolute or both, "
                f"got {measure!r}"
            )
        absolute = settings.get("overlap_threshold", ConflictDetector.OVERLAP_DISTANCE_THRESHOLD)
        relative = settings.get(
            "relative_overlap_threshold", ConflictDetector.OVERLAP_RATIO_THRESHOLD
        )
        if overlap_threshold is not None:
            absolute = overlap_threshold
        if relative_threshold is not None:
            relative = relative_threshold

        abs_arg = absolute if measure in ("absolute", "both") else None
        rel_arg = relative if measure in ("relative", "both") else None

        conflicts.extend(ConflictDetector._check_temporal_consistency(calibrations))
        conflicts.extend(
            ConflictDetector._check_distribution_overlap(
                calibrations,
                overlap_threshold=abs_arg,
                relative_threshold=rel_arg,
            )
        )
        conflicts.extend(ConflictDetector._check_monophyly_conflicts(calibrations))

        return conflicts

    @staticmethod
    def _compare_parent_child(parent: CalibrationPoint, child: CalibrationPoint) -> List[Conflict]:
        """Compare one calibration that must be older than another.

        Args:
            parent: The clade that nests ``child`` (or the root).
            child: The nested clade.

        Returns:
            Zero, one or two conflicts describing hard and soft violations.
        """
        conflicts: List[Conflict] = []
        stats_a = ConflictDetector._get_stats(parent)
        stats_b = ConflictDetector._get_stats(child)
        if not stats_a or not stats_b:
            # Reporting nothing here let a pair that was never *evaluated* look
            # identical to a pair that was evaluated and found consistent, so
            # `one_on_x` calibrations -- whose quantiles are unbounded --
            # produced a green "No conflicts detected.".
            unclear = [c.name for c, s in ((parent, stats_a), (child, stats_b)) if not s]
            return [
                Conflict(
                    conflict_type="unchecked",
                    severity="info",
                    calibration_a=parent.name,
                    calibration_b=child.name,
                    description=(
                        f"Not checked: {', '.join(unclear)} uses a distribution "
                        f"whose mean or quantiles are undefined, so this pair "
                        f"could not be compared."
                    ),
                    recommendation=(
                        "Treat the absence of a conflict below as covering only "
                        "the pairs that were compared, or switch to a prior with "
                        "finite quantiles."
                    ),
                )
            ]

        if stats_a["hpd_upper"] < stats_b["hpd_lower"]:
            conflicts.append(
                Conflict(
                    conflict_type="temporal",
                    severity="error",
                    calibration_a=parent.name,
                    calibration_b=child.name,
                    description=(
                        f"Temporal inconsistency: '{parent.name}' (parent clade) "
                        f"has HPD upper {stats_a['hpd_upper']:.2f}, younger than "
                        f"'{child.name}' (nested clade) HPD lower "
                        f"{stats_b['hpd_lower']:.2f}. "
                        f"Parent clade MRCA must be older than nested clade MRCA."
                    ),
                    recommendation=(
                        f"Adjust calibration distributions so that '{parent.name}' "
                        f"spans older times than '{child.name}'."
                    ),
                )
            )
        elif stats_a["mean"] < stats_b["mean"]:
            conflicts.append(
                Conflict(
                    conflict_type="temporal",
                    severity="warning",
                    calibration_a=parent.name,
                    calibration_b=child.name,
                    description=(
                        f"Potential temporal issue: '{parent.name}' (parent clade) "
                        f"mean {stats_a['mean']:.2f} is younger than "
                        f"'{child.name}' (nested clade) mean "
                        f"{stats_b['mean']:.2f}. Parent should generally be "
                        f"older, though HPD intervals may overlap."
                    ),
                    recommendation=("Review calibration distributions for consistency."),
                )
            )
        return conflicts

    @staticmethod
    def _check_temporal_consistency(
        calibrations: List[CalibrationPoint],
    ) -> List[Conflict]:
        """Check temporal consistency between nested clades.

        If clade A contains clade B (A.taxa ⊃ B.taxa), then A's MRCA time
        should be >= B's MRCA time. A root calibration (``taxa is None``) nests
        every other clade, so it is compared as the parent of each of them:
        without those comparisons, "root younger than its own descendant" — the
        most severe violation available — would be the one thing the detector
        could never report.
        """
        conflicts: List[Conflict] = []

        for i, cal_a in enumerate(calibrations):
            for j, cal_b in enumerate(calibrations):
                if i >= j:
                    continue

                # Check taxon set containment, with the root as universal parent
                if cal_a.taxa is None or cal_b.taxa is None:
                    if cal_a.taxa is None and cal_b.taxa is None:
                        continue  # two roots: nothing to order
                    root_cal = cal_a if cal_a.taxa is None else cal_b
                    other = cal_b if cal_a.taxa is None else cal_a
                    if other.tipsonly:
                        continue  # tip dates are not clade ages
                    conflicts.extend(ConflictDetector._compare_parent_child(root_cal, other))
                    continue

                set_a = set(cal_a.taxa)
                set_b = set(cal_b.taxa)

                # If A contains B (A is parent/ancestor clade)
                if set_b.issubset(set_a) and set_b != set_a:
                    conflicts.extend(ConflictDetector._compare_parent_child(cal_a, cal_b))

                # If B contains A (B is parent/ancestor clade)
                elif set_a.issubset(set_b) and set_a != set_b:
                    conflicts.extend(ConflictDetector._compare_parent_child(cal_b, cal_a))

        return conflicts

    @staticmethod
    def _check_distribution_overlap(
        calibrations: List[CalibrationPoint],
        overlap_threshold: Optional[float] = None,
        relative_threshold: Optional[float] = None,
    ) -> List[Conflict]:
        """Check for excessive distribution overlap between calibrations.

        Two criteria are available:

        * absolute — the 1-Wasserstein distance below a fixed number of time
          units (opt-in through ``overlap_measure``, ``OVERLAP_DISTANCE_THRESHOLD``);
        * relative — that distance divided by the mean width of the two
          priors' 95% intervals, which is dimensionless.

        The absolute criterion is scale-dependent: W between two priors grows
        with node depth, so at a fixed threshold shallow nodes are flagged as
        redundant almost automatically while deep ones never are, and the alarm
        rate ends up measuring phylogenetic position rather than information
        redundancy. The relative criterion is therefore
        the default; both can be reported at once.

        Args:
            calibrations: Calibration points to compare pairwise.
            overlap_threshold: Absolute threshold in time units, or None to
                skip the absolute criterion.
            relative_threshold: Dimensionless threshold, or None to skip it.

        Returns:
            List of detected conflicts.
        """
        conflicts: List[Conflict] = []

        for i, cal_a in enumerate(calibrations):
            for j, cal_b in enumerate(calibrations):
                if i >= j:
                    continue

                if not cal_a.distribution or not cal_b.distribution:
                    continue

                # Skip pairs involving a root calibration (taxa=None): the root
                # clade nests every other clade, which is handled by the
                # temporal consistency check instead.
                if cal_a.taxa is None or cal_b.taxa is None:
                    continue

                # Skip nested clades (already checked in temporal consistency)
                set_a = set(cal_a.taxa)
                set_b = set(cal_b.taxa)
                if set_a.issubset(set_b) or set_b.issubset(set_a):
                    continue

                # Compute Wasserstein distance
                try:
                    distance = compute_wasserstein_distance(
                        cal_a.distribution.type,
                        cal_a.distribution.parameters,
                        cal_a.distribution.offset,
                        cal_b.distribution.type,
                        cal_b.distribution.parameters,
                        cal_b.distribution.offset,
                    )
                except ValueError:
                    # Unsupported distribution (e.g. OneOnX): overlap undefined
                    continue

                stats_a = ConflictDetector._get_stats(cal_a)
                stats_b = ConflictDetector._get_stats(cal_b)
                pooled_scale = None
                if stats_a and stats_b:
                    width_a = stats_a["hpd_upper"] - stats_a["hpd_lower"]
                    width_b = stats_b["hpd_upper"] - stats_b["hpd_lower"]
                    pooled_scale = (width_a + width_b) / 2.0
                ratio = distance / pooled_scale if pooled_scale and pooled_scale > 0 else None

                flagged = False
                reasons: List[str] = []
                if (
                    relative_threshold is not None
                    and ratio is not None
                    and ratio < relative_threshold
                ):
                    flagged = True
                    reasons.append(
                        f"relative Wasserstein distance = {ratio:.3f} of the pooled "
                        f"95%-interval width (threshold {relative_threshold})"
                    )
                if overlap_threshold is not None and distance < overlap_threshold:
                    flagged = True
                    reasons.append(
                        f"Wasserstein distance = {distance:.3f} below the absolute "
                        f"threshold {overlap_threshold}"
                    )

                if flagged:
                    conflicts.append(
                        Conflict(
                            conflict_type="distribution_overlap",
                            severity="warning",
                            calibration_a=cal_a.name,
                            calibration_b=cal_b.name,
                            description=(
                                f"Distributions of '{cal_a.name}' and '{cal_b.name}' "
                                f"are very similar ({'; '.join(reasons)}). "
                                f"If these calibrate different nodes, consider whether they "
                                f"provide independent information."
                            ),
                            recommendation=(
                                "Consider whether overlapping calibrations provide "
                                "redundant information or if they could be merged."
                            ),
                        )
                    )

        return conflicts

    @staticmethod
    def _check_monophyly_conflicts(
        calibrations: List[CalibrationPoint],
    ) -> List[Conflict]:
        """Check for monophyly conflicts (overlapping but non-nested taxon sets).

        Two calibrations with intersecting but non-nested taxon sets cannot
        both be monophyletic on the same tree.
        """
        conflicts: List[Conflict] = []

        for i, cal_a in enumerate(calibrations):
            for j, cal_b in enumerate(calibrations):
                if i >= j:
                    continue

                if not cal_a.taxa or not cal_b.taxa:
                    continue
                if not cal_a.monophyletic or not cal_b.monophyletic:
                    continue

                set_a = set(cal_a.taxa)
                set_b = set(cal_b.taxa)

                intersection = set_a & set_b
                if intersection and not set_a.issubset(set_b) and not set_b.issubset(set_a):
                    conflicts.append(
                        Conflict(
                            conflict_type="monophyly",
                            severity="error",
                            calibration_a=cal_a.name,
                            calibration_b=cal_b.name,
                            description=(
                                f"Monophyly conflict: '{cal_a.name}' and '{cal_b.name}' have "
                                f"overlapping but non-nested taxon sets "
                                f"(shared: {', '.join(sorted(intersection))}). "
                                f"Both cannot be monophyletic on the same tree."
                            ),
                            recommendation=(
                                "Ensure calibration taxon sets are either disjoint or nested, "
                                "or relax monophyly constraint for one of the calibrations."
                            ),
                        )
                    )

        return conflicts

    @staticmethod
    def _get_stats(cal: CalibrationPoint) -> Optional[Dict[str, float]]:
        """Get distribution statistics for a calibration point."""
        if not cal.distribution:
            return None
        try:
            return compute_distribution_stats(
                cal.distribution.type,
                cal.distribution.parameters,
                cal.distribution.offset,
            )
        except ValueError:
            # Unsupported distribution type: skip time-consistency checking
            # for this calibration pair rather than aborting the analysis.
            return None


# ============================================================================
# Sensitivity Analyzer
# ============================================================================


class SensitivityAnalyzer:
    """Analyze the sensitivity of posterior to calibration priors."""

    @staticmethod
    def generate_sensitivity_xmls(
        config: BEASTConfig,
        output_dir: str,
    ) -> List[str]:
        """Generate sensitivity analysis XML files.

        Generates:
        1. Full calibration + prior-only sampling XML
        2. Leave-one-out: remove each calibration + prior-only sampling XML

        Args:
            config: BEAST configuration.
            output_dir: Output directory for XML files.

        Returns:
            List of generated XML file paths.
        """
        os.makedirs(output_dir, exist_ok=True)
        xml_files: List[str] = []

        # 1. Full calibration + prior-only
        xml_path = SensitivityAnalyzer._generate_prior_only_xml(
            config, output_dir, suffix="00_full_calibration"
        )
        xml_files.append(xml_path)

        # 2. Leave-one-out
        for i, cal in enumerate(config.calibrations):
            reduced_config = SensitivityAnalyzer._remove_calibration(config, i)
            xml_path = SensitivityAnalyzer._generate_prior_only_xml(
                reduced_config,
                output_dir,
                # Index prefix keeps names unique even when calibration names
                # collapse to the same sanitized identifier (e.g. "a b"/"a-b")
                suffix=f"{i + 1:02d}_without_{cal.name}",
            )
            xml_files.append(xml_path)

        if len(set(xml_files)) != len(xml_files):
            raise ValueError(
                "Sensitivity analysis generated duplicate file names; "
                "calibration names must differ after normalization"
            )
        return xml_files

    @staticmethod
    def _remove_calibration(
        config: BEASTConfig,
        index: int,
    ) -> BEASTConfig:
        """Create a copy of config with one calibration removed.

        Args:
            config: Original config.
            index: Index of calibration to remove.

        Returns:
            New config with the calibration removed.
        """
        new_config = copy.deepcopy(config)
        new_config.calibrations.pop(index)
        return new_config

    @staticmethod
    def _generate_prior_only_xml(
        config: BEASTConfig,
        output_dir: str,
        suffix: str = "",
    ) -> str:
        """Generate a prior-only sampling XML.

        Args:
            config: BEAST configuration (modified to sample from prior).
            output_dir: Output directory.
            suffix: Suffix for the filename.

        Returns:
            Path to the generated XML file.
        """
        # Create a copy with sample_from_prior=True
        prior_config = copy.deepcopy(config)
        prior_config.mcmc.sample_from_prior = True
        prior_config.mcmc.chain_length = min(config.mcmc.chain_length, 1000000)

        # Generate XML
        from .xml_writer import XMLWriter

        writer = XMLWriter(prior_config)
        xml_str = writer.generate_xml()

        # Write to file
        from .utils import sanitize_id

        def _safe_name(value: str) -> str:
            # sanitize_id already maps path separators to "_"; collapsing dot
            # sequences as well makes the traversal-safety explicit.
            return sanitize_id(value).replace("..", "_").strip("._")

        analysis_name = _safe_name(config.metadata.get("analysis_name", "analysis"))
        filename = f"{analysis_name}_sensitivity_{_safe_name(suffix)}.xml"
        filepath = os.path.join(output_dir, filename)
        Path(filepath).write_text(xml_str, encoding="utf-8")

        return filepath


# ============================================================================
# Visualization Generator
# ============================================================================


class VisualizationGenerator:
    """Generate visualizations of calibration distributions."""

    @staticmethod
    def plot_calibration_distributions(
        calibrations: List[CalibrationPoint],
        output_file: str,
    ) -> str:
        """Generate a density plot of all calibration distributions.

        Args:
            calibrations: List of calibration points.
            output_file: Output file path (PNG/PDF/SVG).

        Returns:
            Path to the generated file.
        """
        import matplotlib

        matplotlib.use("Agg")  # Non-interactive backend
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))

        colors = plt.cm.tab10(np.linspace(0, 1, max(len(calibrations), 1)))

        all_x_min = float("inf")
        all_x_max = float("-inf")

        for idx, cal in enumerate(calibrations):
            if not cal.distribution:
                continue

            try:
                stats = compute_distribution_stats(
                    cal.distribution.type,
                    cal.distribution.parameters,
                    cal.distribution.offset,
                )
            except ValueError:
                # Undefined statistics (e.g. OneOnX): cannot plot a range
                continue

            # Determine x range
            x_min = max(0, stats["hpd_lower"] - 0.5 * (stats["hpd_upper"] - stats["hpd_lower"]))
            x_max = stats["hpd_upper"] + 0.5 * (stats["hpd_upper"] - stats["hpd_lower"])

            all_x_min = min(all_x_min, x_min)
            all_x_max = max(all_x_max, x_max)

            x = np.linspace(x_min, x_max, 500)
            density = VisualizationGenerator._compute_density(cal.distribution, x)

            if density is None:
                continue

            color = colors[idx % len(colors)]
            ax.plot(x, density, label=cal.name, color=color, linewidth=2)
            ax.fill_between(x, density, alpha=0.15, color=color)

            # Mark mean
            ax.axvline(stats["mean"], color=color, linestyle="--", alpha=0.5, linewidth=1)

        ax.set_xlabel("Time (million years)", fontsize=12)
        ax.set_ylabel("Density", fontsize=12)
        ax.set_title("Calibration Prior Distributions", fontsize=14)
        ax.legend(fontsize=10, loc="best")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        fig.savefig(output_file, dpi=150, bbox_inches="tight")
        plt.close(fig)

        return output_file

    @staticmethod
    def _compute_density(
        dist: DistributionConfig,
        x: np.ndarray,
    ) -> Optional[np.ndarray]:
        """Compute probability density at points x for a distribution.

        Args:
            dist: Distribution configuration.
            x: Array of x values.

        Returns:
            Array of density values, or None if the density is undefined
            (discrete or improper distributions).
        """
        from .utils import _frozen_distribution

        try:
            dist_obj = _frozen_distribution(dist.type, dist.parameters, 0.0)
        except ValueError:
            return None

        if dist.type.lower() == "poisson":
            # Discrete distribution: no smooth density to draw
            return None

        return dist_obj.pdf(x - dist.offset)


# ============================================================================
# Report Generator
# ============================================================================


class ReportGenerator:
    """Generate HTML diagnostic reports."""

    @staticmethod
    def generate_report(
        config: BEASTConfig,
        conflicts: List[Conflict],
        summaries: List[CalibrationSummary],
        sensitivity_xmls: List[str],
        visualization_file: Optional[str],
        output_file: str,
    ) -> str:
        """Generate an HTML diagnostic report.

        Args:
            config: BEAST configuration.
            conflicts: List of detected conflicts.
            summaries: List of calibration summaries.
            sensitivity_xmls: List of sensitivity analysis XML paths.
            visualization_file: Path to visualization image.
            output_file: Output HTML file path.

        Returns:
            Path to the generated report.
        """
        # Build calibration summary table
        summary_rows = []
        for s in summaries:
            taxa_str = ", ".join(s.taxa) if s.taxa else "All (root)"
            params_str = ", ".join(f"{k}={v}" for k, v in s.distribution_params.items())
            if "warning" in s.stats:
                stats_cells = f'<td colspan="2">⚠ {html.escape(s.stats["warning"])}</td>'
            else:
                interval_kind = s.stats.get("interval_type", "hpd")
                kind_label = "95% HPD" if interval_kind == "hpd" else "95% central"
                stats_cells = (
                    f"<td>{s.stats.get('mean', 0):.2f}</td>"
                    f"<td>{kind_label}: [{s.stats.get('hpd_lower', 0):.2f}, "
                    f"{s.stats.get('hpd_upper', 0):.2f}]</td>"
                )
            summary_rows.append(f"""
                <tr>
                    <td>{html.escape(s.name)}</td>
                    <td>{html.escape(taxa_str)}</td>
                    <td>{s.taxa_count}</td>
                    <td>{html.escape(s.distribution_type)}</td>
                    <td>{html.escape(params_str)}</td>
                    {stats_cells}
                    <td>{html.escape(s.provenance_source)}</td>
                    <td>{html.escape(s.calibration_type)}</td>
                </tr>
            """)

        # Build conflict table
        conflict_rows = []
        for c in conflicts:
            severity_class = "conflict-error" if c.severity == "error" else "conflict-warning"
            conflict_rows.append(f"""
                <tr class="{severity_class}">
                    <td>{html.escape(c.conflict_type)}</td>
                    <td>{html.escape(c.severity)}</td>
                    <td>{html.escape(c.calibration_a)}</td>
                    <td>{html.escape(c.calibration_b)}</td>
                    <td>{html.escape(c.description)}</td>
                    <td>{html.escape(c.recommendation)}</td>
                </tr>
            """)

        # Build sensitivity XML list
        sensitivity_list = ""
        if sensitivity_xmls:
            items = [f"<li>{html.escape(os.path.basename(x))}</li>" for x in sensitivity_xmls]
            sensitivity_list = f"<ul>{''.join(items)}</ul>"

        # Visualization section (embedded as base64 so the report is
        # self-contained and can be moved without breaking the image)
        viz_section = ""
        if visualization_file and os.path.exists(visualization_file):
            import base64

            with open(visualization_file, "rb") as img_f:
                img_b64 = base64.b64encode(img_f.read()).decode("ascii")
            viz_section = f"""
                <div class="section">
                    <h2>5. Calibration Distribution Visualization</h2>
                    <img src="data:image/png;base64,{img_b64}"
                         alt="Calibration Distributions" style="max-width: 100%;">
                </div>
            """

        n_errors = sum(1 for c in conflicts if c.severity == "error")
        n_warnings = sum(1 for c in conflicts if c.severity == "warning")

        conflicts_table = (
            "<table><thead><tr><th>Type</th><th>Severity</th>"
            "<th>Calibration A</th><th>Calibration B</th>"
            "<th>Description</th><th>Recommendation</th></tr></thead>"
            f"<tbody>{''.join(conflict_rows)}</tbody></table>"
        )

        unchecked_pairs = [c for c in conflicts if c.conflict_type == "unchecked"]
        unchecked_note = (
            "<p class='stat-warn'><b>"
            f"{len(unchecked_pairs)} calibration pair(s) could NOT be checked: "
            "a prior with undefined quantiles has no comparable interval. The "
            "absence of a conflict above covers only the pairs that were "
            "compared.</b></p>"
            if unchecked_pairs
            else ""
        )

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Beast2Py Diagnostic Report</title>
    <style>
        body {{
          font-family: 'Segoe UI', Arial, sans-serif;
          margin: 20px; background: #f8f9fa; color: #333;
        }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #2c3e50; margin-top: 30px; }}
        .section {{ background: white; padding: 20px; margin: 15px 0; border-radius: 8px;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: left; }}
        th {{ background: #3498db; color: white; }}
        tr:nth-child(even) {{ background: #f2f2f2; }}
        .conflict-error {{ background: #fee !important; }}
        .conflict-warning {{ background: #fffde7 !important; }}
        .summary {{ display: flex; gap: 20px; margin: 15px 0; }}
        .stat-card {{ background: white; padding: 15px; border-radius: 8px; flex: 1;
                      text-align: center; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .stat-number {{ font-size: 28px; font-weight: bold; }}
        .stat-label {{ color: #666; font-size: 14px; }}
        .stat-errors {{ color: #e74c3c; }}
        .stat-warnings {{ color: #f39c12; }}
        .stat-ok {{ color: #27ae60; }}
    </style>
</head>
<body>
    <h1>Beast2Py Diagnostic Report</h1>

    <div class="section">
        <h2>1. Analysis Overview</h2>
        <p><strong>Analysis:</strong> {html.escape(config.metadata.get('analysis_name', 'N/A'))}</p>
        <p><strong>Description:</strong>
           {html.escape(config.metadata.get('analysis_description', 'N/A'))}</p>
        <p><strong>Partitions:</strong> {len(config.partitions)}</p>
        <p><strong>Calibration points:</strong> {len(config.calibrations)}</p>
        <p><strong>Tree prior:</strong> {config.tree_prior_type.value}</p>
    </div>

    <div class="summary">
        <div class="stat-card">
            <div class="stat-number {'stat-errors' if n_errors else 'stat-ok'}">{n_errors}</div>
            <div class="stat-label">Errors</div>
        </div>
        <div class="stat-card">
            <div class="stat-number
                 {'stat-warnings' if n_warnings else 'stat-ok'}">
              {n_warnings}</div>
            <div class="stat-label">Warnings</div>
        </div>
        <div class="stat-card">
            <div class="stat-number">{len(config.calibrations)}</div>
            <div class="stat-label">Calibration Points</div>
        </div>
    </div>

    <div class="section">
        <h2>2. Calibration Points Summary</h2>
        <table>
            <thead>
                <tr>
                    <th>Name</th><th>Taxa</th><th># Taxa</th>
                    <th>Distribution</th><th>Parameters</th>
                    <th>Mean</th><th>95% HPD</th>
                    <th>Source</th><th>Type</th>
                </tr>
            </thead>
            <tbody>
                {''.join(summary_rows)}
            </tbody>
        </table>
    </div>

    <div class="section">
        <h2>3. Conflict Detection Results</h2>
        {f'<p class="stat-ok">✓ No conflicts detected.</p>' if not conflicts else ''}
        {unchecked_note}
        {conflicts_table if conflicts else ''}
    </div>

    <div class="section">
        <h2>4. Prior Sensitivity Analysis</h2>
        <p>The following XML files have been generated for prior sensitivity analysis.
        Run each with BEAST2 (short chains are sufficient, e.g., 1M generations)
        and compare the resulting tree height posteriors.</p>
        {sensitivity_list if sensitivity_xmls else '<p>No sensitivity analysis XMLs generated.</p>'}
    </div>

    {viz_section}

    <div class="section">
        <p><small>Generated by Beast2Py v{__version__}</small></p>
    </div>
</body>
</html>"""

        Path(output_file).write_text(html_content, encoding="utf-8")

        return output_file


# ============================================================================
# Diagnostic Engine (orchestrator)
# ============================================================================


class DiagnosticEngine:
    """Orchestrate the full diagnostic workflow."""

    @staticmethod
    def run(
        config: BEASTConfig,
        output_dir: str,
        run_sensitivity: bool = True,
        run_visualization: bool = True,
    ) -> DiagnosticReport:
        """Run the complete diagnostic workflow.

        Args:
            config: BEAST configuration.
            output_dir: Output directory for reports and files.
            run_sensitivity: Whether to run sensitivity analysis.
            run_visualization: Whether to generate visualizations.

        Returns:
            DiagnosticReport with all results.
        """
        os.makedirs(output_dir, exist_ok=True)

        # 1. Build calibration summaries
        summaries: List[CalibrationSummary] = []
        for cal in config.calibrations:
            stats: Dict[str, Any] = {}
            stats_warning = ""
            if cal.distribution:
                try:
                    stats = compute_distribution_stats(
                        cal.distribution.type,
                        cal.distribution.parameters,
                        cal.distribution.offset,
                    )
                except ValueError as e:
                    stats_warning = str(e)

            summaries.append(
                CalibrationSummary(
                    name=cal.name,
                    taxa=cal.taxa,
                    taxa_count=len(cal.taxa) if cal.taxa else 0,
                    distribution_type=cal.distribution.type if cal.distribution else "none",
                    distribution_params=cal.distribution.parameters if cal.distribution else {},
                    offset=cal.distribution.offset if cal.distribution else 0.0,
                    stats=stats,
                    provenance_source=cal.provenance.source if cal.provenance else "",
                    provenance_reference=cal.provenance.reference if cal.provenance else "",
                    calibration_type=cal.provenance.calibration_type if cal.provenance else "",
                )
            )
            if stats_warning:
                summaries[-1].stats = {"warning": stats_warning}

        # 2. Detect conflicts
        conflicts = ConflictDetector.detect_conflicts(
            config.calibrations, settings=getattr(config, "diagnostics", None)
        )

        # 3. Sensitivity analysis
        sensitivity_xmls: List[str] = []
        if run_sensitivity:
            sens_dir = os.path.join(output_dir, "sensitivity")
            sensitivity_xmls = SensitivityAnalyzer.generate_sensitivity_xmls(config, sens_dir)

        # 4. Visualization
        viz_file: Optional[str] = None
        viz_error: Optional[str] = None
        if run_visualization and config.calibrations:
            viz_file = os.path.join(output_dir, "calibration_distributions.png")
            try:
                VisualizationGenerator.plot_calibration_distributions(config.calibrations, viz_file)
            except Exception as e:
                # Do not abort the diagnostics, but surface the reason
                # instead of silently producing a report without a figure.
                viz_error = f"Visualization failed: {e}"
                viz_file = None

        return DiagnosticReport(
            config=config,
            conflicts=conflicts,
            summaries=summaries,
            sensitivity_xmls=sensitivity_xmls,
            visualization_file=viz_file,
            visualization_error=viz_error,
        )

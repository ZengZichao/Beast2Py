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
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from . import __version__
from .models import BEASTConfig, CalibrationPoint, DistributionConfig
from .utils import compute_distribution_stats, compute_wasserstein_distance

# ============================================================================
# Gate 2 capability boundary
# ============================================================================

# Gate 2 reads the YAML, not a trace, so its reach has to be stated as a list
# rather than implied by a green tick.  These two tuples are the single source
# of that statement: the CLI, the HTML report and the manuals all render them.

#: The structural properties Gate 2 can decide from the configuration alone.
DETECTABLE_CHECKS: Tuple[str, ...] = (
    "an incompatible or duplicated pair of calibrations on the same node",
    "a nested pair whose parent prior is wholly younger than its child prior",
    "a set of monophyly declarations that no rooted tree can satisfy",
    "a semantically illegal combination (root + originate, tipsonly + originate, "
    "an out-of-support or non-finite distribution parameter)",
)

#: What Gate 2 cannot see, and which therefore cannot be ruled out by Gate 2.
NOT_DETECTABLE_CHECKS: Tuple[str, ...] = (
    "a distorted effective marginal prior: the tree prior, the topological "
    "constraints and the remaining calibrations can truncate or reshape a "
    "specified prior without any structural contradiction existing",
    "posterior dominance: deciding it needs the likelihood, i.e. matched "
    "posterior traces, which Beast2Py neither runs nor reads",
    "MCMC non-convergence, poor mixing or an inadequate ESS",
    "too little data to inform a node, which is a property of the alignment and "
    "the model rather than of the configuration",
)

#: One-name summaries of the entries above, for the single-line CLI notice.
#: A test pins that the two tuples stay the same length and order.
NOT_DETECTABLE_LABELS: Tuple[str, ...] = (
    "distortion of the effective marginal prior",
    "posterior dominance",
    "MCMC non-convergence",
    "insufficient data to inform a node",
)

#: The four distributions Gate 2's vocabulary is deliberately restricted to the
#: first of, and which a prior-only run is needed to see the second of.
PRIOR_LAYERS: Tuple[Tuple[str, str], ...] = (
    (
        "specified marginal prior",
        "the distribution written in the YAML, which is all Gate 2 reads",
    ),
    (
        "effective marginal prior",
        "the marginal a node actually samples once the tree prior, the "
        "constraints and the other calibrations act on it; observable only by "
        "running a prior-only chain",
    ),
    (
        "joint effective prior",
        "the simultaneous distribution of two or more node ages; a single trace "
        "column gives a marginal, so the joint needs joint samples (per-column "
        "medians of one .log), not a pile of marginals",
    ),
    (
        "posterior",
        "the distribution after the likelihood is applied; out of scope for "
        "Beast2Py, which produces no estimate of it",
    ),
)


def gate2_scope_lines() -> List[str]:
    """Render the Gate 2 boundary as plain text.

    The CLI, the ``diagnose`` output and the manuals all print this, so the
    wording of what a green Gate 2 does not mean cannot drift apart between
    them.
    """
    lines = ["Gate 2 can detect (static structure only):"]
    lines += [f"  - {item}" for item in DETECTABLE_CHECKS]
    lines.append("Gate 2 cannot detect, so a passing Gate 2 rules none of these out:")
    lines += [f"  - {item}" for item in NOT_DETECTABLE_CHECKS]
    return lines


def gate2_scope_html() -> str:
    """The same boundary, for the HTML report."""
    detectable = "".join(f"<li>{html.escape(item)}</li>" for item in DETECTABLE_CHECKS)
    undetectable = "".join(f"<li>{html.escape(item)}</li>" for item in NOT_DETECTABLE_CHECKS)
    layers = "".join(
        f"<dt>{html.escape(name)}</dt><dd>{html.escape(gloss)}</dd>" for name, gloss in PRIOR_LAYERS
    )
    return (
        "<h3>What this gate decides, and what it cannot</h3>"
        f"<p>Gate 2 reads the configuration, not a trace. It checks the "
        f"<b>specified marginal priors</b> and the taxon sets they attach to. "
        f"The four distributions below are distinct; only the first is visible "
        f"to this gate.</p>"
        f"<dl>{layers}</dl>"
        f"<p><b>Detected here:</b></p><ul>{detectable}</ul>"
        f"<p><b>Not detectable here — a passing Gate 2 does not rule these "
        f"out:</b></p><ul>{undetectable}</ul>"
        "<p>Pairs that could not be compared are counted separately as "
        "<code>unchecked</code> and are never reported as passes.</p>"
    )


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
    # mean/median/hpd_lower/hpd_upper are floats; ``warning`` is a string the
    # report shows in their place when the distribution has no defined quantiles.
    stats: Dict[str, Any]
    provenance_source: str = ""
    provenance_reference: str = ""
    calibration_type: str = ""


@dataclass
class DiagnosticReport:
    """Complete diagnostic report."""

    config: BEASTConfig
    conflicts: List[Conflict]
    summaries: List[CalibrationSummary]
    assessments: List[PairAssessment] = field(default_factory=list)
    sensitivity_xmls: List[str] = field(default_factory=list)
    visualization_file: Optional[str] = None
    report_path: Optional[str] = None
    visualization_error: Optional[str] = None


# ============================================================================
# Conflict Detector
# ============================================================================


@dataclass
class PairAssessment:
    """Structural assessment of one calibration pair.

    ``relation`` names the taxon-set relation; ``node_a``/``node_b`` name the
    node each calibration resolves to (root, crown, stem or tip).  ``status``
    is ``checked`` when a structural rule was applied, ``unchecked`` when the
    available information could not decide the relation, and
    ``not_applicable`` when no ordering or topology rule applies (disjoint
    clades, tip dates).
    """

    calibration_a: str
    calibration_b: str
    relation: str
    node_a: str
    node_b: str
    status: str
    conflicts: List["Conflict"] = field(default_factory=list)


class ConflictDetector:
    """Detect structural conflicts between calibration points.

    Only three families of conflict are diagnosed, each with a topological or
    temporal justification:

    * ``same_node`` -- two calibrations constrain the same node to incompatible
      distributions (or repeat it).  Crown-crown and stem-stem pairs on the
      same taxon set are the same node; crown-stem pairs are parent/child.
    * ``temporal`` -- a strictly nesting pair whose parent prior is wholly
      younger than its child prior.
    * ``topology`` -- two intersecting, non-nested taxon sets both declared
      monophyletic, which no rooted tree can satisfy.

    Similar marginal densities of *independent* clades are not a conflict.
    They are reported, only when the user opts in with an explicit threshold,
    as ``prior_similarity`` with severity ``info``.  Such a report carries no
    merge recommendation and never blocks Gate 2.

    :data:`DETECTABLE_CHECKS` and :data:`NOT_DETECTABLE_CHECKS` state the
    boundary of that decision in full; nothing else is inside Gate 2.
    """

    # No similarity threshold ships as a default.  The report is off unless the
    # configuration both enables it and supplies a number, because no reference
    # set has calibrated one: see ``_normalise_settings``.

    @staticmethod
    def _normalise_settings(
        settings: Optional[Dict[str, Any]],
        overlap_threshold: Optional[float],
        relative_threshold: Optional[float],
    ) -> Dict[str, Any]:
        cfg = dict(settings or {})
        enabled = bool(cfg.get("prior_similarity", False))
        measure = str(cfg.get("similarity_measure", "relative")).lower()
        threshold = cfg.get("similarity_threshold")

        # Explicit function arguments are an opt-in for callers that do not go
        # through the YAML parser (library users, older notebooks).
        if overlap_threshold is not None:
            enabled = True
            measure = "absolute"
            threshold = overlap_threshold
        if relative_threshold is not None:
            enabled = True
            measure = "relative"
            threshold = relative_threshold

        if measure not in ("relative", "absolute", "both"):
            raise ValueError(
                "diagnostics.similarity_measure must be relative, absolute or both, "
                f"got {measure!r}"
            )
        if enabled and threshold is None:
            raise ValueError(
                "prior similarity was requested without a threshold. Beast2Py ships "
                "no calibrated default; pass diagnostics.similarity_threshold."
            )
        if threshold is not None:
            threshold = float(threshold)
            if threshold <= 0:
                raise ValueError("similarity threshold must be > 0")
        return {
            "prior_similarity": enabled,
            "similarity_measure": measure,
            "similarity_threshold": threshold,
        }

    @staticmethod
    def assess(
        calibrations: List[CalibrationPoint],
        overlap_threshold: Optional[float] = None,
        relative_threshold: Optional[float] = None,
        settings: Optional[Dict[str, Any]] = None,
        all_taxa: Optional[List[str]] = None,
    ) -> List[PairAssessment]:
        """Assess every calibration pair once.

        Args:
            calibrations: Calibration points.
            overlap_threshold: Absolute similarity threshold in time units.
                Passing it opts in to the descriptive report.
            relative_threshold: Dimensionless similarity threshold.
            settings: The configuration's ``diagnostics`` block.
            all_taxa: Every taxon in the analysis, when known.  A non-null
                taxon set that contains all of them resolves to the root, so
                ``use_originate`` on it is invalid and is reported as unchecked
                here (the parser rejects it when it has the same information).

        Returns:
            One :class:`PairAssessment` per unordered pair.
        """
        opts = ConflictDetector._normalise_settings(settings, overlap_threshold, relative_threshold)
        all_taxa_set = set(all_taxa) if all_taxa is not None else None
        assessments: List[PairAssessment] = []

        for i, cal_a in enumerate(calibrations):
            for j, cal_b in enumerate(calibrations):
                if i >= j:
                    continue
                assessments.append(ConflictDetector._assess_pair(cal_a, cal_b, opts, all_taxa_set))
        return assessments

    @staticmethod
    def detect_conflicts(
        calibrations: List[CalibrationPoint],
        overlap_threshold: Optional[float] = None,
        relative_threshold: Optional[float] = None,
        settings: Optional[Dict[str, Any]] = None,
        all_taxa: Optional[List[str]] = None,
    ) -> List[Conflict]:
        """Return the flattened list of conflicts from :meth:`assess`."""
        conflicts: List[Conflict] = []
        for assessment in ConflictDetector.assess(
            calibrations,
            overlap_threshold=overlap_threshold,
            relative_threshold=relative_threshold,
            settings=settings,
            all_taxa=all_taxa,
        ):
            conflicts.extend(assessment.conflicts)
        return conflicts

    @staticmethod
    def summarize(assessments: List[PairAssessment]) -> Dict[str, int]:
        """Count pair statuses and conflict severities for reports."""
        summary = {
            "pairs": len(assessments),
            "checked": 0,
            "unchecked": 0,
            "not_applicable": 0,
            "errors": 0,
            "warnings": 0,
            "info": 0,
            "prior_similarity": 0,
        }
        for assessment in assessments:
            summary[assessment.status] = summary.get(assessment.status, 0) + 1
            for conflict in assessment.conflicts:
                summary[conflict.severity] = summary.get(conflict.severity, 0) + 1
                if conflict.conflict_type == "prior_similarity":
                    summary["prior_similarity"] += 1
        return summary

    # -- pair assessment ----------------------------------------------------

    @staticmethod
    def _node_kind(cal: CalibrationPoint, all_taxa: Optional[set]) -> str:
        if cal.tipsonly:
            return "tip"
        if cal.taxa is None:
            return "root"
        if all_taxa is not None and set(cal.taxa) == all_taxa:
            # A non-null taxon list that spans the whole ingroup still resolves
            # to the root: MRCA(all taxa) = root and the root has no stem.
            return "root"
        return "stem" if cal.use_originate else "crown"

    @staticmethod
    def _contains(
        outer: Optional[Set[str]],
        inner: Optional[Set[str]],
    ) -> bool:
        """Whether ``outer`` strictly contains ``inner``, False when either is absent.

        ``nested`` and ``intersecting`` are only ever reported for a pair that
        carries two taxon sets, so the ``None`` arms do not arise in those
        branches; spelling the guard out keeps the classifier's contract in one
        place instead of relying on a caller remembering it.
        """
        if outer is None or inner is None:
            return False
        return inner < outer

    @staticmethod
    def _assess_pair(
        cal_a: CalibrationPoint,
        cal_b: CalibrationPoint,
        opts: Dict[str, Any],
        all_taxa: Optional[set],
    ) -> PairAssessment:
        node_a = ConflictDetector._node_kind(cal_a, all_taxa)
        node_b = ConflictDetector._node_kind(cal_b, all_taxa)
        set_a = set(cal_a.taxa) if cal_a.taxa is not None else None
        set_b = set(cal_b.taxa) if cal_b.taxa is not None else None

        assessment = PairAssessment(
            calibration_a=cal_a.name,
            calibration_b=cal_b.name,
            relation="",
            node_a=node_a,
            node_b=node_b,
            status="not_applicable",
        )

        # Tip dates are observations, not clade-age constraints.  They cannot
        # be ordered against a clade MRCA without a dated tree, so no
        # structural rule applies.
        if node_a == "tip" or node_b == "tip":
            assessment.relation = "tip"
            return assessment

        both_root = node_a == "root" and node_b == "root"
        same_set = both_root or (set_a is not None and set_b is not None and set_a == set_b)

        if same_set:
            if node_a == node_b or both_root:
                assessment.relation = "same_node"
                assessment.status = "checked"
                assessment.conflicts.extend(ConflictDetector._check_same_node(cal_a, cal_b, node_a))
            elif {node_a, node_b} == {"crown", "stem"}:
                # Same taxon set, one crown and one stem: stem(parent) is
                # strictly older than crown(child).  This is a determinate
                # parent/child relation, not a duplicate.
                parent, child = (cal_a, cal_b) if node_a == "stem" else (cal_b, cal_a)
                assessment.relation = "same_set_crown_stem"
                assessment.status = "checked"
                assessment.conflicts.extend(
                    ConflictDetector._compare_parent_child(
                        parent,
                        child,
                        relation=(
                            "the stem node is the parent of the crown node for the "
                            "same taxon set"
                        ),
                    )
                )
            else:
                assessment.relation = "same_set_unknown"
                assessment.status = "unchecked"
                assessment.conflicts.append(
                    ConflictDetector._unchecked(
                        cal_a,
                        cal_b,
                        f"the same taxon set resolves to nodes ({node_a}, {node_b}) "
                        f"whose ordering cannot be decided from the configuration",
                    )
                )
            ConflictDetector._maybe_similarity(cal_a, cal_b, opts, assessment)
            if any(c.conflict_type == "unchecked" for c in assessment.conflicts):
                assessment.status = "unchecked"
            return assessment

        # One or both calibrations are the root, or the sets are nested,
        # intersecting or disjoint.
        relation = ConflictDetector._taxon_relation(set_a, set_b, node_a, node_b)
        assessment.relation = relation

        if relation == "root_vs_clade":
            parent, child = (cal_a, cal_b) if node_a == "root" else (cal_b, cal_a)
            assessment.status = "checked"
            assessment.conflicts.extend(
                ConflictDetector._compare_parent_child(
                    parent, child, relation="the root is an ancestor of every clade"
                )
            )
        elif relation == "nested":
            nested = ConflictDetector._contains(set_a, set_b)
            parent, child = (cal_a, cal_b) if nested else (cal_b, cal_a)
            parent_node, child_node = (node_a, node_b) if nested else (node_b, node_a)
            if parent_node == "crown" and child_node == "stem":
                # C(parent) >= C(child) and S(child) > C(child); the stem of
                # the child may sit above or below the crown of the parent.
                assessment.status = "unchecked"
                assessment.conflicts.append(
                    ConflictDetector._unchecked(
                        parent,
                        child,
                        "the nested child is calibrated on its stem while the "
                        "parent is calibrated on its crown, so their relative "
                        "ages are not determined by taxon containment",
                    )
                )
            else:
                assessment.status = "checked"
                assessment.conflicts.extend(
                    ConflictDetector._compare_parent_child(
                        parent,
                        child,
                        relation="the parent taxon set strictly contains the child",
                    )
                )
        elif relation == "intersecting":
            if cal_a.monophyletic and cal_b.monophyletic:
                assessment.status = "checked"
                assessment.conflicts.append(
                    ConflictDetector._topology_conflict(
                        cal_a, cal_b, (set_a or set()) & (set_b or set())
                    )
                )
            else:
                assessment.status = "unchecked"
                assessment.conflicts.append(
                    ConflictDetector._unchecked(
                        cal_a,
                        cal_b,
                        "the taxon sets intersect but neither contains the other, "
                        "and at least one calibration does not assert monophyly; "
                        "the topology constraint cannot be evaluated",
                    )
                )
        else:  # disjoint
            assessment.status = "not_applicable"

        ConflictDetector._maybe_similarity(cal_a, cal_b, opts, assessment)
        if any(c.conflict_type == "unchecked" for c in assessment.conflicts):
            assessment.status = "unchecked"
        return assessment

    @staticmethod
    def _taxon_relation(
        set_a: Optional[set],
        set_b: Optional[set],
        node_a: str,
        node_b: str,
    ) -> str:
        if set_a is None or set_b is None:
            # Only one side can be the root here (both-root was handled as the
            # same node); the root nests every clade.
            return "root_vs_clade"
        if set_a < set_b or set_b < set_a:
            return "nested"
        if set_a & set_b:
            return "intersecting"
        return "disjoint"

    # -- rule implementations ----------------------------------------------

    @staticmethod
    def _get_stats(cal: CalibrationPoint) -> Optional[Dict[str, float]]:
        if not cal.distribution:
            return None
        try:
            return compute_distribution_stats(
                cal.distribution.type,
                cal.distribution.parameters,
                cal.distribution.offset,
            )
        except ValueError:
            return None

    @staticmethod
    def _signature(cal: CalibrationPoint) -> Any:
        dist = cal.distribution
        if dist is None:
            # No distribution to compare: two such points on one node are
            # indistinguishable, which is what a duplicate reading means here.
            return ("<none>", (), 0.0)
        return (
            dist.type.lower(),
            tuple(sorted((str(k), str(v)) for k, v in (dist.parameters or {}).items())),
            float(dist.offset or 0.0),
        )

    @staticmethod
    def _unchecked(
        cal_a: CalibrationPoint,
        cal_b: CalibrationPoint,
        reason: str,
    ) -> Conflict:
        return Conflict(
            conflict_type="unchecked",
            severity="info",
            calibration_a=cal_a.name,
            calibration_b=cal_b.name,
            description=f"Not checked: {reason}.",
            recommendation=(
                "No pass is implied for this pair. Supply a proper distribution "
                "with finite quantiles and a determinate node relation, or verify "
                "the pair manually."
            ),
        )

    @staticmethod
    def _check_same_node(
        cal_a: CalibrationPoint,
        cal_b: CalibrationPoint,
        node_kind: str,
    ) -> List[Conflict]:
        stats_a = ConflictDetector._get_stats(cal_a)
        stats_b = ConflictDetector._get_stats(cal_b)
        if not stats_a or not stats_b:
            unclear = [c.name for c, s in ((cal_a, stats_a), (cal_b, stats_b)) if not s]
            return [
                ConflictDetector._unchecked(
                    cal_a,
                    cal_b,
                    f"{', '.join(unclear)} uses a distribution whose mean or "
                    f"quantiles are undefined, so the two constraints on the "
                    f"same {node_kind} node cannot be compared",
                )
            ]

        if ConflictDetector._signature(cal_a) == ConflictDetector._signature(cal_b):
            return [
                Conflict(
                    conflict_type="duplicate",
                    severity="warning",
                    calibration_a=cal_a.name,
                    calibration_b=cal_b.name,
                    description=(
                        f"'{cal_a.name}' and '{cal_b.name}' declare the same "
                        f"{node_kind} node (identical taxon set) with an identical "
                        f"distribution. The prior is counted twice."
                    ),
                    recommendation=(
                        "Keep one calibration, or make the two constraints "
                        "intentionally different and document why the node is "
                        "constrained twice."
                    ),
                )
            ]

        disjoint = (
            stats_a["hpd_upper"] < stats_b["hpd_lower"]
            or stats_b["hpd_upper"] < stats_a["hpd_lower"]
        )
        if disjoint:
            return [
                Conflict(
                    conflict_type="same_node",
                    severity="error",
                    calibration_a=cal_a.name,
                    calibration_b=cal_b.name,
                    description=(
                        f"'{cal_a.name}' and '{cal_b.name}' constrain the same "
                        f"{node_kind} node (identical taxon set) to distributions "
                        f"with disjoint 95% intervals "
                        f"([{stats_a['hpd_lower']:.2f}, {stats_a['hpd_upper']:.2f}] "
                        f"vs [{stats_b['hpd_lower']:.2f}, {stats_b['hpd_upper']:.2f}]). "
                        f"At most one of the two can describe the node's fossil "
                        f"evidence."
                    ),
                    recommendation=(
                        "Decide which fossil/estimate applies to this node and "
                        "remove or replace the other calibration. This is a "
                        "same-node contradiction, not a similarity between "
                        "independent clades."
                    ),
                )
            ]

        return [
            Conflict(
                conflict_type="same_node",
                severity="warning",
                calibration_a=cal_a.name,
                calibration_b=cal_b.name,
                description=(
                    f"'{cal_a.name}' and '{cal_b.name}' constrain the same "
                    f"{node_kind} node (identical taxon set) with two different "
                    f"distributions. Their 95% intervals overlap, so this is not a "
                    f"hard contradiction, but the node prior is their product, "
                    f"which is usually not what either entry describes."
                ),
                recommendation=(
                    "Check that two constraints on the same node are intended; "
                    "otherwise keep one."
                ),
            )
        ]

    @staticmethod
    def _compare_parent_child(
        parent: CalibrationPoint,
        child: CalibrationPoint,
        relation: str,
    ) -> List[Conflict]:
        stats_parent = ConflictDetector._get_stats(parent)
        stats_child = ConflictDetector._get_stats(child)
        if not stats_parent or not stats_child:
            unclear = [c.name for c, s in ((parent, stats_parent), (child, stats_child)) if not s]
            return [
                ConflictDetector._unchecked(
                    parent,
                    child,
                    f"{', '.join(unclear)} uses a distribution whose mean or "
                    f"quantiles are undefined, so the temporal ordering "
                    f"({relation}) cannot be compared",
                )
            ]

        if stats_parent["hpd_upper"] < stats_child["hpd_lower"]:
            return [
                Conflict(
                    conflict_type="temporal",
                    severity="error",
                    calibration_a=parent.name,
                    calibration_b=child.name,
                    description=(
                        f"Temporal inconsistency: '{parent.name}' (parent) has "
                        f"95% upper {stats_parent['hpd_upper']:.2f}, younger than "
                        f"'{child.name}' (child) 95% lower "
                        f"{stats_child['hpd_lower']:.2f}; {relation}. The parent "
                        f"MRCA cannot be younger than the nested child MRCA."
                    ),
                    recommendation=(
                        f"Adjust the priors so that '{parent.name}' spans older "
                        f"times than '{child.name}', or correct the taxon sets."
                    ),
                )
            ]
        # Overlapping intervals are not a contradiction: a soft parent prior
        # with a long tail can have a younger mean than its child's prior while
        # remaining perfectly compatible.  Only a wholly younger parent
        # interval is reported.
        return []

    @staticmethod
    def _topology_conflict(
        cal_a: CalibrationPoint,
        cal_b: CalibrationPoint,
        shared: set,
    ) -> Conflict:
        return Conflict(
            conflict_type="topology",
            severity="error",
            calibration_a=cal_a.name,
            calibration_b=cal_b.name,
            description=(
                f"Topology conflict: '{cal_a.name}' and '{cal_b.name}' have "
                f"intersecting but non-nested taxon sets (shared: "
                f"{', '.join(sorted(shared))}) and both assert monophyly. Two "
                f"monophyletic clades in a rooted tree must be nested or disjoint, "
                f"so no tree can satisfy both."
            ),
            recommendation=(
                "Make the taxon sets nested or disjoint, or set monophyletic: "
                "false on one calibration if it is not a clade constraint."
            ),
        )

    # -- descriptive similarity --------------------------------------------

    @staticmethod
    def _maybe_similarity(
        cal_a: CalibrationPoint,
        cal_b: CalibrationPoint,
        opts: Dict[str, Any],
        assessment: PairAssessment,
    ) -> None:
        if not opts.get("prior_similarity"):
            return
        if assessment.relation == "same_node":
            return
        if not cal_a.distribution or not cal_b.distribution:
            return
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
            if assessment.status != "checked":
                assessment.status = "unchecked"
                assessment.conflicts.append(
                    ConflictDetector._unchecked(
                        cal_a,
                        cal_b,
                        "prior similarity could not be computed because at least "
                        "one distribution has no finite quantiles",
                    )
                )
            return

        stats_a = ConflictDetector._get_stats(cal_a)
        stats_b = ConflictDetector._get_stats(cal_b)
        ratio = None
        if stats_a and stats_b:
            width_a = stats_a["hpd_upper"] - stats_a["hpd_lower"]
            width_b = stats_b["hpd_upper"] - stats_b["hpd_lower"]
            pooled = (width_a + width_b) / 2.0
            if pooled > 0:
                ratio = distance / pooled

        measure = opts["similarity_measure"]
        threshold = opts["similarity_threshold"]
        reasons: List[str] = []
        flagged = False
        if measure in ("relative", "both") and ratio is not None:
            reasons.append(
                f"relative Wasserstein distance = {ratio:.3f} of the pooled "
                f"95%-interval width (threshold {threshold})"
            )
            if ratio < threshold:
                flagged = True
        if measure in ("absolute", "both"):
            reasons.append(
                f"Wasserstein distance = {distance:.3f} time units " f"(threshold {threshold})"
            )
            if distance < threshold:
                flagged = True

        if flagged:
            assessment.conflicts.append(
                Conflict(
                    conflict_type="prior_similarity",
                    severity="info",
                    calibration_a=cal_a.name,
                    calibration_b=cal_b.name,
                    description=(
                        f"Prior densities of '{cal_a.name}' and '{cal_b.name}' are "
                        f"similar ({'; '.join(reasons)}). This is a descriptive "
                        f"summary of two specified marginal priors, not evidence "
                        f"that fossil information is redundant. Independent clades "
                        f"may legitimately share a similar age prior."
                    ),
                    recommendation=(
                        "None. This entry is not a conflict, does not block Gate 2, "
                        "and does not justify merging or removing a calibration."
                    ),
                )
            )

    # -- retained helper for external users --------------------------------

    @staticmethod
    def _check_monophyly_conflicts(
        calibrations: List[CalibrationPoint],
    ) -> List[Conflict]:
        """Backward-compatible wrapper returning only topology conflicts."""
        return [
            conflict
            for conflict in ConflictDetector.detect_conflicts(calibrations)
            if conflict.conflict_type == "topology"
        ]


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
        # Create a copy with sample_from_prior=True.  The chain length is the
        # configured one: a hard-coded 1M cap produced too few samples for a
        # 200-ESS convergence criterion.  Each scenario gets a distinct seed so
        # the chains are independent and their $(seed) log names cannot collide.
        prior_config = copy.deepcopy(config)
        prior_config.mcmc.sample_from_prior = True
        base_seed = int(getattr(prior_config.mcmc, "seed", 1) or 1)
        if suffix.startswith("00_"):
            prior_config.mcmc.seed = base_seed
        else:
            try:
                index = int(suffix.split("_", 1)[0])
            except ValueError:
                index = 1
            prior_config.mcmc.seed = base_seed + index

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
        import matplotlib
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))

        colors = matplotlib.colormaps["tab10"](np.linspace(0, 1, max(len(calibrations), 1)))

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
        assessments: Optional[List[PairAssessment]] = None,
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
            severity_class = {
                "error": "conflict-error",
                "warning": "conflict-warning",
            }.get(c.severity, "conflict-info")
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
        stats = (
            ConflictDetector.summarize(assessments)
            if assessments is not None
            else {
                "pairs": 0,
                "checked": 0,
                "unchecked": sum(1 for c in conflicts if c.conflict_type == "unchecked"),
                "not_applicable": 0,
                "prior_similarity": 0,
            }
        )

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
            "a prior with undefined quantiles, an undecidable crown/stem relation, "
            "or a non-monophyletic intersection has no comparable structural "
            "prediction. The absence of a structural conflict covers only the "
            "pairs that were compared; it is not a pass for these pairs.</b></p>"
            if unchecked_pairs
            else ""
        )
        similarity_note = (
            "<p class='stat-warn'><b>Prior similarity is descriptive only.</b> "
            "Similar marginal priors of independent clades are not evidence of "
            "redundant fossil information, are not a conflict, and must not be "
            "used to merge or remove a calibration.</p>"
            if any(c.conflict_type == "prior_similarity" for c in conflicts)
            else ""
        )
        pass_note = (
            '<p class="stat-ok">\u2713 No structural conflicts among the checked pairs.</p>'
            if not n_errors and not n_warnings
            else ""
        )
        capability_block = gate2_scope_html()

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
        .conflict-info {{ background: #eef6ff !important; }}
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
            <div class="stat-label">Structural errors</div>
        </div>
        <div class="stat-card">
            <div class="stat-number
                 {'stat-warnings' if n_warnings else 'stat-ok'}">
              {n_warnings}</div>
            <div class="stat-label">Warnings</div>
        </div>
        <div class="stat-card">
            <div class="stat-number {'stat-warnings' if stats.get('unchecked') else 'stat-ok'}">
              {stats.get('unchecked', 0)}</div>
            <div class="stat-label">Unchecked pairs</div>
        </div>
        <div class="stat-card">
            <div class="stat-number">{stats.get('checked', 0)}</div>
            <div class="stat-label">Checked pairs</div>
        </div>
        <div class="stat-card">
            <div class="stat-number">{stats.get('not_applicable', 0)}</div>
            <div class="stat-label">No applicable relation</div>
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
        <h2>3. Structural Conflict Detection Results</h2>
        <p><b>{stats.get('pairs', 0)} pairs assessed:</b>
           {stats.get('checked', 0)} checked,
           {stats.get('unchecked', 0)} unchecked,
           {stats.get('not_applicable', 0)} with no applicable structural relation.
           "Checked" means the pair was compared, not that the whole model is correct.</p>
        {pass_note}
        {unchecked_note}
        {similarity_note}
        {capability_block}
        {conflicts_table if conflicts else ''}
    </div>

    <div class="section">
        <h2>4. Prior-Only Sensitivity Design</h2>
        <p>The following XML files record a prior-only analysis design. Running
        them samples the specified and effective marginal priors; it does not
        sample or diagnose the posterior. Inspect ESS and mixing, and report the
        chain length, burn-in and convergence criteria used.</p>
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

        # 2. Detect conflicts.  The full taxon list lets a calibration whose
        # non-null taxon list spans the whole ingroup resolve to the root
        # rather than to a nonexistent stem.
        all_taxa: List[str] = []
        for alignment in getattr(config, "all_alignments", []) or []:
            all_taxa.extend(alignment.taxa_names)
        assessments = ConflictDetector.assess(
            config.calibrations,
            settings=getattr(config, "diagnostics", None),
            all_taxa=sorted(set(all_taxa)) or None,
        )
        conflicts = [conflict for assessment in assessments for conflict in assessment.conflicts]

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
            assessments=assessments,
            sensitivity_xmls=sensitivity_xmls,
            visualization_file=viz_file,
            visualization_error=viz_error,
        )

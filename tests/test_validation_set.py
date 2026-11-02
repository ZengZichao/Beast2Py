"""Synthetic standard-answer validation set for the calibration diagnostics.

Each scenario in ``tests/fixtures/diagnostics_validation.json`` has an
expected status and conflict count. All boundary ages and priors are synthetic
values chosen only to exercise the ordering rules; they are not taken from any
publication or empirical dataset.
"""

import json
from pathlib import Path

import pytest

from beast2py.diagnostics import ConflictDetector
from beast2py.models import CalibrationPoint, DistributionConfig

FIXTURE = Path(__file__).parent / "fixtures" / "diagnostics_validation.json"


def _load():
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def _cal(raw):
    return CalibrationPoint(
        name=raw["name"],
        taxa=raw.get("taxa"),
        monophyletic=raw.get("monophyletic", True),
        use_originate=raw.get("use_originate", False),
        tipsonly=raw.get("tipsonly", False),
        distribution=DistributionConfig(
            type=raw["distribution"]["type"],
            parameters=raw["distribution"].get("parameters", {}),
            offset=raw["distribution"].get("offset", 0.0),
        ),
    )


def _ids():
    return [s["id"] for s in _load()["scenarios"]]


@pytest.mark.parametrize("scenario_id", _ids())
def test_validation_scenario(scenario_id):
    scenario = next(s for s in _load()["scenarios"] if s["id"] == scenario_id)
    calibrations = [_cal(c) for c in scenario["calibrations"]]
    assessments = ConflictDetector.assess(calibrations)
    conflicts = [c for a in assessments for c in a.conflicts]
    expected = scenario["expected"]

    assert sum(1 for c in conflicts if c.severity == "error") == expected["errors"], scenario[
        "description"
    ]
    assert sum(1 for c in conflicts if c.severity == "warning") == expected.get(
        "warnings", 0
    ), scenario["description"]
    for conflict_type in ("topology", "temporal", "same_node"):
        assert sum(1 for c in conflicts if c.conflict_type == conflict_type) == expected.get(
            conflict_type, 0
        ), scenario["description"]
    assert sum(1 for c in conflicts if c.conflict_type == "unchecked") == expected.get(
        "unchecked", 0
    ), scenario["description"]
    assert assessments[0].status == expected["status"], scenario["description"]
    assert not any(
        conflict.conflict_type == "prior_similarity" for conflict in conflicts
    ), "prior similarity must be off unless explicitly requested"


def test_independent_similar_priors_do_not_block_gate_2():
    """Regression for the old redundancy rule: similarity is informational."""
    scenario = next(s for s in _load()["scenarios"] if s["id"] == "independent_similar")
    calibrations = [_cal(c) for c in scenario["calibrations"]]
    conflicts = ConflictDetector.detect_conflicts(
        calibrations,
        settings={
            "prior_similarity": True,
            "similarity_measure": "relative",
            "similarity_threshold": 0.15,
        },
    )
    assert [c for c in conflicts if c.severity == "error"] == []
    similarity = [c for c in conflicts if c.conflict_type == "prior_similarity"]
    assert len(similarity) == 1
    assert similarity[0].severity == "info"
    assert "not evidence" in similarity[0].description
    assert "not a conflict" in similarity[0].recommendation


def test_validation_set_covers_every_required_answer():
    ids = set(_ids())
    assert {
        "independent_similar",
        "independent_different",
        "same_node_mutually_exclusive",
        "duplicate_same_node",
        "parent_younger_than_child",
        "intersecting_non_nested_monophyletic",
        "intersecting_non_monophyletic",
        "crown_stem_same_set_violation",
        "nested_parent_crown_child_stem",
        "stem_stem_nested",
        "root_vs_clade",
        "one_on_x_unchecked",
        "tipsonly_not_a_clade_age",
    } <= ids

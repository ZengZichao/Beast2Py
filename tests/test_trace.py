"""Tests for the BEAST2 trace reader and the prior-only summary arithmetic.

All logs are synthetic BEAST2-shaped tables written by the test helpers. No
archived run, published number, or external dataset is required.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from beast2py.trace import (
    TraceFormatError,
    combine_traces,
    effective_sample_size,
    hpd_interval,
    node_age_column,
    prior_shift,
    read_trace,
    summarise_column,
)


def write_log(
    path: Path,
    columns,
    rows,
    *,
    header_comment="# BEAST2 2.7.7 output\n#\n",
    id_name="Sample",
):
    """Write a BEAST2-shaped tab-separated log."""
    names = [id_name] + list(columns)
    lines = [header_comment.rstrip("\n")]
    lines.append("\t".join(names))
    for index, row in enumerate(rows):
        lines.append("\t".join([str((index + 1) * 100)] + [repr(float(v)) for v in row]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# --- parsing ----------------------------------------------------------------


def test_read_trace_uses_the_header_not_the_column_order(tmp_path):
    path = write_log(tmp_path / "a.log", ["treeHeight", "mrca.age(X.taxonset)"], [[1.0, 2.0]])
    trace = read_trace(path)
    assert trace.columns == ["Sample", "treeHeight", "mrca.age(X.taxonset)"]
    assert trace.numeric_columns == ["treeHeight", "mrca.age(X.taxonset)"]
    assert list(trace.column("mrca.age(X.taxonset)")) == [2.0]


def test_read_trace_accepts_both_sample_id_spellings(tmp_path):
    write_log(tmp_path / "rate.log", ["treeHeight"], [[1.0]], id_name="SampleRate")
    write_log(tmp_path / "plain.log", ["treeHeight"], [[1.0]], id_name="Sample")
    assert read_trace(tmp_path / "rate.log").id_column == "SampleRate"
    assert read_trace(tmp_path / "plain.log").id_column == "Sample"


def test_missing_column_is_an_error_not_a_zero(tmp_path):
    path = write_log(tmp_path / "a.log", ["treeHeight"], [[1.0]])
    with pytest.raises(TraceFormatError, match="does not log the required column"):
        read_trace(path, required_columns=["mrca.age(Gone.taxonset)"])


def test_non_numeric_cell_is_a_read_failure(tmp_path):
    path = tmp_path / "a.log"
    path.write_text("# c\nSample\ttreeHeight\n100\t1.0\n200\tNaN-ish\n", encoding="utf-8")
    with pytest.raises(TraceFormatError, match="non-numeric value"):
        read_trace(path)


def test_ragged_row_is_a_read_failure(tmp_path):
    path = tmp_path / "a.log"
    path.write_text("# c\nSample\ttreeHeight\ntreePrior\n100\t1.0\n200\n", encoding="utf-8")
    with pytest.raises(TraceFormatError, match=r"field\(s\) but the header names"):
        read_trace(path)


def test_duplicate_sample_ids_are_rejected_as_concatenation(tmp_path):
    path = tmp_path / "a.log"
    path.write_text("# c\nSample\ttreeHeight\n100\t1.0\n200\t2.0\n100\t3.0\n", encoding="utf-8")
    with pytest.raises(TraceFormatError, match="concatenated"):
        read_trace(path)


def test_header_without_a_sample_index_fails(tmp_path):
    path = tmp_path / "a.log"
    path.write_text("# c\ntreeHeight\ttreePrior\n1.0\t2.0\n", encoding="utf-8")
    with pytest.raises(TraceFormatError, match=r"has no 'SampleRate or Sample'"):
        read_trace(path)


def test_empty_and_header_only_files_fail(tmp_path):
    empty = tmp_path / "empty.log"
    empty.write_text("# nothing here\n", encoding="utf-8")
    with pytest.raises(TraceFormatError, match="no header line"):
        read_trace(empty)
    header_only = tmp_path / "head.log"
    header_only.write_text("Sample\ttreeHeight\n", encoding="utf-8")
    with pytest.raises(TraceFormatError, match="no samples"):
        read_trace(header_only)


# --- burn-in ----------------------------------------------------------------


def test_burn_in_by_count_and_by_fraction_agree_at_the_boundary(tmp_path):
    rows = [[float(index)] for index in range(100)]
    path = write_log(tmp_path / "a.log", ["treeHeight"], rows)
    by_count = read_trace(path, burn_in=10)
    by_fraction = read_trace(path, burn_in_fraction=0.1)
    assert len(by_count) == len(by_fraction) == 90
    assert by_count.burn_in_removed == by_fraction.burn_in_removed == 10
    assert by_count.n_samples_total == 100
    # floor semantics: 101 samples, 10% -> 10 dropped, matching BEAST2's own
    # post-processing rather than a rounded estimate.
    odd = write_log(tmp_path / "b.log", ["treeHeight"], rows + [[999.0]])
    assert read_trace(odd, burn_in_fraction=0.1).burn_in_removed == 10


def test_burn_in_choices_are_mutually_exclusive(tmp_path):
    path = write_log(tmp_path / "a.log", ["treeHeight"], [[1.0]])
    with pytest.raises(TraceFormatError, match="not both"):
        read_trace(path, burn_in=1, burn_in_fraction=0.1)


def test_burn_in_that_empties_the_trace_is_refused(tmp_path):
    path = write_log(tmp_path / "a.log", ["treeHeight"], [[1.0], [2.0]])
    with pytest.raises(TraceFormatError, match="removes all"):
        read_trace(path, burn_in=5)
    with pytest.raises(TraceFormatError, match=r"must be in \[0, 1\)"):
        read_trace(path, burn_in_fraction=1.0)
    with pytest.raises(TraceFormatError, match="must be >= 0"):
        read_trace(path, burn_in=-1)


# --- combining --------------------------------------------------------------


def test_combine_requires_the_same_columns(tmp_path):
    first = read_trace(write_log(tmp_path / "a.log", ["treeHeight"], [[1.0], [2.0]]))
    second = read_trace(write_log(tmp_path / "b.log", ["treePrior"], [[1.0], [2.0]]))
    with pytest.raises(TraceFormatError, match="must log the same columns"):
        combine_traces([first, second])


def test_combine_refuses_unequal_chain_lengths(tmp_path):
    first = read_trace(write_log(tmp_path / "a.log", ["treeHeight"], [[1.0], [2.0], [3.0]]))
    second = read_trace(write_log(tmp_path / "b.log", ["treeHeight"], [[1.0], [2.0]]))
    with pytest.raises(TraceFormatError, match="unequal chain lengths"):
        combine_traces([first, second])


def test_combine_pools_equal_independent_runs(tmp_path):
    first = read_trace(write_log(tmp_path / "a.log", ["treeHeight"], [[1.0], [2.0]]))
    second = read_trace(write_log(tmp_path / "b.log", ["treeHeight"], [[3.0], [4.0]]))
    pooled = combine_traces([first, second])
    assert list(pooled.column("treeHeight")) == [1.0, 2.0, 3.0, 4.0]
    assert pooled.n_samples_total == 4


# --- statistics -------------------------------------------------------------


def test_ess_is_reported_as_not_estimable_for_degenerate_input():
    assert math.isnan(effective_sample_size([1.0] * 50))
    assert math.isnan(effective_sample_size([1.0, 2.0]))


def test_ess_of_independent_samples_is_close_to_n():
    rng = np.random.default_rng(20261028)
    values = rng.normal(size=4000)
    ess = effective_sample_size(values)
    # The estimator is allowed to be optimistic, not wrong by a factor: an
    # uncorrelated stream must land in the same order as its sample count.
    assert 3000 < ess <= 4000


def test_ess_of_a_autocorrelated_stream_is_below_its_length():
    rng = np.random.default_rng(7)
    values = np.empty(4000)
    values[0] = rng.normal()
    for index in range(1, values.size):
        values[index] = 0.98 * values[index - 1] + 0.2 * rng.normal()
    assert effective_sample_size(values) < values.size / 5


def test_hpd_is_the_shortest_covering_interval():
    values = list(range(1000))
    lower, upper = hpd_interval(values, mass=0.95)
    assert (lower, upper) == (0, 950)
    assert math.isnan(hpd_interval([])[0])
    with pytest.raises(ValueError, match="mass must be"):
        hpd_interval([1.0, 2.0], mass=1.5)


def test_summarise_column_reports_the_row_a_figure_needs(tmp_path):
    rows = [[float(index)] for index in range(1, 21)]
    trace = read_trace(write_log(tmp_path / "a.log", ["treeHeight"], rows))
    summary = summarise_column(trace, "treeHeight", ess_minimum=5.0)
    assert summary.n == 20
    assert summary.median == pytest.approx(10.5)
    # The criterion only labels the row; a low-ESS row is still reported rather
    # than dropped, which is what keeps a not-converged scenario visible.
    assert summary.converged is (summary.ess >= 5.0)
    assert summarise_column(trace, "treeHeight").converged is None


# --- node name resolution ---------------------------------------------------


def test_node_age_column_prefers_the_exact_name_then_the_taxonset_spelling(tmp_path):
    trace = read_trace(
        write_log(
            tmp_path / "a.log",
            ["treeHeight", "mrca.age(A.taxonset)", "mrca.age(B)"],
            [[1.0, 2.0, 3.0]],
        )
    )
    assert node_age_column(trace, "A") == "mrca.age(A.taxonset)"
    assert node_age_column(trace, "B") == "mrca.age(B)"
    assert node_age_column(trace, "treeHeight") == "treeHeight"


def test_ambiguous_and_absent_node_names_are_errors(tmp_path):
    trace = read_trace(
        write_log(
            tmp_path / "a.log",
            ["mrca.age(C.taxonset)", "median.age(C)"],
            [[1.0, 2.0]],
        )
    )
    with pytest.raises(TraceFormatError, match="ambiguous"):
        node_age_column(trace, "C")
    with pytest.raises(TraceFormatError, match="no node-age column"):
        node_age_column(trace, "Missing")


# --- leave-one-out comparison ----------------------------------------------


def test_prior_shift_reports_median_difference_and_width_ratio(tmp_path):
    reference = read_trace(
        write_log(tmp_path / "ref.log", ["treeHeight"], [[float(index)] for index in range(50)])
    )
    scenario = read_trace(
        write_log(
            tmp_path / "scn.log", ["treeHeight"], [[float(index) + 10.0] for index in range(50)]
        )
    )
    shift = prior_shift(reference, scenario, ["treeHeight"], scenario_name="without_X")[0]
    assert shift.median_shift == pytest.approx(10.0)
    assert shift.interval_width_ratio == pytest.approx(1.0)
    assert shift.reference_median == pytest.approx(24.5)


def test_zero_reference_width_gives_no_ratio_instead_of_infinity(tmp_path):
    reference = read_trace(write_log(tmp_path / "ref.log", ["treeHeight"], [[5.0]] * 20))
    scenario = read_trace(
        write_log(tmp_path / "scn.log", ["treeHeight"], [[float(index)] for index in range(20)])
    )
    shift = prior_shift(reference, scenario, ["treeHeight"])[0]
    assert math.isnan(shift.interval_width_ratio)


# --- the CLI that exposes it -------------------------------------------------


class TestPriorSummaryCLI:
    def test_prints_one_row_per_column_and_labels_it_prior_only(self, tmp_path, capsys):
        from beast2py.main import main

        log = write_log(
            tmp_path / "full.log",
            ["treeHeight", "mrca.age(A.taxonset)"],
            [[float(i), float(i) + 5.0] for i in range(40)],
        )
        assert main(["prior-summary", "--trace", str(log)]) == 0
        out = capsys.readouterr().out
        assert out.splitlines()[0].startswith("scenario\tparameter\tmedian")
        assert "treeHeight" in out and "A" in out or "mrca.age(A.taxonset)" in out
        assert "prior-only marginals" in out
        assert "posterior-dominance claim" in out

    def test_reference_adds_the_shift_columns(self, tmp_path, capsys):
        from beast2py.main import main

        reference = write_log(tmp_path / "ref.log", ["treeHeight"], [[float(i)] for i in range(40)])
        scenario = write_log(
            tmp_path / "scn.log", ["treeHeight"], [[float(i) + 3.0] for i in range(40)]
        )
        ret = main(
            [
                "prior-summary",
                "--trace",
                str(scenario),
                "--reference",
                str(reference),
                "--column",
                "treeHeight",
                "--label",
                "without_X",
                "--burn-in-fraction",
                "0.0",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        header = out.splitlines()[0].split("\t")
        assert "median_shift" in header and "reference_median" in header
        row = dict(zip(header, out.splitlines()[1].split("\t")))
        assert row["scenario"] == "without_X"
        assert float(row["median"]) == pytest.approx(22.5)
        assert float(row["reference_median"]) == pytest.approx(19.5)
        assert float(row["median_shift"]) == pytest.approx(3.0)

    def test_a_column_the_scenario_does_not_log_is_a_usage_error(self, tmp_path, capsys):
        """A leave-one-out run has no column for the calibration it removed.

        Reading zero there would put a fake number in a figure, so the command
        exits 2 instead.
        """
        from beast2py.main import main

        log = write_log(tmp_path / "a.log", ["treeHeight"], [[1.0], [2.0]])
        ret = main(["prior-summary", "--trace", str(log), "--column", "Gone"])
        assert ret == 2
        # _error() writes to stderr, like every other refusal in the CLI.
        assert "no node-age column for 'Gone'" in capsys.readouterr().err

    def test_burn_in_options_are_mutually_exclusive_at_the_parser_level(self, tmp_path):
        log = write_log(tmp_path / "a.log", ["treeHeight"], [[1.0]])
        with pytest.raises(SystemExit):
            from beast2py.main import build_parser

            build_parser().parse_args(
                [
                    "prior-summary",
                    "--trace",
                    str(log),
                    "--burn-in",
                    "1",
                    "--burn-in-fraction",
                    "0.1",
                ]
            )

    def test_low_ess_rows_are_reported_not_dropped(self, tmp_path, capsys):
        from beast2py.main import main

        # A monotone ramp is almost perfectly autocorrelated, so its ESS is far
        # below the criterion: the row must still be printed, with a warning.
        log = write_log(tmp_path / "a.log", ["treeHeight"], [[float(i)] for i in range(60)])
        ret = main(
            [
                "prior-summary",
                "--trace",
                str(log),
                "--column",
                "treeHeight",
                "--burn-in-fraction",
                "0.0",
                "--ess-min",
                "1000",
            ]
        )
        captured = capsys.readouterr()
        assert ret == 0
        assert "treeHeight" in captured.out.splitlines()[1]
        assert "below the ESS criterion" in captured.err

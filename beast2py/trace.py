"""Read BEAST2 trace logs and summarise what a prior-only run actually sampled.

Beast2Py generates XML and stops; it runs no MCMC and estimates no posterior.
What it can legitimately do with a BEAST2 output file is read the marginal
columns of a chain that was told to ``sampleFromPrior`` and report the median,
the 95% interval, the effective sample size and the shift each calibration
produces in the **effective marginal prior**.  That is the arithmetic behind the
leave-one-out figure, so it lives here, under tests, rather than in a one-off
script next to the run directories.

Column access is by name throughout: a BEAST2 logger writes ``SampleRate`` plus
one column per traced parameter, and the order depends on the XML, not on this
package.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "TraceFormatError",
    "Trace",
    "read_trace",
    "combine_traces",
    "node_age_column",
    "effective_sample_size",
    "hpd_interval",
    "summarise_column",
    "prior_shift",
]

#: The column BEAST2 names for the sample index.  A ``FileLogFrame``/``Logger``
#: with ``log=true`` writes ``SampleRate``; the plain ``<log>``` used by the MCMC
#: and by ``beast -sampleFromPrior`` writes ``Sample``.  Both are accepted, and a
#: file naming neither is a read failure rather than a mystery offset column.
SAMPLE_ID_COLUMNS = ("SampleRate", "Sample")

#: Kept for callers that speak of the "sample rate column"; it is the header
#: spelling, not a promise about which one a given file uses.
SAMPLE_RATE_COLUMN = SAMPLE_ID_COLUMNS[0]


class TraceFormatError(ValueError):
    """The file is not a trace this reader can turn into numbers.

    Raised instead of silently dropping rows: a scenario whose trace is
    truncated, mis-ordered or partially non-numeric must fail, because the
    summary would otherwise look exactly like a successful read.
    """


def _sample_id_column(columns: Sequence[str], path: Optional[Path] = None) -> str:
    """Which header field carries the sample index."""
    for name in SAMPLE_ID_COLUMNS:
        if name in columns:
            return name
    raise TraceFormatError(
        f"{_where(path)} header has no {' or '.join(SAMPLE_ID_COLUMNS)!r} column; "
        f"found {', '.join(list(columns)[:12])}"
    )


@dataclass
class Trace:
    """One parsed trace, with burn-in already applied.

    Attributes:
        columns: Column names in file order, sample index first when present.
        sample_ids: The sample-index values kept after burn-in.
        data: Mapping of column name to the kept samples.
        n_samples_total: Rows before burn-in.
        burn_in_removed: Rows dropped as burn-in.
        source: Where the trace was read from, for the report.
        id_column: Which header field was used as the sample index.
    """

    columns: List[str]
    sample_ids: np.ndarray
    data: Dict[str, np.ndarray]
    n_samples_total: int
    burn_in_removed: int
    source: Optional[str] = None
    id_column: str = SAMPLE_ID_COLUMNS[0]

    def __len__(self) -> int:
        return int(self.sample_ids.shape[0])

    def column(self, name: str) -> np.ndarray:
        """The kept samples of one column, addressed by its header name."""
        try:
            return self.data[name]
        except KeyError:
            raise TraceFormatError(
                f"column {name!r} is not in the trace"
                + (f" ({self.source})" if self.source else "")
                + f"; the file holds {len(self.columns)} column(s): "
                + ", ".join(self.columns[:12])
                + (" ..." if len(self.columns) > 12 else "")
            ) from None

    @property
    def numeric_columns(self) -> List[str]:
        return [name for name in self.columns if name != self.id_column]


def _split_header(line: str, path: Optional[Path]) -> List[str]:
    """Split a BEAST2 header on tabs, falling back to whitespace.

    BEAST2 writes tabs.  A log re-saved through an editor that expanded them is
    still one header per column, so whitespace splitting is a fallback rather
    than a guess: it is only reached when the line contains no tab at all.
    """
    if "\t" in line:
        return [name.strip() for name in line.split("\t") if name.strip()]
    names = [name for name in line.split() if name]
    if len(names) < 2:
        raise TraceFormatError(
            f"{_where(path)} header line names a single column ({names[0] if names else 'empty'});"
            " a BEAST2 trace log has a tab-separated header"
        )
    return names


def _where(path: Optional[Path]) -> str:
    return str(path) if path else "in-memory trace"


def read_trace(
    path: Path,
    burn_in: Optional[float] = None,
    burn_in_fraction: Optional[float] = None,
    required_columns: Sequence[str] = (),
    source: Optional[str] = None,
) -> Trace:
    """Parse one BEAST2 trace log.

    Args:
        path: The ``.log`` file.
        burn_in: Number of leading samples to drop.  Mutually exclusive with
            ``burn_in_fraction``.
        burn_in_fraction: Fraction of the logged samples to drop, e.g. ``0.1``.
        required_columns: Columns the caller needs; absence is an error rather
            than a silently empty summary.
        source: Label used in error messages in place of the file path.

    Raises:
        TraceFormatError: no header, no data rows, a ragged row, a non-numeric
            value, a duplicate sample id, a bad burn-in, or a missing required
            column.
    """
    if burn_in is not None and burn_in_fraction is not None:
        raise TraceFormatError(
            "pass either burn_in (a sample count) or burn_in_fraction (a "
            "fraction), not both: that is two different answers to one question"
        )
    path = Path(path)
    header: Optional[List[str]] = None
    rows: List[List[float]] = []
    sample_ids: List[int] = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n").rstrip("\r")
            stripped = line.strip()
            # BEAST2 puts a comment block above the header and blank lines
            # between blocks; both are structure, not data.
            if not stripped or stripped.startswith("#"):
                continue
            if header is None:
                header = _split_header(line, path)
                id_column = _sample_id_column(header, path)
                rate_index = header.index(id_column)
                continue
            fields = line.split("\t") if "\t" in line else stripped.split()
            if len(fields) != len(header):
                raise TraceFormatError(
                    f"{_where(path)} line {line_number} has {len(fields)} "
                    f"field(s) but the header names {len(header)}: a truncated "
                    "or re-saved trace log is not a partial result, it is a "
                    "read failure"
                )
            try:
                values = [float(field) for field in fields]
            except ValueError as error:
                bad = next(
                    (f for f in fields if _not_numeric(f)),
                    "?",
                )
                raise TraceFormatError(
                    f"{_where(path)} line {line_number} holds a non-numeric "
                    f"value {bad!r} in column "
                    f"{header[min(fields.index(bad), len(header) - 1)] if bad != '?' else '?'}: "
                    f"({error})"
                ) from error
            sample_ids.append(int(values[rate_index]))
            rows.append(values)

    if header is None:
        raise TraceFormatError(f"{_where(path)} contains no header line")
    if not rows:
        raise TraceFormatError(f"{_where(path)} contains a header but no samples")

    duplicates = _duplicate_sample_ids(sample_ids)
    if duplicates:
        shown = ", ".join(str(value) for value in duplicates[:5])
        raise TraceFormatError(
            f"{_where(path)} repeats sample id(s) {shown}; two chains have been "
            "concatenated into one file. Keep them separate and combine them "
            "with combine_traces(), which checks that they are the same length"
        )

    matrix = np.asarray(rows, dtype=float)
    rate_index = header.index(id_column)
    data = {name: matrix[:, index].copy() for index, name in enumerate(header)}
    trace = Trace(
        columns=list(header),
        sample_ids=matrix[:, rate_index].astype(int),
        data=data,
        n_samples_total=len(rows),
        burn_in_removed=0,
        source=source or str(path),
        id_column=id_column,
    )
    trace = _apply_burn_in(trace, burn_in, burn_in_fraction, path)

    missing = [name for name in required_columns if name not in trace.data]
    if missing:
        raise TraceFormatError(
            f"{_where(path)} does not log the required column(s): "
            + ", ".join(missing)
            + ". Either the XML has no logger for them or the calibration was "
            "removed in this scenario; do not read a zero where there was no column."
        )
    return trace


def _not_numeric(field: str) -> bool:
    try:
        float(field)
    except ValueError:
        return True
    return False


def _duplicate_sample_ids(sample_ids: List[int]) -> List[int]:
    seen = set()
    duplicates = []
    for value in sample_ids:
        if value in seen and value not in duplicates:
            duplicates.append(value)
        seen.add(value)
    return duplicates


def _apply_burn_in(
    trace: Trace,
    burn_in: Optional[float],
    burn_in_fraction: Optional[float],
    path: Optional[Path],
) -> Trace:
    """Drop leading samples by count or fraction, in place on a copy."""
    if burn_in is None and burn_in_fraction is None:
        return trace
    if burn_in is not None:
        removed = int(burn_in)
        if removed < 0:
            raise TraceFormatError(f"{_where(path)} burn_in must be >= 0")
    else:
        # read_trace() has already returned when neither knob was given, so
        # the only way here is a fraction; 0.0 is its own no-op value.
        fraction = float(burn_in_fraction or 0.0)
        if not 0.0 <= fraction < 1.0:
            raise TraceFormatError(
                f"{_where(path)} burn_in_fraction must be in [0, 1), got {fraction}"
            )
        # Same rule BEAST2's own post-processing uses: the fraction counts the
        # logged samples, so the boundary is a floor, not a round.
        removed = int(math.floor(fraction * trace.n_samples_total))
    if removed >= trace.n_samples_total:
        raise TraceFormatError(
            f"{_where(path)} burn-in removes all {trace.n_samples_total} "
            "sample(s); shorten it or lengthen the chain"
        )
    kept = slice(removed, None)
    return Trace(
        columns=list(trace.columns),
        sample_ids=trace.sample_ids[kept],
        data={name: values[kept] for name, values in trace.data.items()},
        n_samples_total=trace.n_samples_total,
        burn_in_removed=removed,
        source=trace.source,
        id_column=trace.id_column,
    )


def combine_traces(traces: Sequence[Trace], *, independent: bool = True) -> Trace:
    """Stack independent runs of the same scenario column-wise identically.

    Two chains of one scenario are only comparable when they logged the same
    columns; and unequal lengths mean one run stopped early, which changes what
    a pooled mean means.  Both are errors, not warnings.

    Args:
        traces: Traces to pool, burn-in already applied to each.
        independent: When ``True`` (the default) duplicate sample ids across
            runs are allowed, because they are separate chains; within one file
            they are a concatenation and ``read_trace`` already refused that.
    """
    if not traces:
        raise TraceFormatError("combine_traces() got no traces to combine")
    first = traces[0]
    for other in traces[1:]:
        if other.columns != first.columns:
            raise TraceFormatError(
                f"{other.source} logs {', '.join(other.columns[:8])} but "
                f"{first.source} logs {', '.join(first.columns[:8])}; two runs "
                "of one scenario must log the same columns before they can be pooled"
            )
        if len(other) != len(first):
            raise TraceFormatError(
                f"{other.source} kept {len(other)} sample(s) and {first.source} "
                f"kept {len(first)}; unequal chain lengths pooled would weight "
                "the longer run more heavily. Truncate deliberately or report "
                "them as separate runs."
            )
    columns = list(first.columns)
    id_column = first.id_column
    stacked = {name: np.concatenate([trace.data[name] for trace in traces]) for name in columns}
    sample_ids = np.concatenate([trace.sample_ids for trace in traces])
    if not independent:
        duplicates = _duplicate_sample_ids(list(sample_ids))
        if duplicates:
            raise TraceFormatError(f"pooled traces repeat sample id(s) {duplicates[:5]}")
    return Trace(
        columns=columns,
        sample_ids=sample_ids,
        data=stacked,
        n_samples_total=sum(trace.n_samples_total for trace in traces),
        burn_in_removed=sum(trace.burn_in_removed for trace in traces),
        source=", ".join(str(trace.source) for trace in traces),
        id_column=id_column,
    )


def node_age_column(trace: Trace, name: str) -> str:
    """Resolve a calibration name to its node-age column, by header name.

    BEAST2 names the logged node by the *taxon set* id, so a calibration
    ``cladeA`` whose taxon set is written as ``cladeA.taxonset`` is logged
    as ``mrca.age(cladeA.taxonset)``, while an inline taxon set may be logged
    as ``mrca.age(cladeA)``.  Matching on the column index instead would
    silently follow whatever order the XML happened to produce, so the search is
    on names: an exact hit first, then the ``.taxonset`` spelling, then any
    single column whose ``mrca.age(...)``/``median.age(...)`` label starts with
    the requested name.  Zero or several candidates are errors.

    Raises:
        TraceFormatError: the name matches no column, or matches more than one.
    """
    if name in trace.numeric_columns:
        return name
    candidates = []
    for column in trace.numeric_columns:
        for prefix in ("mrca.age(", "median.age("):
            if column.startswith(prefix) and column.endswith(")"):
                label = column[len(prefix) : -1]
                if label == name or label.split(".")[0] == name:
                    candidates.append(column)
    exact = [column for column in candidates if column == f"mrca.age({name})"]
    if exact:
        return exact[0]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise TraceFormatError(
            f"no node-age column for {name!r} in {trace.source}; "
            f"the trace logs {len(trace.numeric_columns)} column(s). A "
            "leave-one-out run does not log the node of the calibration it "
            "removed (its MRCAPrior is gone from the XML), so that column has to "
            "come from the .trees file, not from this log"
        )
    raise TraceFormatError(
        f"{name!r} is ambiguous in {trace.source}: "
        + ", ".join(sorted(candidates))
        + ". Use the full column name instead of the calibration name."
    )


def effective_sample_size(values: Iterable[float]) -> float:
    """Effective sample size by the initial positive sequence estimator.

    ``ESS = n / (1 + 2 * sum of the positive autocovariances)`` with the sum
    truncated at the first non-positive autocorrelation: the estimator is stable
    on the positively-correlated node-age columns a prior-only run produces, and
    it errs on the optimistic side exactly as the published BEAST2-trace
    literature says it does.  The returned value is what the run's ESS criterion
    is judged against; it is not a statement about convergence.

    Returns ``nan`` for a trace too short to estimate, and ``inf`` is never
    returned for a constant column: a constant has zero variance, and reporting
    an infinite sample size for it would be a lie about information.
    """
    x = np.asarray(list(values), dtype=float)
    n = x.shape[0]
    if n < 10:
        return float("nan")
    x = x - x.mean()
    variance = x.var(ddof=1)
    if variance == 0:
        return float("nan")
    autocovariance = np.correlate(x, x, mode="full")[n - 1 :] / (variance * n)
    total = 1.0
    for k in range(1, n):
        if autocovariance[k] <= 0:
            break
        total += 2.0 * autocovariance[k]
    return float(n / total)


def hpd_interval(values: Iterable[float], mass: float = 0.95) -> Tuple[float, float]:
    """Shortest interval holding ``mass`` of the samples, ``(lower, upper)``.

    The equal-tailed interval is what most readers assume from "95% HPD"; this
    is the maximum-a-posteriori-density interval, computed as the shortest
    covering window over the sorted samples, which is the definition BEAST2's
    own log summing uses.  ``(nan, nan)`` for no samples.
    """
    if not 0.0 < mass < 1.0:
        raise ValueError(f"mass must be in (0, 1), got {mass}")
    x = np.sort(np.asarray(list(values), dtype=float))
    n = x.shape[0]
    if n == 0:
        return (float("nan"), float("nan"))
    keep = max(1, int(math.floor(mass * n)))
    widths = x[keep:] - x[: n - keep]
    if widths.size == 0:
        return (float(x[0]), float(x[-1]))
    start = int(np.argmin(widths))
    return (float(x[start]), float(x[start + keep]))


@dataclass
class ColumnSummary:
    """One column of one scenario, as a figure or table reports it."""

    column: str
    median: float
    lower: float
    upper: float
    ess: float
    n: int
    converged: Optional[bool] = None
    scenario: str = ""

    def as_row(self) -> Dict[str, object]:
        return {
            "scenario": self.scenario,
            "parameter": self.column,
            "median": self.median,
            "hpd_lower": self.lower,
            "hpd_upper": self.upper,
            "ess": self.ess,
            "n": self.n,
            "converged": self.converged,
        }


def summarise_column(
    trace: Trace,
    column: str,
    ess_minimum: Optional[float] = None,
    scenario: str = "",
) -> ColumnSummary:
    """Median, 95% interval and ESS for one column.

    ``ess_minimum`` only labels the row; deciding a run converged is the
    caller's criterion, and a row below it is still reported rather than hidden.
    """
    values = trace.column(column)
    lower, upper = hpd_interval(values)
    ess = effective_sample_size(values)
    return ColumnSummary(
        column=column,
        median=float(np.median(values)),
        lower=lower,
        upper=upper,
        ess=ess,
        n=int(values.shape[0]),
        converged=None if ess_minimum is None else bool(ess >= ess_minimum),
        scenario=scenario,
    )


@dataclass
class PriorShift:
    """How far one scenario's effective prior moved from the reference."""

    column: str
    reference_median: float
    scenario_median: float
    median_shift: float
    reference_lower: float
    reference_upper: float
    scenario_lower: float
    scenario_upper: float
    interval_width_ratio: float
    scenario_ess: float
    scenario: str = ""

    def as_row(self) -> Dict[str, object]:
        return {
            "scenario": self.scenario,
            "parameter": self.column,
            "reference_median": self.reference_median,
            "median": self.scenario_median,
            "median_shift": self.median_shift,
            "reference_hpd_lower": self.reference_lower,
            "reference_hpd_upper": self.reference_upper,
            "hpd_lower": self.scenario_lower,
            "hpd_upper": self.scenario_upper,
            "interval_width_ratio": self.interval_width_ratio,
            "ess": self.scenario_ess,
        }


def prior_shift(
    reference: Trace,
    scenario: Trace,
    columns: Sequence[str],
    *,
    scenario_name: str = "",
) -> List[PriorShift]:
    """Per-column shift of a leave-one-out scenario against the full run.

    These are differences between prior-only marginals.  They say what removing
    a calibration did to the effective marginal prior; they say nothing about a
    posterior, and the caller must not relabel them as such.
    """
    shifts: List[PriorShift] = []
    for column in columns:
        reference_values = reference.column(column)
        scenario_values = scenario.column(column)
        reference_lower, reference_upper = hpd_interval(reference_values)
        scenario_lower, scenario_upper = hpd_interval(scenario_values)
        reference_width = reference_upper - reference_lower
        scenario_width = scenario_upper - scenario_lower
        # A zero reference width would make the ratio meaningless rather than
        # infinite: report nan and let the table say so.
        ratio = scenario_width / reference_width if reference_width > 0 else float("nan")
        shifts.append(
            PriorShift(
                column=column,
                reference_median=float(np.median(reference_values)),
                scenario_median=float(np.median(scenario_values)),
                median_shift=float(np.median(scenario_values) - np.median(reference_values)),
                reference_lower=reference_lower,
                reference_upper=reference_upper,
                scenario_lower=scenario_lower,
                scenario_upper=scenario_upper,
                interval_width_ratio=float(ratio),
                scenario_ess=effective_sample_size(scenario_values),
                scenario=scenario_name,
            )
        )
    return shifts

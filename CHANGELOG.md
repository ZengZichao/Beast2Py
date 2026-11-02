# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Beast2Py is pre-1.0. Releases are published to PyPI as
[`beast2py`](https://pypi.org/project/beast2py/) -- install a release with
`pip install beast2py` -- and each release is archived on Zenodo -- see
`CITATION.cff` and `.zenodo.json` for the concept DOI
10.5281/zenodo.23088756, which always resolves to the latest archived version,
with the per-version DOIs recorded in `CITATION.cff`. It is developed in the
open at https://github.com/ZengZichao/Beast2Py. Version numbers follow SemVer,
and a change to generated BEAST2 XML or to a validation gate is treated as a
change to observed behaviour even if no API changes.

## [0.2.0] - 2026-10-27

### Changed

- Calibration diagnostics are now structural and crown/stem-aware. The three
  blocking rules are same-node contradictions, nested temporal order and
  intersecting non-nested monophyly. Root + originate, tipsonly + originate and
  originate on a whole-ingroup taxon set are parse-time errors.
- Distribution similarity between independent clades is no longer a conflict.
  It is an opt-in descriptive report (`diagnostics.prior_similarity`) with an
  explicit, uncalibrated threshold; it carries no merge recommendation and
  never blocks Gate 2. Undecidable pairs are reported as `unchecked`.
- Beta (without offset), Poisson and OneOnX now raise distribution-usage
  warnings; `diagnostics.strict_distributions: true` turns them into errors.
  Beta + offset remains legal.
- All formal examples use explicitly synthetic calibration values and
  provenance cards. No literature-derived calibration values, DOIs, fossil
  ages, or third-party sequence data are redistributed with the repository.
- Sensitivity XMLs use the configured chain length and distinct seeds instead
  of a 1M-generation cap.
- Coverage: 36/36 declared capabilities have an external Gate 3 record; fixtures cover
  the 8 calibration distributions, 4 amino-acid matrices and
  coalescent_constant; a synthetic non-zero tip-date fixture covers the dated
  tip operator.
- The prior-only sensitivity and trace-summary workflow is exercised with
  synthetic inputs; no empirical dataset or archived run is bundled.

### Fixed

- ChiSquare `df` is emitted as an integer attribute, as BEAST2's
  IntegerParameter requires.
- CalibratedYule CalibrationPoint now writes `forParent` when use_originate is
  true.
- Birth-death supports `conditional_on_root`.
- The distribution-usage warning no longer treats a soft minimum prior as an
  error when its upper tail extends beyond the fossil's age bracket.

## [Unreleased]

### Added

- `beast2py prior-summary`, and `beast2py/trace.py` behind it: a tested reader for
  BEAST2 trace logs. It addresses columns by header name (`Sample` or
  `SampleRate`, and `mrca.age(<taxonset>)` resolved from a calibration name),
  applies a burn-in by count or fraction, reports median / 95% highest-density
  interval / initial-positive-sequence ESS / sample count per column, and with
  `--reference` the shift of a leave-one-out scenario against the full
  calibration set including the interval width ratio. It refuses a missing
  column, a non-numeric cell, a ragged row, duplicated sample ids (a
  concatenated file) and unequal chain lengths when pooling, and it labels a
  low-ESS row rather than dropping it. The archived leave-one-out table is
  recomputed from the raw logs by `tests/test_trace.py`, so every figure number
  is reproducible from the traces. This is prior-only arithmetic: it estimates no
  posterior, and both its help text and its output say so.
- Gate 2 now states its own boundary. `DETECTABLE_CHECKS`,
  `NOT_DETECTABLE_CHECKS`, `NOT_DETECTABLE_LABELS` and `PRIOR_LAYERS` in
  `diagnostics.py` are the single source for what the static gate decides and
  what it cannot see (effective-marginal-prior distortion, posterior dominance,
  MCMC non-convergence, insufficient data), rendered into the `generate`/`quick`
  Gate 2 line, the full `diagnose` listing and the HTML report, and mirrored in
  both READMEs and both user manuals with the four-term vocabulary (specified
  marginal prior / effective marginal prior / joint effective prior / posterior).
- `tests/test_disclosed_regressions.py`: a registry of every defect disclosed in
  the release notes, each mapped to the tests that guard it, validated against
  the collected suite so a renamed or deleted regression case fails the registry.
- `.github/actions/setup-beast2`: the BEAST2 runtime install the CI integration
  job had inline, now shared with the release path.

### Changed

- The Release workflow gates a tag on the suite, not on version strings alone:
  `build` now needs a `test` job that installs BEAST2, runs the whole suite, and
  fails the release if any Gate 3 test skipped for want of a runtime.
- Documentation that described behaviour the code does not have: the English
  manuals promised a warning when only the prior means of a nested pair disagree
  (the detector reports an error only when the parent's interval is wholly
  younger, and deliberately not otherwise), and both manuals told readers to run
  a short chain and compare "tree-height posteriors" for the leave-one-out
  design. Those chains sample from the prior; the wording now says so and asks
  for a stated ESS criterion instead of an assumed-adequate chain length.
- The similarity report ships no threshold at all: the old
  `OVERLAP_RATIO_THRESHOLD = 0.15` / `OVERLAP_DISTANCE_THRESHOLD = 2.0` class
  constants are gone rather than kept as "historical names", and the manual and
  example configuration no longer show 0.15 as if it were a value with meaning.

### Fixed

- The CI lint gate was red before any of this: the tree carried findings that
  `black --check .` and `flake8 beast2py tests` reject (unused imports, two
  locals computed and never read, lines over the 100-character limit), and the
  mypy ratchet sat at `MYPY_BASELINE: 20` while `mypy` reported 32 errors, so
  every run failed and the ratchet guarded nothing. The backlog is cleared --
  including `Optional` taxon sets compared with `<` and `&`, a stale `Optional`
  return on `_parse_calibration` (which always returns or raises), a `str`
  parameter rebound to `Path`, `plt.cm.tab10` outside the stubs and a
  heterogeneous parameter-range table -- and the baseline is now 0, with a test
  that fails if it is raised again.
- `use_originate` on a `CalibratedYule` point now has tests: `forParent="true"`
  is asserted in the emitted XML, the crown/stem truth table covers all sixteen
  combinations of node reading against taxon-set relation, and a stem reading of
  the root is refused by both builders.
- The READMEs recomputed their test counts from the collected suite, and a guard
  now fails if a documented count stops matching `pytest --collect-only`.
- The provenance placeholder check that lived only in the submission bundle's
  `verify.py` is now part of the suite: every calibration in `examples/` must
  carry a card with an allowed `source`, a `calibration_type`, notes long enough
  to explain how the boundary became the prior, and a DOI-bearing reference; no
  formal example may contain `TBD`, `placeholder` or `no source consulted`.
- Zenodo failed to archive v0.1.2 ("Extra metadata load failed"): the audit
  batch had moved `.zenodo.json`'s licence to the compound SPDX expression
  `MIT AND LGPL-2.1-only`, but Zenodo's licence vocabulary lists single
  licences only -- the expression is fine in text fields, not as the licence
  id. The archive metadata declares the package's own MIT licence in
  `license.id` again (the value v0.1.1 archived successfully with) and carries
  the full compound in `notes` as well as `description`; the release-integrity
  guard pins that arrangement and documents the constraint instead of
  rejecting the bare id. PyPI and the GitHub release are unaffected. The
  v0.1.2 record was published directly under the concept record as
  [10.5281/zenodo.23452749](https://doi.org/10.5281/zenodo.23452749) and is
  recorded in `CITATION.cff`.

## [0.1.2] - 2026-10-18

Documentation-and-reporting release: no change to generated BEAST2 XML
semantics or validation rules. The CLI now prints the four-check scheme the
documentation describes, and `validate --beast2` fails closed the way
`generate --beast2-validate` already did. 0.1.1 became the first PyPI release,
published 2026-10-18. Tag and release
[v0.1.2](https://github.com/ZengZichao/Beast2Py/releases/tag/v0.1.2), published
2026-10-18.

### Changed

- Gate numbering (2026-10-08, audit comment LANG-012): the CLI now prints the four
  numbered checks the documentation describes. The parse-time validation of the
  input -- the calibration provenance rules and the alignment character
  validation raised by `ConfigParser.parse` -- is labelled
  `Gate 0 (parse-time, before the 3 generation gates)`, on both the passing and
  the failing path of `generate` and `quick`, while the three generation gates
  keep their existing `Gate 1/3`, `Gate 2/3` and `Gate 3/3` lines: the `/3`
  denominator counts the generation gates, and the Gate 0 label says so.
  `validate` keeps its two-level vocabulary and now names the mapping on every
  line (`Level 1/2 passed (Gate 1, structural checks): …`,
  `Level 2/2 passed (Gate 3): …`, `Level 2/2 not run (Gate 3, …)`), so a
  `validate` log and a `generate` log can be read against the same numbering.
  Labels and wording only: no exit code, no flag name and no gate decision
  changed, and a Gate 0 failure is still a `Configuration error` that exits 1
  and writes nothing.
- Documentation: the manuals, READMEs and architecture guides state the
  four-check scheme, and the release-environment note now records that
  MODEL_SELECTION -- which NS 1.2.0's `version.xml` declares as a dependency --
  was **not** installed there while Gate 3 still passed, because BEAST2's
  `XMLParser` resolves the classes a file names and does not enforce add-on
  package dependencies at parse time. That is a statement about the parser, not
  a guarantee about a run.
- `--allow-unvalidated` help text on both `generate` and `validate` now states
  the symmetry explicitly: one flag name, one refusal (exit 2 when a requested
  BEAST2 check could not run), one opt-in downgrade.

### Fixed

- `tests/test_release_integrity.py` imported `tomllib`, which is stdlib only from
  Python 3.11, so the licence-classifiers guard errored on the 3.10 legs of the CI
  matrix the first time this branch ran there. The 3.11+ path still reads the parsed
  `pyproject.toml`; on 3.10 the guard falls back to a line-anchored scan that a real
  classifier entry -- a quoted string alone on its line -- trips and the explanatory
  comment does not.
- `validate --beast2` refuses with exit code 2 and a `Level 2/2 not run` report
  when BEAST2 is requested but cannot be reached, which is the same fail-closed
  behaviour `generate --beast2-validate` already had; `--allow-unvalidated` is
  the explicit structural-only opt-in that restores exit 0. Regression tests now
  pin the refusal from both the in-process CLI and a stripped environment with
  an empty `HOME`, a missing `JAVA_HOME` and an unusable `BEAST2_JAR`.
- CI: the BEAST2 `integration` job now also fetches the MODEL_SELECTION add-on
  that NS 1.2.0 declares as a dependency, best-effort like the BDSKY and NS
  fetches, with a comment explaining that the parse gate does not enforce it.
  The matrix (two OSes x Python 3.10/3.11/3.12) and the black/flake8 gate are
  unchanged.

### Added

- PyPI publication is live: 0.1.1 is available as
  [`beast2py`](https://pypi.org/project/beast2py/), uploaded by the
  repository's Release workflow through PyPI Trusted Publishing. The READMEs
  carry a PyPI badge and name `pip install beast2py` as the primary install
  path; the roadmap entry and the CHANGELOG header that said "not published on
  PyPI" are gone, CONTRIBUTING.md documents the now-enabled pipeline, and
  SECURITY.md notes that the published 0.1.x versions pin users, so an
  advisory must name the earliest affected version.
- Tests for the printed gate numbering: Gate 0 is asserted to precede the three
  generation gates in order, and both parse-time failure modes -- a provenance
  refusal and an illegal-character refusal -- are asserted to be reported as a
  Gate 0 condition that exits 1 without writing and without running a generation
  gate (`tests/test_cli.py`).
- Zenodo archiving is live: the concept DOI is
  [10.5281/zenodo.23088756](https://doi.org/10.5281/zenodo.23088756) and v0.1.1 is
  archived as [10.5281/zenodo.23088757](https://doi.org/10.5281/zenodo.23088757).
  Both DOIs are recorded in `CITATION.cff`; the READMEs (en/zh) carry a concept-DOI
  badge and updated citation guidance, and `CONTRIBUTING.md` documents the
  per-release archiving checks.

## [0.1.1] - 2026-10-02

Maintenance release: no change to generated BEAST2 XML semantics, validation
rules, or the Python API. Tag and release
[v0.1.1](https://github.com/ZengZichao/Beast2Py/releases/tag/v0.1.1), published
2026-10-02.

### Added

- A structured bug report issue form, a pull request template, and a changelog.
- A security disclosure policy (`SECURITY.md`) and a `CODEOWNERS` file naming the
  sole maintainer.
- The tag-to-release path is now fully automatic: when a tag push finds no
  existing GitHub release, the Release workflow creates one and attaches the
  built sdist and wheel.

### Changed

- CI gate configuration: the lint gate, the test matrix, the BEAST2 integration
  job that installs a real BEAST2 runtime, and the packaging job that builds a
  wheel and checks the bundled validator and console entry point.
- Repository security settings for dependency updates and secret scanning.
- Packaging metadata migrated to standards-based fields. The vendored BEAST2
  validation helper continues to ship as package data, and its LGPL-2.1 notice
  and licence text continue to ship with it.
- Documentation: contributor and security pages added. Documentation remains
  bilingual, so every `docs/*_en.md` change is mirrored in `docs/*_zh.md`.
- CI moved off the deprecated Node 20 runtime: the pinned `actions/checkout`,
  `actions/setup-python`, `actions/setup-java` and `codecov/codecov-action`
  are bumped to their current major versions across all three workflows.
  Major bumps stay excluded from Dependabot, so this was done deliberately
  and verified by the CI gate.
- The opt-in PyPI publish job authenticates via Trusted Publishing (OIDC)
  instead of a stored `PYPI_API_TOKEN`: enabling publication now means
  registering the repository as a trusted publisher on pypi.org and setting
  the `PYPI_PUBLISH` variable, with no long-lived token to rotate. The
  publish action is pinned to a commit SHA.
- Version strings across docs, example configurations and docstrings updated
  to this release.

### Fixed

- `ConfigParser` hardcoded `"0.1.0"` as the default `metadata.tool_version`,
  bypassing `beast2py.__version__` (the declared single source of truth), so
  fingerprints kept a stale version suffix after any release bump. The default
  now reads `__version__`, and the committed basic-example XML, fingerprint
  sidecar and methods paragraph were regenerated with it.

## [0.1.0] - 2026-09-30

First public release. Tag and release
[v0.1.0](https://github.com/ZengZichao/Beast2Py/releases/tag/v0.1.0), published
2026-09-30.

### Added

- **YAML-driven BEAST2 XML generation.** Declarative YAML configurations produce
  BEAST2 XML through the `generate` and `quick` commands, the `beast2py`
  console entry point, and the `Beast2Py` Python API. Nothing is written to disk
  unless every validation gate that ran has passed, and an existing output file
  is refused without `--force`.
- **Three-gate validation.** Structural XML checks (well-formedness, id/idref
  consistency, duplicate ids, initial values inside their bounds, sequence
  sanity, a `Prior` on every estimated parameter, an operator on every estimated
  state node, and no operator on a fixed parameter); calibration-prior conflict
  detection; and an optional headless BEAST2 `XMLParser` check behind
  `--beast2-validate` (`--beast2` for `validate`), which runs on a machine with
  BEAST2 2.7.x and JDK 17 or newer and reports a missing BEAST2 installation as
  unavailable rather than as a pass.
- **Calibration conflict detection.** Temporal consistency between nested clades
  including root calibrations, distribution overlap, and monophyly conflicts,
  reported by severity through the `diagnose` command and the API.
- **Prior sensitivity analysis.** Automatic generation of prior-only XML files so
  the marginal effect of each calibration point can be inspected.
- **Analysis fingerprints.** A deterministic identifier of the form
  `B2P-{sha256(config)[:12]}-{version}` with no date component, plus a separate
  data digest over the alignment content, so the same configuration and
  alignment yield the same bytes on any machine at any time.
- **Generated LaTeX methods descriptions.** A publication-ready methods paragraph
  produced from the configuration.
- **Pipeline generation.** Snakemake and Nextflow workflow files emitted from
  the same configuration.
- **Calibration provenance.** Each calibration point records its fossil source,
  reference DOI, and whether it is a hard or soft bound.
- **Diagnostics report.** An HTML report with density plots for the calibration
  distributions.
- **Bilingual documentation.** `README.md` / `README_zh.md` and five manuals in
  `docs/`, each with an English and a Chinese version, plus 19 example YAML
  configurations and example alignments under `examples/`.
- **Test suite.** Unit, semantic and release-integrity tests, plus 21 BEAST2
  integration tests that generate XML for every example configuration and check
  it against the real BEAST 2.7.8 parser. The integration tests skip, visibly
  under `pytest -rs`, when BEAST2 or JDK 17 is unavailable.
- **CI.** A GitHub Actions `CI` workflow running lint and tests on
  ubuntu-latest and macos-latest across Python 3.10, 3.11 and 3.12, an
  `integration` job that installs a real BEAST2 runtime, and a `packaging` job
  that builds a wheel and asserts the bundled validator and console entry point
  are present.
- **Machine-readable citation metadata** in `CITATION.cff` and a `.zenodo.json`
  describing the software for a later Zenodo submission.

### Known limitations at 0.1.0

- Verification is against BEAST2 v2.7.8 only. Other 2.7.x releases are expected
  to work but are untested; BEAST 2.6.x is not supported.
- `config_bd_skyline.yaml` and `config_nested_sampling.yaml` need the BDSKY and
  nested-sampling add-on packages. A green CI run does not by itself show that
  these two configurations cleared the BEAST2 gate, because a failed add-on
  download downgrades to a skip. Run `pytest tests -rs` to see which integration
  cases executed.
- No validation gate can tell you whether a model answers your biological
  question. A model that initialises cleanly may still be the wrong model.

[Unreleased]: https://github.com/ZengZichao/Beast2Py/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/ZengZichao/Beast2Py/releases/tag/v0.1.1
[0.1.0]: https://github.com/ZengZichao/Beast2Py/releases/tag/v0.1.0

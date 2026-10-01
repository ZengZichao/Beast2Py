# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Beast2Py is pre-1.0. It is not published on PyPI and has no Zenodo DOI; install
it from a clone of this repository. Version numbers follow SemVer, and a change
to generated BEAST2 XML or to a validation gate is treated as a change to
observed behaviour even if no API changes.

## [Unreleased]

Maintenance only. No user-facing behaviour change is recorded here yet.

### Added

- A structured bug report issue form, a pull request template, and a changelog.
- A security disclosure policy (`SECURITY.md`) and a `CODEOWNERS` file naming the
  sole maintainer.

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

### Fixed

- None yet.

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

[Unreleased]: https://github.com/ZengZichao/Beast2Py/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/ZengZichao/Beast2Py/releases/tag/v0.1.0

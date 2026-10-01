# Contributing to Beast2Py

Thanks for reading this, and thank you if you are opening a pull request. Beast2Py
is used in real phylogenetic studies, so the bar for a change is that it is
correct, tested and honestly described.

## Before you start: this is a one-person project

Beast2Py is maintained by a single researcher. Two practical consequences:

- Responses are slow. Review usually takes days, sometimes weeks, and happens
  around a full-time job. An open PR is not ignored; it is queued.
- Small, well-scoped pull requests are far more likely to land than large ones. A
  focused change can be reviewed carefully; a sweeping refactor cannot, and will
  usually be asked to be split.

To discuss an approach before writing code, open an issue first. It costs little
and avoids a week of work in a direction that cannot be taken.

## Development setup

```bash
git clone https://github.com/ZengZichao/Beast2Py.git
cd Beast2Py
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
pip install -e ".[dev]"
```

The `dev` extra installs pytest, pytest-cov, flake8, black, mypy and types-PyYAML, all
pinned to exact versions. The pins are deliberate: these tools decide whether the CI gate
passes, so a silent minor bump must not be able to turn a green build red. Python 3.10
or newer is required. Everything except the optional BEAST2 validation gate is
pure Python, so a plain virtualenv is enough. To confirm the install, run
`beast2py --version` and `beast2py list-models`; without the console script,
`python3 -m beast2py.main` is equivalent.

## Before you open a pull request

These mirror what CI enforces. `main` is protected, and the `ci-gate` check has to pass
before anything can merge, so a red local run is a red pull request.

1. Formatting and lint — both hard gates in CI. flake8 reads its settings from
   `.flake8`, so pass no command-line flags; the 100-column limit used to live only in
   the CI command and could not be reproduced locally.

   ```bash
   black --check .
   flake8 beast2py tests
   ```

2. The test suite. Keep `-rs`, which prints the reason for every skipped test.
   To match CI exactly, add coverage:

   ```bash
   pytest tests/ -v -rs --cov=beast2py
   ```

3. The workflow files, if you touched anything under `.github/workflows/`:

   ```bash
   python .github/validate_workflows.py
   ```

   This catches workflows that are valid YAML but that GitHub will not accept — most
   easily two jobs sharing a display name, which GitHub reports as a run with zero jobs
   under the file's path rather than under the workflow name.

4. Types, if you touched `beast2py/`. mypy is advisory in CI and gated by a ratchet: it
   fails only if the error count rises above the current baseline of 20, recorded as
   `MYPY_BASELINE` in `.github/workflows/ci.yml`. So running it locally and keeping the
   count flat is enough; you are not expected to fix the existing backlog:

   ```bash
   mypy
   ```

5. The BEAST2 integration tests, if you have BEAST2 2.7.x and JDK 17 or newer.
   The 21 tests in `tests/test_beast2_validation.py` generate XML for all 19
   example configurations and check it against the real BEAST2 `XMLParser`.
   Without a BEAST2 install they skip rather than fail, so such a run proves
   nothing about the third gate:

   ```bash
   pytest tests/test_beast2_validation.py -v -rs
   ```

## Running the third validation gate locally

Gate 3 parses and initialises the generated XML with the real BEAST2 parser. It
needs `BEAST.base.jar` from a BEAST2 2.7.x installation and a JDK 17 or newer;
BEAST 2.7 class files will not load on an older JDK.
`beast2py/beast2_validate.sh` finds the jar in this order:

1. `$BEAST2_JAR`, if that variable is set;
2. the newest `~/.beast/2.7/BEAST.base/*/lib/BEAST.base.jar`;
3. `tools/lib/BEAST.base.jar` or `lib/BEAST.base.jar` next to the script.

A normal BEAST2 install is therefore picked up with no configuration, and an
unusual one is handled by setting the variable explicitly, for example
`export BEAST2_JAR=/path/to/BEAST.base.jar`. The same gate is reachable from the
CLI: `--beast2-validate` on `generate` and `quick`, `--beast2` on `validate`.

The script searches a list of candidate JDKs rather than trusting `java` on
`PATH`, because a macOS machine often has a stale Oracle 1.8 shim on `PATH` next
to a perfectly usable unlinked `openjdk@17`. Set `JAVA_HOME` to force a
particular JDK, and `BEAST2_VALIDATE_TRACE=1` to see which `java` was chosen.
Exit codes are 0 (all files valid), 1 (at least one invalid) and 2 (setup error:
jar or JDK missing or too old). If the jar cannot be found, generation stops
with exit code 2 and writes nothing; `--allow-unvalidated` overrides that, but
XML produced that way must not be described as BEAST2-verified. To run the script
alone: `bash beast2py/beast2_validate.sh examples/output_basic.xml`. The
`beast2_validate.sh` at the repository root is a thin wrapper that execs it.

## Important: the validator ships as precompiled bytecode

`beast2py/tools/` contains a Java source file **and committed `.class` files**:

- `beast2py/tools/classes/Beast2Validator.class` — compiled from
  `beast2py/tools/Beast2Validator.java` in this repository
- `beast2py/tools/classes/beast/pkgmgmt/*.class` — compiled from the BEAST 2
  v2.7.8 package-management sources, redistributed unmodified
- `beast2py/tools/launcher.jar` — a thin `Main-Class` wrapper

`beast2py/beast2_validate.sh` runs the committed `Beast2Validator.class` when it
is present, and only falls back to JDK single-file source mode when it is
absent. **If you edit `Beast2Validator.java` and do not recompile and commit the
class file, CI keeps validating with the stale bytecode and your change silently
does nothing.** Those same `.class` files are declared as package data in
`pyproject.toml` and `MANIFEST.in`, so a stale class file also ships in the wheel.
Recompile with a JDK 17+ and the BEAST2 base jar on the class path:

```bash
export BEAST2_JAR=~/.beast/2.7/BEAST.base/2.7.8/lib/BEAST.base.jar
javac -cp "$BEAST2_JAR" -d beast2py/tools/classes beast2py/tools/Beast2Validator.java
```

Then commit `beast2py/tools/classes/Beast2Validator.class` in the same commit as
the source change. A pull request that changes the Java source without the class
file will be sent back. Do not recompile, reformat or replace the
`beast/pkgmgmt/*.class` files: they are verbatim upstream bytecode under
LGPL-2.1, and `beast2py/tools/NOTICE` records their provenance, their unmodified
status and how to obtain the corresponding source.

## Pull request process

- Branch from `main`; do not push to `main`.
- One logical change per pull request. Unrelated cleanups, reformatting and
  behaviour changes in the same PR make it much harder to review.
- Commit messages follow the Conventional Commits style already used by
  Dependabot here: a type prefix, then a short imperative summary, for example
  `fix(calibration): honour lower bound on lognormal offset`.
- Reference the issue in the body with `Closes #<n>`, and fill in the pull
  request template. In particular, paste the `pytest tests/ -rs` summary line
  and state whether the BEAST2 integration job actually ran or skipped.
- Add a test for any behaviour change. A change to XML generation should be
  checked against the real BEAST2 parser when you can run it, not only against a
  golden file.
- Update `CHANGELOG.md` under `## [Unreleased]` for anything a user would notice.

## Documentation is bilingual

Every `docs/*_en.md` file has a `docs/*_zh.md` counterpart, and the same is true
of `README.md` and `README_zh.md`. A change to an English document without the
matching Chinese change is an incomplete change and will be asked for. The
manuals also carry a strict-validation section listing every input the parser
rejects and the error message it produces, which `tests/test_release_integrity.py`
checks in part; if you change a rejection message, a documented command or flag,
or a count of tests or example configurations, update the documentation and those
counts in the same pull request.

When you document a limitation, state it as precisely as the code does: the
difference between "verified against BEAST2 v2.7.8" and "expected to work on
other 2.7.x releases", and between a check that ran and one that was skipped, is
deliberate throughout this project.

## Licensing

Beast2Py itself is MIT (`LICENSE`), and your contribution is offered under the
same terms. The vendored BEAST2 package-management classes under
`beast2py/tools/` are LGPL-2.1-only, and the package is distributed as
`MIT AND LGPL-2.1-only`. Because of that:

- Do not add GPL-only or otherwise incompatible code.
- Do not add a new third-party binary, jar or `.class` file to the repository or
  the wheel without opening an issue first. New vendored material needs a
  provenance note, a licence decision and updates to `beast2py/tools/NOTICE` and
  `MANIFEST.in`, and that is a conversation, not a pull request.
- Do not copy code from BEAST2 into the MIT-licensed Python modules.

## A green CI run is not proof that the add-on cases ran

The CI workflow has five jobs. `lint` runs black, flake8, mypy and the workflow
validator. The `test` matrix covers ubuntu-latest and macos-latest on Python 3.10, 3.11
and 3.12. An `integration` job installs a real BEAST2 runtime. A `packaging` job builds a
wheel and asserts that the bundled validator, the console entry point and the licence
metadata still work. Finally `ci-gate` aggregates the four of them and is the single check
that branch protection requires on `main`, so that changing the matrix never strands a pull
request waiting on a context name that no longer exists.

In the `integration` job, the two examples that need an add-on package —
`config_bd_skyline.yaml` (BDSKY) and `config_nested_sampling.yaml` (nested sampling) —
depend on a `beast -get` that the workflow deliberately treats as a warning rather than an
error. If that download fails, the two cases **skip** and the job is still green. Check the
integration log for the `-rs` skip reasons, or reproduce locally with
`pytest tests/test_beast2_validation.py -v -rs`.

Note also that the `integration` job installs a BEAST2 **2.7.7** runtime and then upgrades
only the `BEAST.base` package to 2.7.8. That upgrade is best-effort, so a green run does not
by itself prove 2.7.8 was the version in use. The job prints the version it resolved.

## Cutting a release

Pushing a `v*` tag runs `.github/workflows/release.yml`, which builds the sdist and wheel,
checks them, and attaches them to the GitHub release. Before tagging, set the version in
these four places — the `verify` job refuses to build if they disagree:

1. `beast2py/__init__.py` — `__version__`, the single source of truth
2. `CITATION.cff` — `version`, and add `date-released`
3. `.zenodo.json` — `version`
4. the git tag itself, `vX.Y.Z`

Then move the `[Unreleased]` section of `CHANGELOG.md` under the new version heading,
commit, and tag:

```bash
git tag -a v0.2.0 -m "v0.2.0"
git push origin main --tags
```

### Archiving a release on Zenodo

**No Zenodo DOI has been minted yet.** A `.zenodo.json` file in the repository does not by
itself create an archive; the repository has to be connected to Zenodo and enabled there
first, and that has not been done. To set it up:

1. Sign in to Zenodo and, under your account settings, enable the GitHub integration and
   then enable this repository.
2. Publish a GitHub release. Zenodo archives the source at that tag and mints both a
   concept DOI and a version-specific DOI.
3. Copy the version-specific DOI badge into `README.md` and `README_zh.md`, and record the
   DOI in `CITATION.cff` under `identifiers`.

Until a DOI exists, cite the release tag as shown in the README.

### Publishing to PyPI

Not enabled. The `pypi` job in `release.yml` is gated on a repository **variable** named
`PYPI_PUBLISH` being set to the string `true`, and it is skipped with a visible notice
otherwise. To enable it:

1. On pypi.org, register this repository as a **trusted publisher** for the project:
   owner `ZengZichao`, repository `Beast2Py`, workflow `release.yml`, environment `pypi`.
2. Add a repository variable `PYPI_PUBLISH` with the value `true`.

Authentication is Trusted Publishing (OIDC): the job requests a short-lived `id-token`
and PyPI validates it against the trusted-publisher entry above. **No API token is
stored in the repository**, so there is no long-lived secret to rotate or leak.

The gate is a variable rather than a token because GitHub rejects an entire workflow
file if `secrets` appears in any `if` expression, at job level or step level. Only
`github`, `needs`, `vars` and `inputs` are available there — a variable keeps the
opt-in visible in the workflow file. The publish action itself
(`pypa/gh-action-pypi-publish`) is pinned to a commit SHA, since it sits on the
release path.

## Security issues

Do not open a public issue for a vulnerability. Follow the disclosure process in
[SECURITY.md](SECURITY.md).

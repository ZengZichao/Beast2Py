# Beast2Py User Manual

- **Version:** 0.1.0
- **Authors:** Zichao Zeng (曾子超) · [ORCID: 0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)
- **Date:** 2026

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Installation](#2-installation)
3. [Two Invocation Modes: CLI and API](#3-two-invocation-modes-cli-and-api)
4. [Quick Start](#4-quick-start)
5. [CLI Commands in Detail](#5-cli-commands-in-detail)
6. [Python API in Detail](#6-python-api-in-detail)
7. [Complete YAML Configuration Reference](#7-complete-yaml-configuration-reference)
8. [Substitution Models](#8-substitution-models)
9. [Clock Models](#9-clock-models)
10. [Tree Priors](#10-tree-priors)
11. [Calibration Distributions](#11-calibration-distributions)
12. [Multi-Partition Analysis](#12-multi-partition-analysis)
13. [Tip Dates](#13-tip-dates)
14. [Calibration Diagnostics](#14-calibration-diagnostics)
15. [Reproducibility Tools](#15-reproducibility-tools)
16. [Validation](#16-validation)
17. [MCMC Settings](#17-mcmc-settings)
18. [Logger Configuration](#18-logger-configuration)
19. [Initialization Strategies](#19-initialization-strategies)
20. [Operator Weights](#20-operator-weights)
21. [Example Configuration Files](#21-example-configuration-files)
22. [Troubleshooting](#22-troubleshooting)
23. [Appendix: Complete Model List](#23-appendix-complete-model-list)

---

## 1. Introduction

Beast2Py is a Python framework for generating, validating, and diagnosing BEAST2 XML configurations for divergence time estimation. It provides comprehensive model support, automated calibration conflict detection, prior sensitivity analysis, and reproducibility tooling. Beast2Py offers both a command-line interface (CLI) and an application programming interface (API), covering the needs of different usage scenarios.

### Key Features

Beast2Py generates a complete BEAST2 XML from a YAML configuration. The supported models are:

- 13 substitution models: 7 nucleotide models and 6 amino acid models
- 4 clock models: strict, UCLN, UCE, RLC
- 8 tree priors: Yule, CalibratedYule, Birth-Death, constant/exponential Coalescent, Bayesian Skyline, EBSP, BD Skyline Serial

The tool supports multi-partition analyses (including codon partitioning) and can link or
unlink models between partitions; all partitions are assembled on a single linked topology
(see [Section 12](#12-multi-partition-analysis)).

On the diagnostics side, Beast2Py provides automated calibration conflict detection (temporal
consistency, distribution overlap, monophyly conflicts) and prior sensitivity analysis
(leave-one-out XML generation). It also supports calibration distribution visualization,
automatic generation of analysis fingerprints and methods descriptions, and Snakemake/Nextflow
pipeline file generation.

Two properties of the pipeline matter for every section below:

- **The parser is fail-closed.** Every key is checked against a schema, every boolean must be an unambiguous boolean, every alignment is checked for equal lengths, unique names and a legal alphabet, and every distribution and tree-prior parameter is range-checked. Anything the generator cannot honour is an error, not a silent default ([Section 22.1](#221-strict-validation-rejected-inputs-and-their-error-messages)).
- **Writing an XML is gated.** `generate` and `quick` run structural checks, calibration conflict detection and (optionally) the real BEAST2 parser+initialiser, and exit non-zero without writing the file if any gate that ran fails ([Section 16](#16-validation)).

### Design Philosophy

Beast2Py adopts a layered architecture that separates the CLI/API entry layer, the configuration management layer, the domain modules, the data model layer, and the utility support layer. Within the domain modules, the diagnostics and reproducibility modules are the core innovations. The model registry pattern (`ModelRegistry`) maps the model names given in the configuration to the corresponding BEAST2 XML spec values, so users can register custom models without modifying core code.

### Relation to Other Tools

| Tool | Language | When to use it |
|---|---|---|
| BEAUti2 (bundled with BEAST2) | GUI | Interactive, exploratory work; v2.7+ also provides an interactive methods-section viewer and a prior-distribution preview chart |
| babette suite | R | End-to-end R workflow: beautier generates the XML, beastier runs BEAST2, tracerer parses logs (ESS/HPD) |
| BEASTling / beast2-xml | Python | Configuration-file / command-line XML generation for specialized domains (linguistic data; FASTA/FASTQ) |
| LPhy + LPhyBeast | Java | Probabilistic model specification language → XML (already targeting BEAST 3) |
| Beast2Py | Python | Declarative YAML + calibration conflict detection + leave-one-out prior-sensitivity XMLs + provenance/fingerprint/methods/pipeline artifacts |

Note that prior-only sampling is also available outside Beast2Py: the BEAST2 command line accepts `beast -sampleFromPrior your.xml` for any existing XML. Beast2Py's contribution is the automated N+1 leave-one-out orchestration, combined with conflict detection and reproducibility artifacts in a single configuration-stage workflow.

---

## 2. Installation

### 2.1 Prerequisites

Beast2Py requires Python 3.10 or higher and the pip package manager. Generation, diagnostics and the structural validator are pure Python. The optional BEAST2 validation gate additionally needs a BEAST2 installation and a **JDK 17 or newer** runtime, because BEAST 2.7.x class files refuse to load on anything older. Everything here has been verified against **BEAST2 v2.7.8**; other 2.7.x releases are expected to work but have not been tested.

### 2.2 Install from source

```bash
# Clone the repository
git clone https://github.com/ZengZichao/Beast2Py.git
cd Beast2Py

# Install in development mode
pip install -e .
```

### 2.3 Dependencies

The following dependencies are installed automatically:

| Package | Minimum version | Purpose |
|---|---|---|
| pyyaml | >=6.0 | YAML configuration parsing |
| biopython | >=1.80 | Sequence file parsing (FASTA/NEXUS) |
| numpy | >=1.20 | Numerical computation and array operations |
| scipy | >=1.8 | Statistical distributions and Wasserstein distance |
| matplotlib | >=3.5 | Calibration distribution visualization |
| jinja2 | >=3.0 | HTML diagnostic report generation |

### 2.4 Optional dependencies

| Package | Purpose |
|---|---|
| lxml | Optional XML parsing backend with faster parsing |
| pytest | Run the test suite |
| pytest-cov | Test coverage measurement |
| black | Code formatting |
| flake8 | Code style checking |
| mypy | Type checking |

Install the development dependencies:

```bash
pip install -e ".[dev]"
```

Install the optional XML backend:

```bash
pip install -e ".[xml]"
```

### 2.5 Verify the installation

```bash
# Check the version
beast2py --version
# Output: Beast2Py v0.1.0

# List supported models
beast2py list-models

# Verify the Python API
python3 -c "import beast2py; print(beast2py.__version__)"
# Output: 0.1.0
```

Without installing the package, every command is reachable from a clone of the
repository root as `python3 -m beast2py.main <subcommand>`; the two forms are the
same entry point (`[project.scripts] beast2py = "beast2py.main:main"`).

### 2.6 Installing BEAST2 (optional, for the BEAST2 validation gate)

BEAST2 can be downloaded from https://www.beast2.org/. Passing `--beast2-validate` (for
`generate` and `quick`) or `--beast2` (for `validate`) enables the third gate; `--beast2-path`
then accepts either a BEAST2 `beast` executable or a `beast2_validate.sh` script. Specify
neither, and the tool locates the bundled `beast2_validate.sh` on its own, trying next to the
installed package, the project root, and `PATH`. The script needs `BEAST.base.jar`, looked up
in this order:

1. the path in `$BEAST2_JAR`, if that variable is set;
2. the newest `~/.beast/2.7/BEAST.base/*/lib/BEAST.base.jar`;
3. `tools/lib/BEAST.base.jar` or `lib/BEAST.base.jar` next to the script.

The script also appends the add-on packages installed under `~/.beast/2.7` (BDSKY, NS and the
like) to the class path, so the two add-on examples validate too.

The script requires JDK 17 or newer, and searches a batch of **candidate paths** for it
instead of trusting `PATH`. The reason is that a macOS box very often ships a stale Oracle
`java` shim on `PATH` while a usable `openjdk@17` sits unlinked in `/opt/homebrew`. The first
candidate whose major version is at least 17 wins, in this order:

1. entries of `BEAST2_JAVA_CANDIDATES` (a colon-separated list), if you set it;
2. `$JAVA_HOME/bin/java`;
3. whatever `/usr/libexec/java_home -v 17+` reports (skippable with `BEAST2_SKIP_JAVA_HOME_TOOL=1`);
4. the Homebrew, `/usr/lib/jvm` and `/Library/Java/JavaVirtualMachines` directories (skippable
   with `BEAST2_SKIP_JDK_PROBES=1`);
5. `java` on `PATH` — last, because it is the least specific of all the candidates.

If none of these reaches 17, the script gives up, exiting **2** and listing the candidates it
rejected:

```bash
bash beast2_validate.sh examples/output_basic.xml
# ... ends with: VALID: examples/output_basic.xml

BEAST2_VALIDATE_TRACE=1 bash beast2_validate.sh examples/output_basic.xml
# ... prints, on stderr: [beast2_validate] using java 17 at /opt/homebrew/opt/openjdk@17/...
```

`JAVA_HOME` remains the way to *force* a particular JDK, and an explicitly JDK-8-only machine
still gets the clear error:

```bash
JAVA_HOME=/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home \
    bash beast2_validate.sh examples/output_basic.xml
```

Exit codes are `0` (all files valid), `1` (at least one invalid) and `2` (setup error: jar or
Java missing/too old). Several files may be validated in one call.

---

## 3. Two Invocation Modes: CLI and API

Beast2Py provides two equivalent invocation modes; choose whichever fits your workflow.

### 3.1 CLI Mode (Command-Line Interface)

CLI mode suits batch analyses, scripted workflows, and pipeline integration in terminal environments. It is invoked through the `beast2py` command and provides 8 subcommands: `generate`, `diagnose`, `quick`, `validate`, `methods`, `pipeline`, `fingerprint`, and `list-models`.

```bash
beast2py generate --config config.yaml --output output.xml
```

### 3.2 API Mode (Python Programming Interface)

API mode suits interactive analysis in Jupyter Notebooks, programmatic invocation from Python scripts, and integration with other bioinformatics tools. It is used by importing the `beast2py` package.

```python
from beast2py.api import Beast2Py

# Create an API instance
b2p = Beast2Py()

# Generate XML
xml_str = b2p.generate_xml("config.yaml", output="output.xml")

# Run diagnostics
report = b2p.diagnose("config.yaml", report_path="report.html")

# Validate XML
result = b2p.validate_xml("output.xml")
```

Module-level convenience functions are also supported:

```python
from beast2py.api import generate_xml, diagnose, validate_xml

xml = generate_xml("config.yaml", output="output.xml")
report = diagnose("config.yaml", report_path="report.html")
result = validate_xml("output.xml")
```

### 3.3 Correspondence between the two modes

| CLI command | API method | Module-level function |
|---|---|---|
| `beast2py generate` | `Beast2Py.generate_xml()` | `generate_xml()` |
| `beast2py quick` | `Beast2Py.quick_generate()` | `quick_generate()` |
| `beast2py diagnose` | `Beast2Py.diagnose()` | `diagnose()` |
| `beast2py validate` | `Beast2Py.validate_xml()` | `validate_xml()` |
| `beast2py methods` | `Beast2Py.generate_methods()` | `generate_methods()` |
| `beast2py pipeline` | `Beast2Py.generate_pipeline()` | `generate_pipeline()` |
| `beast2py fingerprint` | `Beast2Py.generate_fingerprint()` | `generate_fingerprint()` |
| `beast2py list-models` | `Beast2Py.list_models()` | `list_models()` |

---

## 4. Quick Start

### 4.1 Basic workflow

A typical Beast2Py workflow consists of the following steps.

1. Prepare the sequence data (FASTA or NEXUS format).
2. Create a YAML configuration file specifying the models, calibrations and MCMC settings.
3. Generate the BEAST2 XML with Beast2Py.
4. Optional: run calibration diagnostics to detect potential problems.
5. Run BEAST2 with the generated XML and analyse the results.

### 4.2 CLI quick start

```bash
# Step 1: Create a YAML configuration file (see Section 7 for the full reference)
# Step 2: Generate the XML (three gates run first; nothing is written if one fails)
beast2py generate --config config.yaml --output output.xml -v

# Step 3: Run diagnostics (optional)
beast2py diagnose --config config.yaml --report report.html --sensitivity --visualize

# Step 4: Validate the XML (optionally with BEAST2; needs JDK 17)
beast2py validate --xml output.xml --beast2

# Step 5: Run BEAST2
beast -overwrite output.xml
```

`generate` refuses to touch an existing `output.xml` and exits 1 unless you pass `--force`.
The reason is that overwriting the output while its sibling `.fingerprint.json` /
`.methods.tex` still survive would leave two analyses sharing one name.

Working example, runnable against a clone of this repository (verified output):

```text
$ beast2py generate -c examples/config_basic.yaml -o /tmp/output.xml -v
[INFO] Reading configuration: examples/config_basic.yaml
  [verbose] Analysis: primates_basic_calibration
  [verbose] Partitions: 1
  [verbose] Calibrations: 1
  [verbose] Tree prior: yule
[INFO] Generating BEAST2 XML...
[OK] Gate 1/3 passed: Beast2Py structural checks (ids, idrefs, priors, operators).
[OK] Gate 2/3 passed: no contradictory calibration conflicts.
[OK] XML written to: /tmp/output.xml
[OK] Done!
```

With `--beast2-validate` a third line appears:
`[OK] Gate 3/3 passed: BEAST2 parsed and initialised the model.`. Without BEAST2 the line
becomes `[ERROR] Gate 3/3 unavailable: no BEAST2 installation was found, …`, the exit code is
**2**, and no file is written. `--allow-unvalidated` downgrades that to a `[WARN]` and writes
the XML anyway, in which case it must not be described as BEAST2-verified.

### 4.3 API quick start

```python
from beast2py.api import Beast2Py

b2p = Beast2Py()

# One-step generation, validation, fingerprint, and methods description
result = b2p.generate(
    "config.yaml",
    output="output.xml",
    diagnose=True,
    fingerprint=True,
    methods=True,
)

print(f"XML: {result['xml_path']}")
print(f"Validation passed: {result['validation'].is_valid}")
print(f"Fingerprint: {result['fingerprint']}")
print(f"Methods: {result['methods_path']}")

# Inspect the diagnostic results
if 'diagnostic_report' in result:
    report = result['diagnostic_report']
    for c in report.conflicts:
        print(f"[{c.severity}] {c.description}")
```

### 4.4 Minimal configuration example

```yaml
metadata:
  analysis_name: "example"

alignments:
  - id: "alignment"
    file: "primates.fasta"
    data_type: "nucleotide"

partitions:
  - id: "alignment"
    site_model:
      substitution_model: {type: "hky"}
    clock_model: {type: "strict"}
    tree: "shared"

tree_prior:
  type: "yule"
  birth_rate: {value: 1.0}

calibrations:
  - name: "root"
    taxa: null
    monophyletic: true
    distribution:
      type: "normal"
      parameters: {mean: 10.0, sigma: 1.0}

mcmc:
  chain_length: 10000000

output:
  file: "output.xml"
```

---

## 5. CLI Commands in Detail

### 5.1 `generate` — Generate BEAST2 XML

Generates a complete BEAST2 XML file from a YAML configuration, optionally running diagnostics and generating a fingerprint and methods description at the same time. `generate` writes the file only after every gate that ran has passed.

```bash
beast2py generate \
    --config config.yaml \
    --output output.xml \
    [--diagnose] \
    [--fingerprint] \
    [--methods] \
    [--beast2-validate] \
    [--beast2-path PATH] \
    [--force] \
    [--verbose]
```

**Parameters:**

| Parameter | Short | Required | Default | Description |
|---|---|---|---|---|
| `--config` | `-c` | Yes | — | YAML configuration file path |
| `--output` | `-o` | Yes | — | Output XML file path |
| `--diagnose` | | No | False | Also run calibration diagnostics |
| `--fingerprint` | | No | False | Generate an analysis fingerprint file (`.fingerprint.json`) |
| `--methods` | | No | False | Generate a LaTeX methods description (`.methods.tex`) |
| `--beast2-validate` | | No | False | Run the third gate: the BEAST2 parser, **calling `initAndValidate()` on every object** (needs BEAST2 and JDK 17) |
| `--beast2-path` | | No | `beast` | Path to a BEAST2 validation script (`.sh`) or the BEAST2 executable; by default the tool locates the bundled `beast2_validate.sh` itself, trying next to the installed package, the project root, and `PATH` |
| `--force` | | No | False | Overwrite `--output` if it already exists (default: refuse and exit 1) |
| `--verbose` | `-v` | No | False | Show verbose output |

**Gates and exit codes.** Gate 1 (Beast2Py structural checks), Gate 2 (calibration conflict
detection) and, when requested, Gate 3 (BEAST2 parse and initialisation) all run *before*
the file is written. A failing gate prints `[ERROR] Gate n/3 failed: …`, returns exit code
1 and writes no output file; the Python API raises `ValueError` instead. If an earlier run
left a `.fingerprint.json` or `.methods.tex` next to the same output name and this call
omits `--fingerprint`/`--methods`, the CLI warns that those files are going stale. See
[Section 16](#16-validation) for what the three gates do and do not guarantee.

**Examples:**

```bash
# Basic generation
beast2py generate -c examples/config_basic.yaml -o /tmp/output.xml

# Generate and also run diagnostics, fingerprint, and methods
beast2py generate -c examples/config_basic.yaml -o /tmp/output.xml \
    --diagnose --fingerprint --methods -v

# Generate and validate with BEAST2 (third gate; needs JDK 17)
JAVA_HOME=/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home \
    beast2py generate -c examples/config_basic.yaml -o /tmp/output.xml --beast2-validate

# Replace an existing output
beast2py generate -c examples/config_basic.yaml -o /tmp/output.xml --force
```

### 5.2 `diagnose` — Run calibration diagnostics

Runs calibration conflict detection, prior sensitivity analysis, and distribution visualization, producing a self-contained HTML report.

```bash
beast2py diagnose \
    --config config.yaml \
    --report report.html \
    [--sensitivity] \
    [--visualize] \
    [--sensitivity-xml-dir DIR] \
    [--verbose]
```

**Parameters:**

| Parameter | Short | Required | Default | Description |
|---|---|---|---|---|
| `--config` | `-c` | Yes | — | YAML configuration file path |
| `--report` | | Yes | — | Output HTML report path |
| `--sensitivity` | | No | False | Generate prior sensitivity analysis XMLs (N+1 files) |
| `--visualize` | | No | False | Generate calibration distribution density plots |
| `--sensitivity-xml-dir` | | No | `<report dir>/diagnostics` | Writes the sensitivity XMLs into the `sensitivity/` subfolder of this directory |
| `--verbose` | `-v` | No | False | Show verbose output |

**Examples:**

```bash
# Full diagnostics
beast2py diagnose -c examples/config_basic.yaml --report /tmp/report.html --sensitivity --visualize

# Conflict detection only
beast2py diagnose -c examples/config_basic.yaml --report /tmp/report.html
```

Note that `generate` already runs the conflict detector as its Gate 2; `diagnose` is
the standalone, richer report (sensitivity XMLs, density plots, HTML).

### 5.3 `quick` — Quick mode

The user needs no full YAML configuration: `quick` builds a single-partition XML straight from
command-line arguments. It suits simple analyses and quick tests.

```bash
beast2py quick \
    --alignment input.fasta \
    --output output.xml \
    [--tree-prior yule] \
    [--subst-model hky] \
    [--clock-model strict] \
    [--gamma-categories 0] \
    [--calibration-yaml cals.yaml] \
    [--chain-length 10000000] \
    [--pre-burnin 0] \
    [--name "analysis_name"] \
    [--force] \
    [--verbose]
```

**Parameters:**

| Parameter | Short | Required | Default | Description |
|---|---|---|---|---|
| `--alignment` | `-a` | Yes | — | Sequence file path (FASTA/NEXUS), read as nucleotide data |
| `--output` | `-o` | Yes | — | Output XML file path |
| `--tree-prior` | | No | `yule` | Tree prior type |
| `--subst-model` | | No | `hky` | Substitution model type |
| `--clock-model` | | No | `strict` | Clock model type |
| `--gamma-categories` | | No | `0` | Number of Gamma rate categories, 0 meaning none; above 0 it also creates a `gamma_shape` parameter starting at 0.5 |
| `--calibration-yaml` | | No | — | Path to a calibration-point YAML file |
| `--chain-length` | | No | `10000000` | MCMC chain length |
| `--pre-burnin` | | No | `0` | MCMC pre-burnin length |
| `--name` | | No | `quick_analysis` | Analysis name |
| `--force` | | No | False | Overwrite `--output` if it already exists |
| `--verbose` | `-v` | No | False | Show verbose output |

`quick` runs the same Gate 1 and Gate 2 as `generate` and honours `--force` in the same way,
but it has **no `--beast2-validate` option**. Validate the result with
`beast2py validate --xml output.xml --beast2` instead.

`quick` reads the alignment as nucleotide data, so `--subst-model` may only name a nucleotide
model; the data-type cross-check refuses an amino-acid model here just as it does in YAML mode.
The same validator parses the calibration YAML, so a quoted `monophyletic: "false"`, reversed
`uniform` bounds or an unknown distribution type are errors here too.

**Calibration YAML file format (for `--calibration-yaml`):**

```yaml
- name: "humanChimpMRCA"
  taxa: ["human", "chimp"]
  monophyletic: true
  distribution:
    type: "normal"
    parameters: {mean: 6.0, sigma: 0.5}
  provenance:
    source: "fossil"
    reference: "DOI:10.1038/nature12323"
    calibration_type: "soft"

- name: "rootCalibration"
  taxa: null
  monophyletic: true
  distribution:
    type: "uniform"
    parameters: {lower: 20.0, upper: 40.0}
  use_originate: true
```

**Examples:**

```bash
# Simplest quick generation
beast2py quick -a examples/primates.fasta -o /tmp/output_quick.xml

# With calibrations and Gamma
beast2py quick -a examples/primates.fasta -o /tmp/output_quick.xml \
    --tree-prior birth_death \
    --subst-model gtr \
    --clock-model ucln \
    --gamma-categories 4 \
    --calibration-yaml examples/calibrations.yaml \
    --chain-length 1000000 \
    --name "my_analysis"
```

The calibration YAML is parsed by the same validator as a full configuration, so a
quoted `monophyletic: "false"`, an inverted `uniform` bound or an unknown
distribution type is rejected here too.

### 5.4 `validate` — Validate XML

Validates a BEAST2 XML file: first Beast2Py's own structural and semantic checks, then, when
requested, the BEAST2 parse and model initialisation checks.

```bash
beast2py validate \
    --xml output.xml \
    [--beast2] \
    [--beast2-path PATH]
```

**Parameters:**

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--xml` | Yes | — | Path to the XML file to validate |
| `--beast2` | No | False | Also run the BEAST2 parse and initialisation check (needs BEAST2 and JDK 17) |
| `--beast2-path` | No | `beast` | Path to the BEAST2 executable or to `beast2_validate.sh` |

**Level 1 (always).** This level checks XML well-formedness together with the following.

- The root element is `<beast>`
- The required elements `<data>`, `<run>` and `<state>` are present
- id/idref/`@ref` references are consistent and no ID is duplicated
- The posterior reference and the loggers are present
- The sequence data is non-empty and of equal length, and the taxon names are unique

The semantic checks are: every estimated parameter carries a `Prior`, no state node has
inverted bounds, and every initial value lies inside its own bounds and matches the declared
dimension. When `storeEvery` disables checkpointing, this level warns instead of failing.

**Level 2 (`--beast2`).** The headless BEAST2 `XMLParser` parses the file *and* calls
`initAndValidate()` on every object. The two levels report with distinct wording
(`Level 1/2 passed: structural checks …` vs
`Level 2/2 passed: BEAST2 parsed and initialised the model.`). If either level fails, the
exit code is 1.

```bash
beast2py validate --xml examples/output_basic.xml
JAVA_HOME=/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home \
    beast2py validate --xml examples/output_basic.xml --beast2
```

### 5.5 `methods` — Generate a methods description

Automatically generates a LaTeX methods paragraph suitable for inclusion in a manuscript.

```bash
beast2py methods \
    --config config.yaml \
    [--output methods.tex]
```

| Parameter | Short | Required | Default | Description |
|---|---|---|---|---|
| `--config` | `-c` | Yes | — | YAML configuration file path |
| `--output` | `-o` | No | Standard output | Output file path |

**Generated content:** substitution model description (with citations), clock model description, tree prior description, calibration point descriptions (with distribution parameters and provenance), MCMC settings description, and the analysis fingerprint identifier.

### 5.6 `pipeline` — Generate pipeline files

Generates Snakemake or Nextflow pipeline files for end-to-end reproducible analysis workflows.

```bash
beast2py pipeline \
    --config config.yaml \
    --type snakemake|nextflow|both \
    [--output-dir .]
```

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--config` | Yes | — | YAML configuration file path |
| `--type` | No | `snakemake` | Pipeline type: `snakemake`, `nextflow`, or `both` |
| `--output-dir` | No | `.` | Output directory |

The pipeline contains the following rules or processes: XML generation (invoking `beast2py generate`), BEAST2 execution (invoking `beast -overwrite`), calibration diagnostics (invoking `beast2py diagnose`), and methods description generation (invoking `beast2py methods`).

### 5.7 `fingerprint` — Generate an analysis fingerprint

Generates the deterministic analysis identifier: a SHA-256 digest of the scientific
configuration plus a content digest of every alignment.

```bash
beast2py fingerprint \
    --config config.yaml \
    [--output fingerprint.json] \
    [--verbose]
```

| Parameter | Short | Required | Default | Description |
|---|---|---|---|---|
| `--config` | `-c` | Yes | — | YAML configuration file path |
| `--output` | `-o` | No | Terminal output only | Output JSON file path |
| `--verbose` | `-v` | No | False | Show detailed fingerprint information |

**Fingerprint format:** `B2P-{first 12 hex digits of the config digest}-{tool version}`,
e.g. `B2P-ba18c26ed61f-0.1.0`. The identifier carries **no date and no time zone**, so
rerunning the same configuration on the same data always reproduces it. Twelve hex
digits are 48 bits of digest: by the birthday paradox, a 50% collision probability would
need about 2²⁴ (≈ 1.7 × 10⁷) distinct analyses.

**Contents of the fingerprint JSON file** (verified output for `examples/config_basic.yaml`):

```json
{
  "fingerprint": "B2P-ba18c26ed61f-0.1.0",
  "config_hash": "ba18c26ed61fcf4b157723d410f13f7534b0533c9c4810f6a3f422caa3ed9d8e",
  "full_hash": "ba18c26ed61fcf4b157723d410f13f7534b0533c9c4810f6a3f422caa3ed9d8e",
  "data_hash": "30be5613d3ff6862",
  "hash_bits": 48,
  "collision_note": "The identifier truncates the configuration digest to 48 bits; a 50% collision probability would require about 2**24 distinct analyses.",
  "tool_version": "0.1.0",
  "beast2_version": "2.7.8",
  "generation_time": "2024-01-15T10:30:00.000000",
  "analysis_name": "primates_basic_calibration"
}
```

`data_hash` condenses the digest of each alignment's content down to the first 16 hex
characters, and travels in the XML as `| Data: {hash16}` inside the fingerprint comment.
`generation_time` appears **only** in this sidecar file; the XML itself carries a
timestamp-free comment (see [Section 15.1](#151-analysis-fingerprint)).

### 5.8 `list-models` — List supported models

Displays all model types supported by Beast2Py.

```bash
beast2py list-models
```

This command takes no arguments; the output includes the complete lists of substitution models, clock models, tree priors, calibration distributions, and data types.

### 5.9 Global options

```bash
beast2py --version       # Show the version number
beast2py --help          # Show help
beast2py <command> --help  # Show help for a specific command
```

---

## 6. Python API in Detail

### 6.1 The `Beast2Py` class

The `Beast2Py` class is the central entry point of the API and encapsulates all functionality.

```python
from beast2py.api import Beast2Py

b2p = Beast2Py(base_dir=".")  # base_dir resolves relative file paths
```

**Constructor parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `base_dir` | `str` | Current working directory | Base directory for resolving relative file paths in the configuration |

### 6.2 `generate_xml` — Generate XML

```python
xml_str = b2p.generate_xml(
    config,           # Configuration file path, dict, or BEASTConfig object
    output=None,      # Optional: output file path
    force=False,      # Optional: overwrite an existing output file
)
```

The `config` parameter accepts three types: a YAML file path (string), a configuration dictionary (`dict`), or a parsed `BEASTConfig` object. Returns the generated XML string.

When `output` is given, `generate_xml()` puts the XML through **the same gates as the CLI** (structural checks, then calibration conflict detection) and writes it only after all of them pass. A failing gate raises `ValueError` and leaves the file untouched; overwriting an existing file needs `force=True`.

This behaviour is deliberate: an earlier version wrote the file and checked afterwards, so the library and the command line disagreed about what "validated" means. Call `XMLWriter(config).generate_xml()` directly if you only want the raw string with no gating.

### 6.3 `quick_generate` — Quick generation

```python
xml_str = b2p.quick_generate(
    alignment="examples/primates.fasta",
    output="output.xml",
    tree_prior="yule",
    subst_model="hky",
    clock_model="strict",
    gamma_categories=0,
    calibration_yaml="examples/calibrations.yaml",  # optional
    chain_length=10_000_000,
    pre_burnin=0,
    name="my_analysis",
    force=False,
)
```

Same gates, same `force` semantics as `generate_xml`.

### 6.4 `diagnose` — Run diagnostics

```python
report = b2p.diagnose(
    config,
    output_dir="diagnostics",
    run_sensitivity=True,
    run_visualization=True,
    report_path="report.html",  # optional
)

# report is a DiagnosticReport object
for conflict in report.conflicts:
    print(f"[{conflict.severity}] {conflict.conflict_type}: {conflict.description}")

for summary in report.summaries:
    print(f"{summary.name}: mean={summary.stats['mean']:.2f}")

print(f"Sensitivity XMLs: {report.sensitivity_xmls}")
print(f"Visualization file: {report.visualization_file}")
```

### 6.5 `detect_conflicts` — Detect conflicts only

```python
conflicts = b2p.detect_conflicts(config)
for c in conflicts:
    print(f"[{c.severity}] {c.calibration_a} <-> {c.calibration_b}: {c.description}")
```

### 6.6 `generate_sensitivity_xmls` — Generate sensitivity XMLs

```python
xml_paths = b2p.generate_sensitivity_xmls(config, output_dir="sensitivity")
print(f"Generated {len(xml_paths)} XML files")
```

### 6.7 `plot_calibrations` — Plot calibration distributions

```python
b2p.plot_calibrations(config, output_file="calibrations.png")
```

PNG, PDF, and SVG formats are supported.

### 6.8 `validate_xml` — Validate XML

```python
# Validate a file
result = b2p.validate_xml("output.xml")

# Validate an XML string
result = b2p.validate_xml(xml_string)

# Also use BEAST2 native validation
result = b2p.validate_xml("output.xml", beast2_validate=True, beast2_path="/opt/beast/bin/beast")

print(f"Valid: {result.is_valid}")
print(f"Errors: {result.errors}")
print(f"Warnings: {result.warnings}")
```

### 6.9 `generate_fingerprint` — Generate a fingerprint

```python
fingerprint = b2p.generate_fingerprint(config, output="fingerprint.json")
print(f"Fingerprint: {fingerprint}")
```

### 6.10 `generate_methods` — Generate a methods description

```python
methods_text = b2p.generate_methods(config, output="methods.tex")
print(methods_text[:200])
```

### 6.11 `generate_pipeline` — Generate pipeline files

```python
# Generate Snakemake
path = b2p.generate_pipeline(config, pipeline_type="snakemake", output_dir="pipeline/")

# Generate Nextflow
path = b2p.generate_pipeline(config, pipeline_type="nextflow", output_dir="pipeline/")

# Generate both
sm_path, nf_path = b2p.generate_pipeline(config, pipeline_type="both", output_dir="pipeline/")
```

### 6.12 `generate` — All-in-one generation

The `generate` method is an enhanced version of `generate_xml` that can perform validation, diagnostics, fingerprinting, and methods generation in one call.

```python
result = b2p.generate(
    config="config.yaml",
    output="output.xml",
    diagnose=True,
    fingerprint=True,
    methods=True,
    beast2_validate=False,
    beast2_path="beast",
    force=False,
)

# The returned dict contains:
# 'xml_path': path to the XML file
# 'validation': ValidationResult object
# 'fingerprint': fingerprint string (if fingerprint=True)
# 'fingerprint_path': path to the fingerprint file (if fingerprint=True)
# 'methods_path': path to the methods file (if methods=True)
# 'diagnostic_report': DiagnosticReport object (if diagnose=True)
# 'beast2_validation': {'success': bool, 'output': str} (if beast2_validate=True)
```

`generate()` delegates the write to `generate_xml()`, so it inherits the gates and the
`force` behaviour: a failed gate raises `ValueError` before anything is written. With
`beast2_validate=True` and no BEAST2 installation present, the third gate reports itself
as *skipped* rather than passed; `XMLValidator.beast2_validation_status()` distinguishes
`"passed"` / `"failed"` / `"unavailable"` so a missing checker is reported as unavailable and
never as a successful validation.

### 6.13 `list_models` — List models

```python
models = b2p.list_models()
print(models["substitution_models"])  # ['blosum62', 'cprev', 'dayhoff', 'gtr', 'hky', ...]
print(models["clock_models"])         # ['rlc', 'strict', 'uce', 'ucln']
print(models["tree_priors"])          # ['bayesian_skyline', 'bd_skyline_serial', ...]
```

### 6.14 `parse_config` — Parse a configuration

```python
config = b2p.parse_config("config.yaml")
print(f"Analysis name: {config.metadata['analysis_name']}")
print(f"Number of partitions: {len(config.partitions)}")
print(f"Number of calibrations: {len(config.calibrations)}")
print(f"Tree prior: {config.tree_prior_type.value}")
```

### 6.15 `read_sequence` — Read sequences

```python
alignment = b2p.read_sequence("primates.fasta", data_type="nucleotide")
print(f"Number of taxa: {alignment.n_taxa}")
print(f"Number of sites: {alignment.n_sites}")
print(f"Taxon names: {alignment.taxa_names}")
```

### 6.16 Module-level convenience functions

All `Beast2Py` class methods have corresponding module-level functions that can be called without creating an instance:

```python
from beast2py.api import (
    generate_xml,
    quick_generate,
    diagnose,
    validate_xml,
    generate_fingerprint,
    generate_methods,
    generate_pipeline,
    list_models,
    parse_config,
    read_sequence,
)

# Call directly
xml = generate_xml("config.yaml", output="output.xml")
result = validate_xml("output.xml")
models = list_models()
```

### 6.17 Direct use of configuration types

The API also exposes all configuration dataclasses, which can be used to build configurations programmatically:

```python
from beast2py import (
    BEASTConfig,
    Partition,
    SiteModelConfig,
    ClockModelConfig,
    ClockModelType,
    TreePriorType,
    MCMCConfig,
    CalibrationPoint,
    DistributionConfig,
    RealParameter,
)

# Build a configuration programmatically
site_model = SiteModelConfig(
    substitution_model={"type": "hky"},
    gamma_categories=4,
    gamma_shape=RealParameter(value=0.5, lower=0.0),
)

clock_model = ClockModelConfig(type=ClockModelType.UCLN)

# ... build Partition, Alignment, and other objects
```

---

## 7. Complete YAML Configuration Reference

### 7.1 Complete configuration template

Every key below is part of the schema; a key not listed for a section is rejected
([Section 22.1](#221-strict-validation-rejected-inputs-and-their-error-messages)).

```yaml
# === Metadata ===
metadata:
  analysis_name: "my_analysis"          # Analysis name (default: "unnamed_analysis")
  analysis_description: "Description"   # Analysis description (optional)
  author: "Your Name"                   # Author (optional)
  date: "2026-08-01"                    # Free-text date (optional)
  beast2_version: "2.7.8"              # BEAST2 version (default 2.7.8)
  tool_version: "0.1.0"                # Beast2Py version (default 0.1.0)
  embed_fingerprint: true              # Accepted here and in `output`; `output` wins

# === Sequence data ===
alignments:
  - id: "gene1"                         # Alignment ID (required)
    file: "gene1.fasta"                 # Sequence file path, relative to the YAML file
    format: "auto"                      # Format: fasta / nexus / auto (default auto)
    data_type: "nucleotide"             # Data type: nucleotide / aminoacid (default nucleotide)

  # FilteredAlignment (codon partitioning). Either declare the source separately and
  # point at it with data_id (shown here, as the examples do), or give `file` and
  # `filter` together with `data_id` naming the source. A filter may never point at
  # its own alignment id, and the source must be declared before the filter.
  - id: "codon3"
    data_id: "gene1"                    # Source alignment id (required with filter)
    filter: "3::3"                      # BEAST2 filter grammar, 1-based sites
    data_type: "nucleotide"

# === Partition configuration ===
partitions:
  - id: "gene1"                         # Partition ID (required, must match an alignment ID)
    site_model:
      substitution_model:
        type: "hky"                     # Model type (see Section 8)
        # Model-specific parameters (see Section 8)
      gamma_categories: 4               # Number of Gamma rate categories (0 = none, default 0)
      gamma_shape:                      # Gamma shape; only allowed when gamma_categories > 1
        value: 0.5
        lower: 0.0
        estimate: true                  # Estimate or not (default true)
      proportion_invariant:             # Proportion of invariant sites (optional)
        value: 0.01
        lower: 0.0
        upper: 1.0
        estimate: true
      linked_to: null                   # Share another partition's site model (optional)
    clock_model:
      type: "strict"                    # Clock type (see Section 9)
      # Clock-specific parameters (see Section 9)
      linked_to: null                   # Link to another partition's clock (optional)
    tree: "shared"                      # Only "shared" is accepted; anything else is refused

# === Tip dates (optional) ===
tip_dates:
  enabled: false                        # Enable or not (default false)
  trait_name: "date-forward"            # date-forward / date-backward (nothing else)
  units: "year"                         # year / month / day / coordinate
  dates:                                # Must cover EVERY taxon of EVERY partition
    taxonA: 2000
    taxonB: 2005
  # file: "tip_dates.csv"              # Alternative/extra source: two-column CSV

# === Tree prior ===
tree_prior:
  type: "yule"                          # Tree prior type (see Section 10)
  # Tree-prior-specific parameters (see Section 10)

# === Calibration points ===
calibrations:
  - name: "humanChimpMRCA"              # Calibration name (required)
    taxa: ["human", "chimp"]            # Taxon list (null = root node)
    monophyletic: true                  # Enforce monophyly (default true)
    distribution:
      type: "normal"                    # Distribution type (see Section 11)
      parameters: {mean: 6.0, sigma: 0.5}  # Distribution parameters (range-checked)
      offset: 0.0                       # Offset (default 0.0)
      hyperpriors:                      # Hyperpriors (optional); keys must be
        mean:                           # parameters of the distribution above
          type: "uniform"
          parameters: {lower: 0.0, upper: 20.0}
    use_originate: false                # Use the stem node (default false)
    tipsonly: false                     # Tip-only calibration (default false)
    provenance:                         # Provenance metadata (optional)
      source: "fossil"                  # Source type: fossil / biogeographic / secondary / other
      reference: "DOI:10.1038/xxxxx"    # Reference DOI or URL
      calibration_type: "soft"          # soft / hard
      original_age_min: 5.5            # Original minimum age
      original_age_max: 7.0            # Original maximum age
      notes: "Fossil description"

# === Calibration method ===
calibration_method: "mrca_prior"        # mrca_prior / calibrated_yule (default mrca_prior)

# === Parameter priors (optional; `priors:` is accepted as a synonym) ===
parameter_priors:
  - parameter: "birthRate"              # Parameter name or alias (see Section 10.5)
    distribution:
      type: "uniform"
      parameters: {lower: 0.0, upper: 1000.0}

# === MCMC settings ===
mcmc:
  type: "standard"                      # standard / nested_sampling (default standard)
  chain_length: 10000000                # Chain length, must be > 0 (default 10000000)
  pre_burnin: 0                         # Pre-burnin, must be < chain_length (default 0)
  store_every: 10000                    # Checkpoint interval in logged samples;
                                        # default chain_length // 10000; -1 disables
  sample_from_prior: false              # Sample from the prior only (default false)
  seed: 1                               # Accepted and range-checked but NOT written
                                        # into the XML: BEAST2's MCMC has no `seed`
                                        # input, the seed is `beast -seed N`
  # Nested Sampling parameters (when type=nested_sampling):
  particle_count: 1                     # Number of particles (default 1)
  sub_chain_length: 10000               # Sub-chain length (default 10000)

# === Loggers ===
loggers:
  trace_log:
    file_name: "output.$(seed).log"     # Log file name (supports the $(seed) variable)
    log_every: 1000                     # Log interval
  tree_log:
    file_name: "output.$(seed).trees"
    log_every: 1000
  screen_log:
    log_every: 10000                    # Screen log interval

# === Initialization ===
initialization:
  tree_type: "random"                   # random / upgma / newick (default random)
  newick_file: null                     # Tree file for the newick type

# === Operator weights (optional overrides) ===
operator_weights:
  tree_scaler: 3.0
  uniform: 30.0
  # ... other operators (free-form mapping; keys are not schema-checked)

# === Calibration diagnostics settings (optional) ===
diagnostics:
  overlap_measure: relative             # relative / absolute / both (default relative)
  overlap_threshold: 2.0                # absolute Wasserstein distance, time units
  relative_overlap_threshold: 0.15      # dimensionless ratio (default 0.15)

# === Output settings ===
output:
  file: "output.xml"                    # Output XML file path (default output.xml)
  embed_fingerprint: true               # Embed the fingerprint comment in the XML (default true)
```

### 7.2 Configuration parameter notes

**RealParameter format:** many parameters accept the `RealParameter` format, which can be a plain scalar or a full dictionary:

```yaml
# Plain scalar
clock_rate: 1.0

# Full dictionary
clock_rate:
  value: 1.0        # Parameter value
  lower: 0.0        # Lower bound (optional)
  upper: 10.0       # Upper bound (optional)
  dimension: 1      # Dimension (default 1)
  estimate: true    # Estimate or not (default true)
```

**Bounds are overrides, not requirements.** Omit `lower`/`upper` and the generator supplies the
bounds implied by the parameter's statistical role (`model.py::PARAM_SUPPORT`): frequencies and
probabilities `[0, 1]`, rates and population sizes `> 0`, counts and group sizes `>= 1`, and
signed quantities such as `growth_rate` unbounded. A `lower`/`upper` you do write replaces the
default. The parser rejects a `value` outside the resulting interval, and a `lower >= upper`.

**`estimate: false` is honoured for every parameter class** — substitution-model rates,
base frequencies, gamma shape, proportion invariant, clock rate, `ucld.mean` and
`ucld.stdev`. The generator writes a fixed parameter inline (`<parameter estimate="false">`),
so it never becomes a state node, never receives an operator and never appears in the trace log.

If you also wrote a `parameter_priors` entry for that parameter, the prior now targets
something that does not exist, and the parser refuses the run. That is the point: delete the
prior.

```yaml
substitution_model:
  type: "hky"
  kappa: {value: 2.0, estimate: false}   # fixed: no kappaScaler, no trace column
```

**Booleans** must be real YAML booleans or one of `true/false/yes/no/on/off/1/0`. The parser
accepts quoted forms and reads them literally. `bool("false")` is `True` in Python, so a quoted
`"false"` used to flip the setting. It now parses as false, and the parser rejects anything
ambiguous such as `"nope"`.

---

## 8. Substitution Models

### 8.1 Nucleotide substitution models

The "Parameters" column lists the **only** keys each model accepts inside
`substitution_model:` (besides `type`; models that use frequencies also take
`frequencies`). The parser rejects any other key — for example `rates` on a model
that does not accept it.

| Model | Config name | Parameters | Description |
|---|---|---|---|
| Jukes-Cantor | `jc69` | none | Equal rates and equal frequencies; spec `JukesCantor` |
| HKY | `hky` | `kappa` | Transition/transversion bias model |
| TN93 | `tn93` | `kappa1`, `kappa2` | Two transition rates |
| GTR | `gtr` | `rateAC`, `rateAG`, `rateAT`, `rateCG`, `rateCT`, `rateGT` (or `rates`) | General time-reversible model |
| SYM | `sym` | the six `rate*` keys (or `rates`) | Symmetric model: six free exchangeabilities with **equal** base frequencies — request them explicitly, see the note below |
| TIM | `tim` | `rateAG`, `rateCT`, `rateTransversions1`, `rateTransversions2` | Transitional model (two equal transversion rates ×2) |
| TVM | `tvm` | `rateAC`, `rateAT`, `rateCG`, `rateGT`, `rateTransitions` | Transversional model (four free transversions + one shared transition) |

**`rates:` is a GTR/SYM convenience only.** It must list **exactly six** values in the
order `rateAC rateAG rateAT rateCG rateCT rateGT`, and it cannot be combined with
explicitly named `rate*` keys, because the expansion would silently overwrite one of them.
The parser no longer pads a short vector with 1.0 and no longer truncates a long one.

TIM and TVM impose rate *equalities*, which a six-free-rate vector cannot express. Writing
`rates` on those two models is therefore an error, and the message names the parameters to
use instead.

```yaml
# GTR via the vector
substitution_model:
  type: "gtr"
  rates: {value: "1.0 4.0 1.0 1.0 4.0 1.0", lower: 0.0}
  frequencies: {mode: "estimated", value: 0.25, dimension: 4}

# SYM via the vector
substitution_model:
  type: "sym"
  rates: {value: "1.0 4.0 1.0 1.0 4.0 1.0", dimension: 6, lower: 0.0}

# TIM with its own parameters (see examples/config_subst_tim.yaml)
substitution_model:
  type: "tim"
  rateAG: {value: 4.0, lower: 0.0}
  rateCT: {value: 4.0, lower: 0.0}
  rateTransversions1: {value: 1.0, lower: 0.0}
  rateTransversions2: {value: 1.0, lower: 0.0}
  frequencies: {mode: "empirical"}

# TVM with its own parameters (see examples/config_subst_models.yaml)
substitution_model:
  type: "tvm"
  rateAC: {value: 1.0, lower: 0.0}
  rateAT: {value: 1.0, lower: 0.0}
  rateCG: {value: 1.0, lower: 0.0}
  rateGT: {value: 1.0, lower: 0.0}
  rateTransitions: {value: 4.0, lower: 0.0}
```

**Equilibrium frequencies (`frequencies:`).** Frequencies are available on the
nucleotide models that need them (`hky`, `tn93`, `gtr`, `sym`, `tim`, `tvm`); `jc69` and
the amino-acid matrices carry their own. Four explicit modes exist:

| `mode` | XML produced | Meaning |
|---|---|---|
| `estimated` (default) | `<frequencies frequencies="@…freqParameter">` + state node | A free probability vector, normalised to sum to 1, sampled and logged |
| `fixed` | `<frequencies><parameter estimate="false" value="…"/></frequencies>` | The generator writes your vector verbatim, outside the state space; it is neither sampled nor logged |
| `uniform` | `<frequencies spec="Frequencies" data="@aln" estimate="false"/>` | BEAST2's uniform-over-characters frequencies |
| `empirical` | `<frequencies spec="Frequencies" data="@aln" estimate="true"/>` | Frequencies that BEAST2 counts from the alignment |

> **Why `mode` exists.** In BEAST2, `<frequencies estimate="false">` means *uniform over
> characters*, and does **not** mean "use the numbers I gave you". An `estimate: false`
> therefore used to discard a user-supplied vector silently.

> Say what you mean now: `mode: fixed` writes your vector, `mode: uniform` asks BEAST2 for
> equal frequencies, `mode: empirical` asks BEAST2 to count them. For backward
> compatibility, the parser reads `estimate: false` on a vector without an explicit `mode`
> as `mode: fixed`. `uniform` and `empirical` need a `data_id` pointing at an alignment.

> A `fixed` or `estimated` vector must be a proper probability vector: every entry is
> greater than 0, the length matches `dimension`, and the entries sum to 1 within a 2%
> tolerance and are then renormalised to sum to exactly 1.

> **SYM caveat.** SYM *is* "GTR with equal equilibrium frequencies", but the generator does
> not impose that constraint for you: with no `frequencies:` block it emits a free
> (estimated) frequency vector, which makes the model effectively GTR. Write
> `frequencies: {mode: uniform}` to get the equal-frequency SYM that
> `examples/config_subst_sym.yaml` uses.

### 8.2 Amino acid substitution models

| Model | Config name | Description |
|---|---|---|
| WAG | `wag` | WAG matrix |
| JTT | `jtt` | JTT matrix |
| Dayhoff | `dayhoff` | Dayhoff matrix |
| BLOSUM62 | `blosum62` | BLOSUM62 matrix |
| CPREV | `cprev` | CPREV matrix |
| MTREV | `mtrev` | MTREV matrix (mitochondrial) |

```yaml
substitution_model: {type: "wag"}
```

These are `EmpiricalSubstitutionModel` instances that carry their own empirical
equilibrium frequencies, so the generator emits no `<frequencies>` element and accepts no
`frequencies:` block.

The parser cross-checks every model three ways: the declared `data_type`, the model family,
and the alphabet actually observed in the alignment. It therefore refuses `wag` on
nucleotide data, `hky` on amino-acid data, and `data_type: aminoacid` on a file that only
contains `A C G T`. The message names the problem, instead of a BEAST2 parse error far from
the cause.

```yaml
# Amino-acid data: the alignment must say so, and only amino-acid models are allowed
alignments:
  - id: "protein"
    file: "examples/aminoacid.fasta"     # 10 taxa x 234 aa
    data_type: "aminoacid"
partitions:
  - id: "protein"
    site_model:
      substitution_model: {type: "jtt"}
      gamma_categories: 4
      proportion_invariant: {value: 0.05, lower: 0.0, upper: 1.0}
```

### 8.3 Site model options

The site model adds rate heterogeneity on top of the substitution model:

```yaml
site_model:
  substitution_model: {type: "hky"}
  gamma_categories: 4          # Number of Gamma rate categories (0 = no Gamma)
  gamma_shape:                 # Gamma shape parameter; requires gamma_categories > 1
    value: 0.5
    lower: 0.0
    estimate: true             # Estimate or fix (default true)
  proportion_invariant:        # Proportion of invariant sites (optional)
    value: 0.01
    lower: 0.0
    upper: 1.0
    estimate: true
```

**Gamma+I combination:** setting both `gamma_categories > 0` and `proportion_invariant` yields the Gamma+I model.

**`gamma_shape` needs more than one category.** BEAST2's `SiteModel` documents the shape
input as "Ignored if gammaCategoryCount 1 or less", so the parser refuses `gamma_shape` when
`gamma_categories <= 1`: an estimated shape there would be a pure prior random
walk contributing nothing to the likelihood, yet still showing up as a trace column that
a reader would report as "the γ shape parameter".

---

## 9. Clock Models

| Model | Config name | Key parameters | Description |
|---|---|---|---|
| Strict clock | `strict` | `clock_rate` | Same rate on all branches |
| UC relaxed lognormal | `ucln` | `ucld_mean`, `ucld_stdev` | Uncorrelated relaxed lognormal clock |
| UC relaxed exponential | `uce` | `ucld_mean` | Uncorrelated relaxed exponential clock |
| Random local clock | `rlc` | (automatic) | Random local clock model |

`ucld_mean` and `ucld_stdev` are accepted **only** for `ucln` and `uce`. Writing either on
a `strict` or `rlc` clock is an error, with
`… only applies to the ucln/uce clocks; the strict clock has no such parameter`.
The UCLN rate distribution is lognormal with `M = 1` fixed and `S = ucld.stdev`;
the UCE rate distribution is exponential with mean 1.

**Configuration examples:**

```yaml
# Strict clock
clock_model:
  type: "strict"
  clock_rate: {value: 1.0, lower: 0.0}

# UCLN relaxed clock — both bounds are honoured, so ucld.stdev cannot drift
clock_model:
  type: "ucln"
  ucld_mean: {value: 1.0, lower: 0.0, upper: 10.0}
  ucld_stdev: {value: 0.333, lower: 0.0, upper: 1.0}

# UCE relaxed exponential clock
clock_model:
  type: "uce"
  ucld_mean: {value: 1.0, lower: 0.0}

# Random local clock (rate and indicator state nodes are created automatically)
clock_model:
  type: "rlc"

# Link to another partition's clock
clock_model:
  linked_to: "gene1"    # Share the clock model of gene1
```

**Notes**

- **The generator honours both bounds.** It used to drop `upper` for `clock_rate`,
  `ucld.mean` and `gamma_shape`, which left `ucld.stdev` unbounded above — the classic UCLN
  standard-deviation drift. Both `lower` and `upper` now reach the XML, and when you give
  neither the generator applies the role defaults (`ucld.mean` and `clock.rate` above 0;
  `ucld.stdev` within `[0, 1]`).
- **A fixed clock rate must not be the only time scale.** The generator now honours
  `clock_rate: {estimate: false}` as well: it writes the parameter inline with
  `estimate="false"`, so that parameter leaves the state space, the operator set and the
  trace. If no calibration accompanies it, absolute time is unidentifiable and the parser
  refuses the run (the message is listed
  in [Section 22.1](#221-strict-validation-rejected-inputs-and-their-error-messages)).
- **`ucld.mean` normalisation.** The generator emits every relaxed-clock XML with
  `normalize="true"`, matching BEAUti's reference layout, so rates and node ages stay
  identifiable.

---

## 10. Tree Priors

| Tree prior | Config name | Key parameters | Description |
|---|---|---|---|
| Yule | `yule` | `birth_rate` | Pure birth process |
| CalibratedYule | `calibrated_yule` | `birth_rate` | Yule model with integrated calibrations |
| Birth-Death | `birth_death` | `birth_rate`, `death_rate`, `sample_probability` (or `sampling_rate`) | Birth-death process (Gernhard 2008 parameterisation) |
| Coalescent constant | `coalescent_constant` | `pop_size` | Constant population coalescent |
| Coalescent exponential | `coalescent_exponential` | `pop_size`, `growth_rate` | Exponential growth coalescent |
| Bayesian Skyline | `bayesian_skyline` | (auto-computed) | Bayesian Skyline Plot |
| EBSP | `ebsp` | (auto-computed) | Extended Bayesian Skyline Plot (built into the BEAST 2.7 core package) |
| BD Skyline Serial | `bd_skyline_serial` | `birth_rate`, `death_rate`, `sampling_rate`, `rho` | Piecewise-constant birth-death skyline with serial sampling (requires the BDSKY add-on package) |

**Admissible domains.** The parser enforces the constraints below. Without them the XML would
be written and would only fail later, inside BEAST2:

- `birth_death` requires `0 <= death_rate < birth_rate` **and** `0 < sample_probability <= 1`.
  `BirthDeathGernhard08Model` is parameterised by `birth - death` (>= 0) and
  `death / birth` (< 1), so a super-critical process (`death_rate >= birth_rate`) has no
  valid state there; the error message points you at `bd_skyline_serial` (BDSKY), which
  estimates birth, death and sampling rates separately and does allow `death >= birth`.
  `birth_rate` and `death_rate` must also be estimated or fixed *together*, since only the
  derived quantities exist as state nodes.
- `bd_skyline_serial` requires `birth_rate >= 0`, `death_rate >= 0`, and `sampling_rate`
  and `rho` within `[0, 1]`.
- `yule`, `calibrated_yule`, `birth_death` need `birth_rate > 0`; the coalescent priors
  need `pop_size > 0`. The `birth_rate` of a Yule-type prior is the prior's own parameter,
  so the parser rejects `estimate: false` on it instead of ignoring it silently.
- `coalescent_exponential`'s `growth_rate` is *signed* (a negative value means a declining
  population) and is therefore left unbounded, with a broad normal default prior.

**Default priors.** When an estimated continuous parameter has no `parameter_priors` entry,
the generator gives it a weakly informative default
(`xml_writer.py::_default_prior_element`): positive rates get a broadly dispersed lognormal
centred on the initial value; `ucld.stdev` gets `Normal(mean 0.4, sigma 0.3334)` (BEAUti's
convention), truncated to the parameter's range; probability and frequency parameters get
`Uniform(lower, upper)`; signed parameters get a broad normal. Discrete state nodes
(`IntegerParameter`, `BooleanParameter`) need no prior.

Without this step, a flat prior on (0, ∞) would be an improper density, and the analysis
could not be reconstructed from the XML.

**Configuration examples:**

```yaml
# Yule model
tree_prior:
  type: "yule"
  birth_rate: {value: 1.0, lower: 0.0}

# Birth-Death model — death_rate must be strictly smaller than birth_rate
tree_prior:
  type: "birth_death"
  birth_rate: {value: 1.0, lower: 0.0}
  death_rate: {value: 0.5, lower: 0.0}
  sample_probability: {value: 0.1}

# Constant population coalescent
tree_prior:
  type: "coalescent_constant"
  pop_size: {value: 10.0, lower: 0.0}

# Exponential growth coalescent (growth_rate may be negative)
tree_prior:
  type: "coalescent_exponential"
  pop_size: {value: 10.0, lower: 0.0}
  growth_rate: {value: 0.1}

# Bayesian Skyline
tree_prior:
  type: "bayesian_skyline"

# CalibratedYule (calibration_method must also be set)
calibration_method: "calibrated_yule"
tree_prior:
  type: "calibrated_yule"
  birth_rate: {value: 1.0, lower: 0.0}
```

**Naming parameters in `parameter_priors`.** A `parameter_priors:` entry is matched against
the state nodes the generator actually built. The alias table below lists, on the left, the
name written in the configuration and, on the right, the id that appears in the XML:

| You write | Matches XML id(s) |
|---|---|
| `birth_rate` | `birthRate`, `birthDiffRate` |
| `death_rate` | `relativeDeathRate` |
| `sample_probability`, `sampling_rate` | `sampleProbability`, `bdSkylineSamplingRate` |
| `pop_size` | `popSize` |
| `pop_sizes` | `bPopSizes`, `populationSizes` |
| `growth_rate` | `growthRate` |
| `clock_rate` | `clockRate` |
| `ucld_mean` / `ucld_stdev` | `ucld.mean` / `ucld.stdev` |
| `gamma_shape` / `proportion_invariant` | `gammaShape` / `proportionInvariant` |
| `kappa` | `kappa` (any partition's `…hky.kappa`) |

A prior whose target cannot be resolved is a **hard error**, and the message lists the state
nodes that are currently estimated, because such a prior would silently never be evaluated.
The same applies in reverse: naming a parameter you fixed with `estimate: false` is an error
too, and the fix is to delete that prior.

---

## 11. Calibration Distributions

### 11.1 Supported distribution types

| Distribution | Config name | Parameters | Constraints | Typical use |
|---|---|---|---|---|
| Normal | `normal` | `mean`, `sigma` | `sigma > 0` | Well-constrained symmetric calibrations |
| LogNormal | `lognormal` | `M`, `S`, `mean_in_real_space` | `S > 0`; `M > 0` only when `mean_in_real_space: true` | Right-skewed calibrations (common for fossils) |
| Uniform | `uniform` | `lower`, `upper` | both required, `lower < upper` | Hard boundary constraints |
| Exponential | `exponential` | `mean` | `mean > 0` | Minimum age constraints |
| Gamma | `gamma` | `alpha`, `beta` | both `> 0` | Flexible skewed distribution |
| Beta | `beta` | `alpha`, `beta` | both `> 0` | Bounded distribution |
| Laplace | `laplace` | `mu`, `scale` | `scale > 0` | Symmetric distribution with heavier tails |
| InverseGamma | `inverse_gamma` | `alpha`, `beta` | both `> 0` | Right-skewed distribution on positive reals |
| OneOnX | `one_on_x` | (none) | — | Scale-invariant prior |
| Poisson | `poisson` | `lambda` | `lambda > 0` | Count data |
| ChiSquare | `chi_square` | `dof` | `dof > 0` (written to BEAST2's `df` input) | Chi-square distribution |

`uniform` is checked twice over: the parameters themselves must satisfy `lower < upper`
(after adding `offset`), because an inverted `Uniform` has zero density everywhere; and
`provenance.original_age_min` may not exceed `original_age_max`. For `lognormal`, `M` is
the **log-scale** mean (median = `exp(M)`) unless `mean_in_real_space` is set, in which
case `M` is a real-space mean and must be positive.

**Configuration examples:**

```yaml
# Normal distribution
distribution:
  type: "normal"
  parameters: {mean: 6.0, sigma: 0.5}

# LogNormal distribution (M is the mean in log space)
distribution:
  type: "lognormal"
  parameters: {M: 1.0, S: 0.5, mean_in_real_space: false}
  offset: 12.0    # The distribution starts at 12.0 Ma

# LogNormal distribution (with mean_in_real_space=true, M is the mean in real space)
distribution:
  type: "lognormal"
  parameters: {M: 15.0, S: 0.5, mean_in_real_space: true}

# Uniform distribution
distribution:
  type: "uniform"
  parameters: {lower: 20.0, upper: 40.0}

# Exponential distribution
distribution:
  type: "exponential"
  parameters: {mean: 5.0}

# Gamma distribution
distribution:
  type: "gamma"
  parameters: {alpha: 2.0, beta: 0.5}
```

### 11.2 Offset support

All distributions support the `offset` attribute, which shifts the origin of the distribution. This is useful when a calibration point has a minimum age constraint:

```yaml
distribution:
  type: "lognormal"
  parameters: {M: 1.0, S: 0.5}
  offset: 12.0    # The distribution starts at 12.0 Ma, i.e. actual age = offset + LogNormal(M, S)
```

### 11.3 Hyperpriors

Calibration distribution parameters can themselves carry prior distributions (hyperpriors):

```yaml
calibrations:
  - name: "humanChimpMRCA"
    taxa: ["human", "chimp"]
    distribution:
      type: "normal"
      parameters: {mean: 6.0, sigma: 0.5}
      hyperpriors:
        mean:
          type: "uniform"
          parameters: {lower: 0.0, upper: 20.0}
        sigma:
          type: "uniform"
          parameters: {lower: 0.0, upper: 5.0}
```

Each hyperprior key must name a parameter of the distribution it sits under — a
`hyperpriors: {mu: …}` block on a `normal` (whose parameters are `mean`/`sigma`)
makes the parser refuse the run, and the hyperprior's own parameters are range-checked
like any calibration distribution. A hyperprior turns that parameter into an estimated state node, which is
why the generator also gives it an operator and a `Prior`.

### 11.4 Calibration provenance metadata

Each calibration point can carry provenance metadata recording the fossil source and calibration type:

```yaml
provenance:
  source: "fossil"              # Source type: fossil / biogeographic / secondary / other
  reference: "DOI:10.1038/xxx"  # Reference DOI or URL
  calibration_type: "soft"      # soft (soft bound) / hard (hard bound)
  original_age_min: 5.5         # Original minimum fossil age
  original_age_max: 7.0         # Original maximum fossil age
  notes: "Fossil description"   # Free-text notes
```

**The record and the prior must not contradict each other.** The parser compares
`original_age_min` / `original_age_max` against the 95% interval of the distribution actually
used (computed with the tool's own `compute_distribution_stats`, offset included):

| Situation | Result |
|---|---|
| `calibration_type: hard` and the 95% interval lies entirely outside the recorded age range | **error** — a hard-bound calibration must really contain the age it claims to enforce |
| `calibration_type: soft` and the same mismatch | **warning** — legitimate: a soft bound deliberately puts mass below the fossil minimum (for example a lognormal whose 2.5% quantile sits at the hard minimum), but the methods text should say so |
| `original_age_min > original_age_max` | **error** |

```text
Calibration 'root': provenance records an original minimum age of 20.0 but the 95%
interval of the normal prior is [4.02, 5.98], i.e. entirely younger than the fossil.
A hard-bound calibration must contain the age it claims to enforce; either change
calibration_type to soft (with a soft-bound distribution) or move the prior.
```

### 11.5 Other calibration options

```yaml
calibrations:
  - name: "rootCalibration"
    taxa: null                    # null means a root calibration
    monophyletic: true            # Enforce the monophyly constraint (default true)
    use_originate: true           # Use the stem (originate) node instead of the MRCA
    tipsonly: false               # Tip-only calibration (default false)
```

### 11.6 Further reading

- [CladeAge](https://github.com/BEAST2-Dev/cladeage) (BEAST2 add-on): generates calibration densities from birth-death rates and fossil first-occurrence ages by simulation. It is an alternative source of calibration priors to hand-picking Normal/LogNormal parameters.
- [FossilCalibrations](https://fossilcalibrations.org): a curated database of fossil calibrations; its DOI and age-range fields map directly onto the `provenance` block above.

---

## 12. Multi-Partition Analysis

### 12.1 Multiple sequence files

Use multiple independent sequence files as separate partitions:

```yaml
alignments:
  - id: "gene1"
    file: "gene1.fasta"
    data_type: "nucleotide"
  - id: "gene2"
    file: "gene2.fasta"
    data_type: "nucleotide"

partitions:
  - id: "gene1"
    site_model:
      substitution_model: {type: "hky"}
    clock_model:
      type: "ucln"
      ucld_mean: {value: 1.0}
      ucld_stdev: {value: 0.333}
    tree: "shared"
  - id: "gene2"
    site_model:
      substitution_model: {type: "gtr"}
    clock_model:
      linked_to: "gene1"    # Share the clock model with gene1
    tree: "shared"
```

### 12.2 Codon partitioning

Use FilteredAlignment to implement codon partitioning over different site subsets of the same sequence file:

```yaml
alignments:
  - id: "alignment"
    file: "coding.fasta"
    data_type: "nucleotide"
  - id: "codon12"
    data_id: "alignment"       # Reference to the source alignment (required)
    filter: "1::3,2::3"        # 1st + 2nd codon positions
    data_type: "nucleotide"
  - id: "codon3"
    data_id: "alignment"
    filter: "3::3"             # 3rd codon positions
    data_type: "nucleotide"

partitions:
  - id: "codon12"
    site_model:
      substitution_model: {type: "hky"}
      gamma_categories: 4
    clock_model:
      type: "ucln"
      ucld_mean: {value: 1.0}
      ucld_stdev: {value: 0.333}
    tree: "shared"
  - id: "codon3"
    site_model:
      substitution_model: {type: "gtr"}
      gamma_categories: 4
    clock_model:
      linked_to: "codon12"    # Share the clock of codon12
    tree: "shared"
```

**Filter expression syntax.** The expression follows BEAST2's own `FilteredAlignment`
grammar and is **1-based**; it is a comma-separated list of:

| Form | Meaning |
|---|---|
| `N` | a single site |
| `from-to` | the inclusive range; an empty `from` means site 1, an empty `to` means the last site |
| `from-to\step` | the same range, taking every `step`-th site |
| `from:to:step` | an iterator; `from:to:` or `from::step` leave the empty parts at their defaults |

So `1::3,2::3` selects positions 1, 2, 4, 5, 7, 8 and so on (codon positions 1 and 2), and
`3::3` selects 3, 6, 9 and so on (codon position 3). The parser range-checks every term
against the source alignment, and reports an error naming the offending term for an
unparsable term, a `step < 1`, a coordinate beyond the alignment length, or a term that
selects nothing.

**`data_id` is required whenever a filter is declared.** The source alignment named by
`data_id` must already have been declared, and may not be the filtered alignment itself: a
filter pointing at its own id would produce `data="@self"` and drop all sequence data.

The pattern above — declare the source once, then point each filtered partition at it with
`data_id` — is what `examples/config_multi_partition.yaml` and
`examples/config_advanced.yaml` use. It is also what puts the unfiltered copy of the data
into the XML exactly once. If a filtered alignment omits `data_id`, the parser reports an
error (`… has a filter but no data_id: it would reference itself and lose all sequence data`).

**Two partitions must not cover the same sites.** The parser does not stop the same FASTA
from being declared twice under two ids, with a `TreeLikelihood` built over each. But the
resulting posterior is the product of two identical likelihoods, which counts the same 898
sites twice, and BEAST2 runs it without complaint. Split model comparisons into separate
runs (one configuration per candidate model, as in
`examples/config_subst_jc69.yaml`, `config_subst_sym.yaml`, `config_subst_tim.yaml`,
`config_subst_models.yaml`) and use disjoint filters for partitioned analyses.

### 12.3 Model linking rules

- **`clock_model.linked_to: "partition_id"`**: the partition shares the clock model of another partition.
- **`site_model.linked_to: "partition_id"`**: the partition shares the site model of another partition (see the `config_advanced.yaml` example). When not specified, each partition has an independent site model.
- A `linked_to` value must name another existing partition — self-links and unknown ids are errors.

### 12.4 Topology: only one linked tree is assembled

**`tree: "shared"` is the only legal value.** Per-partition or unlinked topologies are
not implemented, so the parser refuses the setting outright instead of parsing, storing and
then silently ignoring it:

```text
Partition 'codon3': tree: 'separate' requests a tree other than the shared one, but
Beast2Py assembles every partition on a single linked topology (per-partition /
unlinked trees are not implemented). Silently falling back to the shared tree would
turn an intended *BEAST-style multi-locus analysis into a fully linked one, so this
is an error. Remove the 'tree' key to analyse the shared topology, or run one
analysis per independent partition.
```

`tree: <some partition id>` is rejected the same way, with the same message. The consequence
to keep in mind: a multi-locus analysis here assumes one shared species/gene tree, and each
partition's relaxed clock indexes that same tree. Unlinked gene trees (the `*BEAST` and
StarDivergence use case) are on the roadmap, not in this release.

---

## 13. Tip Dates

Tip dates serve analyses whose samples were taken at different times, such as infectious
disease studies.

```yaml
tip_dates:
  enabled: true
  trait_name: "date-forward"    # date-forward or date-backward, nothing else
  units: "year"                # year, month, day, or coordinate
  dates:                        # one entry per taxon — coverage must be complete
    sample_2020: 2020
    sample_2019: 2019
    sample_2018: 2018
  # Alternative/extra source: a two-column taxon,date CSV
  # file: "tip_dates.csv"
```

**Rules the parser now enforces**

- **Complete coverage.** `dates` must cover every taxon of every partition. BEAST2 treats an
  undated tip as sampled at time zero, so a partial mapping quietly turned the missing
  lineages into contemporaneous samples. The parser now reports exactly which taxa are
  missing (`tip_dates.dates covers 9 of 12 taxa; these are missing: …`), and it refuses taxa
  that are absent from the alignments.
- **`trait_name` ∈ {`date-forward`, `date-backward`}.** `date-forward` counts forward from the
  oldest sampling point (0 is the oldest sample); `date-backward` counts back from the most
  recent one (0 is the youngest sample). Any other value still yields valid XML, but BEAST2
  does not read it as a date, so the analysis runs as contemporaneous sampling. The parser
  rejects `date_backward` (underscore) and lists the valid names.
- **`units` ∈ {`year`, `month`, `day`, `coordinate`}.**
- **Dates must be numeric.** The parser reports an error for a CSV row whose second column is
  not a number, except a header line.
- **Keep one time axis.** `units`, the tip values and the calibrations must share one time
  axis. `units: year` with dates 0–20 next to a root calibration measured in Ma puts the two
  six orders of magnitude apart. Express everything in years (and scale the clock rate, e.g.
  2 × 10⁻³ subs/site/Myr = 2 × 10⁻⁶ subs/site/year), or keep everything in Ma.

**What gets emitted.** Enabling tip dates makes the generator attach the `TraitSet` to the
tree *and* emit a `TipDatesRandomWalker` operator over the dated taxa, with no `tipsonly`
calibration involved. Previously that operator was gated on a `tipsonly` calibration, so the
common BEAUti-style configuration that treats "dates are a trait, not a calibration" produced
a `TraitSet` no operator ever moved, and the run returned a tree that read as dated but was
not.

The operator covers only tips whose date is non-zero. A fully contemporaneous sample
therefore emits no `TipDatesRandomWalker`, because there is no date to move — that is the
case of `examples/config_tip_dates.yaml`, whose 12 primates are all extant. The behaviour is
verified against BEAST 2.7.8: a variant of that configuration with one tip dated 1.5 Ma
earlier parses and initialises.

**`trait_name` semantics:** `date-forward` means dates increase with sample age from the
oldest sample; `date-backward` means the value is the distance back from the most recent
sample, so 0 = the youngest tip (this is what BEAUti's "Tip dating" tab produces).

---

## 14. Calibration Diagnostics

### 14.1 Overview

Calibration diagnostics is the core innovation of Beast2Py and consists of three sub-modules: the conflict detector (`ConflictDetector`), the sensitivity analyzer (`SensitivityAnalyzer`), and the visualization generator (`VisualizationGenerator`).

### 14.2 Conflict detection

The conflict detector analyzes pairs of calibration points for three classes of problems:

**Temporal consistency conflicts:** for nested clades (taxon set A ⊂ B), it checks whether the parent clade (B) is calibrated older than the child clade (A). It compares the 95% HPD intervals: the detector reports an error when the parent's HPD upper bound is younger than the child's HPD lower bound, and a warning when only the means disagree while the HPD intervals still overlap.

**Distribution overlap warnings:** it computes the 1-Wasserstein distance (Earth Mover's Distance) between calibration distributions. For a non-nested pair whose distance is too small, it reports a warning that points at possibly redundant information. **The default criterion is the dimensionless one** — the Wasserstein distance divided by the mean width of the two priors' 95% intervals, warned below `relative_overlap_threshold` (default 0.15). The absolute criterion (distance below `overlap_threshold`, default 2.0 time units) is opt-in, because a fixed distance is scale-dependent: it grows with node depth, so one threshold flags shallow nodes almost automatically and never reaches the deep ones. See [Section 14.9](#149-the-diagnostics-block) for how to switch and tune this.

**Monophyly conflicts:** it detects taxon sets that intersect but are not nested. When both calibrations enforce monophyly, that topology cannot exist, and the detector reports an error.

### 14.3 Sensitivity analysis

The sensitivity analyzer generates N+1 prior-only sampling XMLs (N is the number of calibration points): one version with the full calibration set, and N leave-one-out versions (removing one calibration point at a time). The user runs these short-chain XMLs with BEAST2 (1 million generations is enough) and compares the tree-height posteriors, which reveals the marginal influence of each calibration point. (Sampling from the prior alone is a built-in BEAST2 capability, available for any XML via `beast -sampleFromPrior`; what the sensitivity analyzer automates is the N+1 leave-one-out design.) The files land in a `sensitivity/` subfolder of the report's `diagnostics/` directory; when you set `--sensitivity-xml-dir`, they land in the `sensitivity/` subfolder of that directory instead.

### 14.4 Visualization

The visualization generator uses matplotlib to create density plots of all calibration distributions, annotating the 95% HPD intervals (shaded regions) and means (dashed lines). It supports PNG, PDF and SVG output.

### 14.5 HTML report

The diagnostic report is generated as a self-contained HTML file containing an analysis overview, a calibration summary table, conflict detection results, distribution visualizations, and sensitivity analysis instructions.

### 14.6 Running diagnostics from the CLI

```bash
beast2py diagnose \
    --config config.yaml \
    --report report.html \
    --sensitivity \
    --visualize \
    --sensitivity-xml-dir sensitivity_xmls/
```

### 14.7 Running diagnostics from the API

```python
from beast2py.api import Beast2Py

b2p = Beast2Py()

# Full diagnostics
report = b2p.diagnose(
    "config.yaml",
    output_dir="diagnostics/",
    run_sensitivity=True,
    run_visualization=True,
    report_path="report.html",
)

# Conflict detection only
conflicts = b2p.detect_conflicts("config.yaml")

# Sensitivity XMLs only
xmls = b2p.generate_sensitivity_xmls("config.yaml", output_dir="sensitivity/")

# Density plots only
b2p.plot_calibrations("config.yaml", output_file="calibrations.pdf")
```

### 14.8 Interpreting conflict detection results

Conflict detection returns a list of `Conflict` objects, each with the following attributes: `conflict_type` (conflict type: temporal, distribution_overlap, monophyly), `severity` (severity: error, warning, info), `calibration_a` and `calibration_b` (the calibration names involved), `description` (detailed description), and `recommendation` (suggested action).

- **error** indicates a problem that must be fixed (e.g. a temporal inconsistency or monophyly conflict).
- **warning** indicates an issue that deserves attention but is not fatal (e.g. distribution overlap or mean inconsistency).
- **info** indicates informational messages.

### 14.9 The `diagnostics:` block

The `diagnostics:` section configures the thresholds behind the distribution-overlap check. The same section is passed to Gate 2 of `generate`:

```yaml
diagnostics:
  overlap_measure: "relative"           # relative (default) | absolute | both
  overlap_threshold: 2.0                # absolute criterion: Wasserstein distance in time units
  relative_overlap_threshold: 0.15      # dimensionless criterion: W / mean 95%-interval width
```

- `overlap_measure` decides which criteria run. `relative` (the default) is the scale-free ratio described above; `absolute` restores the previous fixed-distance test; `both` reports a warning when either criterion trips. The parser rejects an unrecognised value with `diagnostics.overlap_measure must be relative, absolute or both, got ...`.
- `overlap_threshold` (default `2.0`) and `relative_overlap_threshold` (default `0.15`) are the two thresholds, and only the criteria named by `overlap_measure` take effect.

Once the dimensionless criterion became the default, a fixed `2.0` stopped meaning the same thing across a tree with widely divergent node ages. A warning now fires when two priors are nearly interchangeable *relative to their own widths*, which is the scale-free notion of redundant calibration information.

---

## 15. Reproducibility Tools

### 15.1 Analysis fingerprint

The analysis fingerprint is a deterministic identifier generated from the SHA-256 hash of the *scientific* configuration, in the date-free format `B2P-{hash12}-{version}`, e.g. `B2P-ba18c26ed61f-0.1.0`. It carries **no date and no time zone**, so the same configuration over the same alignment data reproduces byte-for-byte on every rerun.

These changes alter the fingerprint, because `models.py::to_dict` includes them in the hash: the sequence content, the tree setting (`Partition.tree`), a crown/stem (`use_originate`) or `tipsonly` flag, a hyperprior, a `parameter_priors` entry, a tip date, or an operator weight. Renaming the output file alone does not.

Besides the identifier, the XML comment records a separate digest of each alignment's content: `<!-- Analysis Fingerprint: B2P-{hash12}-{version} | Data: {hash16} -->`, where `Data:` is 16 hex characters taken from the per-alignment sequence digests. The fingerprint travels as an XML comment and can also be written to a JSON sidecar file. The wall-clock `generation_time` appears **only** in that sidecar, never in the identifier or in the XML. See [Section 5.7](#57-fingerprint--generate-an-analysis-fingerprint) for the verified JSON contents.

```bash
# CLI
beast2py fingerprint --config config.yaml --output fingerprint.json --verbose

# Or during generation
beast2py generate --config config.yaml --output output.xml --fingerprint
```

```python
# API
fp = b2p.generate_fingerprint("config.yaml", output="fingerprint.json")
print(f"Analysis fingerprint: {fp}")
```

### 15.2 Methods description

Automatically generates a LaTeX methods paragraph describing the substitution model (with citations), clock model, tree prior, calibration points (with distribution parameters and provenance), MCMC settings, and the analysis fingerprint. It can be inserted directly into a manuscript.

```bash
beast2py methods --config config.yaml --output methods.tex
```

```python
methods_text = b2p.generate_methods("config.yaml", output="methods.tex")
```

### 15.3 Pipeline file generation

Generates Snakemake or Nextflow pipeline files containing complete rules/processes for XML generation, BEAST2 execution, diagnostics, and methods description.

```bash
# Snakemake
beast2py pipeline --config config.yaml --type snakemake --output-dir pipeline/

# Nextflow
beast2py pipeline --config config.yaml --type nextflow --output-dir pipeline/

# Both
beast2py pipeline --config config.yaml --type both --output-dir pipeline/
```

```python
# API
path = b2p.generate_pipeline("config.yaml", pipeline_type="snakemake", output_dir="pipeline/")
sm_path, nf_path = b2p.generate_pipeline("config.yaml", pipeline_type="both", output_dir="pipeline/")
```

The generated Snakemake file is `Snakefile` and the Nextflow file is `main.nf`.

---

## 16. Validation

### 16.1 Structural validation

Structural validation (`validator.py::XMLValidator`) runs two tiers in one call.

**Structural.** This tier checks XML well-formedness and the following:

- The root element is `<beast>`
- The required elements `<data>` and `<run>` are present
- id/idref/`@ref` references are consistent: every idref names a defined id, and no id repeats
- The `<state>` element is present
- The posterior distribution reference and the logger elements are present
- The namespace attribute contains `beast.base`

**Semantic** (`validator.py::_check_semantic`) covers the checks that make an analysis
*correct*, not merely parseable:

- Every estimated state node carries a `Prior`. A flat prior on an unbounded support is an
  improper density and cannot be reconstructed from the XML; bounded discrete
  `IntegerParameter`/`BooleanParameter` nodes are exempt.
- Every estimated state node carries an operator, and no operator touches a fixed parameter.
- No state node has inverted bounds (`lower < upper`).
- Every initial value lies inside its bounds and matches the declared dimension.
- Each `<data>` block has non-empty, equal-length sequences with unique taxon names.

A `storeEvery` that disables checkpointing is reported as a warning rather than an error.

```bash
beast2py validate --xml output.xml
```

```python
result = b2p.validate_xml("output.xml")
print(f"Valid: {result.is_valid}")
print(f"Errors: {result.errors}")
print(f"Warnings: {result.warnings}")
```

### 16.2 BEAST2 native validation

`--beast2` validates through BEAST2's own `XMLParser`, which is the stricter check. Beast2Py
ships a headless validation tool (`beast2py/tools/Beast2Validator.java`) that runs without a
JavaFX environment. It parses the file **and** calls `initAndValidate()` on every object, so
what is checked is whether the model can **complete initialisation**, not merely whether it
parses. A misspelled optional attribute, a never-firing operator or an illegal constraint all
fail here instead of printing `VALID`.

The Java helper, its `launcher.jar` and `classes/` live inside the package and ship with the
wheel as `package-data`. `validator.py` locates `beast2_validate.sh` next to the installed
package, in the project root, on `PATH`, or wherever `--beast2-path` points.

```bash
beast2py validate --xml output.xml --beast2 --beast2-path /path/to/beast
```

```python
result = b2p.validate_xml("output.xml", beast2_validate=True, beast2_path="/path/to/beast")
```

This is the same code path as Gate 3 of `generate`: it needs a BEAST2 installation and a
**JDK 17 or newer**, and `--beast2-path` may point at a specific installation or script. The
bundled script looks for that JDK among the candidates listed in
[Section 2.6](#26-installing-beast2-optional-for-the-beast2-validation-gate) instead of
trusting `PATH` alone, because an older `java` on `PATH` should not shadow a usable
`openjdk@17`. If no candidate qualifies, the script exits 2 with a clear message.

`beast2_validation_status` distinguishes three outcomes: `passed`, `failed` and `unavailable`.
A missing BEAST2 installation is therefore reported as *unavailable*, never as a pass. Because
you asked for the check explicitly, generation stops with exit code 2 unless you pass
`--allow-unvalidated`.

### 16.3 Validation status

All **19** example configuration files pass all three gates: structural checks, calibration-conflict detection, and the BEAST2 `parseFile` and `initAndValidate` check. The third gate uses the native BEAST2 **v2.7.8** parser in headless mode, so no JavaFX is required. Two of the files depend on add-on packages — `config_bd_skyline.yaml` and `config_nested_sampling.yaml` — and were verified against the packages installed under `~/.beast/2.7/`, namely **BDSKY 1.5.1** and **NS 1.2.0**.

The test suite comprises **235** collected tests: **214** unit/semantic tests and **21** BEAST2 integration tests. The integration tests generate XML from all 19 example configurations plus the `quick` path and the committed `output_basic.xml`, and check each with the real BEAST 2.7.8 parser and model initialisation. When BEAST2 or JDK 17 is absent, those tests simply skip instead of failing. All are passing.

---

## 17. MCMC Settings

### 17.1 Standard MCMC

```yaml
mcmc:
  type: "standard"
  chain_length: 10000000      # Chain length (total number of iterations)
  pre_burnin: 0               # Pre-burnin iterations (must be < chain_length)
  store_every: 1000           # Checkpoint interval in *logged samples* (see below)
  sample_from_prior: false    # Sample from the prior only (for sensitivity analyses)
  seed: 1                     # RNG seed for `beast -seed` (NOT written into the XML; see below)
```

**`store_every` counts logged samples, not generations.** BEAST2's `MCMC` tests `(sampleNr + 1) % storeEvery == 0` (`MCMC.java:490`). With the default `logEvery`, `store_every` is therefore the number of *saved samples* between `.xml.state` checkpoints, that is, how many resume points an interrupted run has. The default **scales with the chain length**, `max(1, chain_length // 10000)`: `1000` for a 10,000,000-generation chain, `100` for a 1,000,000 one. Every run thus gets roughly ten checkpoints instead of only the last. An explicit `-1` is still legal but disables checkpointing, and the validator then warns `storeEvery='-1' ... cannot be resumed`.

**The generator deliberately keeps `mcmc.seed` out of the XML.** BEAST2's `MCMC` has no `seed` *input*, so writing one makes the model unparseable (`has no input with name seed`). The RNG seed is a command-line argument (`beast -seed N`). Beast2Py keeps `mcmc.seed` in the configuration for the CLI and the launcher to consume. Logger file names may still carry the `$(seed)` macro, which BEAST2 expands at run time (see [Section 18](#18-logger-configuration)).

### 17.2 Nested Sampling MCMC

Nested Sampling is used for model comparison and marginal likelihood estimation and requires the BEAST2 NS add-on package.

```yaml
mcmc:
  type: "nested_sampling"
  chain_length: 10000000      # Length of each sub-chain
  particle_count: 1           # Number of particles
  sub_chain_length: 10000     # Sub-chain length
```

### 17.3 Prior-only sampling

Setting `sample_from_prior` to `true` generates an XML that samples from the prior only (sequence data unused). This is normally set automatically by the sensitivity analyzer, but it can also be used manually:

```yaml
mcmc:
  type: "standard"
  chain_length: 1000000       # A short chain suffices
  sample_from_prior: true
```

**Tip:** the same effect is available for any existing BEAST2 XML without editing the file, through the BEAST2 command line: `beast -sampleFromPrior output.xml`. Beast2Py's added value is generating the full set of N+1 leave-one-out sensitivity XMLs automatically (see Section 14.3).

---

## 18. Logger Configuration

Loggers control the frequency and format of the output files written during a BEAST2 run.

```yaml
loggers:
  trace_log:                    # Parameter trace log (.log file)
    file_name: "output.$(seed).log"
    log_every: 1000             # Log every 1000 generations
  tree_log:                     # Tree log (.trees file)
    file_name: "output.$(seed).trees"
    log_every: 1000
  screen_log:                   # Screen output
    log_every: 10000            # Print every 10000 generations
```

**The `$(seed)` variable:** BEAST2 replaces `$(seed)` in a file name with the random seed value, so repeated runs do not overwrite each other.

---

## 19. Initialization Strategies

Beast2Py supports three initial tree generation strategies:

```yaml
initialization:
  tree_type: "random"     # Random tree (default)
  # tree_type: "upgma"    # UPGMA clustering tree
  # tree_type: "newick"   # User-provided Newick tree
  newick_file: null        # Tree file path for the newick type
```

- **`random`**: generates a random topology, suitable for most analyses.
- **`upgma`**: builds a UPGMA clustering tree from sequence distances, which gives a better starting point.
- **`newick`**: reads a user-provided Newick-format tree file, for analyses that need a specific starting topology.

---

## 20. Operator Weights

Users may override the default MCMC operator weights. The default tree operator weights are: `tree_scaler` 3.0, `uniform` 30.0, `subtree_slide` 15.0, `exchange_narrow` 15.0, `exchange_wide` 3.0, `wilson_balding` 1.0, and `up_down` 3.0.

```yaml
operator_weights:
  tree_scaler: 3.0        # Tree scaling (default 3.0)
  uniform: 30.0           # Uniform operator (default 30.0)
  subtree_slide: 15.0     # Subtree slide (default 15.0)
  exchange_narrow: 15.0   # Narrow exchange (default 15.0)
  exchange_wide: 3.0      # Wide exchange (default 3.0)
  wilson_balding: 1.0     # Wilson-Balding operator (default 1.0)
  # ... other operators
```

This feature is intended for advanced users. Normally you should leave the default weights alone.

---

## 21. Example Configuration Files

The `examples/` directory contains **19** complete example configuration files. Counted
honestly, they exercise **all 4 clock models** (strict, UCLN, UCE, RLC),
**7 of the 8 tree priors** (every one except `coalescent_constant`), and
**9 of the 13 substitution models**: all 7 nucleotide models (JC69, HKY, TN93, GTR, SYM,
TIM, TVM) plus two amino-acid matrices (WAG, JTT).

The four unexemplified models are the amino-acid Dayhoff, BLOSUM62, CPREV and MTREV.
They do **not** put "every model in one file": substitution-model comparison is split
across `config_subst_jc69/sym/tim/models.yaml`, one model per run. Several likelihoods on
one shared topology give a posterior self-product, not a comparison (see Section 12.2).

| File | Scenario | Tree prior | Substitution model | Clock model | Calibrations | Special features |
|---|---|---|---|---|---|---|
| `config_basic.yaml` | Basic analysis | Yule | HKY+Γ4 | Strict | 1 | Minimal configuration; explicit `parameter_priors` + provenance; source of `output_basic.*` |
| `config_advanced.yaml` | Advanced multi-calibration | Birth-Death | GTR+Γ4+I / linked site model | UCLN + linked | 3 | Codon partitions + site *and* clock linking + hyperpriors + `diagnostics` block + stem calibration |
| `config_subst_models.yaml` | Substitution model (TVM) | Yule | TVM+Γ4 | Strict | 1 | One model per run: 4 transversion + 1 transition rate |
| `config_subst_jc69.yaml` | Substitution model (JC69) | Yule | JC69+Γ4 | Strict | 1 | Equal rates/frequencies baseline model |
| `config_subst_sym.yaml` | Substitution model (SYM) | Birth-Death | SYM+Γ4 | UCLN | 2 | Six free exchangeabilities, `frequencies.mode: uniform` |
| `config_subst_tim.yaml` | Substitution model (TIM) | Yule | TIM+Γ4 | Strict | 1 | `rateAG`/`rateCT` + two transversion rates, `frequencies.mode: empirical` |
| `config_aa_model.yaml` | Amino acid model | Birth-Death | WAG+Γ4 | UCLN | 2 | Amino acid data + Laplace/LogNormal calibrations |
| `config_aa_jtt.yaml` | Amino acid model (JTT) | Yule | JTT+Γ4+I | Strict | 2 | Gamma+I on protein data |
| `config_tip_dates.yaml` | Tip dates | Birth-Death | HKY | UCLN | 1 | `TraitSet` tip dates (all 12 taxa dated) |
| `config_calibrated_yule.yaml` | CalibratedYule | CalibratedYule | HKY+Γ4 | Strict | 3 | `calibration_method` switching + stem calibration |
| `config_multi_partition.yaml` | Multi-partition | Birth-Death | HKY+Γ4 / GTR+Γ4 | UCLN + linked | 2 | `FilteredAlignment` codon partitioning + clock linking, one linked topology |
| `config_relaxed_clock.yaml` | Relaxed lognormal clock | Birth-Death | GTR+Γ4 | UCLN | 2 | Bounded `ucld_mean` / `ucld_stdev` |
| `config_relaxed_exponential.yaml` | Relaxed exponential clock | Yule | GTR+Γ4 | UCE | 2 | UCE relaxed clock |
| `config_rlc.yaml` | Random local clock | Yule | HKY+Γ4 | RLC | 1 | Indicator + rate state nodes |
| `config_coalescent_exponential.yaml` | Exponential-growth coalescent | Coalescent (exponential) | HKY+Γ4 | Strict | 1 | Signed `growth_rate` |
| `config_skyline.yaml` | Bayesian Skyline | Bayesian Skyline | TN93+Γ4 | Strict | 1 | Skyline population history |
| `config_ebsp.yaml` | EBSP | EBSP | HKY+Γ4 | Strict | 1 | Extended Bayesian Skyline Plot (BEAST 2.7 core package) |
| `config_bd_skyline.yaml` | BD Skyline Serial | BD Skyline Serial | HKY+Γ4 | Strict | 1 | Requires the bdsky add-on (verified against BDSKY 1.5.1) |
| `config_nested_sampling.yaml` | Nested Sampling | Yule | HKY+Γ4 | Strict | 1 | NS MCMC (requires the NS add-on, verified against NS 1.2.0) |

The counts in the table are the numbers of `calibrations:` entries actually present in each file, and the model strings are what the files declare per partition (`/` separates partitions).

In addition, `calibrations.yaml` is a standalone calibration-point file for the `quick`
subcommand. The example data are `primates.fasta` (12 taxa × 898 bp) and `aminoacid.fasta`
(10 taxa × 234 aa), and the example output artifacts are `output_basic.xml`,
`output_basic.fingerprint.json` and `output_basic.methods.tex`.

**Running the examples:**

```bash
# Run the basic example
beast2py generate -c examples/config_basic.yaml -o output_basic.xml -v

# Run the advanced example (with diagnostics)
beast2py generate -c examples/config_advanced.yaml -o output_advanced.xml --diagnose --fingerprint --methods -v

# Run the multi-partition example
beast2py generate -c examples/config_multi_partition.yaml -o output_mpart.xml -v

# Quick mode example
beast2py quick -a examples/primates.fasta -o output_quick.xml \
    --tree-prior yule --subst-model hky --clock-model strict \
    --calibration-yaml examples/calibrations.yaml
```

```python
# Run the examples via the API
from beast2py.api import Beast2Py

b2p = Beast2Py(base_dir="examples/")

# Basic example
xml = b2p.generate_xml("config_basic.yaml", output="output_basic.xml")

# Advanced example (with all features)
result = b2p.generate(
    "config_advanced.yaml",
    output="output_advanced.xml",
    diagnose=True,
    fingerprint=True,
    methods=True,
)
```

---

## 22. Troubleshooting

### 22.1 Strict validation, rejected inputs and their error messages

Beast2Py is **fail-closed**: when the generator cannot honour a configuration, it raises `ConfigError` *before* writing any file (CLI: `[ERROR] Configuration error: …`, exit 1; Python API: an exception). The messages below are the actual strings the current code emits.

**Schema and unknown keys.** The parser checks every section against a key whitelist (`config.py::KNOWN_KEYS`). It refuses a key outside that section's schema instead of silently defaulting it:

```text
mcmc: unknown key(s) chainlength_typo. Allowed keys for this section: chain_length,
particle_count, pre_burnin, sample_from_prior, seed, store_every, sub_chain_length,
type. Beast2Py refuses unknown keys rather than silently applying a default, because
a typo in a model or MCMC setting changes the analysis without changing the XML's
validity.
```

**Booleans.** The parser accepts only a real YAML boolean or one of `true/false/yes/no/on/off/1/0`, and rejects anything ambiguous. A quoted `"false"` is still accepted and means *false*:

```text
'alignment.hky.kappa': 'maybe' is not a boolean; use true/false (unquoted in YAML) or
one of ['0', '1', 'false', 'no', 'off', 'on', 'true', 'yes']
```

**Gamma shape without categories** (aligns with `SiteModel.java:64`, "Ignored if gammaCategoryCount 1 or less"):

```text
Partition 'alignment': gamma_shape was supplied but gamma_categories is 1. SiteModel.java
documents the shape input as 'Ignored if gammaCategoryCount 1 or less', so an estimated
shape would be a pure prior random walk that contributes nothing to the likelihood while
still appearing as a column in the trace.
```

**Fixed clock with no calibration** (absolute time is unidentifiable):

```text
The clock rate is fixed and no calibration is supplied, so absolute time is not
identifiable: node heights can only be estimated in expected-substitutions-per-site
units. Either estimate the clock rate or add at least one absolute (e.g. fossil)
calibration.
```

**Unlinked and per-partition trees** (`tree: separate` or a partition id):

```text
Partition 'alignment': tree: 'separate' requests a tree other than the shared one, but
Beast2Py assembles every partition on a single linked topology (per-partition /
unlinked trees are not implemented). ... Remove the 'tree' key to analyse the shared
topology, or run one analysis per independent partition.
```

**Tip-date coverage** (every taxon must be dated; unknown/missing taxa refused):

```text
tip_dates.dates covers 2 of 12 taxa; these are missing: Gorilla, Homo_sapiens, ...
BEAST2 treats an undated tip as sampled at time zero, which silently turns those
lineages into contemporaneous samples. Supply a date for every tip, or disable tip_dates.
```

**MCMC length rules** (`chain_length > 0`, `pre_burnin >= 0` and `< chain_length`, `store_every` positive or `-1`):

```text
mcmc.pre_burnin (5000) must be smaller than chain_length (100); otherwise the chain
never samples.
```

**Filter grammar and range** (`data_id` required for a filtered alignment, 1-based ranges checked against the source length):

```text
Alignment 'codon3': 'filter' requires 'data_id' naming the source alignment; without it
the filter would reference itself
```

**Frequency vector must be proper** (length matches dimension, entries > 0, sums to 1 within 2 % then renormalised):

```text
alignment.frequencies: frequencies sum to 1.2, not 1.0. Supply a proper probability
vector (a rounding tolerance of 0.02 is accepted and the vector is then renormalised).
```

**Inverted calibration bound:**

```text
Calibration 'root': uniform lower (40.0) must be strictly below upper (20.0); as given,
the density is zero on the whole support.
```

**Birth–death domain** (`death_rate < birth_rate`; use BDSKY for `death >= birth`):

```text
tree_prior: death_rate (0.8) must be smaller than birth_rate (0.5) for the
BirthDeathGernhard08Model, whose parameters are birth - death (>= 0) and death / birth
(< 1). A super-critical process (death >= birth) needs the BDSKY add-on
(tree_prior.type: bd_skyline_serial), which parameterises birth, death and sampling
rates separately.
```

**Unresolvable `parameter_priors` name** (a prior that matches no estimated node is an error, not a silent no-op):

```text
parameter_priors: 'totally_bogus_name' does not match any estimated parameter of this
analysis. Estimated state nodes are: alignment.clockRate, alignment.freqParameter,
alignment.hky.kappa, birthRate. A prior whose target does not exist is silently never
evaluated.
```

### 22.2 Common errors

**`Configuration file not found`**

Make sure the YAML file path is correct, either absolute or relative to the configuration file's own directory. Beast2Py resolves sequence file paths relative to the directory that contains the YAML configuration.

**`Calibration taxon not found in any alignment`**

Check that the taxon names in the calibrations match the sequence headers exactly. A FASTA header takes the first word after `>` as the taxon name, and the match is case-sensitive.

**`Partition ID not matching alignment`**

Make sure the partition ID matches an alignment ID. If you use `FilteredAlignment` (codon partitioning), the partition ID should match the filtered alignment's ID, not the base alignment's.

**`Duplicate partition ID`**

Every partition must have a unique ID. Check the YAML configuration for duplicate partition IDs.

**`Invalid substitution model`, `Invalid clock model`, `Invalid tree prior type`**

Use `beast2py list-models` to see the supported model types and make sure the names used in the configuration are correct.

**`Namespace does not contain 'beast.base'`**

This warning indicates that the generated XML uses the BEAST2 v2.6.x namespace format. Beast2Py generates the v2.7.x format. If you use BEAST2 v2.6.x, upgrade to v2.7.0 or newer.

**`BEAST2 executable not found`**

Make sure BEAST2 is installed and the `beast` command is on `PATH`, or specify the full path with `--beast2-path`.

**`Could not find package: bdsky / NS` (missing add-on package)**

The BD Skyline Serial and Nested Sampling examples require the [bdsky](https://github.com/BEAST2-Dev/bdsky) and [NS](https://github.com/BEAST2-Dev/nested-sampling) add-on packages. Install them with the BEAST2 Package Manager GUI, or from the command line (where `launcher.jar` is the one shipped with BEAST2 in its `lib` directory):

```bash
java -cp launcher.jar beast.util.PackageManager -add bdsky
java -cp launcher.jar beast.util.PackageManager -add NS
```

Confirm the installed versions with the Package Manager or `beast -version`; the R package [mauricer](https://github.com/ropensci/mauricer) provides equivalent helpers (`is_beast2_pkg_installed()`, `get_beast2_pkg_names()`).

### 22.3 Common diagnostics issues

**Temporal consistency errors**

The parent clade is calibrated younger than the child clade. Check the distribution parameters of the nested calibration points and make sure the parent clade's HPD interval covers an older time range than the child clade's.

**Distribution overlap warnings**

By the default (relative) criterion, the detector takes the 1-Wasserstein distance between two non-nested calibrations, divides it by the mean width of their 95% intervals, and warns when that ratio falls below `relative_overlap_threshold` (default 0.15). In other words, the two priors are nearly interchangeable relative to their own widths. Check whether these calibrations provide redundant information, and whether they could be merged. See [Section 14.9](#149-the-diagnostics-block) to switch to (or add) the absolute `overlap_threshold` criterion.

**Monophyly conflicts**

Two taxon sets that intersect but are not nested both require monophyly. Make sure the calibration taxon sets are either completely disjoint or completely nested, or relax the monophyly constraint of one calibration (set `monophyletic: false`).

### 22.4 API usage issues

**`ImportError: No module named 'beast2py'`**

Make sure you installed the package with `pip install -e .`, and that the current working directory or the Python path contains the project directory.

**`FileNotFoundError` during API calls**

If the configuration uses relative sequence file paths, make sure the `base_dir` parameter points to the correct directory: `Beast2Py(base_dir="examples/")`.

### 22.5 Getting help

```bash
# General help
beast2py --help

# Command-specific help
beast2py generate --help
beast2py diagnose --help
beast2py quick --help
beast2py validate --help
beast2py methods --help
beast2py pipeline --help
beast2py fingerprint --help
beast2py list-models --help
```

```python
# Python help
import beast2py
help(beast2py.api)
help(beast2py.api.Beast2Py)
```

### 22.6 Reporting issues

Please report bugs and feature requests at https://github.com/ZengZichao/Beast2Py/issues. When submitting an issue, please include: the Beast2Py version (`beast2py --version`), the Python version, the operating system, the configuration file contents (if applicable), and the full error message.

---

## 23. Appendix: Complete Model List

### 23.1 Substitution models (13)

**Nucleotide models (7):** `jc69`, `hky`, `tn93`, `gtr`, `sym`, `tim`, `tvm`

**Amino acid models (6):** `wag`, `jtt`, `dayhoff`, `blosum62`, `cprev`, `mtrev`

### 23.2 Clock models (4)

`strict`, `ucln`, `uce`, `rlc`

### 23.3 Tree priors (8)

`yule`, `calibrated_yule`, `birth_death`, `coalescent_constant`, `coalescent_exponential`, `bayesian_skyline`, `ebsp`, `bd_skyline_serial`

### 23.4 Calibration distributions (11)

`normal`, `lognormal`, `uniform`, `exponential`, `gamma`, `beta`, `laplace`, `inverse_gamma`, `one_on_x`, `poisson`, `chi_square`

### 23.5 Data types (2)

`nucleotide`, `aminoacid`

### 23.6 MCMC types (2)

`standard`, `nested_sampling`

### 23.7 Initialization strategies (3)

`random`, `upgma`, `newick`

### 23.8 Calibration methods (2)

`mrca_prior`, `calibrated_yule`

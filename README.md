# Beast2Py

**A Python framework for reproducible divergence time estimation with automated calibration prior specification, validation, and diagnostics.**
[中文文档](README_zh.md) | [English README](README_en.md)

[![License: MIT AND LGPL-2.1-only](https://img.shields.io/badge/License-MIT%20AND%20LGPL--2.1--only-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![BEAST2 2.7.8](https://img.shields.io/badge/BEAST2-2.7.8-orange.svg)](https://www.beast2.org/)

---

## Overview

Beast2Py is a Python framework for **standardizing, validating, and diagnosing node calibration priors in BEAST2 divergence time estimation**, and generating reproducible BEAST2 analysis configurations. It provides both a command-line interface (CLI) and a Python application programming interface (API), enabling researchers to generate complete BEAST2 XML files from human-readable YAML configurations, automatically detect calibration conflicts, perform prior sensitivity analysis, and produce reproducibility artifacts.

### Core Innovations

> Comparison baseline: BEAUti2 v2.7.8 (BEAST.app package) and beautier v2.6.15.

| Feature | Description |
|---------|-------------|
| **Calibration conflict detection** | Automatic detection of temporal inconsistencies, distribution overlaps, and monophyly conflicts between calibration points |
| **Prior sensitivity analysis** | Automatic generation of "prior-only sampling" XMLs to evaluate the marginal impact of calibration priors |
| **Calibration provenance metadata** | Each calibration point records fossil source, reference DOI, and calibration type (hard/soft bound) |
| **Reproducibility framework** | YAML configuration versioning + analysis fingerprint + LaTeX methods description + Snakemake/Nextflow pipeline generation |

### Three-Tier Feature Set

- **Tier 1 (beautier parity):** All features supported by the R package beautier
- **Tier 2 (beautier unsupported):** all 13 entries of beautier's "Missing features/unsupported" list, counted down to the second level (the list has 10 top-level bullets; the *Clock models* bullet contributes its 2 sub-items and the *Tree priors* bullet its 3, so 8 + 2 + 3 = 13): distribution offset, two or more DNA alignments, two or more site/clock/tree models, two or more MRCA priors, shared site/clock/tree models, amino-acid alignments, hyper parameters, relaxed exponential clock, random local clock, Calibrated Yule, coalescent EBSP, Birth Death Skyline Serial, and initialization. Beyond those 13, Beast2Py also adds SYM, TIM and TVM nucleotide models
- **Tier 3 (unique integration):** seven diagnostic and reproducibility capabilities uniquely integrated at the configuration stage. A few of these exist in isolated form elsewhere (e.g., BEAUti2 v2.7+ ships an interactive methods-section viewer and a prior-distribution preview, and `beast -sampleFromPrior` runs prior-only sampling on any existing XML), but no other tool offers them as a single, scriptable, configuration-stage workflow

## Supported Features

- **Sequence data:** FASTA / NEXUS import (nucleotide and amino acid)
- **Substitution models:** JC69, HKY, TN93, GTR, SYM, TIM, TVM + amino acid models (WAG, JTT, Dayhoff, Blosum62, CPREV, MTREV)
- **Site models:** Gamma rate heterogeneity, Invariant sites, Gamma+I
- **Clock models:** Strict clock, UC Relaxed LogNormal (UCLN), UC Relaxed Exponential (UCE), Random Local Clock (RLC)
- **Tree priors:** Yule, CalibratedYule, Birth-Death, Coalescent (Constant/Exponential), Bayesian Skyline, EBSP; BD Skyline Serial (requires the [bdsky](https://github.com/BEAST2-Dev/bdsky) add-on package)
- **Calibration distributions:** Normal, LogNormal, Uniform, Exponential, Gamma, Beta, Laplace, InverseGamma, OneOnX, Poisson, ChiSquare
- **Multi-partition:** Several alignments (or codon-position `FilteredAlignment` subsets of one alignment) with independent or shared site and clock models (`linked_to`). All partitions are assembled on **one linked topology**: the parser refuses `tree: separate` and any per-partition tree id outright, and says why, because unlinked gene trees (`` `*BEAST` `` / StarDivergence-style analyses) are not implemented
- **Multi MRCA priors:** Unlimited calibration points
- **Tip dates:** A `TraitSet` supports serially sampled data (`TipDatesRandomWalker` is emitted as soon as tip dates are enabled, without needing a `tipsonly` calibration)
- **MCMC:** Standard MCMC + Nested Sampling MCMC (nested sampling requires the [NS](https://github.com/BEAST2-Dev/nested-sampling) add-on package)
- **Initialization:** RandomTree, UPGMA, User Newick
- **Hyperpriors:** Priors on the parameters of a calibration distribution
- **Strict configuration validation:** the parser rejects unknown keys, quoted booleans, ragged or duplicated alignments, mislabelled data types, out-of-domain priors and unidentifiable time scales at parse time, instead of silently applying a default

## Installation

### Prerequisites

- Python 3.10 or higher
- BEAST2 v2.7.8, verified against this version (the generated XML targets the BEAST 2.7 `beast.base.*` namespace; other 2.7.x releases are expected to work but have not been tested). BEAST 2.7.x class files require **JDK 17 or newer**
- JDK 17 or newer, only for the optional BEAST2 validation gate (`--beast2-validate` / `validate --beast2`); everything else is pure Python

### Install from source

```bash
git clone https://github.com/ZengZichao/Beast2Py.git
cd Beast2Py
pip install -e .
```

### Dependencies

```
pyyaml >= 6.0
biopython >= 1.80
numpy >= 1.20
scipy >= 1.8
matplotlib >= 3.5
defusedxml >= 0.7
```

## Quick Start

### CLI Usage

Every command below has been run from a clone of the repository root. After
`pip install -e .` they are available as `beast2py …`; without installing, replace
`beast2py` with `python3 -m beast2py.main` (identical behaviour).

```bash
# Generate BEAST2 XML from a YAML configuration
beast2py generate --config examples/config_basic.yaml --output output_basic.xml -v

# Run calibration diagnostics (conflict detection + sensitivity + visualization)
beast2py diagnose --config examples/config_basic.yaml --report report.html --sensitivity --visualize

# Quick mode for a simple single-partition analysis
beast2py quick \
    --alignment examples/primates.fasta \
    --tree-prior yule \
    --subst-model hky \
    --clock-model strict \
    --calibration-yaml examples/calibrations.yaml \
    --chain-length 1000000 \
    --name quick_analysis \
    --output output_quick.xml

# Validate generated XML (structural checks; add --beast2 for the BEAST2 gate)
beast2py validate --xml examples/output_basic.xml

# Generate LaTeX methods description
beast2py methods --config examples/config_basic.yaml --output methods.tex

# Generate analysis fingerprint
beast2py fingerprint --config examples/config_basic.yaml --output fingerprint.json

# Generate Snakemake/Nextflow pipeline
beast2py pipeline --config examples/config_basic.yaml --type snakemake --output-dir pipeline

# List all supported models
beast2py list-models
```

`generate` and `quick` refuse to overwrite an existing `--output` unless `--force`
is given, and they write nothing unless every validation gate that ran has passed
(see [Validation](#validation)).

### API Usage

```python
from beast2py.api import Beast2Py

# Create an API instance
b2p = Beast2Py()

# Generate XML from a YAML config file. When `output` is given, the same
# validation gates as the CLI run first, and the file is written only if they
# pass; an existing file is refused unless force=True.
xml_str = b2p.generate_xml("examples/config_basic.yaml", output="output_basic.xml")

# Run calibration diagnostics
report = b2p.diagnose("examples/config_basic.yaml", output_dir="diagnostics/", report_path="report.html")
for c in report.conflicts:
    print(f"[{c.severity}] {c.description}")

# Validate generated XML
result = b2p.validate_xml("output_basic.xml")
print(f"Valid: {result.is_valid}")

# Generate LaTeX methods description
methods = b2p.generate_methods("examples/config_basic.yaml", output="methods.tex")

# Generate analysis fingerprint
fingerprint = b2p.generate_fingerprint("examples/config_basic.yaml")
print(f"Analysis fingerprint: {fingerprint}")

# List all supported models
models = b2p.list_models()
print(models["substitution_models"])
```

`generate_xml()`, `quick_generate()` and the all-in-one `generate()` all take
`force=False`; a failed gate raises `ValueError` and leaves the output untouched.

## Documentation

Comprehensive documentation is available in the `docs/` directory:

- **[User Manual](docs/user_manual_en.md)** — Complete handbook covering installation, the CLI/API, the configuration reference, diagnostics, and troubleshooting
- **[Architecture Guide](docs/architecture_en.md)** — Detailed description of the layered architecture and module design
- **[Tutorial](docs/tutorial_en.md)** — Step-by-step tutorial with example configurations
- **[Calibration Guide](docs/calibration_guide_en.md)** — Guide to specifying calibration priors with provenance metadata
- **[XML Format Guide](docs/xml_format_en.md)** — Description of the generated BEAST2 XML format

Each manual also carries a **strict validation** section listing every input the parser now rejects and the error message it produces.

## Example Configurations

The `examples/` directory holds **19** complete YAML configurations.
`tests/test_beast2_validation.py` generates XML for all 19 and checks each one
against the real BEAST 2.7.8 `XMLParser`, which parses *and* initialises the model.
Their actual coverage is:

- **Substitution models:** 9 of the 13 supported models are exemplified — all seven nucleotide models (JC69, HKY, TN93, GTR, SYM, TIM, TVM) plus two of the six amino-acid matrices (WAG, JTT). `dayhoff`, `blosum62`, `cprev` and `mtrev` are supported but have no example
- **Clock models:** all four (strict, UCLN, UCE, RLC)
- **Tree priors:** 7 of the 8 supported priors — `coalescent_constant` has no example of its own
- **Calibration distributions:** 3 of the 11 supported distributions appear (lognormal, normal, uniform)

| File | Scenario | Tree prior | Substitution model | Clock model | Calibrations | Special features |
|---|---|---|---|---|---|---|
| `config_basic.yaml` | Basic analysis | Yule | HKY+Γ4 | Strict | 1 | explicit `parameter_priors`, provenance; source of `output_basic.*` |
| `config_advanced.yaml` | Advanced multi-calibration | Birth-Death | GTR+Γ4+I / HKY | UCLN + linked | 3 | codon partitions + site/clock linking + hyperpriors + `diagnostics` block + stem calibration |
| `config_subst_jc69.yaml` | Model comparison, JC69 | Yule | JC69+Γ4 | Strict | 1 | one analysis per candidate model |
| `config_subst_sym.yaml` | Model comparison, SYM | Birth-Death | SYM+Γ4 | UCLN | 2 | `rates` vector + `frequencies.mode: uniform` |
| `config_subst_tim.yaml` | Model comparison, TIM | Yule | TIM+Γ4 | Strict | 1 | named rates + `frequencies.mode: empirical` |
| `config_subst_models.yaml` | Model comparison, TVM | Yule | TVM+Γ4 | Strict | 1 | five named rate parameters |
| `config_aa_model.yaml` | Amino acid model | Birth-Death | WAG+Γ4 | UCLN | 2 | `aminoacid.fasta` + offset lognormal calibrations |
| `config_aa_jtt.yaml` | Amino acid model | Yule | JTT+Γ4+I | Strict | 2 | Gamma+I on protein data |
| `config_tip_dates.yaml` | Tip dates | Birth-Death | HKY | UCLN | 1 | `TraitSet`, `date-backward`, one date per taxon |
| `config_calibrated_yule.yaml` | CalibratedYule | CalibratedYule | HKY+Γ4 | Strict | 3 | `calibration_method` switching + stem calibration |
| `config_multi_partition.yaml` | Multi-partition | Birth-Death | HKY+Γ4 / GTR+Γ4 | UCLN + linked | 2 | `FilteredAlignment` codon partitioning, clock linking, one linked topology |
| `config_relaxed_clock.yaml` | Relaxed lognormal clock | Birth-Death | GTR+Γ4 | UCLN | 2 | bounded `ucld_mean` / `ucld_stdev` |
| `config_relaxed_exponential.yaml` | Relaxed exponential clock | Yule | GTR+Γ4 | UCE | 2 | UCE relaxed clock |
| `config_rlc.yaml` | Random local clock | Yule | HKY+Γ4 | RLC | 1 | indicator + rate state nodes |
| `config_skyline.yaml` | Bayesian Skyline | Bayesian Skyline | TN93+Γ4 | Strict | 1 | skyline population function |
| `config_coalescent_exponential.yaml` | Coalescent, exponential growth | Coalescent (exponential) | HKY+Γ4 | Strict | 1 | signed `growth_rate` |
| `config_ebsp.yaml` | EBSP | EBSP | HKY+Γ4 | Strict | 1 | Extended Bayesian Skyline Plot (BEAST 2.7 core package) |
| `config_bd_skyline.yaml` | BD Skyline Serial | BD Skyline Serial | HKY+Γ4 | Strict | 1 | requires the [bdsky](https://github.com/BEAST2-Dev/bdsky) add-on package |
| `config_nested_sampling.yaml` | Nested Sampling | Yule | HKY+Γ4 | Strict | 1 | NS MCMC (requires the [NS](https://github.com/BEAST2-Dev/nested-sampling) add-on package) |

In addition, `calibrations.yaml` provides a calibration-only YAML for the `quick` subcommand, and the example alignments `primates.fasta` (12 taxa × 898 bp, converted from the BEAST2 example alignment `Primates.nex`) and `aminoacid.fasta` (10 taxa × 234 aa) are included.

## Repository Structure

- `beast2py/` — Source code of the Beast2Py Python package (CLI, API, configuration parsing, XML generation, diagnostics, reproducibility, and validation modules)
- `examples/` — 19 example configurations, example alignments (FASTA), and example output artifacts
- `tests/` — Test suite (352 collected tests: 331 unit, semantic and release-integrity tests + 21 BEAST2 integration tests)
- `beast2py/tools/` — Headless BEAST2 validation helper bundled *inside* the package (`Beast2Validator.java` source plus the precompiled `classes/` and `launcher.jar`), shipped as `package-data` by `pyproject.toml`; the repository root keeps only a thin `beast2_validate.sh` wrapper
- `docs/` — User manual, architecture guide, tutorial, calibration guide, and XML format guide (Chinese and English)
- `scripts/check_figure_export.py` — Post-export gate for the three draw.io figures, which have no Python source to redraw them: it reads the delivered PDF and PNG back and fails on an export older than its `.drawio` source, a span under 7.92 pt, a font outside Helvetica, a rule under 0.71 pt, or a raster under 600 ppi at the page's own width (`pip install -e .[figures]`)
- `scripts/check_submission_placeholders.py` — Pre-upload gate listing every unresolved archive identifier (Zenodo DOI, Dryad DOI and reviewer URL, TreeBASE accession) still present in the manuscript or cover letter, each with its paragraph or page; exits non-zero while any remain (`pip install -e .[figures]` enables the PDF pass)
- `beast2_validate.sh` — Launcher script for headless BEAST2 validation (needs `BEAST.base.jar` and JDK 17 or newer; it auto-detects `~/.beast/2.7/BEAST.base/*/lib/BEAST.base.jar`, or set `BEAST2_JAR`)
- `pyproject.toml` — Package build and dependency configuration

## Validation

`generate` and `quick` run three gates in sequence. If any gate that ran reports a
failure, the command exits non-zero and **writes no output file**. An existing output file
needs `--force`, otherwise it is refused, and the CLI warns when a sibling
`.fingerprint.json` / `.methods.tex` would be left stale. The three gates are:

1. **Structural checks** — XML well-formedness and required elements, id/idref consistency,
   duplicate IDs, initial values inside their bounds, sequence sanity; plus confirmation
   that every estimated parameter carries a `Prior`, that every estimated state node carries
   an operator, and that no operator touches a fixed parameter
2. **Calibration conflict detection** — temporal consistency between nested
   clades (root calibrations included), distribution overlap, monophyly conflicts
3. **BEAST2 check (optional, `--beast2-validate`)** — the headless BEAST2 `XMLParser` calls
   `initAndValidate()` on every object after parsing, so what is checked is whether the
   model can **complete initialisation**, not merely whether it parses

**What validation does not prove:** none of the three gates can tell you whether
the model answers your biological question. A model that BEAST2 happily
initialises may still be the wrong model.

Gate 3 only runs where BEAST2 and JDK 17 are installed. If you ask for
`--beast2-validate` and BEAST2 cannot be found, generation stops with exit code 2 and
writes nothing, so a CI job that reads only the exit code never treats "BEAST2 missing" as
a pass. `--allow-unvalidated` restores the write-anyway behaviour, in which case the XML
must not be described as BEAST2-verified. Everything has been verified against
BEAST2 **v2.7.8** only.

The test suite comprises 352 collected tests (331 unit/semantic/release + 21 BEAST2
integration). The integration tests generate XML from all 19 example
configurations plus the `quick` path and the committed `output_basic.xml`, and
check each with the real BEAST 2.7.8 parser and model initialisation. When the
BEAST2 jars or JDK 17 are missing, those tests simply skip instead of failing.

## Reproducibility

Beast2Py generates the following reproducibility artifacts:

- **Analysis fingerprint** — deterministic identifier `B2P-{sha256(config)[:12]}-{version}` (e.g. `B2P-ba18c26ed61f-0.1.0`). It contains **no date component**, so the same configuration plus the same alignment data yields the same bytes and the same identifier in any time zone and at any rerun time. Renaming the output file does not change the fingerprint; changing the sequence content, the tree setting, a crown/stem flag, a hyperprior, a parameter prior, a tip date or an operator weight does
- **Data digest** — besides the identifier, the XML comment carries a separate digest: 16 hex characters over the alignment *content*. The wall-clock generation time appears only in the `.fingerprint.json` sidecar file, never in the XML
- **LaTeX methods description** — Publication-ready methods paragraph
- **Snakemake/Nextflow pipeline** — End-to-end reproducible workflow files
- **Calibration provenance** — Fossil source, reference DOI, and calibration type for each calibration point

## Roadmap

- Publish to PyPI (the `beast2py` package name has been verified as available)
- Support CladeAge-style calibrations (birth-death simulation-based calibration densities mapped to `cladeage.math.distributions.FossilPrior` XML)
- Automatic ESS/HPD parsing of prior-only logs with threshold alerts (Tracer-style ACT/ESS algorithms)
- BEAST2 run orchestration (pybeast-style isolated run directories, seed management, resume snapshots)
- Add automatic BEAST2 version checking
- BEAST2 v2.6.x namespace compatibility
- pirouette-style end-to-end validation (known tree → simulated alignments → injected conflicting calibrations → detector sensitivity check)
- Track the BEAST 3 / LPhy ecosystem transition and adapt accordingly

## Related Tools

Adjacent open-source tools in the BEAST2 ecosystem (see also the BEAST2 blog post "10 ways to generate BEAST XML"):

| Tool | Language | Scope |
|---|---|---|
| BEAUti2 (bundled with BEAST2) | Java/JavaFX | Graphical XML editor; v2.7+ includes an interactive methods-section viewer and a prior-distribution preview |
| [beautier](https://github.com/ropensci/beautier) / [babette](https://github.com/ropensci/babette) | R | Programmatic XML generation; the babette suite (beastier runner, mauricer package manager, tracerer ESS/HPD diagnostics) provides an end-to-end R workflow |
| [BEASTling](https://github.com/lmaurits/BEASTling) | Python | INI config → BEAST2 XML (linguistic/CLDF data; includes calibration-vs-monophyly constraint checks) |
| [beast2-xml (acorg)](https://github.com/acorg/beast2-xml) | Python | Command-line XML generator built on BEAUti templates (FASTA/FASTQ input) |
| [BEASTmasteR](https://github.com/nmatzke/BEASTmasteR) | R | NEXUS + Excel workbook → XML for fossil tip-dating / FBD |
| [LPhy](https://github.com/LinguaPhylo/linguaPhylo) / [LPhyBeast](https://github.com/LinguaPhylo/LPhyBeast) | Java | Probabilistic model specification language → BEAST XML (already targeting BEAST 3) |
| [pybeast](https://github.com/Wytamma/pybeast) / [beastiary](https://github.com/Wytamma/beastiary) | Python | BEAST2 run wrapper (run directories/seeds/resume) and live MCMC trace monitoring (ESS/HPD) |

Beast2Py's differentiator is that it integrates the following into a single configuration-stage pipeline: declarative YAML configuration, three-way calibration conflict detection (temporal consistency, distribution overlap, monophyly), automated leave-one-out prior-sensitivity XML generation, and calibration provenance (DOI), analysis fingerprint, methods description and pipeline generation.

## Third-party licences

The compiled classes shipped under `beast2py/tools/classes/beast/pkgmgmt/` come from
[BEAST 2](https://github.com/CompEvol/beast2) v2.7.8 and are licensed under the **GNU LGPL v2.1**; the rest of
this package is MIT. The licence text is in `beast2py/tools/LICENSE.BEAST2-LGPL-2.1.txt`, and
`beast2py/tools/NOTICE` records the provenance, whether anything was modified, and how to obtain the full
corresponding source. `beast2py/tools/Beast2Validator.java` ships its source with the package, as LGPL-2.1
requires for that compilation unit.

## License

The package is MIT; the bundled BEAST2 package-management classes are LGPL-2.1 (see `beast2py/tools/NOTICE`). The SPDX identifier is `MIT AND LGPL-2.1-only`. See [LICENSE](LICENSE) for details.

## Citation

If you use Beast2Py in your research, please cite:

> Beast2Py: A Python framework for reproducible divergence time estimation with automated calibration prior specification, validation, and diagnostics. (manuscript in preparation)

## Authors

Author information will be added upon official release.

## Funding

Funding information will be added upon official release.

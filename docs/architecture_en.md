# Architecture

## Overview

Beast2Py follows a layered architecture designed for extensibility and maintainability. The system is divided into five layers: CLI entry, configuration management, domain modules, data models, and utility support.

## Module Responsibilities

| Module | File | Responsibility |
|------|------|------|
| CLI Entry | `main.py` | Command-line argument parsing, the three validation gates, workflow orchestration |
| Python API | `api.py` | `Beast2Py` facade — the same pipeline as the CLI, callable from Python |
| Configuration | `config.py` | YAML config parsing, strict key/boolean/domain validation, default filling |
| Data Models | `models.py` | All dataclass definitions |
| Sequence Reading | `sequence.py` | FASTA/NEXUS file reading via biopython |
| Partition Management | `partition.py` | Multi-partition model link/unlink management |
| Calibration Assembly | `calibration.py` | MRCAPrior, TaxonSet, distribution XML building |
| Model Assembly | `model.py` | Tree prior, substitution model, clock model, operators, state, init XML building |
| XML Generation | `xml_writer.py` | Complete XML assembly using ElementTree, operator/state-node invariants |
| Calibration Diagnostics | `diagnostics.py` | Conflict detection, sensitivity analysis, visualization (core innovation) |
| Reproducibility | `reproducibility.py` | Analysis fingerprint, methods description, pipeline generation (core innovation) |
| Validator | `validator.py` | XML structural and semantic validation, id/idref consistency, BEAST2 gate |
| Utility Functions | `utils.py` | ID generation, XML serialization, distribution calculations |
| Model Registry | `registry.py` | Model name to XML spec mapping, extensible registration |
| BEAST2 Validation Launcher | `beast2_validate.sh` | Headless BEAST2 validation entry point shipped *inside* the package (the repository root keeps only a thin wrapper of the same name) |
| Java Validation Helper | `tools/` | `Beast2Validator.java` plus the precompiled `classes/` and `launcher.jar`, shipped as package data |

## Data Flow

```
1. User Input
   ├── FASTA/NEXUS sequence files (multiple supported)
   └── YAML configuration file (with calibrations + provenance)
        │
        ▼
2. Configuration Parsing (config.py) — strict, fail-closed, before anything is built
   ├── Sequence reading (sequence.py) → List[Alignment]
   ├── Partition management (partition.py) → List[Partition] (alignment + model config + link/unlink)
   └── Rejected here: unknown keys, ambiguous booleans, out-of-domain priors, a bad
       FilteredAlignment expression, `tree: separate`, and a provenance range that a
       *hard* bound contradicts (a *soft* bound only warns)
   → BEASTConfig object (partitions + calibrations + MCMC + diagnostics + provenance)
        │
        ▼
3. Model Assembly (model.py)
   → Tree prior XML + Substitution model XML + Clock model XML per partition,
     one shared linked topology; every `estimate: false` parameter is written
     *inline* instead of entering the state space
        │
        ▼
4. Calibration Assembly (calibration.py)
   → MRCAPrior XML elements + TaxonSet elements (+ TraitSet for tip dates)
        │
        ▼
5. XML Generation (xml_writer.py) — in memory only
   → Complete XML using xml.etree.ElementTree; the writer also asserts that every
     estimated state node has an operator and that no operator touches a fixed one
        │
        ▼
6. Validation Gates (main.py::_validate_and_write) — nothing is written unless every
   gate that ran passes; any failure exits non-zero
   ├── Gate 1 structural + semantic checks (validator.py)
   ├── Gate 2 calibration conflict detection (diagnostics.py) ← Core Innovation
   │     temporal consistency, distribution overlap, monophyly conflicts
   └── Gate 3 optional BEAST2 parse + initAndValidate (beast2_validate.sh + tools/)
        │
        ▼
7. Output (only after the gates pass)
   ├── BEAST2 XML file
   └── Reproducibility packaging (reproducibility.py) ← Core Innovation
       ├── Analysis fingerprint (config hash + data hash + tool version, date-free)
       ├── Methods description (LaTeX fragment)
       └── Snakemake/Nextflow pipeline files (separate `pipeline` command)

Standalone commands on the same building blocks:
   diagnose (diagnostics.py) → conflict report (HTML/text), N+1 prior-only sampling
   XMLs, calibration distribution density plots
   methods / fingerprint / pipeline → the reproducibility artifacts on their own
```

## Model Registry Pattern

The `ModelRegistry` class in `registry.py` implements the Registry Pattern, mapping the model names given in the configuration to the corresponding BEAST2 XML spec values. This brings three benefits:

- **Extensibility**: Users can register custom models without modifying core code
- **Consistency**: All model specs are defined in one place
- **Validation**: Unknown model names are detected at parse time

## Key Design Decisions

1. **XML Generation via ElementTree**: Uses `xml.etree.ElementTree` to build XML DOM trees, ensuring proper escaping and format consistency

2. **YAML Configuration**: Declarative YAML configuration files for reproducibility and version control

3. **Fail-Closed Validation**: The parser checks every configuration section against a key whitelist and refuses ambiguous booleans, out-of-domain priors and inputs the generator cannot honour, rather than silently filling in a default. The three validation gates run *before* any file is written, and any failure exits non-zero. `--force` unlocks only the one thing that is safe to unlock: replacing an existing output file

4. **Estimated vs. Fixed Parameters**: The generator emits a parameter marked `estimate: false` as an inline `<parameter estimate="false">` child of the element that uses it, so that parameter never enters `<state>`, never receives an operator, and produces no trace column (`model.py::_state_node`). The writer then enforces both directions of that invariant: an estimated state node without an operator, or an operator attached to a fixed parameter, is an error

5. **Checkpoint Interval Scales With the Chain**: `mcmc.store_every` defaults to `max(1, chain_length // 10000)` (`config.py`), so any run gets roughly ten checkpoints instead of only one at the very end. The generator deliberately does not serialise the RNG `seed` into the XML: BEAST2's `MCMC` has no `seed` input, so the seed can only be passed as the command-line argument `beast -seed N`

6. **Partition ID Namespacing**: All parameter IDs are prefixed with partition ID (e.g., `gene1.hky.kappa`) for multi-partition support

7. **One Linked Topology**: Every partition is assembled on a single shared tree. The parser refuses `tree: separate` and per-partition tree ids outright: if it silently fell back to linked gene trees, an analysis intended as *BEAST-style would become a fully linked one

8. **Scale-Free Diagnostics**: The distribution-overlap criterion is dimensionless by default — the 1-Wasserstein distance between two priors, divided by the mean width of their 95% intervals, with a threshold of 0.15. The previous absolute criterion (2.0 time units) remains selectable through `diagnostics.overlap_measure`. The check skips distributions that have no defined statistics (OneOnX)

9. **Self-Contained BEAST2 Gate**: The headless validator lives inside the package (`beast2py/beast2_validate.sh`, `beast2py/tools/`), so a plain `pip install .` ships a working gate 3. The script searches a batch of candidate paths for a JDK 17 or newer (`BEAST2_JAVA_CANDIDATES`, `JAVA_HOME`, `/usr/libexec/java_home`, the Homebrew and system JVM directories, and only then `java` on `PATH`), and honours `BEAST2_JAR`, `BEAST2_SKIP_JAVA_HOME_TOOL`, `BEAST2_SKIP_JDK_PROBES` and `BEAST2_VALIDATE_TRACE`

10. **Date-Free Fingerprint**: The analysis fingerprint is `B2P-{hash12}-{version} | Data: {hash16}`, that is, the configuration digest plus a digest over the alignment sequences. The same configuration and the same data therefore produce the same identifier on every machine and on every day; the wall-clock generation time lives only in the JSON sidecar file

11. **BEAST2 v2.7.x Namespace**: Uses `beast.base.*` package names, verified against BEAST 2.7.8 — whose class files require a JDK 17 or newer at run time

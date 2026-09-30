# Tutorial

## 1. Installation

```bash
# Clone the repository
git clone https://github.com/ZengZichao/Beast2Py.git
cd Beast2Py

# Install
pip install -e .
```

Without installing, every command below runs from the repository root as
`python3 -m beast2py.main <subcommand>` — the same entry point.

## 2. Basic Usage: Single Partition Analysis

### Step 1: Prepare your sequence data

Create a FASTA file with your sequences:

```
>human
ATGGCAAATCATT...
>chimp
ATGGCAAATCATT...
>gorilla
ATGGCAAATCATT...
```

### Step 2: Create a YAML configuration

Create `config.yaml`:

```yaml
metadata:
  analysis_name: "my_analysis"
  beast2_version: "2.7.8"

alignments:
  - id: "alignment"
    file: "input.fasta"
    data_type: "nucleotide"

partitions:
  - id: "alignment"
    site_model:
      substitution_model:
        type: "hky"
      gamma_categories: 4
      gamma_shape: {value: 0.5, lower: 0.0}
    clock_model:
      type: "strict"
    tree: "shared"

tree_prior:
  type: "yule"
  birth_rate: {value: 1.0, lower: 0.0}

calibrations:
  - name: "rootCalibration"
    taxa: null
    monophyletic: true
    distribution:
      type: "normal"
      parameters: {mean: 10.0, sigma: 1.0}
    provenance:
      source: "fossil"
      reference: "DOI:10.xxxx"
      calibration_type: "soft"

mcmc:
  chain_length: 10000000

output:
  file: "output.xml"
  embed_fingerprint: true
```

Two defaults are worth noting.

- This example omits `mcmc.store_every`. The generator derives it from the chain length as
  `max(1, chain_length // 10000)` (`1000` for this 10,000,000-generation chain), so one run
  gets about ten `.xml.state` checkpoints instead of only the last one.
- The generator inlines a parameter written `estimate: false` into the element that uses it,
  so that parameter does not enter `<state>`, does not receive an operator, and does not
  appear in the trace.

### Step 3: Generate XML

```bash
beast2py generate --config config.yaml --output output.xml
```

`generate` runs three validation gates before it writes anything: Beast2Py's structural and
semantic checks, calibration conflict detection, and, with `--beast2-validate`, the BEAST2
parser calling `initAndValidate()` on every object. It writes the file only after every gate
that ran has passed, and any failure exits non-zero. Without `--force`, `generate` does not
overwrite an existing `output.xml`. Add `--diagnose`, `--fingerprint` and `--methods` to
produce the report and the reproducibility artifacts in the same call.

### Step 4: Run calibration diagnostics

```bash
beast2py diagnose --config config.yaml --report report.html --sensitivity --visualize
```

### Step 5: Run BEAST2

```bash
beast -overwrite output.xml
```

The RNG seed is a BEAST2 command-line argument, not part of the XML (BEAST2's `MCMC` has
no `seed` input, so writing one makes the model unparseable):

```bash
beast -overwrite -seed 1 output.xml
```

Every `store_every` logged samples, BEAST2 writes an `output.xml.state` checkpoint, which is
what makes an interrupted run resumable. With the default above, a 10,000,000-generation
chain gets about ten checkpoints. Pass `mcmc.store_every: -1` to switch checkpointing off;
the validator then warns that the chain cannot be resumed.

## 3. Quick Mode

For simple single-partition analyses, use the `quick` command:

```bash
beast2py quick \
    --alignment input.fasta \
    --tree-prior yule \
    --subst-model hky \
    --clock-model strict \
    --calibration-yaml calibrations.yaml \
    --chain-length 10000000 \
    --output quick_output.xml
```

The output is named `quick_output.xml`, not §2's `output.xml`: neither command replaces an
existing output unless you pass `--force`, so reusing the name would stop this run with a
refusal rather than mix two analyses under one file.

`quick` also accepts `--gamma-categories N` (default 0, meaning no Gamma), `--pre-burnin N`,
`--name NAME`, `--force` and `-v/--verbose`. It runs the same three gates, and the same
`store_every` default applies.

The calibrations YAML file:

```yaml
- name: "humanChimpMRCA"
  taxa: ["human", "chimp"]
  monophyletic: true
  distribution:
    type: "normal"
    parameters: {mean: 6.0, sigma: 0.5}
```

## 4. Multi-Partition Analysis

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
      linked_to: "gene1"    # Share clock with gene1
    tree: "shared"
```

`tree: "shared"` is the only legal value. Unlinked, per-partition topologies are not
implemented, so the parser refuses `tree: separate` (or a partition id in that slot)
outright instead of silently linking the trees. Otherwise a multi-locus analysis intended
as *BEAST-style would become a fully linked one. Until per-partition trees land, run one
analysis per independent partition. `linked_to` must name another existing partition;
self-links and unknown ids are errors.

## 5. Codon Partitioning

```yaml
alignments:
  - id: "alignment"
    file: "coding.fasta"
    data_type: "nucleotide"
  - id: "codon12"
    data_id: "alignment"
    filter: "1::3,2::3"
    data_type: "nucleotide"
  - id: "codon3"
    data_id: "alignment"
    filter: "3::3"
    data_type: "nucleotide"

partitions:
  - id: "codon12"
    site_model: {substitution_model: {type: "hky"}}
    clock_model: {type: "strict"}
    tree: "shared"
  - id: "codon3"
    site_model: {substitution_model: {type: "gtr"}}
    clock_model: {linked_to: "codon12"}
    tree: "shared"
```

`data_id` is required whenever a `filter` is declared, and it must name an already-declared
alignment other than the filtered alignment itself.

Filter expressions follow BEAST2's own `FilteredAlignment.parseFilterSpec` grammar, and the
positions are **1-based** (counted from 1). An expression is a comma-separated list whose
terms may be:

- a single site (`N`)
- an inclusive range (`from-to`, where an empty end means the last site)
- a range stepping every k-th site (`from-to\step`)
- an iterator (`from:to:step`; the two examples above use `from::step`)

So `1::3,2::3` keeps positions 1, 2, 4, 5, 7, 8 and so on (the first two codon positions),
and `3::3` keeps 3, 6, 9 and so on (the third position).

The parser range-checks every term against the source alignment. It reports an error naming
the offending term for an unparsable expression, a step below 1, a coordinate beyond the
alignment length, or a term that selects no site at all.

## 6. Generating Methods Description

```bash
beast2py methods --config config.yaml --output methods.tex
```

This generates a LaTeX paragraph describing the analysis configuration, suitable for inclusion in a manuscript. Omit `--output` to print it instead.

## 7. Generating Pipeline Files

```bash
# Snakemake pipeline
beast2py pipeline --config config.yaml --type snakemake --output-dir pipeline

# Nextflow pipeline
beast2py pipeline --config config.yaml --type nextflow --output-dir pipeline

# Generate both
beast2py pipeline --config config.yaml --type both --output-dir pipeline
```

`--output-dir` defaults to the current directory; `--type` defaults to `snakemake`.

## 8. Validation

```bash
# Structural validation
beast2py validate --xml output.xml

# Also validate with BEAST2
beast2py validate --xml output.xml --beast2
```

`--beast2` runs the headless BEAST2 parser. It parses the file *and* calls
`initAndValidate()` on every object, so the check is about whether the model can be
initialised, not merely whether it is well-formed. This is the same code path as Gate 3 of
`generate` (`--beast2-validate`), and it needs BEAST2 plus a **JDK 17 or newer**.

The bundled `beast2_validate.sh` (shipped inside the package, with only a thin wrapper at
the repository root) looks for `BEAST.base.jar` in this order:

1. `$BEAST2_JAR`
2. the newest `~/.beast/2.7/BEAST.base/*/lib/BEAST.base.jar`
3. `tools/lib/` or `lib/` next to the script

The script also appends any add-on packages installed under `~/.beast/2.7` to the class path.

The script then searches a batch of candidate paths for a JDK 17 or newer instead of
trusting `PATH` — an ancient Oracle `java` shim on `PATH` cannot shadow a usable
`openjdk@17`. The search order is: `BEAST2_JAVA_CANDIDATES` (a colon-separated list),
`JAVA_HOME`, `/usr/libexec/java_home`, the Homebrew and system JVM directories, and only
then `java` on `PATH`.

Two switches remove the corresponding sources: `BEAST2_SKIP_JAVA_HOME_TOOL` drops
`/usr/libexec/java_home`, and `BEAST2_SKIP_JDK_PROBES` drops the Homebrew and system JVM
directories. Setting `BEAST2_VALIDATE_TRACE` prints the `java` that was chosen.

Exit codes: `0` means every file is valid, `1` means at least one is not, `2` means a setup
error. The script validates any number of files in one call:

```bash
bash beast2_validate.sh output.xml examples/output_basic.xml
```

## 9. Listing Supported Models

```bash
beast2py list-models
```

"""Test fixtures and data for Beast2Py tests."""

import os
from pathlib import Path

# BEAST 2.7.x needs JDK 17+. If JAVA_HOME is unset, fall back to the Homebrew
# JDK 17 locations that beast2_validate.sh also probes, so the integration
# tests run deterministically on macOS. An explicit JAVA_HOME is respected
# as-is (the validation script reports a clear error if it is too old).
_HOMEBREW_JDK17_HOMES = (
    "/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home",
    "/opt/homebrew/opt/openjdk/libexec/openjdk.jdk/Contents/Home",
    "/usr/local/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home",
    "/usr/local/opt/openjdk/libexec/openjdk.jdk/Contents/Home",
)

if not os.environ.get("JAVA_HOME"):
    for _home in _HOMEBREW_JDK17_HOMES:
        if (Path(_home) / "bin" / "java").is_file():
            os.environ["JAVA_HOME"] = _home
            break

# Path to the examples directory
EXAMPLES_DIR = Path(__file__).parent.parent / "examples"

# Minimal FASTA content for testing.
#
# These sequences carry real variation (four variable sites, three distinct
# patterns). The previous fixture was four *identical* 32 bp sequences, i.e.
# zero informative sites, so every code path that depends on site patterns,
# state counts or frequency estimation was exercised on a degenerate sample
#
TEST_FASTA = """>taxonA
ACGTACGTACGTACGTTTTGCAAGGCATCGATCGGG
>taxonB
ACGTACGTACGTAAGTTTTGCAAGGCATCGATCGGA
>taxonC
ACGAACGTACGTAAGTTTTGCAAGGCATCGATCGGG
>taxonD
ACGAACGTACGTAAGTTTTGCAAGGCATCGATCGGA
"""

# Minimal YAML config for testing
TEST_CONFIG_YAML = """
metadata:
  analysis_name: "test_analysis"
  analysis_description: "Test analysis"
  author: "Test"
  date: "2026-08-01"

alignments:
  - id: "alignment"
    file: "test.fasta"
    format: "auto"
    data_type: "nucleotide"

partitions:
  - id: "alignment"
    site_model:
      substitution_model:
        type: "hky"
      gamma_categories: 0
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

mcmc:
  type: "standard"
  chain_length: 100000
  pre_burnin: 0

output:
  file: "test_output.xml"
  embed_fingerprint: true
"""

# YAML config with multiple calibrations (for conflict detection tests)
TEST_MULTI_CAL_YAML = """
metadata:
  analysis_name: "test_multi_cal"
  analysis_description: "Test with multiple calibrations"

alignments:
  - id: "alignment"
    file: "test.fasta"
    format: "auto"
    data_type: "nucleotide"

partitions:
  - id: "alignment"
    site_model:
      substitution_model:
        type: "hky"
    clock_model:
      type: "strict"
    tree: "shared"

tree_prior:
  type: "yule"

calibrations:
  - name: "calA"
    taxa: ["taxonA", "taxonB"]
    monophyletic: true
    distribution:
      type: "normal"
      parameters: {mean: 5.0, sigma: 0.5}
  - name: "calB"
    taxa: ["taxonA", "taxonB", "taxonC"]
    monophyletic: true
    distribution:
      type: "normal"
      parameters: {mean: 10.0, sigma: 1.0}
  - name: "rootCal"
    taxa: null
    monophyletic: true
    distribution:
      type: "uniform"
      parameters: {lower: 15.0, upper: 30.0}

mcmc:
  type: "standard"
  chain_length: 100000

output:
  file: "test_multi_cal.xml"
"""

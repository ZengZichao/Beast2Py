# Calibration Prior Best Practices Guide

English | [中文](calibration_guide_zh.md)

## Introduction

Node calibration is a critical step in divergence time estimation. Poorly chosen calibrations can lead to biased posterior estimates or "prior dominance" where the prior overwhelms the data. This guide provides best practices for specifying calibration priors in Beast2Py.

## Key Principles

### 1. Use Appropriate Distributions

- **Normal**: For well-constrained calibrations with symmetric uncertainty
- **LogNormal**: For calibrations with right-skewed uncertainty (common in fossil data)
- **Uniform**: For hard boundary constraints only
- **Exponential**: For minimum-age constraints with no upper bound

Beast2Py accepts eleven distribution types in total. The four above are the ones that usually
make biological sense for a node age. The others are `gamma`, `beta`, `laplace`,
`inverse_gamma`, `poisson` and `chi_square` (parameterisable shapes), and `one_on_x`
(`f(x) = 1/x`). The last one is an **improper** density: it has no mean and no 95% interval,
so Beast2Py reports its statistics as "undefined" and excludes it from the overlap check.
Use `one_on_x` as a deliberately flat prior only on a parameter that something else bounds.

### 2. Soft vs. Hard Bounds

- **Soft bounds** (recommended): Allow the posterior to extend beyond the prior range, letting data inform the estimate
- **Hard bounds**: Constrain the posterior strictly; use only when fossil evidence is definitive

Each calibration point records this choice in `provenance.calibration_type` (a free string; anything other than `"hard"` is treated as a soft bound, and omitting it defaults to a **soft** bound). Beast2Py also *enforces* the claim you record: at parse time it compares the prior's 95% interval with `provenance.original_age_min` / `original_age_max`, and

- a calibration declared **hard** whose prior lies entirely outside the recorded fossil range is an **error** — the run stops before any file is written;
- a calibration declared **soft** in the same situation is a **warning**, because a soft bound is *meant* to put most of its mass above the fossil minimum (a lognormal whose 2.5% quantile sits at the hard minimum is the BEAST2 convention).

### 3. Document Provenance

Always specify the `provenance` section for each calibration point:

```yaml
calibrations:
  - name: "homininaMRCA"
    taxa: ["Homo_sapiens", "Pan"]
    distribution:
      type: "lognormal"
      parameters: {M: 0.2, S: 0.35}
      offset: 5.0
    provenance:
      source: "other"            # fossil / biogeographic / secondary / other
      reference: "Synthetic example data; no publication source"
      calibration_type: "soft"    # soft / hard
      original_age_min: 5.0
      original_age_max: 6.0
      notes: "Synthetic example calibration; values are illustrative only"
```

`provenance` accepts exactly `source`, `reference`, `calibration_type`, `original_age_min`, `original_age_max` and `notes`; the parser refuses any other key at parse time. `original_age_min` must not exceed `original_age_max`.

### 4. Check for Conflicts

Use the `diagnose` command to detect:

- **Temporal inconsistency**: A nested pair whose parent prior is *wholly* younger than its child prior. The detector reports an error only when the parent's 95% upper bound is younger than the child's 95% lower bound. Overlapping intervals are not reported, including a parent whose mean is younger than its child's: a soft parent prior with a long tail can have a younger mean and still be perfectly compatible. A root calibration (`taxa: null`) takes part as the parent of every other clade
- **Prior similarity (descriptive, off by default)**: Two independent clades may legitimately have similar age priors. If you enable `diagnostics.prior_similarity` with an explicit `similarity_threshold`, the tool reports `prior_similarity` entries with severity `info`. They carry no merge recommendation and never block Gate 2. `OneOnX` has no defined mean or interval, so a pair involving it is reported as unchecked.
- **Monophyly conflict**: Overlapping but non-nested taxon sets

The same conflict detection runs automatically as Gate 2 of `generate` and `quick`: when it reports an error-severity conflict, the gate writes nothing and the command exits non-zero.

### 5. Run Sensitivity Analysis

Always assess the marginal influence of each calibration point:

```bash
beast2py diagnose --config config.yaml --report report.html --sensitivity
```

This generates N+1 prior-only sampling XMLs (N is the number of calibrations): one with the full calibration set, plus N leave-one-out variants. They land in a `sensitivity/` subfolder of `<report dir>/diagnostics`, or of the directory given by `--sensitivity-xml-dir`. Run each XML with BEAST2, then compare the tree height and the calibrated node ages between the full run and each leave-one-out run. These chains sample from the prior, so the comparison is about the **effective marginal prior** each calibration set induces; it says nothing about a posterior, because no likelihood is involved. Choose the chain length, sampling interval and burn-in yourself, decide an ESS criterion before running, and report all three with the result — a prior-only run is not automatically cheap or automatically converged.

## Common Pitfalls

### 1. Prior Dominance

**Problem**: The calibration prior overwhelms the sequence data, making the posterior indistinguishable from the prior.

**Solution**: Check the prior-only sampling result against the posterior. If they are identical, the prior is too informative. The same effect is visible without an extra run when the model is over-constrained: the trace simply tracks the prior HPD.

### 2. Conflicting Calibrations

**Problem**: Two calibrations specify incompatible time ranges for nested clades.

**Solution**: Run `diagnose` to detect temporal conflicts, then adjust calibrations so parent clade ranges are older than child clade ranges. Note that `generate` refuses to write anything when Gate 2 reports an error-severity conflict, so a pipeline that exits 0 has already passed this check.

### 3. Overly Narrow Priors

**Problem**: Using a very narrow Normal distribution (a very small `sigma`) can artificially constrain the posterior.

**Solution**: Use wider distributions or LogNormal distributions for fossil calibrations, which typically have asymmetric uncertainty.

### 4. A Provenance Record That Disagrees With the Prior

**Problem**: The recorded fossil range and the prior actually used are two different ages. BEAST2 only ever sees the distribution, so this kind of recording error stays invisible to it.

**Solution**: Fill in `original_age_min` / `original_age_max` honestly. Beast2Py compares the two recorded ranges with the prior's 95% interval at parse time: a hard-bound disagreement is an error, a soft-bound one is a warning that the two should be explained together in the methods text.

## References

- Heled & Drummond (2012) Syst Biol 61: 138-149, DOI:10.1093/sysbio/syr087 — Calibrated Yule and MRCA calibration methods
- Rieux & Balloux (2016) Mol Ecol 25: 1911-1924, DOI:10.1111/mec.13586 — Calibration best practices
- Ho & Phillips (2009) Syst Biol 58: 367-380, DOI:10.1093/sysbio/syp035 — Calibrating the molecular clock
- Parham et al. (2012) Syst Biol 61: 346-359, DOI:10.1093/sysbio/syr107 — Best practices for justifying fossil calibrations

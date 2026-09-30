# BEAST2 XML Format Guide

## Overview

Beast2Py generates BEAST2 v2.7.x compatible XML files using `xml.etree.ElementTree`,
verified against BEAST 2.7.8. This document describes the XML structure the tool produces.
The generator assembles the XML in memory and writes it to disk only after the validation
gates (structural and semantic checks, calibration conflict detection, and — on request —
the BEAST2 parser plus `initAndValidate()`) have all passed. The files on disk therefore
always have the structure described here.

## XML Structure

```xml
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast version='2.0' namespace='...'>
    <!-- 1. Sequence alignment data -->
    <data id="alignment" dataType="nucleotide">
        <sequence taxon="human">ATGGCAAAT...</sequence>
        ...
    </data>

    <!-- Optional: FilteredAlignment for codon partitioning -->
    <data data="@alignment" dataType="nucleotide" filter="3::3" id="thirdSites"
          spec="beast.base.evolution.alignment.FilteredAlignment"/>

    <!-- 2. TaxonSets for calibration -->
    <!-- Implementation convention: the first occurrence of a taxon defines it
         with id=, later occurrences reference it with idref=. Writing every
         entry as idref makes BEAST2 fail with "Could not find object
         associated with idref ..." when no top-level taxon definition exists. -->
    <taxonset spec='TaxonSet' id='humanChimpMRCA.taxonset'>
        <taxon spec='Taxon' id='human'/><taxon spec='Taxon' idref='chimp'/>
    </taxonset>

    <!-- 3. Posterior distribution -->
    <distribution id="posterior" spec="CompoundDistribution">
        <distribution id="prior" spec="CompoundDistribution">
            <!-- Tree prior -->
            <!-- MRCA calibration priors -->
            <!-- Parameter priors -->
        </distribution>
        <distribution id="likelihood" spec="CompoundDistribution">
            <!-- TreeLikelihood for each partition -->
        </distribution>
    </distribution>

    <!-- 4. MCMC run configuration -->
    <run id="mcmc" spec="MCMC" chainLength="10000000" preBurnin="0" storeEvery="1000">
        <init spec="beast.base.evolution.tree.coalescent.RandomTree" id="randomTree"
              initial="@tree" taxa="@alignment" rootHeight="10.275">...</init>
        <state id="state" storeEvery="1000">...</state>
        <distribution idref="posterior"/>
        <operator>...</operator>
        <logger id="tracelog" .../>
        <logger id="treelog" .../>
        <logger id="screenlog" .../>
    </run>
</beast>
```

### The `<state>` element holds only estimated parameters

A parameter the user marks `estimate: false` never enters `<state>`. The generator writes it
*inline* as a child of the element that consumes it, so it gets no operator and produces no
trace column:

```xml
<!-- estimated: a state node referenced from the substitution model -->
<parameter id="alignment.hky.kappa" name="stateNode" value="2" estimate="true" lower="0"/>
<substModel spec="HKY" id="alignment.hky">
    <parameter idref="alignment.hky.kappa" name="kappa"/>
</substModel>

<!-- fixed: an inline parameter, no state node, no operator, not traced -->
<substModel spec="HKY" id="alignment.hky">
    <parameter id="alignment.hky.kappa" name="kappa" value="2" estimate="false"
               dimension="1" lower="0"/>
</substModel>
```

The writer guards both halves of that invariant and raises an error rather than emitting a
silently wrong model: every estimated state node must carry an operator (otherwise BEAST2
would only warn and then pin the parameter to its start value), and an operator may only
target an estimated parameter.

### Attributes of the `<run>` element

`chainLength`, `preBurnin` and `storeEvery` are written on `<run>` (and `storeEvery` also
on `<state>`). `storeEvery` counts **logged samples**, not generations, and defaults to
`max(1, chain_length // 10000)` unless `mcmc.store_every` says otherwise. There is
deliberately **no `seed` attribute**: BEAST2's `MCMC` has no `seed` input, and writing one
makes the XML unparseable — the RNG seed is the command-line argument `beast -seed N`. The
`$(seed)` macro may still appear in logger file names, where BEAST2 expands it at run time.

For `mcmc.type: nested_sampling` the `<run>` element switches to the add-on's own engine:

```xml
<run id="mcmc" spec="nestedsampling.gss.NS" chainLength="10000000" preBurnin="0"
     storeEvery="10000" particleCount="1" subChainLength="10000">
```

## Namespace

For BEAST2 v2.7.x, the namespace uses `beast.base.*` package names:

```
beast.base.evolution.alignment:beast.base.evolution.speciation:beast.pkgmgmt:beast.base.core:...
```

## Multi-Partition Support

The generator writes a separate `<data>` element and a separate `<TreeLikelihood>` element for each partition:

```xml
<data id='gene1' dataType='nucleotide'>...</data>
<data id='gene2' dataType='nucleotide'>...</data>

<distribution id="likelihood" spec="CompoundDistribution">
    <distribution data='@gene1' id="treeLikelihood.gene1" spec="TreeLikelihood" tree='@tree'>
        <siteModel id="gene1.siteModel" spec="SiteModel">...</siteModel>
        <branchRateModel id="gene1.clockModel"
                         spec="beast.base.evolution.branchratemodel.UCRelaxedClockModel"
                         normalize="true">...</branchRateModel>
    </distribution>
    <distribution data='@gene2' id="treeLikelihood.gene2" spec="TreeLikelihood" tree='@tree'>
        <siteModel idref="gene1.siteModel"/>           <!-- Shared site model -->
        <branchRateModel idref="gene1.clockModel"/>    <!-- Shared clock -->
    </distribution>
</distribution>
```

The generator emits a shared clock exactly once, with its own `id`, and every other
partition points at it with `idref` (a shared site model works the same way). All
partitions are assembled on **one** `tree='@tree'`: per-partition topologies are not
implemented, so the parser refuses `tree: separate` outright instead of silently linking
the trees. An XML file with several `<tree>` state nodes is therefore never produced.

## Calibration Distributions

All distributions support the `offset` attribute (default 0.0). The `spec` column below
gives the class name; what is written into the XML is the fully qualified form, e.g.
`beast.base.inference.distribution.LogNormalDistributionModel`.

| Distribution | spec | Parameters |
|------|------|------|
| Normal | `Normal` | mean, sigma |
| LogNormal | `LogNormalDistributionModel` | M, S, meanInRealSpace |
| Uniform | `Uniform` | lower, upper |
| Exponential | `Exponential` | mean |
| Gamma | `Gamma` | alpha, beta |
| Beta | `Beta` | alpha, beta |
| OneOnX | `OneOnX` | (none) |
| Laplace | `LaplaceDistribution` | mu, scale |
| InverseGamma | `InverseGamma` | alpha, beta |
| Poisson | `Poisson` | lambda |
| ChiSquare | `ChiSquare` | df (config key: `dof`) |

`OneOnX` (`f(x) = 1/x`) is an improper density: its mean and 95% interval are undefined,
so the diagnostics report them as such and exclude it from the distribution-overlap check.

## Tree Priors

| Tree Prior | spec |
|------|------|
| Yule | `YuleModel` |
| CalibratedYule | `CalibratedYuleModel` |
| Birth-Death | `BirthDeathGernhard08Model` |
| Coalescent (Constant) | `Coalescent` + `ConstantPopulation` |
| Coalescent (Exponential) | `Coalescent` + `ExponentialGrowth` |
| Bayesian Skyline | `BayesianSkyline` |
| EBSP | `Coalescent` + `CompoundPopulationFunction` (BEAST 2.7 base package) |
| BD Skyline Serial | `bdsky.evolution.speciation.BirthDeathSkylineModel` (requires the bdsky add-on) |

## Clock Models

| Clock Model | spec |
|------|------|
| Strict Clock | `StrictClockModel` |
| UC Relaxed LogNormal | `UCRelaxedClockModel` |
| UC Relaxed Exponential | `UCRelaxedClockModel` |
| Random Local Clock | `RandomLocalClockModel` |

UCLN and UCE share the `UCRelaxedClockModel` spec and differ only in the branch-rate
distribution: UCLN carries a `LogNormalDistributionModel` with `M = 1` and
`S = ucld.stdev`, UCE an `Exponential` with mean `ucld.mean`. The generator emits every
relaxed clock with `normalize="true"`, matching BEAUti's reference layout, so rates and
node ages stay identifiable. The random local clock adds two state nodes of its own,
`rates` and `indicators`.

## Analysis Fingerprint

The tool embeds an analysis fingerprint as an XML comment:

```xml
<!-- Analysis Fingerprint: B2P-04ab3d09277b-0.1.0 | Data: 30be5613d3ff6862 -->
```

Format: `B2P-{config_hash_12chars}-{tool_version} | Data: {data_hash_16chars}`

`config_hash` is the first 12 hex characters of the SHA-256 digest of the
serialised configuration, and `data_hash` is the first 16 hex characters of the
SHA-256 digest of the alignment *sequences*.

The identifier deliberately carries no calendar date, no time zone and not the
output file name, so the same configuration and the same alignments yield the
same fingerprint on every machine and on every day. Changing any model component,
calibration, MCMC setting or sequence changes the fingerprint with it.

The generation time is recorded separately in the JSON sidecar file and never
becomes part of the identifier. The manuscript's Supplementary Section S4
documents the digest inputs and the truncation length.

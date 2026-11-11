# Declared-capability coverage matrix (v0.2.0)

This matrix covers the **36 declared model capabilities** that Beast2Py
advertises: 4 clock models, 8 tree priors, 13 substitution models and 11
calibration distributions. Every row has

* a YAML fixture that actually uses the capability,
* the XML class the generator must emit,
* optional external Gate 3 evidence (parse + `initAndValidate`), and
* an automated test or an explicit reason why the check is not a test.

Gate 3 evidence is produced by running the fixtures through BEAST2 and is kept
outside this repository. The repository ships only synthetic fixtures and the
automated checks in `tests/test_coverage_fixtures.py`; no Gate 3 log or XML is
bundled.

> Scope note: "36/36 covered" means *native parse and model-initialisation
> coverage*. It is not a claim of biological correctness, MCMC convergence or
> posterior accuracy.

## Clock models (4)

| # | Capability | Fixture | XML class | Gate 3 evidence (external) | Test |
|---|---|---|---|---|---|
| 1 | Strict clock | `examples/config_basic.yaml` | `beast.base.evolution.branchratemodel.StrictClockModel` | `logs/examples/config_basic.log` | `test_every_shipped_example_generates_and_validates` |
| 2 | Uncorrelated lognormal (UCLN) | `examples/config_advanced.yaml` | `beast.base.evolution.branchratemodel.UCRelaxedClockModel` (lognormal) | `logs/examples/config_advanced.log` | same |
| 3 | Uncorrelated exponential (UCE) | `examples/config_relaxed_exponential.yaml` | `beast.base.evolution.branchratemodel.UCRelaxedClockModel` (exponential) | `logs/examples/config_relaxed_exponential.log` | same |
| 4 | Random local clock (RLC) | `examples/config_rlc.yaml` | `beast.base.evolution.branchratemodel.RandomLocalClockModel` | `logs/examples/config_rlc.log` | same |

## Tree priors (8)

| # | Capability | Fixture | XML class | Gate 3 evidence (external) | Test |
|---|---|---|---|---|---|
| 5 | Yule | `examples/config_basic.yaml` | `beast.base.evolution.speciation.YuleModel` | `logs/examples/config_basic.log` | `test_every_shipped_example_generates_and_validates` |
| 6 | Calibrated Yule | `examples/config_calibrated_yule.yaml` | `beast.base.evolution.speciation.CalibratedYuleModel` | `logs/examples/config_calibrated_yule.log` | same |
| 7 | Birth-death (Gernhard 2008) | `examples/config_advanced.yaml` | `beast.base.evolution.speciation.BirthDeathGernhard08Model` | `logs/examples/config_advanced.log` | same |
| 8 | Coalescent, constant size | `tests/fixtures/coverage/config_coalescent_constant.yaml` | `Coalescent` + `ConstantPopulation` | `logs/config_coalescent_constant.log` | `test_coalescent_constant_fixture_initialises_a_population_model` |
| 9 | Coalescent, exponential growth | `examples/config_coalescent_exponential.yaml` | `Coalescent` + `ExponentialGrowth` | `logs/examples/config_coalescent_exponential.log` | `test_every_shipped_example_generates_and_validates` |
| 10 | Bayesian skyline | `examples/config_skyline.yaml` | `BayesianSkyline` | `logs/examples/config_skyline.log` | same |
| 11 | Extended Bayesian skyline (EBSP) | `examples/config_ebsp.yaml` | `Coalescent` + `CompoundPopulationFunction` | `logs/examples/config_ebsp.log` | same |
| 12 | Birth-death skyline, serial (BDSKY add-on 1.5.1) | `examples/config_bd_skyline.yaml` | `bdsky.evolution.speciation.BirthDeathSkylineModel` | `logs/examples/config_bd_skyline.log` | same; add-on version recorded in the archived `coverage/README.md`, which lives outside this repository like the rest of the Gate 3 evidence |

## Substitution models (13)

| # | Capability | Fixture | XML class | Gate 3 evidence (external) | Test |
|---|---|---|---|---|---|
| 13 | JC69 | `examples/config_subst_jc69.yaml` | `JukesCantor` | `logs/examples/config_subst_jc69.log` | `test_every_shipped_example_generates_and_validates` |
| 14 | HKY | `examples/config_basic.yaml` | `HKY` | `logs/examples/config_basic.log` | same |
| 15 | TN93 | `examples/config_skyline.yaml` | `TN93` | `logs/examples/config_skyline.log` | same |
| 16 | GTR | `examples/config_advanced.yaml` | `GTR` | `logs/examples/config_advanced.log` | same |
| 17 | SYM | `examples/config_subst_sym.yaml` | `SYM` | `logs/examples/config_subst_sym.log` | same |
| 18 | TIM | `examples/config_subst_tim.yaml` | `TIM` | `logs/examples/config_subst_tim.log` | same |
| 19 | TVM | `examples/config_subst_models.yaml` | `TVM` | `logs/examples/config_subst_models.log` | same |
| 20 | WAG | `examples/config_aa_model.yaml` | `WAG` | `logs/examples/config_aa_model.log` | same |
| 21 | JTT | `examples/config_aa_jtt.yaml` | `JTT` | `logs/examples/config_aa_jtt.log` | same |
| 22 | Dayhoff | `tests/fixtures/coverage/config_aa_dayhoff.yaml` | `Dayhoff` | `logs/config_aa_dayhoff.log` | `test_aminoacid_matrix_fixture_emits_expected_class` |
| 23 | BLOSUM62 | `tests/fixtures/coverage/config_aa_blosum62.yaml` | `Blosum62` | `logs/config_aa_blosum62.log` | same |
| 24 | CPREV | `tests/fixtures/coverage/config_aa_cprev.yaml` | `CPREV` | `logs/config_aa_cprev.log` | same |
| 25 | MTREV | `tests/fixtures/coverage/config_aa_mtrev.yaml` | `MTREV` | `logs/config_aa_mtrev.log` | same |

## Calibration distributions (11)

| # | Capability | Fixture | XML class | Gate 3 evidence (external) | Test |
|---|---|---|---|---|---|
| 26 | Normal | `tests/fixtures/coverage/config_dist_normal.yaml` | `beast.base.inference.distribution.Normal` | `logs/config_dist_normal.log` | `test_distribution_fixture_emits_expected_class` |
| 27 | Lognormal | `examples/config_basic.yaml` | `beast.base.inference.distribution.LogNormalDistributionModel` | `logs/examples/config_basic.log` | same |
| 28 | Uniform | `examples/config_relaxed_clock.yaml` (Crown Catarrhini) | `beast.base.inference.distribution.Uniform` | `logs/examples/config_relaxed_clock.log` | same |
| 29 | Exponential | `examples/config_calibrated_yule.yaml` | `beast.base.inference.distribution.Exponential` | `logs/examples/config_calibrated_yule.log` | same |
| 30 | Gamma | `tests/fixtures/coverage/config_dist_gamma.yaml` | `beast.base.inference.distribution.Gamma` | `logs/config_dist_gamma.log` | `test_distribution_fixture_emits_expected_class` |
| 31 | Beta | `tests/fixtures/coverage/config_dist_beta.yaml` | `beast.base.inference.distribution.Beta` | `logs/config_dist_beta.log` | same; `test_beta_offset_is_legal_and_beta_without_offset_warns` |
| 31a | Beta + offset (legal exception) | `tests/fixtures/coverage/config_dist_beta_offset.yaml` | `beast.base.inference.distribution.Beta` | `logs/config_dist_beta_offset.log` | same (legal exception, no warning) |
| 32 | OneOnX | `tests/fixtures/coverage/config_dist_one_on_x.yaml` | `beast.base.inference.distribution.OneOnX` | `logs/config_dist_one_on_x.log` | same + `test_one_on_x_and_poisson_warn_but_parse`; diagnostics report unchecked |
| 33 | Laplace | `tests/fixtures/coverage/config_dist_laplace.yaml` | `beast.base.inference.distribution.LaplaceDistribution` | `logs/config_dist_laplace.log` | `test_distribution_fixture_emits_expected_class` |
| 34 | Inverse gamma | `tests/fixtures/coverage/config_dist_inverse_gamma.yaml` | `beast.base.inference.distribution.InverseGamma` | `logs/config_dist_inverse_gamma.log` | same |
| 35 | Poisson | `tests/fixtures/coverage/config_dist_poisson.yaml` | `beast.base.inference.distribution.Poisson` | `logs/config_dist_poisson.log` | same + discrete warning test |
| 36 | Chi-square | `tests/fixtures/coverage/config_dist_chi_square.yaml` | `beast.base.inference.distribution.ChiSquare` | `logs/config_dist_chi_square.log` | `test_chi_square_df_is_an_integer_attribute` |

## Notes and explicit limits

* **Beta and Beta+offset are separate rows** because they have different
  support and different diagnostic behaviour. The 11 *distribution families*
  are rows 26-36; Beta and Beta+offset share one family but are two use cases.
  Counting families gives 11, counting fixtures gives 12.
* **OneOnX** is parsed and initialised by BEAST2, but its mean and quantiles
  are undefined. The diagnostics report such pairs as `unchecked`, never as
  passed.
* **Poisson** is a discrete distribution. It is accepted by BEAST2, but the
  configuration parser warns (and `strict_distributions: true` rejects) when
  it is used as a continuous node-age prior.
* **BDSKY** requires add-on package BDSKY 1.5.1; the Gate 3 run used
  the project-local copy at `tools/beast_home/.beast/2.7/BDSKY/`. The core
  capabilities use BEAST 2.7.7 with its bundled OpenJDK 17.0.3.
* Gate 3 proves that BEAST2 can instantiate the model. It does not prove that
  the MCMC mixes, that the priors are biologically appropriate, or that the
  analysis converges.

"""Tests for the ModelBuilder module."""

from beast2py.models import (
    Alignment,
    BEASTConfig,
    ClockModelConfig,
    ClockModelType,
    DataType,
    InitTreeType,
    InitializationConfig,
    MCMCConfig,
    Partition,
    RealParameter,
    Sequence,
    SiteModelConfig,
    TreePriorType,
)
from beast2py.model import ModelBuilder
from beast2py.partition import PartitionManager


def _make_config(
    tree_prior_type=TreePriorType.YULE,
    tree_prior_params=None,
    clock_type=ClockModelType.STRICT,
    subst_model="hky",
):
    aln = Alignment(
        id="alignment",
        data_type=DataType.NUCLEOTIDE,
        sequences=[Sequence(taxon=f"t{i}", sequence="ACGT") for i in range(4)],
    )
    sm = SiteModelConfig(substitution_model={"type": subst_model})
    cm = ClockModelConfig(type=clock_type)
    p = Partition(id="alignment", alignment=aln, site_model=sm, clock_model=cm)
    pm = PartitionManager([p])

    return (
        BEASTConfig(
            metadata={"analysis_name": "test"},
            partitions=[p],
            tree_prior_type=tree_prior_type,
            tree_prior_params=tree_prior_params or {"birth_rate": {"value": 1.0}},
            calibrations=[],
            mcmc=MCMCConfig(chain_length=100000),
            initialization=InitializationConfig(tree_type=InitTreeType.RANDOM),
        ),
        pm,
    )


def test_build_substitution_model_hky():
    sm = SiteModelConfig(substitution_model={"type": "hky", "kappa": {"value": 2.0}})
    elem, state_nodes = ModelBuilder.build_substitution_model("alignment", sm)
    assert elem.tag == "substModel"
    assert elem.get("spec") == "HKY"
    assert len(state_nodes) >= 1  # kappa + frequencies


def test_build_substitution_model_gtr():
    sm = SiteModelConfig(
        substitution_model={
            "type": "gtr",
            "rates": {"value": "1.0 1.0 1.0 1.0 1.0 1.0", "dimension": 6},
        }
    )
    elem, state_nodes = ModelBuilder.build_substitution_model("alignment", sm)
    assert elem.get("spec") == "GTR"
    assert len(state_nodes) >= 1


def test_build_site_model():
    sm = SiteModelConfig(
        substitution_model={"type": "hky"},
        gamma_categories=4,
        gamma_shape=RealParameter(value=0.5, lower=0.0),
    )
    elem, state_nodes = ModelBuilder.build_site_model("alignment", sm)
    assert elem.tag == "siteModel"
    assert elem.get("gammaCategoryCount") == "4"


def test_build_clock_model_strict():
    cm = ClockModelConfig(type=ClockModelType.STRICT)
    elem, state_nodes = ModelBuilder.build_clock_model("alignment", cm, "tree", 4)
    assert elem.tag == "branchRateModel"
    assert "StrictClockModel" in elem.get("spec")


def test_build_clock_model_ucln():
    cm = ClockModelConfig(
        type=ClockModelType.UCLN,
        ucld_mean=RealParameter(value=1.0),
        ucld_stdev=RealParameter(value=0.333),
    )
    elem, state_nodes = ModelBuilder.build_clock_model("alignment", cm, "tree", 4)
    assert "UCRelaxedClockModel" in elem.get("spec")
    assert any("ucld.stdev" in sn[0] for sn in state_nodes)


def test_build_clock_model_rlc():
    cm = ClockModelConfig(type=ClockModelType.RLC)
    elem, state_nodes = ModelBuilder.build_clock_model("alignment", cm, "tree", 4)
    assert "RandomLocalClockModel" in elem.get("spec")
    assert any("indicators" in sn[0] for sn in state_nodes)


def test_build_tree_prior_yule():
    config, pm = _make_config(TreePriorType.YULE, {"birth_rate": {"value": 1.0}})
    elem, state_nodes, prior_ids = ModelBuilder.build_tree_prior(config, "tree", pm)
    assert elem.tag == "distribution"
    assert elem.get("spec") == "beast.base.evolution.speciation.YuleModel"
    assert "birthRate" in state_nodes[0][0]


def test_build_tree_prior_birth_death():
    config, pm = _make_config(
        TreePriorType.BIRTH_DEATH,
        {"birth_rate": {"value": 1.0}, "death_rate": {"value": 0.5}},
    )
    elem, state_nodes, _ = ModelBuilder.build_tree_prior(config, "tree", pm)
    assert "BirthDeathGernhard08Model" in elem.get("spec")


def test_build_tree_prior_coalescent_constant():
    config, pm = _make_config(
        TreePriorType.COALESCENT_CONSTANT,
        {"pop_size": {"value": 10.0, "lower": 0.0}},
    )
    elem, state_nodes, _ = ModelBuilder.build_tree_prior(config, "tree", pm)
    assert elem.get("spec") == "Coalescent"


def test_build_tree_prior_calibrated_yule():
    config, pm = _make_config(
        TreePriorType.CALIBRATED_YULE,
        {"birth_rate": {"value": 1.0}},
    )
    elem, state_nodes, _ = ModelBuilder.build_tree_prior(config, "tree", pm)
    assert elem.get("spec") == "beast.base.evolution.speciation.CalibratedYuleModel"


def test_build_tree_prior_bayesian_skyline():
    config, pm = _make_config(TreePriorType.BAYESIAN_SKYLINE)
    elem, state_nodes, _ = ModelBuilder.build_tree_prior(config, "tree", pm)
    assert elem.get("spec") == "BayesianSkyline"


def test_build_tree_prior_ebsp():
    config, pm = _make_config(TreePriorType.EBSP)
    elem, state_nodes, _ = ModelBuilder.build_tree_prior(config, "tree", pm)
    # EBSP is a PopulationFunction, so it is wrapped in a Coalescent prior
    assert elem.get("spec") == "Coalescent"
    cpf = elem.find("populationModel")
    assert cpf is not None
    assert cpf.get("spec").endswith("CompoundPopulationFunction")
    assert cpf.find("itree") is not None


def test_build_tree_prior_bd_skyline_serial():
    config, pm = _make_config(
        TreePriorType.BD_SKYLINE_SERIAL,
        {"birth_rate": {"value": 1.0}, "death_rate": {"value": 0.5}},
    )
    elem, state_nodes, _ = ModelBuilder.build_tree_prior(config, "tree", pm)
    # Uses the bdsky add-on package class
    assert elem.get("spec") == "bdsky.evolution.speciation.BirthDeathSkylineModel"
    assert elem.get("contemp") == "false"
    assert elem.get("conditionOnRoot") == "true"


def test_build_init_random():
    config, pm = _make_config()
    elem = ModelBuilder.build_init(config, pm, "tree")
    assert elem.tag == "init"
    assert "RandomTree" in elem.get("spec")


def test_build_init_upgma():
    config, pm = _make_config()
    config.initialization.tree_type = InitTreeType.UPGMA
    elem = ModelBuilder.build_init(config, pm, "tree")
    assert "ClusterTree" in elem.get("spec")
    assert elem.get("clusterType") == "upgma"


def test_build_operators_basic():
    config, pm = _make_config()
    state_ids = ["birthRate"]
    operators = ModelBuilder.build_operators(config, pm, state_ids)
    assert len(operators) > 0
    # Should have tree operators
    specs = [op.get("spec") for op in operators]
    assert "ScaleOperator" in specs
    assert "beast.base.evolution.operator.Uniform" in specs
    assert "SubtreeSlide" in specs

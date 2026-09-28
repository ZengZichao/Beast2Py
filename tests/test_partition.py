"""Tests for the PartitionManager module."""

import pytest
from beast2py.models import (
    Alignment,
    ClockModelConfig,
    ClockModelType,
    DataType,
    Partition,
    Sequence,
    SiteModelConfig,
)
from beast2py.partition import PartitionManager


def _make_alignment(aln_id="alignment", taxa=None):
    taxa = taxa or ["taxonA", "taxonB", "taxonC"]
    return Alignment(
        id=aln_id,
        data_type=DataType.NUCLEOTIDE,
        sequences=[Sequence(taxon=t, sequence="ACGT") for t in taxa],
    )


def _make_partition(pid="alignment", linked_to=None, taxa=None):
    aln = _make_alignment(pid, taxa=taxa)
    sm = SiteModelConfig(substitution_model={"type": "hky"})
    cm = ClockModelConfig(type=ClockModelType.STRICT, linked_to=linked_to)
    return Partition(id=pid, alignment=aln, site_model=sm, clock_model=cm)


def test_partition_manager_basic():
    p = _make_partition()
    pm = PartitionManager([p])
    assert pm.n_partitions == 1
    assert pm.partition_ids == ["alignment"]


def test_partition_manager_shared_tree():
    p1 = _make_partition("gene1")
    p2 = _make_partition("gene2")
    pm = PartitionManager([p1, p2])
    assert pm.get_tree_id() == "tree"


def test_partition_manager_link_validation():
    p1 = _make_partition("gene1")
    p2 = _make_partition("gene2", linked_to="gene1")
    pm = PartitionManager([p1, p2])
    assert pm.get_clock_ref("gene2") == "gene1.clockModel"
    assert pm.get_clock_ref("gene1") == "gene1.clockModel"


def test_partition_manager_invalid_link():
    p1 = _make_partition("gene1", linked_to="nonexistent")
    with pytest.raises(ValueError, match="non-existent"):
        PartitionManager([p1])


def test_partition_manager_n_taxa():
    p = _make_partition("alignment", taxa=["A", "B", "C", "D"])
    pm = PartitionManager([p])
    assert pm.get_n_taxa() == 4


def test_partition_manager_first_alignment_id():
    p1 = _make_partition("gene1")
    p2 = _make_partition("gene2")
    pm = PartitionManager([p1, p2])
    assert pm.get_first_alignment_id() == "gene1"

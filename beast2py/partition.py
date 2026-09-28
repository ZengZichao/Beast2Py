"""Partition management for multi-partition BEAST2 analyses.

Manages multiple sequence alignments with independent or shared models
(site model, clock model, tree) via link/unlink mechanisms.
"""

from __future__ import annotations

from typing import Dict, List

from .models import Partition


class PartitionManager:
    """Manage multi-partition model link/unlink relationships.

    In BEAST2, multiple partitions can share models through id/idref references.
    This class tracks which models are shared and provides helpers for
    generating unique model IDs.
    """

    def __init__(self, partitions: List[Partition]):
        """Initialize the partition manager.

        Args:
            partitions: List of Partition objects.
        """
        self.partitions: Dict[str, Partition] = {p.id: p for p in partitions}
        self._validate_links()

    def _validate_links(self) -> None:
        """Validate that all linked_to references are valid.

        Links must target an existing partition, must not be self-referential,
        and the target itself must not link further (chained or cyclic links
        would leave every affected partition referencing a model that is never
        generated).
        """
        for p in self.partitions.values():
            for kind, target, getter in (
                ("clock model", p.clock_model.linked_to, lambda q: q.clock_model.linked_to),
                ("site model", p.site_model.linked_to, lambda q: q.site_model.linked_to),
            ):
                if not target:
                    continue
                if target == p.id:
                    raise ValueError(f"Partition '{p.id}' links its {kind} to itself")
                if target not in self.partitions:
                    raise ValueError(
                        f"Partition '{p.id}' links {kind} to " f"non-existent partition '{target}'"
                    )
                nxt = self.partitions[target]
                if getter(nxt):
                    raise ValueError(
                        f"Partition '{p.id}' links its {kind} to '{target}', "
                        f"which itself links to '{getter(nxt)}'; chained or "
                        f"cyclic links are not supported"
                    )

    @property
    def n_partitions(self) -> int:
        """Return the number of partitions."""
        return len(self.partitions)

    @property
    def partition_ids(self) -> List[str]:
        """Return list of partition IDs."""
        return list(self.partitions.keys())

    def get_partition(self, partition_id: str) -> Partition:
        """Get a partition by ID.

        Args:
            partition_id: Partition identifier.

        Returns:
            The Partition object.

        Raises:
            KeyError: If partition_id is not found.
        """
        return self.partitions[partition_id]

    def get_site_ref(self, partition_id: str) -> str:
        """Get the site model ID/reference for a partition.

        Args:
            partition_id: Partition identifier.

        Returns:
            The site model ID (either own or referenced).
        """
        p = self.partitions[partition_id]
        if p.site_model.linked_to:
            return f"{p.site_model.linked_to}.siteModel"
        return f"{partition_id}.siteModel"

    def get_clock_ref(self, partition_id: str) -> str:
        """Get the clock model ID/reference for a partition.

        Args:
            partition_id: Partition identifier.

        Returns:
            The clock model ID (either own or referenced).
        """
        p = self.partitions[partition_id]
        if p.clock_model.linked_to:
            return f"{p.clock_model.linked_to}.clockModel"
        return f"{partition_id}.clockModel"

    def get_first_alignment_id(self) -> str:
        """Get the alignment id of the first partition.

        Returns:
            The *alignment* id (``partition.alignment.id``). This used to
            return the partition id, which only worked because the parser
            forces the two to be equal; constructing partitions with different
            ids through the Python API produced a duplicated ``<data>`` element
            and an ``@`` reference to an alignment that never existed
            .
        """
        first = next(iter(self.partitions.values()))
        return first.alignment.id

    def alignment_ids(self) -> List[str]:
        """Return the distinct alignment ids in partition order.

        Returns:
            Ordered list of unique alignment ids.
        """
        seen: List[str] = []
        for part in self.partitions.values():
            if part.alignment.id not in seen:
                seen.append(part.alignment.id)
        return seen

    def get_tree_id(self) -> str:
        """Get the tree ID.

        Beast2Py assembles every partition on one shared, linked topology, so
        this is always ``"tree"``. ``Partition.tree`` values other than
        ``"shared"`` are rejected by the parser (see
        ``ConfigParser._validate_tree_ref``) instead of being silently ignored
        as they were before was fixed.

        Returns:
            The tree ID string.
        """
        return "tree"

    def get_n_taxa(self) -> int:
        """Get the number of taxa (from the first partition).

        Returns:
            Number of taxa.
        """
        first = list(self.partitions.values())[0]
        return first.alignment.n_taxa

"""Sequence file reading for Beast2Py.

Reads FASTA and NEXUS sequence files using biopython, returning Alignment objects.
Supports both nucleotide and amino acid sequences.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Optional, Set

from .models import Alignment, DataType, Sequence

# IUPAC ambiguity codes and BEAST2 missing-data tokens are legal; anything else
# in an alignment is a data error that used to reach <sequence> unchecked
# . Note that '?' is BEAST2's wildcard token, '.' and '-' are gaps.
NUCLEOTIDE_ALPHABET: Set[str] = set("ACGTUacgtu" "RYSWKMBDHVN" "ryswnkmbdhvn" "-.?*")
AMINOACID_ALPHABET: Set[str] = set("ACDEFGHIKLMNPQRSTVWY" "acdefghiklmnpqrstvwy" "BZUOX*" "-.?")


class AlignmentError(ValueError):
    """Raised when an alignment violates basic comparability invariants."""


def validate_alignment(alignment: Alignment) -> None:
    """Check an alignment for the invariants BEAST2 silently assumes.

    Rejects (rather than warns about): zero sequences, duplicated taxon names,
    ragged sequence lengths, empty or fully-missing sequences, and characters
    outside the declared alphabet. unequal lengths (10/14/14/14)
    and a repeated ``>Homo_sapiens`` header used to be written out unchecked,
    while ``n_sites`` reported only the first sequence's length.

    Args:
        alignment: The parsed alignment to check.

    Raises:
        AlignmentError: If any invariant is violated.
    """
    seqs = alignment.sequences
    if not seqs:
        raise AlignmentError(f"Alignment '{alignment.id}' contains no sequences")

    # Duplicated taxon names
    seen: dict[str, int] = {}
    for idx, seq in enumerate(seqs, start=1):
        if not seq.taxon:
            raise AlignmentError(f"Alignment '{alignment.id}': sequence #{idx} has an empty name")
        seen.setdefault(seq.taxon, 0)
        seen[seq.taxon] += 1
    dupes = sorted(name for name, count in seen.items() if count > 1)
    if dupes:
        raise AlignmentError(
            f"Alignment '{alignment.id}': duplicate taxon name(s): {', '.join(dupes)}. "
            f"BEAST2 requires taxon names to be unique within an alignment."
        )

    # Ragged lengths
    lengths = {len(seq.sequence) for seq in seqs}
    if len(lengths) > 1:
        detail = ", ".join(f"{seq.taxon}={len(seq.sequence)}" for seq in seqs)
        raise AlignmentError(
            f"Alignment '{alignment.id}': sequences are not all the same length "
            f"(sizes: {detail}). Re-align the data before generating XML."
        )
    n_sites = lengths.pop() if lengths else 0
    if n_sites == 0:
        raise AlignmentError(f"Alignment '{alignment.id}': sequences have zero length")

    alphabet = (
        NUCLEOTIDE_ALPHABET if alignment.data_type == DataType.NUCLEOTIDE else AMINOACID_ALPHABET
    )
    observed: Set[str] = set()
    for seq in seqs:
        bad = sorted({ch for ch in seq.sequence if ch not in alphabet})
        if bad:
            raise AlignmentError(
                f"Alignment '{alignment.id}': taxon '{seq.taxon}' contains characters "
                f"illegal for {alignment.data_type.value} data: {''.join(bad)}. "
                f"Check data_type or clean the alignment."
            )
        observed |= {ch for ch in seq.sequence if ch not in "-.?"}
    if not observed:
        raise AlignmentError(
            f"Alignment '{alignment.id}': every site is missing data; nothing can be"
            f" used to compute a likelihood"
        )
    alignment.observed_characters = "".join(sorted(observed))


def observed_characters(alignment: Alignment) -> str:
    """Return the non-gap characters present in an alignment, cached.

    Args:
        alignment: The alignment to inspect.

    Returns:
        A sorted string of observed characters.
    """
    cached = getattr(alignment, "observed_characters", None)
    if cached is not None:
        return str(cached)
    observed: Set[str] = set()
    for seq in alignment.sequences:
        observed |= {ch for ch in seq.sequence if ch not in "-.?"}
    return "".join(sorted(observed))


def detect_data_type(sequences) -> str:
    """Infer "nucleotide" or "aminoacid" from the residues actually present.

    The `quick` subcommand had no way to say which kind of data it was given
    and always read the file as nucleotide, so an amino-acid alignment with an
    amino-acid model could not be expressed through it at all (2).

    Args:
        sequences: Iterable of objects exposing ``sequence``.

    Returns:
        "nucleotide" when every residue is a nucleotide/alphabet placeholder,
        otherwise "aminoacid".
    """
    residues = set()
    for seq in sequences:
        residues.update(str(seq.sequence or ""))
    # The nucleotide alphabet is a proper subset of the amino-acid one, so it
    # has to be tested first or proteins would never be reported.
    usable = residues - set(" \t\r\n")
    if not usable or usable <= NUCLEOTIDE_ALPHABET:
        return "nucleotide"
    return "aminoacid"


class SequenceReader:
    """Read sequence data from FASTA / NEXUS files."""

    @staticmethod
    def read(
        file_path: str,
        format: str = "auto",
        data_type: str = "nucleotide",
        filter: Optional[str] = None,
        alignment_id: str = "alignment",
    ) -> Alignment:
        """Read a sequence file and return an Alignment object.

        Args:
            file_path: Path to the sequence file.
            format: File format ("fasta", "nexus", or "auto" for auto-detection).
            data_type: "nucleotide", "aminoacid", or "auto" to infer it from
                the residues actually present.
            filter: Optional FilteredAlignment filter expression.
            alignment_id: ID for the alignment.

        Returns:
            An Alignment object.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the format cannot be determined or parsed.
        """
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"Sequence file not found: {file_path}")

        # Auto-detect format
        if format == "auto":
            ext = path.suffix.lower()
            if ext in (".fas", ".fasta", ".fa", ".fna", ".ffn"):
                format = "fasta"
            elif ext in (".nex", ".nexus", ".nxs"):
                format = "nexus"
            else:
                # Try to detect from file content
                with open(file_path, "r", encoding="utf-8") as f:
                    first_line = f.read(100).strip()
                if first_line.startswith("#NEXUS") or first_line.startswith("#nexus"):
                    format = "nexus"
                else:
                    format = "fasta"

        if data_type == "auto":
            # Peek at the residues so a caller that does not want to declare the
            # kind of data can still route an amino-acid alignment to an
            # amino-acid model (2).
            peek = (
                SequenceReader._read_fasta(str(path))
                if format == "fasta"
                else SequenceReader._read_nexus(str(path))
            )
            data_type = detect_data_type(peek)
        dt = DataType(data_type)

        if format == "fasta":
            sequences = SequenceReader._read_fasta(file_path)
        elif format == "nexus":
            sequences = SequenceReader._read_nexus(file_path)
        else:
            raise ValueError(f"Unsupported sequence format: {format}")

        return Alignment(
            id=alignment_id,
            data_type=dt,
            sequences=sequences,
            filter=filter,
        )

    @staticmethod
    def _read_fasta(file_path: str) -> List[Sequence]:
        """Parse a FASTA file.

        Uses biopython's SeqIO for robust parsing.

        Args:
            file_path: Path to the FASTA file.

        Returns:
            List of Sequence objects.
        """
        try:
            from Bio import SeqIO
        except ImportError:
            # Fallback: simple parser
            return SequenceReader._read_fasta_simple(file_path)

        sequences = []
        headers = []
        for record in SeqIO.parse(file_path, "fasta"):
            sequences.append(
                Sequence(
                    taxon=str(record.id),
                    sequence=str(record.seq).upper(),
                )
            )
            headers.append(str(getattr(record, "description", "") or "").strip()
                           or str(record.id))
        SequenceReader._check_fasta_headers(sequences, headers, file_path)
        return sequences

    @staticmethod
    def _check_fasta_headers(
        sequences: List[Sequence], headers: List[str], file_path: str
    ) -> None:
        """Report FASTA headers that a taxon name was truncated from.

        The FASTA id is the first whitespace-delimited token of the header and
        the rest is a description; Biopython's parser does the same, and a taxon
        name containing spaces cannot be emitted into Newick without quoting.
        The name is therefore kept short, but a collision *caused by that
        truncation* must name the real headers rather than accusing the author
        of duplicating names that were never duplicated (3).

        Args:
            sequences: Parsed sequences.
            headers: The full header text each sequence came from.
            file_path: Source file, for messages.

        Raises:
            AlignmentError: If distinct headers collapse to one taxon name.
        """
        by_name: Dict[str, List[str]] = {}
        for seq, header in zip(sequences, headers):
            by_name.setdefault(seq.taxon, []).append(header)
        collapsed = {
            name: sorted(set(hdrs))
            for name, hdrs in by_name.items()
            if len(set(hdrs)) > 1
        }
        if collapsed:
            detail = "; ".join(
                f"{name!r} <- {', '.join(repr(h) for h in hdrs)}"
                for name, hdrs in sorted(collapsed.items())
            )
            raise AlignmentError(
                f"FASTA headers differ but collapse to the same taxon name once the "
                f"description after the first whitespace is dropped: {detail} "
                f"({file_path}). Rename the records so their first word is already "
                f"unique; BEAST2 requires unique taxon names."
            )
        truncated = [
            (seq.taxon, header)
            for seq, header in zip(sequences, headers)
            if header and header != seq.taxon
        ]
        if truncated:
            shown = ", ".join(f"{t!r} (from {h!r})" for t, h in truncated[:5])
            more = f" (+{len(truncated) - 5} more)" if len(truncated) > 5 else ""
            warnings.warn(
                f"{len(truncated)} FASTA header(s) carry a description that is not "
                f"part of the taxon name; names were truncated to their first word: "
                f"{shown}{more}. Logs and trees will use the shortened names.",
                stacklevel=3,
            )

    @staticmethod
    def _read_fasta_simple(file_path: str) -> List[Sequence]:
        """Simple FASTA parser without biopython dependency.

        Args:
            file_path: Path to the FASTA file.

        Returns:
            List of Sequence objects.
        """
        sequences = []
        # The FASTA id is the first whitespace-delimited token of the header and
        # the remainder is a description; Biopython's own parser does the same,
        # and a taxon name containing spaces cannot be written into Newick
        # without quoting. The name is therefore kept as-is, but the original
        # header is remembered so that a collision caused by *truncation* can
        # name the real headers instead of blaming the author's data
        # (3).
        source_headers: List[str] = []
        current_taxon = None
        current_header = None
        current_seq_parts: List[str] = []

        def flush() -> None:
            if current_taxon is not None:
                sequences.append(
                    Sequence(
                        taxon=current_taxon,
                        sequence="".join(current_seq_parts).upper(),
                    )
                )
                source_headers.append(current_header or current_taxon)

        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    flush()
                    current_header = line[1:].strip()
                    tokens = current_header.split()
                    if not tokens:
                        raise AlignmentError(
                            f"FASTA record with an empty header: '{line}' in {file_path}"
                        )
                    current_taxon = tokens[0]
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line)

            flush()

        SequenceReader._check_fasta_headers(sequences, source_headers, file_path)
        return sequences

    @staticmethod
    def _read_nexus(file_path: str) -> List[Sequence]:
        """Parse a NEXUS file.

        Uses biopython's AlignIO or Nexus module for robust parsing.

        Args:
            file_path: Path to the NEXUS file.

        Returns:
            List of Sequence objects.
        """
        try:
            from Bio import AlignIO
        except ImportError:
            return SequenceReader._read_nexus_simple(file_path)

        sequences = []
        try:
            alignment = AlignIO.read(file_path, "nexus")
            for record in alignment:
                sequences.append(
                    Sequence(
                        taxon=str(record.id),
                        sequence=str(record.seq).upper(),
                    )
                )
        except Exception:
            # Try alternative parsing with Nexus module
            try:
                from Bio.Nexus import Nexus

                nex = Nexus.Nexus()
                nex.read(file_path)
                for taxon_name, seq_data in nex.matrix.items():
                    sequences.append(
                        Sequence(
                            taxon=taxon_name,
                            sequence=str(seq_data).upper(),
                        )
                    )
            except Exception as e:
                raise ValueError(f"Failed to parse NEXUS file {file_path}: {e}")

        return sequences

    @staticmethod
    def _read_nexus_simple(file_path: str) -> List[Sequence]:
        """Simple NEXUS parser without biopython dependency.

        Handles basic non-interleaved NEXUS data blocks.

        Args:
            file_path: Path to the NEXUS file.

        Returns:
            List of Sequence objects.
        """
        sequences = []
        in_data_block = False
        in_matrix = False
        seq_map: dict[str, List[str]] = {}
        taxon_order: List[str] = []

        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                stripped = line.strip()

                if not in_data_block:
                    if stripped.lower().startswith("begin data;") or stripped.lower().startswith(
                        "begin characters;"
                    ):
                        in_data_block = True
                    continue

                if stripped.lower() == "end;":
                    break

                if stripped.lower().startswith("matrix"):
                    in_matrix = True
                    continue

                if in_matrix and stripped and stripped != ";":
                    # Parse taxon sequence line
                    parts = stripped.split(None, 1)
                    if len(parts) >= 2:
                        taxon = parts[0]
                        seq = parts[1].strip().rstrip(";")
                        if taxon not in seq_map:
                            seq_map[taxon] = []
                            taxon_order.append(taxon)
                        seq_map[taxon].append(seq)

        for taxon in taxon_order:
            sequences.append(
                Sequence(
                    taxon=taxon,
                    sequence="".join(seq_map[taxon]).upper(),
                )
            )

        return sequences

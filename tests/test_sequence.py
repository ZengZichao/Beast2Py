"""Tests for the sequence module."""

import tempfile
from pathlib import Path
import pytest

from beast2py.sequence import SequenceReader
from beast2py.models import DataType

from .conftest import TEST_FASTA


class TestSequenceReader:
    """Test the SequenceReader class."""

    def setup_method(self, method):
        """Create a temp directory for each test."""
        self.tmpdir = Path(tempfile.mkdtemp())
        self.fasta_path = self.tmpdir / "test.fasta"
        self.fasta_path.write_text(TEST_FASTA)

    def test_read_fasta(self):
        """Test reading a FASTA file."""
        alignment = SequenceReader.read(
            self.fasta_path,
            format="auto",
            data_type="nucleotide",
            alignment_id="test_aln",
        )

        assert alignment.id == "test_aln"
        assert alignment.data_type == DataType.NUCLEOTIDE
        assert alignment.n_taxa == 4
        assert alignment.n_sites == 36
        assert "taxonA" in alignment.taxa_names
        assert "taxonB" in alignment.taxa_names
        assert "taxonC" in alignment.taxa_names
        assert "taxonD" in alignment.taxa_names

    def test_read_fasta_explicit_format(self):
        """Test reading a FASTA file with explicit format."""
        alignment = SequenceReader.read(
            self.fasta_path,
            format="fasta",
            data_type="nucleotide",
        )
        assert alignment.n_taxa == 4

    def test_read_nonexistent_file(self):
        """Test that reading a non-existent file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            SequenceReader.read("/nonexistent/file.fasta")

    def test_read_empty_fasta(self):
        """Test reading an empty FASTA file."""
        empty_path = self.tmpdir / "empty.fasta"
        empty_path.write_text("")

        alignment = SequenceReader.read(
            empty_path,
            format="fasta",
            data_type="nucleotide",
        )
        assert alignment.n_taxa == 0

    def test_read_fasta_uppercase(self):
        """Test that sequences are converted to uppercase."""
        lowercase_fasta = ">test_taxon\nacgtacgt\n"
        lower_path = self.tmpdir / "lower.fasta"
        lower_path.write_text(lowercase_fasta)

        alignment = SequenceReader.read(
            lower_path,
            format="fasta",
            data_type="nucleotide",
        )
        assert alignment.sequences[0].sequence == "ACGTACGT"


import gzip

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fasta_split import split_fasta
from genomics_toolkit.operations.fasta_stats import read_fasta


@pytest.fixture
def sample_fasta(tmp_path):
    fasta = tmp_path / "sample.fasta"
    fasta.write_text(
        ">seq1\nACGT\n"
        ">seq2\nGGCC\n"
        ">seq3\nATAT\n"
        ">seq4\nTTAA\n"
        ">seq5\nCCCC\n"
    )
    return fasta


def test_split_into_chunks(sample_fasta, tmp_path):
    output_dir = tmp_path / "chunks"

    result = split_fasta(
        sample_fasta,
        output_dir,
        records_per_file=2,
    )

    files = sorted(output_dir.glob("*.fasta"))

    assert len(files) == 3
    assert result["total_records"] == 5
    assert result["file_count"] == 3

    assert list(read_fasta(files[0])) == [
        ("seq1", "ACGT"),
        ("seq2", "GGCC"),
    ]

    assert list(read_fasta(files[1])) == [
        ("seq3", "ATAT"),
        ("seq4", "TTAA"),
    ]

    assert list(read_fasta(files[2])) == [
        ("seq5", "CCCC"),
    ]


def test_gzipped_output(sample_fasta, tmp_path):
    output_dir = tmp_path / "chunks"

    split_fasta(
        sample_fasta,
        output_dir,
        records_per_file=2,
        gzip_output=True,
    )

    files = sorted(output_dir.glob("*.fasta.gz"))

    assert len(files) == 3

    with gzip.open(files[0], "rt", encoding="utf-8") as handle:
        assert ">seq1" in handle.read()


def test_gzipped_input(sample_fasta, tmp_path):
    compressed = tmp_path / "sample.fasta.gz"

    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(sample_fasta.read_text())

    output_dir = tmp_path / "chunks"

    result = split_fasta(
        compressed,
        output_dir,
        records_per_file=2,
    )

    assert result["total_records"] == 5


def test_invalid_chunk_size(sample_fasta, tmp_path):
    with pytest.raises(ValueError, match="at least 1"):
        split_fasta(
            sample_fasta,
            tmp_path / "chunks",
            records_per_file=0,
        )


def test_existing_output_directory_rejected(sample_fasta, tmp_path):
    output_dir = tmp_path / "chunks"
    output_dir.mkdir()

    with pytest.raises(FileExistsError):
        split_fasta(sample_fasta, output_dir)


def test_empty_fasta(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("")

    output_dir = tmp_path / "chunks"
    result = split_fasta(fasta, output_dir)

    assert result["total_records"] == 0
    assert result["file_count"] == 0
    assert output_dir.exists()
    assert list(output_dir.iterdir()) == []


def test_malformed_fasta_cleans_up(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text("ACGT\n>seq1\nGGCC\n")

    output_dir = tmp_path / "chunks"

    with pytest.raises(ValueError):
        split_fasta(fasta, output_dir)

    assert not output_dir.exists()
    assert list(tmp_path.glob(".chunks.*")) == []


def test_cli_split(sample_fasta, tmp_path):
    output_dir = tmp_path / "chunks"

    exit_code = main([
        "fasta-split",
        "--input", str(sample_fasta),
        "--output-dir", str(output_dir),
        "--records-per-file", "2",
    ])

    assert exit_code == 0
    assert len(list(output_dir.glob("*.fasta"))) == 3

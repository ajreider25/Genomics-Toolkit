
import gzip

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fastq.fastq_convert import (
    fastq_to_fasta,
)
from genomics_toolkit.operations.fasta.fasta_stats import read_fasta


@pytest.fixture
def sample_fastq(tmp_path):
    path = tmp_path / "reads.fastq"
    path.write_text(
        "@read1\nACGTACGT\n+\nIIIIIIII\n"
        "@read2\nGGCCNNAT\n+\n55555555\n",
        encoding="utf-8",
    )
    return path


def test_basic_conversion(sample_fastq, tmp_path):
    output = tmp_path / "reads.fasta"

    result = fastq_to_fasta(sample_fastq, output)

    assert result["record_count"] == 2
    assert result["total_bases"] == 16
    assert list(read_fasta(output)) == [
        ("read1", "ACGTACGT"),
        ("read2", "GGCCNNAT"),
    ]


def test_quality_scores_removed(sample_fastq, tmp_path):
    output = tmp_path / "reads.fasta"

    fastq_to_fasta(sample_fastq, output)

    content = output.read_text()

    assert "IIIIIIII" not in content
    assert "55555555" not in content
    assert ">read1" in content


def test_line_wrapping(sample_fastq, tmp_path):
    output = tmp_path / "reads.fasta"

    fastq_to_fasta(sample_fastq, output, line_width=4)

    assert output.read_text() == (
        ">read1\n"
        "ACGT\n"
        "ACGT\n"
        ">read2\n"
        "GGCC\n"
        "NNAT\n"
    )


def test_gzip_input(sample_fastq, tmp_path):
    compressed = tmp_path / "reads.fastq.gz"

    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(sample_fastq.read_text())

    output = tmp_path / "reads.fasta"

    result = fastq_to_fasta(compressed, output)

    assert result["record_count"] == 2


def test_gzip_output(sample_fastq, tmp_path):
    output = tmp_path / "reads.fasta.gz"

    fastq_to_fasta(sample_fastq, output)

    with gzip.open(output, "rt", encoding="utf-8") as handle:
        content = handle.read()

    assert ">read1" in content
    assert ">read2" in content


def test_empty_input(tmp_path):
    source = tmp_path / "empty.fastq"
    output = tmp_path / "empty.fasta"

    source.write_text("")

    result = fastq_to_fasta(source, output)

    assert result["record_count"] == 0
    assert result["total_bases"] == 0
    assert output.read_text() == ""


def test_invalid_line_width(sample_fastq, tmp_path):
    with pytest.raises(ValueError):
        fastq_to_fasta(
            sample_fastq,
            tmp_path / "output.fasta",
            line_width=0,
        )


def test_same_input_output(sample_fastq):
    with pytest.raises(ValueError):
        fastq_to_fasta(sample_fastq, sample_fastq)


def test_missing_input(tmp_path):
    with pytest.raises(FileNotFoundError):
        fastq_to_fasta(
            tmp_path / "missing.fastq",
            tmp_path / "output.fasta",
        )


def test_malformed_input_no_partial_output(tmp_path):
    source = tmp_path / "bad.fastq"
    output = tmp_path / "output.fasta"

    source.write_text(
        "@good\nACGT\n+\nIIII\n"
        "@bad\nACGT\n+\nIII\n"
    )

    with pytest.raises(ValueError):
        fastq_to_fasta(source, output)

    assert not output.exists()


def test_cli_conversion(sample_fastq, tmp_path):
    output = tmp_path / "converted.fasta"

    exit_code = main([
        "fastq-convert",
        "--input", str(sample_fastq),
        "--output", str(output),
    ])

    assert exit_code == 0
    assert list(read_fasta(output)) == [
        ("read1", "ACGTACGT"),
        ("read2", "GGCCNNAT"),
    ]

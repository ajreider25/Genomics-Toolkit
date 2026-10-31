import gzip

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fastq.fastq_trim import (
    trim_read,
    trim_fastq,
)
from genomics_toolkit.operations.fastq.fastq_stats import read_fastq


@pytest.fixture
def sample_fastq(tmp_path):
    path = tmp_path / "reads.fastq"
    path.write_text(
        "@read1\nACGTACGT\n+\nIIIIII!!\n"
        "@read2\nGGCC\n+\nIIII\n"
        "@read3\nATAT\n+\n!!!!\n",
        encoding="utf-8",
    )
    return path


def test_trim_read():
    sequence, quality = trim_read("ACGTAC", "IIII!!", 20)

    assert sequence == "ACGT"
    assert quality == "IIII"


def test_trim_read_no_trimming():
    assert trim_read("ACGT", "IIII", 20) == ("ACGT", "IIII")


def test_trim_read_all_low_quality():
    assert trim_read("ACGT", "!!!!", 20) == ("", "")


def test_trim_read_threshold_inclusive():
    assert trim_read("ACGT", "5555", 20) == ("ACGT", "5555")


def test_trim_read_invalid_quality():
    with pytest.raises(ValueError):
        trim_read("ACGT", "IIII", 100)


def test_trim_read_mismatched_lengths():
    with pytest.raises(ValueError):
        trim_read("ACGT", "III", 20)


def test_trim_fastq_basic(sample_fastq, tmp_path):
    output = tmp_path / "trimmed.fastq"

    result = trim_fastq(
        sample_fastq,
        output,
        min_quality=20,
        min_length=1,
    )

    assert result["total_reads"] == 3
    assert result["kept_reads"] == 2
    assert result["trimmed_reads"] == 2
    assert result["discarded_reads"] == 1
    assert result["bases_removed"] == 6

    records = list(read_fastq(output))

    assert records == [
        ("read1", "ACGTAC", "IIIIII"),
        ("read2", "GGCC", "IIII"),
    ]


def test_trim_fastq_min_length(sample_fastq, tmp_path):
    output = tmp_path / "trimmed.fastq"

    result = trim_fastq(
        sample_fastq,
        output,
        min_quality=20,
        min_length=7,
    )

    assert result["kept_reads"] == 0
    assert result["discarded_reads"] == 3


def test_trim_fastq_gzip_output(sample_fastq, tmp_path):
    output = tmp_path / "trimmed.fastq.gz"

    trim_fastq(
        sample_fastq,
        output,
        min_quality=20,
        min_length=1,
    )

    with gzip.open(output, "rt", encoding="utf-8") as handle:
        content = handle.read()

    assert "@read1" in content
    assert "@read2" in content
    assert "@read3" not in content


def test_trim_fastq_gzip_input(sample_fastq, tmp_path):
    compressed = tmp_path / "reads.fastq.gz"

    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(sample_fastq.read_text())

    output = tmp_path / "trimmed.fastq"

    result = trim_fastq(
        compressed,
        output,
        min_quality=20,
        min_length=1,
    )

    assert result["kept_reads"] == 2


def test_trim_fastq_same_input_output(sample_fastq):
    with pytest.raises(ValueError, match="different"):
        trim_fastq(sample_fastq, sample_fastq)


def test_trim_fastq_invalid_min_length(sample_fastq, tmp_path):
    with pytest.raises(ValueError):
        trim_fastq(
            sample_fastq,
            tmp_path / "output.fastq",
            min_length=-1,
        )


def test_trim_fastq_malformed_input(tmp_path):
    source = tmp_path / "bad.fastq"
    output = tmp_path / "trimmed.fastq"

    source.write_text(
        "@good\nACGT\n+\nIIII\n"
        "@bad\nACGT\n+\nIII\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError):
        trim_fastq(source, output)

    assert not output.exists()


def test_trim_fastq_empty_file(tmp_path):
    source = tmp_path / "empty.fastq"
    output = tmp_path / "trimmed.fastq"

    source.write_text("")

    result = trim_fastq(source, output)

    assert result["total_reads"] == 0
    assert result["kept_reads"] == 0
    assert output.read_text() == ""


def test_cli_fastq_trim(sample_fastq, tmp_path):
    output = tmp_path / "trimmed.fastq"

    exit_code = main([
        "fastq-trim",
        "--input", str(sample_fastq),
        "--output", str(output),
        "--min-quality", "20",
        "--min-length", "1",
    ])

    assert exit_code == 0
    assert len(list(read_fastq(output))) == 2
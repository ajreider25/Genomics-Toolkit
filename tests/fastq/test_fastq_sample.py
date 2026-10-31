
import gzip

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fastq.fastq_sample import sample_fastq
from genomics_toolkit.operations.fastq.fastq_stats import read_fastq


@pytest.fixture
def sample_reads(tmp_path):
    path = tmp_path / "reads.fastq"

    with path.open("w", encoding="utf-8") as handle:
        for index in range(20):
            handle.write(
                f"@read{index}\nACGT\n+\nIIII\n"
            )

    return path


def test_exact_count(sample_reads, tmp_path):
    output = tmp_path / "sampled.fastq"

    result = sample_fastq(sample_reads, output, count=5)

    assert result["total_reads"] == 20
    assert result["sampled_reads"] == 5
    assert result["discarded_reads"] == 15
    assert len(list(read_fastq(output))) == 5


def test_count_larger_than_input(sample_reads, tmp_path):
    output = tmp_path / "sampled.fastq"

    result = sample_fastq(sample_reads, output, count=100)

    assert result["sampled_reads"] == 20


def test_count_zero(sample_reads, tmp_path):
    output = tmp_path / "sampled.fastq"

    result = sample_fastq(sample_reads, output, count=0)

    assert result["sampled_reads"] == 0
    assert output.read_text() == ""


def test_fraction_zero(sample_reads, tmp_path):
    output = tmp_path / "sampled.fastq"

    result = sample_fastq(sample_reads, output, fraction=0)

    assert result["sampled_reads"] == 0


def test_fraction_one(sample_reads, tmp_path):
    output = tmp_path / "sampled.fastq"

    result = sample_fastq(sample_reads, output, fraction=1)

    assert result["sampled_reads"] == 20


def test_fraction_sampling(sample_reads, tmp_path):
    output = tmp_path / "sampled.fastq"

    result = sample_fastq(
        sample_reads,
        output,
        fraction=0.5,
        seed=42,
    )

    assert 0 < result["sampled_reads"] < 20


def test_reproducibility(sample_reads, tmp_path):
    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"

    sample_fastq(sample_reads, first, count=7, seed=123)
    sample_fastq(sample_reads, second, count=7, seed=123)

    assert first.read_bytes() == second.read_bytes()


def test_preserves_input_order(sample_reads, tmp_path):
    output = tmp_path / "sampled.fastq"

    sample_fastq(sample_reads, output, count=8, seed=42)

    identifiers = [
        int(read_id.replace("read", ""))
        for read_id, _, _ in read_fastq(output)
    ]

    assert identifiers == sorted(identifiers)


def test_gzip_input_output(sample_reads, tmp_path):
    compressed = tmp_path / "reads.fastq.gz"
    output = tmp_path / "sampled.fastq.gz"

    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(sample_reads.read_text())

    result = sample_fastq(compressed, output, count=5)

    assert result["sampled_reads"] == 5
    assert len(list(read_fastq(output))) == 5


def test_invalid_fraction(sample_reads, tmp_path):
    with pytest.raises(ValueError):
        sample_fastq(
            sample_reads,
            tmp_path / "sampled.fastq",
            fraction=1.5,
        )


def test_invalid_count(sample_reads, tmp_path):
    with pytest.raises(ValueError):
        sample_fastq(
            sample_reads,
            tmp_path / "sampled.fastq",
            count=-1,
        )


def test_both_sampling_methods(sample_reads, tmp_path):
    with pytest.raises(ValueError):
        sample_fastq(
            sample_reads,
            tmp_path / "sampled.fastq",
            count=5,
            fraction=0.5,
        )


def test_no_sampling_method(sample_reads, tmp_path):
    with pytest.raises(ValueError):
        sample_fastq(
            sample_reads,
            tmp_path / "sampled.fastq",
        )


def test_same_input_output(sample_reads):
    with pytest.raises(ValueError):
        sample_fastq(sample_reads, sample_reads, count=5)


def test_malformed_input_leaves_no_output(tmp_path):
    source = tmp_path / "bad.fastq"
    output = tmp_path / "sampled.fastq"

    source.write_text(
        "@good\nACGT\n+\nIIII\n"
        "@bad\nACGT\n+\nIII\n"
    )

    with pytest.raises(ValueError):
        sample_fastq(source, output, count=1)

    assert not output.exists()


def test_empty_input(tmp_path):
    source = tmp_path / "empty.fastq"
    output = tmp_path / "sampled.fastq"

    source.write_text("")

    result = sample_fastq(source, output, count=5)

    assert result["total_reads"] == 0
    assert result["sampled_reads"] == 0
    assert output.read_text() == ""


def test_cli_sample(sample_reads, tmp_path):
    output = tmp_path / "sampled.fastq"

    exit_code = main([
        "fastq-sample",
        "--input", str(sample_reads),
        "--output", str(output),
        "--count", "5",
        "--seed", "123",
    ])

    assert exit_code == 0
    assert len(list(read_fastq(output))) == 5

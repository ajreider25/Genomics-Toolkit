
import gzip

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fastq.fastq_filter import filter_fastq
from genomics_toolkit.operations.fastq.fastq_stats import read_fastq


@pytest.fixture
def sample_fastq(tmp_path):
    path = tmp_path / "reads.fastq"
    path.write_text(
        "@high\nACGTACGT\n+\nIIIIIIII\n"
        "@low\nACGT\n+\n!!!!\n"
        "@ambiguous\nNNNNACGT\n+\n????????\n",
        encoding="utf-8",
    )
    return path


def test_filter_by_min_length(sample_fastq, tmp_path):
    output = tmp_path / "filtered.fastq"

    result = filter_fastq(sample_fastq, output, min_length=8)

    assert result["total_reads"] == 3
    assert result["kept_reads"] == 2
    assert result["filtered_reads"] == 1


def test_filter_by_max_length(sample_fastq, tmp_path):
    output = tmp_path / "filtered.fastq"

    result = filter_fastq(sample_fastq, output, max_length=4)

    assert result["kept_reads"] == 1
    assert [r[0] for r in read_fastq(output)] == ["low"]


def test_filter_by_quality(sample_fastq, tmp_path):
    output = tmp_path / "filtered.fastq"

    result = filter_fastq(sample_fastq, output, min_quality=30)

    assert result["kept_reads"] == 2
    assert [r[0] for r in read_fastq(output)] == [
        "high", "ambiguous"
    ]


def test_filter_by_n_percent(sample_fastq, tmp_path):
    output = tmp_path / "filtered.fastq"

    result = filter_fastq(
        sample_fastq, output, max_n_percent=25
    )

    assert result["kept_reads"] == 2


def test_combined_filters(sample_fastq, tmp_path):
    output = tmp_path / "filtered.fastq"

    result = filter_fastq(
        sample_fastq,
        output,
        min_length=8,
        min_quality=30,
        max_n_percent=25,
    )

    assert result["kept_reads"] == 1
    assert [r[0] for r in read_fastq(output)] == ["high"]


def test_no_filters_keeps_all_reads(sample_fastq, tmp_path):
    output = tmp_path / "filtered.fastq"

    result = filter_fastq(sample_fastq, output)

    assert result["kept_reads"] == 3


def test_gzip_output(sample_fastq, tmp_path):
    output = tmp_path / "filtered.fastq.gz"

    filter_fastq(sample_fastq, output, min_length=8)

    with gzip.open(output, "rt", encoding="utf-8") as handle:
        contents = handle.read()

    assert "@high" in contents
    assert "@low" not in contents


def test_invalid_length_range(sample_fastq, tmp_path):
    with pytest.raises(ValueError):
        filter_fastq(
            sample_fastq,
            tmp_path / "filtered.fastq",
            min_length=100,
            max_length=10,
        )


def test_same_input_output(sample_fastq):
    with pytest.raises(ValueError):
        filter_fastq(sample_fastq, sample_fastq)


def test_malformed_input_leaves_no_output(tmp_path):
    source = tmp_path / "bad.fastq"
    output = tmp_path / "filtered.fastq"

    source.write_text(
        "@good\nACGT\n+\nIIII\n"
        "@bad\nACGT\n+\nIII\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError):
        filter_fastq(source, output)

    assert not output.exists()


def test_cli_filter(sample_fastq, tmp_path):
    output = tmp_path / "filtered.fastq"

    exit_code = main([
        "fastq-filter",
        "--input", str(sample_fastq),
        "--output", str(output),
        "--min-quality", "30",
    ])

    assert exit_code == 0
    assert output.exists()
    assert len(list(read_fastq(output))) == 2

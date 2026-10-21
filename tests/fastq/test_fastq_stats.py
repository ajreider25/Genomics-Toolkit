
import gzip
import json

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fastq.fastq_stats import (
    read_fastq,
    calculate_fastq_stats,
    write_fastq_stats,
)


@pytest.fixture
def sample_fastq(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "@read1\nACGT\n+\nIIII\n"
        "@read2\nGGNN\n+\n5555\n"
    )
    return fastq


def test_read_fastq(sample_fastq):
    records = list(read_fastq(sample_fastq))

    assert records == [
        ("read1", "ACGT", "IIII"),
        ("read2", "GGNN", "5555"),
    ]


def test_basic_stats(sample_fastq):
    stats = calculate_fastq_stats(sample_fastq)

    assert stats["read_count"] == 2
    assert stats["total_bases"] == 8
    assert stats["min_read_length"] == 4
    assert stats["max_read_length"] == 4
    assert stats["mean_read_length"] == 4.0

    # ACGT + GGNN: 4 GC bases among 6 called bases.
    assert stats["gc_percent"] == pytest.approx(66.6667, rel=1e-4)


def test_quality_stats(sample_fastq):
    stats = calculate_fastq_stats(sample_fastq)

    # 'I' = Q40, '5' = Q20 under Phred+33.
    assert stats["mean_quality"] == 30.0
    assert stats["q20_percent"] == 100.0
    assert stats["q30_percent"] == 50.0


def test_per_base_quality(sample_fastq):
    stats = calculate_fastq_stats(sample_fastq)

    assert len(stats["per_base_quality"]) == 4

    for position in stats["per_base_quality"]:
        assert position["mean_quality"] == 30.0
        assert position["read_count"] == 2


def test_variable_read_lengths(tmp_path):
    fastq = tmp_path / "variable.fastq"
    fastq.write_text(
        "@read1\nAC\n+\nII\n"
        "@read2\nACGT\n+\n5555\n"
    )

    stats = calculate_fastq_stats(fastq)

    assert stats["read_count"] == 2
    assert stats["min_read_length"] == 2
    assert stats["max_read_length"] == 4
    assert stats["mean_read_length"] == 3.0

    assert [
        position["read_count"]
        for position in stats["per_base_quality"]
    ] == [2, 2, 1, 1]


def test_empty_fastq(tmp_path):
    fastq = tmp_path / "empty.fastq"
    fastq.write_text("")

    stats = calculate_fastq_stats(fastq)

    assert stats["read_count"] == 0
    assert stats["total_bases"] == 0
    assert stats["mean_quality"] == 0.0
    assert stats["per_base_quality"] == []


def test_gzipped_fastq(sample_fastq, tmp_path):
    compressed = tmp_path / "reads.fastq.gz"

    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(sample_fastq.read_text())

    stats = calculate_fastq_stats(compressed)

    assert stats["read_count"] == 2
    assert stats["total_bases"] == 8


def test_invalid_quality_length(tmp_path):
    fastq = tmp_path / "invalid.fastq"
    fastq.write_text("@read1\nACGT\n+\nIII\n")

    with pytest.raises(ValueError, match="lengths do not match"):
        calculate_fastq_stats(fastq)


def test_invalid_base(tmp_path):
    fastq = tmp_path / "invalid.fastq"
    fastq.write_text("@read1\nACGX\n+\nIIII\n")

    with pytest.raises(ValueError, match="invalid nucleotide"):
        calculate_fastq_stats(fastq)


def test_json_output(sample_fastq, tmp_path):
    stats = calculate_fastq_stats(sample_fastq)
    output = tmp_path / "stats.json"

    write_fastq_stats(stats, output, input_path=sample_fastq)

    saved = json.loads(output.read_text())

    assert saved["read_count"] == 2
    assert saved["mean_quality"] == 30.0


def test_cannot_overwrite_input(sample_fastq):
    stats = calculate_fastq_stats(sample_fastq)

    with pytest.raises(ValueError):
        write_fastq_stats(
            stats,
            sample_fastq,
            input_path=sample_fastq,
        )


def test_cli_fastq_stats(sample_fastq, tmp_path):
    output = tmp_path / "stats.json"

    exit_code = main([
        "fastq-stats",
        "--input", str(sample_fastq),
        "--output", str(output),
    ])

    assert exit_code == 0
    assert output.exists()

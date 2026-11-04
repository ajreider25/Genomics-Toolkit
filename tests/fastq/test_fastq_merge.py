
import gzip

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fastq.fastq_merge import merge_fastq
from genomics_toolkit.operations.fastq.fastq_stats import read_fastq


@pytest.fixture
def input_files(tmp_path):
    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"

    first.write_text(
        "@read1\nACGT\n+\nIIII\n"
        "@read2\nGGCC\n+\n5555\n",
        encoding="utf-8",
    )

    second.write_text(
        "@read3\nATAT\n+\n????\n",
        encoding="utf-8",
    )

    return first, second


def test_basic_merge(input_files, tmp_path):
    output = tmp_path / "merged.fastq"

    result = merge_fastq(list(input_files), output)

    assert result["input_files"] == 2
    assert result["total_reads"] == 3
    assert result["total_bases"] == 12

    assert list(read_fastq(output)) == [
        ("read1", "ACGT", "IIII"),
        ("read2", "GGCC", "5555"),
        ("read3", "ATAT", "????"),
    ]


def test_preserves_input_order(input_files, tmp_path):
    output = tmp_path / "merged.fastq"

    merge_fastq(list(reversed(input_files)), output)

    identifiers = [
        read_id for read_id, _, _ in read_fastq(output)
    ]

    assert identifiers == ["read3", "read1", "read2"]


def test_gzip_input(input_files, tmp_path):
    first, second = input_files
    compressed = tmp_path / "first.fastq.gz"

    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(first.read_text())

    output = tmp_path / "merged.fastq"

    result = merge_fastq([compressed, second], output)

    assert result["total_reads"] == 3


def test_gzip_output(input_files, tmp_path):
    output = tmp_path / "merged.fastq.gz"

    merge_fastq(list(input_files), output)

    with gzip.open(output, "rt", encoding="utf-8") as handle:
        content = handle.read()

    assert "@read1" in content
    assert "@read3" in content


def test_duplicate_identifiers_allowed(tmp_path):
    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    output = tmp_path / "merged.fastq"

    first.write_text("@same\nACGT\n+\nIIII\n")
    second.write_text("@same\nGGCC\n+\n5555\n")

    result = merge_fastq([first, second], output)

    assert result["total_reads"] == 2


def test_empty_input_list(tmp_path):
    with pytest.raises(ValueError):
        merge_fastq([], tmp_path / "merged.fastq")


def test_missing_input(input_files, tmp_path):
    with pytest.raises(FileNotFoundError):
        merge_fastq(
            [input_files[0], tmp_path / "missing.fastq"],
            tmp_path / "merged.fastq",
        )


def test_output_matches_input(input_files):
    first, second = input_files

    with pytest.raises(ValueError):
        merge_fastq([first, second], first)


def test_empty_files(tmp_path):
    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    output = tmp_path / "merged.fastq"

    first.write_text("")
    second.write_text("")

    result = merge_fastq([first, second], output)

    assert result["total_reads"] == 0
    assert output.read_text() == ""


def test_malformed_input_leaves_no_output(tmp_path):
    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    output = tmp_path / "merged.fastq"

    first.write_text("@good\nACGT\n+\nIIII\n")
    second.write_text("@bad\nACGT\n+\nIII\n")

    with pytest.raises(ValueError):
        merge_fastq([first, second], output)

    assert not output.exists()


def test_cli_merge(input_files, tmp_path):
    output = tmp_path / "merged.fastq"

    exit_code = main([
        "fastq-merge",
        "--inputs",
        str(input_files[0]),
        str(input_files[1]),
        "--output",
        str(output),
    ])

    assert exit_code == 0
    assert len(list(read_fastq(output))) == 3

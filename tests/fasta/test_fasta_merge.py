
import gzip

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fasta.fasta_merge import merge_fasta
from genomics_toolkit.operations.fasta.fasta_stats import read_fasta


@pytest.fixture
def sample_files(tmp_path):
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"

    first.write_text(
        ">seq1\nACGT\n"
        ">seq2\nGGCC\n"
    )

    second.write_text(
        ">seq3\nATAT\n"
        ">seq4\nTTAA\n"
    )

    return first, second


def test_merge_preserves_order(sample_files, tmp_path):
    output = tmp_path / "merged.fasta"

    result = merge_fasta(sample_files, output)

    assert result["input_files"] == 2
    assert result["total_records"] == 4
    assert result["total_bases"] == 16
    assert result["duplicate_ids"] == 0

    assert list(read_fasta(output)) == [
        ("seq1", "ACGT"),
        ("seq2", "GGCC"),
        ("seq3", "ATAT"),
        ("seq4", "TTAA"),
    ]


def test_duplicate_ids_rejected(sample_files, tmp_path):
    first, second = sample_files
    second.write_text(">seq1\nAAAA\n")

    output = tmp_path / "merged.fasta"

    with pytest.raises(ValueError, match="Duplicate"):
        merge_fasta([first, second], output)

    assert not output.exists()


def test_duplicates_allowed(sample_files, tmp_path):
    first, second = sample_files
    second.write_text(">seq1\nAAAA\n")

    output = tmp_path / "merged.fasta"

    result = merge_fasta(
        [first, second],
        output,
        allow_duplicates=True,
    )

    assert result["duplicate_ids"] == 1
    assert result["total_records"] == 3

    assert list(read_fasta(output))[-1] == (
        "seq1",
        "AAAA",
    )


def test_gzipped_input(sample_files, tmp_path):
    first, second = sample_files
    compressed = tmp_path / "second.fasta.gz"

    with gzip.open(
        compressed, "wt", encoding="utf-8"
    ) as handle:
        handle.write(second.read_text())

    output = tmp_path / "merged.fasta"

    result = merge_fasta(
        [first, compressed],
        output,
    )

    assert result["total_records"] == 4


def test_gzipped_output(sample_files, tmp_path):
    output = tmp_path / "merged.fasta.gz"

    merge_fasta(sample_files, output)

    with gzip.open(
        output, "rt", encoding="utf-8"
    ) as handle:
        content = handle.read()

    assert ">seq1" in content
    assert ">seq4" in content


def test_empty_input_list(tmp_path):
    with pytest.raises(ValueError, match="At least one"):
        merge_fasta([], tmp_path / "merged.fasta")


def test_same_input_output_rejected(sample_files):
    first, second = sample_files

    with pytest.raises(ValueError):
        merge_fasta([first, second], first)


def test_existing_output_preserved_on_failure(
    sample_files, tmp_path
):
    first, second = sample_files
    second.write_text(">seq1\nAAAA\n")

    output = tmp_path / "merged.fasta"
    output.write_text("original content\n")

    with pytest.raises(ValueError):
        merge_fasta([first, second], output)

    assert output.read_text() == "original content\n"


def test_malformed_input_cleans_up(tmp_path):
    valid = tmp_path / "valid.fasta"
    invalid = tmp_path / "invalid.fasta"
    output = tmp_path / "merged.fasta"

    valid.write_text(">seq1\nACGT\n")
    invalid.write_text("ACGT\n>seq2\nGGCC\n")

    with pytest.raises(ValueError):
        merge_fasta([valid, invalid], output)

    assert not output.exists()
    assert list(
        tmp_path.glob(f".{output.name}.*.tmp")
    ) == []


def test_cli_merge(sample_files, tmp_path):
    first, second = sample_files
    output = tmp_path / "merged.fasta"

    exit_code = main([
        "fasta-merge",
        "--inputs",
        str(first),
        str(second),
        "--output",
        str(output),
    ])

    assert exit_code == 0
    assert len(list(read_fasta(output))) == 4

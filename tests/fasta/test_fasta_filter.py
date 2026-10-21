
import gzip

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fasta.fasta_filter import filter_fasta
from genomics_toolkit.operations.fasta.fasta_stats import read_fasta


@pytest.fixture
def sample_fasta(tmp_path):
    fasta = tmp_path / "sample.fasta"
    fasta.write_text(
        ">seq1\nACGT\n"
        ">seq2\nACGTACGT\n"
        ">seq3\nGGCCGGCCGGCC\n"
    )
    return fasta


def test_filter_by_min_length(sample_fasta, tmp_path):
    output = tmp_path / "filtered.fasta"

    result = filter_fasta(
        sample_fasta,
        output,
        min_length=8,
    )

    records = list(read_fasta(output))

    assert [record[0] for record in records] == ["seq2", "seq3"]
    assert result["retained_records"] == 2


def test_filter_by_max_length(sample_fasta, tmp_path):
    output = tmp_path / "filtered.fasta"

    filter_fasta(sample_fasta, output, max_length=4)

    records = list(read_fasta(output))
    assert [record[0] for record in records] == ["seq1"]


def test_filter_by_ids(sample_fasta, tmp_path):
    output = tmp_path / "filtered.fasta"

    filter_fasta(
        sample_fasta,
        output,
        ids={"seq1", "seq3"},
    )

    records = list(read_fasta(output))
    assert [record[0] for record in records] == ["seq1", "seq3"]


def test_combined_filters(sample_fasta, tmp_path):
    output = tmp_path / "filtered.fasta"

    filter_fasta(
        sample_fasta,
        output,
        min_length=8,
        max_length=12,
        ids={"seq1", "seq2"},
    )

    records = list(read_fasta(output))
    assert [record[0] for record in records] == ["seq2"]


def test_same_input_and_output_rejected(sample_fasta):
    with pytest.raises(ValueError, match="same file"):
        filter_fasta(sample_fasta, sample_fasta)


def test_invalid_length_range(sample_fasta, tmp_path):
    with pytest.raises(ValueError, match="Minimum length"):
        filter_fasta(
            sample_fasta,
            tmp_path / "output.fasta",
            min_length=20,
            max_length=10,
        )


def test_gzipped_input(sample_fasta, tmp_path):
    compressed = tmp_path / "sample.fasta.gz"

    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(sample_fasta.read_text())

    output = tmp_path / "filtered.fasta"

    result = filter_fasta(compressed, output, min_length=8)

    assert result["retained_records"] == 2


def test_gzipped_output(sample_fasta, tmp_path):
    output = tmp_path / "filtered.fasta.gz"

    filter_fasta(sample_fasta, output, min_length=8)

    with gzip.open(output, "rt", encoding="utf-8") as handle:
        content = handle.read()

    assert ">seq2" in content
    assert ">seq3" in content
    assert ">seq1" not in content


def test_cli_filter(sample_fasta, tmp_path):
    output = tmp_path / "filtered.fasta"

    exit_code = main([
        "fasta-filter",
        "--input", str(sample_fasta),
        "--output", str(output),
        "--min-length", "8",
        "--ids", "seq1,seq2",
    ])

    assert exit_code == 0
    assert list(read_fasta(output)) == [("seq2", "ACGTACGT")]


def test_malformed_input_preserves_existing_output(tmp_path):
    fasta = tmp_path / "malformed.fasta"
    output = tmp_path / "filtered.fasta"

    fasta.write_text(
        ">seq1\nACGT\n"
        ">seq2\nGGCC\n"
        "ACGTZ\n"
    )

    output.write_text("original content\n")

    # The existing FASTA parser validates structure,
    # not nucleotide alphabet. Use a structural error.
    fasta.write_text("ACGT\n>seq1\nGGCC\n")

    with pytest.raises(ValueError):
        filter_fasta(fasta, output)

    assert output.read_text() == "original content\n"


def test_corrupted_gzip_preserves_existing_output(tmp_path):
    fasta = tmp_path / "corrupted.fasta.gz"
    output = tmp_path / "filtered.fasta"

    fasta.write_bytes(b"not a gzip file")
    output.write_text("original content\n")

    with pytest.raises((OSError, EOFError)):
        filter_fasta(fasta, output)

    assert output.read_text() == "original content\n"


def test_failed_filter_does_not_create_output(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    output = tmp_path / "filtered.fasta"

    fasta.write_text("ACGT\n")

    with pytest.raises(ValueError):
        filter_fasta(fasta, output)

    assert not output.exists()


def test_empty_input_creates_empty_output(tmp_path):
    fasta = tmp_path / "empty.fasta"
    output = tmp_path / "filtered.fasta"

    fasta.write_text("")

    result = filter_fasta(fasta, output)

    assert output.exists()
    assert output.read_text() == ""
    assert result["total_records"] == 0
    assert result["retained_records"] == 0


def test_no_matching_sequences(tmp_path):
    fasta = tmp_path / "sample.fasta"
    output = tmp_path / "filtered.fasta"

    fasta.write_text(">seq1\nACGT\n")

    result = filter_fasta(
        fasta,
        output,
        min_length=100,
    )

    assert output.read_text() == ""
    assert result["retained_records"] == 0
    assert result["removed_records"] == 1


def test_negative_min_length_rejected(sample_fasta, tmp_path):
    output = tmp_path / "filtered.fasta"

    with pytest.raises(ValueError, match="negative"):
        filter_fasta(
            sample_fasta,
            output,
            min_length=-1,
        )

    assert not output.exists()


def test_output_wrapped_at_80_bases(tmp_path):
    fasta = tmp_path / "long.fasta"
    output = tmp_path / "filtered.fasta"

    fasta.write_text(">long\n" + "A" * 165 + "\n")

    filter_fasta(fasta, output)

    lines = output.read_text().splitlines()

    assert lines[0] == ">long"
    assert [len(line) for line in lines[1:]] == [80, 80, 5]


def test_atomic_gzipped_output(tmp_path):
    fasta = tmp_path / "sample.fasta"
    output = tmp_path / "filtered.fasta.gz"

    fasta.write_text(">seq1\nACGT\n")

    result = filter_fasta(fasta, output)

    with gzip.open(output, "rt", encoding="utf-8") as handle:
        assert handle.read() == ">seq1\nACGT\n"

    assert result["retained_records"] == 1

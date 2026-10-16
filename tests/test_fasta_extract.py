
import gzip

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fasta_extract import (
    extract_fasta,
    read_ids_file,
)
from genomics_toolkit.operations.fasta_stats import read_fasta


@pytest.fixture
def sample_fasta(tmp_path):
    fasta = tmp_path / "sample.fasta"
    fasta.write_text(
        ">seq1\nACGT\n"
        ">seq2\nGGCC\n"
        ">seq3\nATAT\n"
    )
    return fasta


def test_read_ids_file(tmp_path):
    ids_file = tmp_path / "ids.txt"
    ids_file.write_text(
        "# Selected sequences\n"
        "seq1\n"
        "\n"
        "seq3\n"
    )

    assert read_ids_file(ids_file) == {"seq1", "seq3"}


def test_duplicate_ids_rejected(tmp_path):
    ids_file = tmp_path / "ids.txt"
    ids_file.write_text("seq1\nseq1\n")

    with pytest.raises(ValueError, match="Duplicate ID"):
        read_ids_file(ids_file)


def test_empty_ids_file_rejected(tmp_path):
    ids_file = tmp_path / "ids.txt"
    ids_file.write_text("# No IDs\n\n")

    with pytest.raises(ValueError, match="no sequence IDs"):
        read_ids_file(ids_file)


def test_extract_selected_sequences(sample_fasta, tmp_path):
    ids_file = tmp_path / "ids.txt"
    output = tmp_path / "extracted.fasta"

    ids_file.write_text("seq3\nseq1\n")

    result = extract_fasta(
        sample_fasta,
        ids_file,
        output,
    )

    assert list(read_fasta(output)) == [
        ("seq1", "ACGT"),
        ("seq3", "ATAT"),
    ]

    assert result["extracted_records"] == 2
    assert result["found_ids"] == 2


def test_missing_ids_reported(sample_fasta, tmp_path):
    ids_file = tmp_path / "ids.txt"
    output = tmp_path / "extracted.fasta"

    ids_file.write_text("seq1\nmissing_seq\n")

    result = extract_fasta(
        sample_fasta,
        ids_file,
        output,
    )

    assert result["missing_ids"] == ["missing_seq"]
    assert result["extracted_records"] == 1


def test_gzipped_input_and_output(sample_fasta, tmp_path):
    compressed = tmp_path / "sample.fasta.gz"
    output = tmp_path / "extracted.fasta.gz"
    ids_file = tmp_path / "ids.txt"

    ids_file.write_text("seq2\n")

    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(sample_fasta.read_text())

    result = extract_fasta(compressed, ids_file, output)

    assert list(read_fasta(output)) == [("seq2", "GGCC")]
    assert result["extracted_records"] == 1


def test_output_cannot_overwrite_input(sample_fasta, tmp_path):
    ids_file = tmp_path / "ids.txt"
    ids_file.write_text("seq1\n")

    original = sample_fasta.read_text()

    with pytest.raises(ValueError, match="same file"):
        extract_fasta(sample_fasta, ids_file, sample_fasta)

    assert sample_fasta.read_text() == original


def test_output_cannot_overwrite_ids_file(sample_fasta, tmp_path):
    ids_file = tmp_path / "ids.txt"
    ids_file.write_text("seq1\n")

    with pytest.raises(ValueError, match="same file"):
        extract_fasta(sample_fasta, ids_file, ids_file)

    assert ids_file.read_text() == "seq1\n"


def test_malformed_fasta_preserves_output(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    ids_file = tmp_path / "ids.txt"
    output = tmp_path / "output.fasta"

    fasta.write_text("ACGT\n>seq1\nGGCC\n")
    ids_file.write_text("seq1\n")
    output.write_text("original content\n")

    with pytest.raises(ValueError):
        extract_fasta(fasta, ids_file, output)

    assert output.read_text() == "original content\n"


def test_duplicate_records_extracted(tmp_path):
    fasta = tmp_path / "duplicate.fasta"
    ids_file = tmp_path / "ids.txt"
    output = tmp_path / "output.fasta"

    fasta.write_text(
        ">seq1\nACGT\n"
        ">seq1\nGGCC\n"
    )
    ids_file.write_text("seq1\n")

    result = extract_fasta(fasta, ids_file, output)

    assert result["extracted_records"] == 2
    assert result["found_ids"] == 1


def test_cli_extract(sample_fasta, tmp_path):
    ids_file = tmp_path / "ids.txt"
    output = tmp_path / "output.fasta"

    ids_file.write_text("seq2\n")

    exit_code = main([
        "fasta-extract",
        "--input", str(sample_fasta),
        "--ids-file", str(ids_file),
        "--output", str(output),
    ])

    assert exit_code == 0
    assert list(read_fasta(output)) == [("seq2", "GGCC")]


def test_corrupted_gzip_preserves_output(tmp_path):
    fasta = tmp_path / "corrupted.fasta.gz"
    ids_file = tmp_path / "ids.txt"
    output = tmp_path / "extracted.fasta"

    fasta.write_bytes(b"not a gzip file")
    ids_file.write_text("seq1\n")
    output.write_text("original content\n")

    with pytest.raises((OSError, EOFError)):
        extract_fasta(fasta, ids_file, output)

    assert output.read_text() == "original content\n"


def test_no_matching_ids(tmp_path):
    fasta = tmp_path / "sample.fasta"
    ids_file = tmp_path / "ids.txt"
    output = tmp_path / "extracted.fasta"

    fasta.write_text(">seq1\nACGT\n")
    ids_file.write_text("missing_seq\n")

    result = extract_fasta(fasta, ids_file, output)

    assert output.exists()
    assert output.read_text() == ""
    assert result["extracted_records"] == 0
    assert result["found_ids"] == 0
    assert result["missing_ids"] == ["missing_seq"]


def test_temporary_files_cleaned_after_failure(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    ids_file = tmp_path / "ids.txt"
    output = tmp_path / "extracted.fasta"

    fasta.write_text("ACGT\n>seq1\nGGCC\n")
    ids_file.write_text("seq1\n")

    with pytest.raises(ValueError):
        extract_fasta(fasta, ids_file, output)

    assert not output.exists()
    assert list(tmp_path.glob(f".{output.name}.*.tmp")) == []


def test_output_wrapped_at_80_bases(tmp_path):
    fasta = tmp_path / "sample.fasta"
    ids_file = tmp_path / "ids.txt"
    output = tmp_path / "extracted.fasta"

    fasta.write_text(">seq1\n" + "A" * 165 + "\n")
    ids_file.write_text("seq1\n")

    extract_fasta(fasta, ids_file, output)

    lines = output.read_text().splitlines()

    assert lines[0] == ">seq1"
    assert [len(line) for line in lines[1:]] == [80, 80, 5]

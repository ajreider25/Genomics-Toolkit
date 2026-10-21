
import gzip

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fasta.fasta_revcomp import (
    reverse_complement,
    write_reverse_complement,
)
from genomics_toolkit.operations.fasta.fasta_stats import read_fasta


def test_reverse_complement_basic():
    assert reverse_complement("ACGTTAGC") == "GCTAACGT"


def test_reverse_complement_lowercase():
    assert reverse_complement("acgt") == "acgt"


def test_reverse_complement_ambiguity_codes():
    assert reverse_complement("RYSWKMBDHVN") == "NBDHVKMWSRY"


def test_reverse_complement_empty_sequence():
    assert reverse_complement("") == ""


def test_reverse_complement_invalid_base():
    with pytest.raises(ValueError, match="Invalid DNA"):
        reverse_complement("ACGTZ")


def test_write_reverse_complement(tmp_path):
    fasta = tmp_path / "input.fasta"
    output = tmp_path / "output.fasta"

    fasta.write_text(
        ">seq1\nACGTTAGC\n"
        ">seq2\nGGCC\n"
    )

    result = write_reverse_complement(fasta, output)

    assert list(read_fasta(output)) == [
        ("seq1", "GCTAACGT"),
        ("seq2", "GGCC"),
    ]

    assert result["sequence_count"] == 2
    assert result["total_bases"] == 12


def test_gzipped_input_and_output(tmp_path):
    fasta = tmp_path / "input.fasta.gz"
    output = tmp_path / "output.fasta.gz"

    with gzip.open(fasta, "wt", encoding="utf-8") as handle:
        handle.write(">seq1\nACGTTAGC\n")

    write_reverse_complement(fasta, output)

    assert list(read_fasta(output)) == [
        ("seq1", "GCTAACGT")
    ]


def test_same_input_output_rejected(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">seq1\nACGT\n")

    with pytest.raises(ValueError, match="same file"):
        write_reverse_complement(fasta, fasta)


def test_invalid_sequence_preserves_output(tmp_path):
    fasta = tmp_path / "input.fasta"
    output = tmp_path / "output.fasta"

    fasta.write_text(
        ">seq1\nACGT\n"
        ">seq2\nACGTZ\n"
    )

    output.write_text("existing content\n")

    with pytest.raises(ValueError, match="Invalid DNA"):
        write_reverse_complement(fasta, output)

    assert output.read_text() == "existing content\n"


def test_cli_reverse_complement(tmp_path):
    fasta = tmp_path / "input.fasta"
    output = tmp_path / "output.fasta"

    fasta.write_text(">seq1\nACGTTAGC\n")

    exit_code = main([
        "fasta-revcomp",
        "--input", str(fasta),
        "--output", str(output),
    ])

    assert exit_code == 0
    assert list(read_fasta(output)) == [
        ("seq1", "GCTAACGT")
    ]


def test_malformed_fasta_preserves_existing_output(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    output = tmp_path / "output.fasta"

    fasta.write_text("ACGT\n>seq1\nGGCC\n")
    output.write_text("original content\n")

    with pytest.raises(ValueError):
        write_reverse_complement(fasta, output)

    assert output.read_text() == "original content\n"


def test_corrupted_gzip_preserves_existing_output(tmp_path):
    fasta = tmp_path / "corrupted.fasta.gz"
    output = tmp_path / "output.fasta"

    fasta.write_bytes(b"not a gzip file")
    output.write_text("original content\n")

    with pytest.raises((OSError, EOFError)):
        write_reverse_complement(fasta, output)

    assert output.read_text() == "original content\n"


def test_invalid_sequence_does_not_create_output(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    output = tmp_path / "output.fasta"

    fasta.write_text(">seq1\nACGTZ\n")

    with pytest.raises(ValueError, match="Invalid DNA"):
        write_reverse_complement(fasta, output)

    assert not output.exists()


def test_empty_fasta(tmp_path):
    fasta = tmp_path / "empty.fasta"
    output = tmp_path / "output.fasta"

    fasta.write_text("")

    result = write_reverse_complement(fasta, output)

    assert result["sequence_count"] == 0
    assert result["total_bases"] == 0
    assert output.read_text() == ""


def test_empty_sequence(tmp_path):
    fasta = tmp_path / "empty_sequence.fasta"
    output = tmp_path / "output.fasta"

    fasta.write_text(">seq1\n")

    result = write_reverse_complement(fasta, output)

    assert result["sequence_count"] == 1
    assert result["total_bases"] == 0
    assert list(read_fasta(output)) == [("seq1", "")]


def test_output_wrapped_at_80_bases(tmp_path):
    fasta = tmp_path / "long.fasta"
    output = tmp_path / "output.fasta"

    fasta.write_text(">seq1\n" + "A" * 165 + "\n")

    write_reverse_complement(fasta, output)

    lines = output.read_text().splitlines()

    assert lines[0] == ">seq1"
    assert [len(line) for line in lines[1:]] == [80, 80, 5]
    assert "".join(lines[1:]) == "T" * 165


def test_temporary_files_cleaned_after_failure(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    output = tmp_path / "output.fasta"

    fasta.write_text(">seq1\nACGTZ\n")

    with pytest.raises(ValueError):
        write_reverse_complement(fasta, output)

    assert list(tmp_path.glob(f".{output.name}.*.tmp")) == []


def test_cli_invalid_sequence_returns_error(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    output = tmp_path / "output.fasta"

    fasta.write_text(">seq1\nACGTZ\n")

    exit_code = main([
        "fasta-revcomp",
        "--input", str(fasta),
        "--output", str(output),
    ])

    assert exit_code == 1
    assert not output.exists()

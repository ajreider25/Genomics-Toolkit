
import csv
import gzip

from genomics_toolkit.operations.fasta_stats import (
    calculate_assembly_stats,
    read_fasta,
    write_assembly_stats,
)
from genomics_toolkit.cli import main


def test_assembly_stats():
    records = [
        ("seq1", "A" * 10),
        ("seq2", "G" * 6),
        ("seq3", "C" * 4),
    ]

    stats = calculate_assembly_stats(records)

    assert stats == {
        "sequence_count": 3,
        "total_bases": 20,
        "min_length": 4,
        "max_length": 10,
        "mean_length": 6.67,
        "n50": 10,
        "l50": 1,
        "gc_percent": 50.0,
    }


def test_empty_assembly():
    stats = calculate_assembly_stats([])

    assert stats["sequence_count"] == 0
    assert stats["total_bases"] == 0
    assert stats["n50"] == 0
    assert stats["l50"] == 0


def test_assembly_with_ambiguous_bases():
    records = [
        ("seq1", "GCNN"),
        ("seq2", "ATNN"),
    ]

    stats = calculate_assembly_stats(records)

    assert stats["total_bases"] == 8
    assert stats["gc_percent"] == 50.0
    assert stats["n50"] == 4
    assert stats["l50"] == 1


def test_assembly_with_empty_sequences():
    records = [
        ("empty", ""),
        ("seq1", "ACGT"),
    ]

    stats = calculate_assembly_stats(records)

    assert stats["sequence_count"] == 2
    assert stats["min_length"] == 0
    assert stats["max_length"] == 4
    assert stats["mean_length"] == 2.0
    assert stats["n50"] == 4
    assert stats["l50"] == 1


def test_write_assembly_stats(tmp_path):
    input_file = tmp_path / "sample.fasta"
    output_file = tmp_path / "results" / "summary.csv"

    input_file.write_text(
        ">seq1\nAAAA\n"
        ">seq2\nGGGG\n"
    )

    write_assembly_stats(input_file, output_file)

    with output_file.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 1
    assert rows[0]["sequence_count"] == "2"
    assert rows[0]["total_bases"] == "8"
    assert rows[0]["n50"] == "4"
    assert rows[0]["l50"] == "1"
    assert rows[0]["gc_percent"] == "50.0"


def test_assembly_stats_gzipped(tmp_path):
    input_file = tmp_path / "sample.fasta.gz"

    with gzip.open(input_file, "wt", encoding="utf-8") as handle:
        handle.write(">seq1\nACGT\n>seq2\nGGCC\n")

    stats = calculate_assembly_stats(read_fasta(input_file))

    assert stats["sequence_count"] == 2
    assert stats["total_bases"] == 8
    assert stats["gc_percent"] == 75.0


def test_cli_with_summary(tmp_path):
    input_file = tmp_path / "sample.fasta"
    output_file = tmp_path / "stats.csv"
    summary_file = tmp_path / "summary.csv"

    input_file.write_text(">seq1\nACGT\n")

    exit_code = main([
        "fasta-stats",
        "--input", str(input_file),
        "--output", str(output_file),
        "--summary", str(summary_file),
    ])

    assert exit_code == 0
    assert output_file.exists()
    assert summary_file.exists()


def test_missing_input_does_not_create_summary(tmp_path):
    input_file = tmp_path / "missing.fasta"
    output_file = tmp_path / "summary.csv"

    try:
        write_assembly_stats(input_file, output_file)
    except FileNotFoundError:
        pass
    else:
        raise AssertionError("Expected FileNotFoundError")

    assert not output_file.exists()

def test_assembly_gc_excludes_ambiguous_bases():
    records = [
        ("seq1", "ACGTACGT"),
        ("seq2", "GGCCNNAT"),
        ("seq3", "ATATATAT"),
    ]

    stats = calculate_assembly_stats(records)

    assert stats["total_bases"] == 24
    assert stats["gc_percent"] == 36.36


import pytest

from genomics_toolkit.operations.fasta_stats import (
    calculate_assembly_stats,
    write_assembly_stats,
)


def test_assembly_with_only_empty_sequences():
    records = [
        ("seq1", ""),
        ("seq2", ""),
    ]

    stats = calculate_assembly_stats(records)

    assert stats["sequence_count"] == 2
    assert stats["total_bases"] == 0
    assert stats["min_length"] == 0
    assert stats["max_length"] == 0
    assert stats["mean_length"] == 0.0
    assert stats["n50"] == 0
    assert stats["l50"] == 0
    assert stats["gc_percent"] == 0.0


def test_invalid_fasta_does_not_create_summary(tmp_path):
    input_file = tmp_path / "invalid.fasta"
    output_file = tmp_path / "summary.csv"

    input_file.write_text("ACGT\n")

    with pytest.raises(ValueError, match="Invalid FASTA"):
        write_assembly_stats(input_file, output_file)

    assert not output_file.exists()


def test_assembly_summary_rejects_input_as_output(tmp_path):
    fasta_file = tmp_path / "sample.fasta"
    original_content = ">seq1\nACGT\n"

    fasta_file.write_text(original_content)

    with pytest.raises(ValueError, match="same file"):
        write_assembly_stats(fasta_file, fasta_file)

    assert fasta_file.read_text() == original_content

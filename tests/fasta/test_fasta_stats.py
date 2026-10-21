import pytest
import csv
from unittest.mock import patch

from genomics_toolkit.workflow_runner import run_workflow

from genomics_toolkit.operations.fasta.fasta_stats import (
    calculate_stats,
    fasta_stats,
    write_fasta_stats,
    read_fasta
)

def test_calculate_stats():
    result = calculate_stats("ACGTNN")

    assert result["length"] == 6
    assert result["gc_percent"] == 50.0
    assert result["at_percent"] == 50.0
    assert result["n_count"] == 2

def test_empty_sequence():
    result = calculate_stats("")

    assert result["length"] == 0
    assert result["gc_percent"] == 0.0
    assert result["at_percent"] == 0.0
    assert result["n_count"] == 0

def test_fasta_stats(tmp_path):
    fasta_file = tmp_path / "sample.fasta"

    fasta_file.write_text(
        ">seq1\nACGT\n"
        ">seq2 description\nGGCCNN\n"
    )

    results = list(fasta_stats(fasta_file))

    assert results == [
        {
            "sequence_id": "seq1",
            "length": 4,
            "gc_percent": 50.0,
            "at_percent": 50.0,
            "n_count": 0,
            "a_count": 1,
            "c_count": 1,
            "g_count": 1,
            "t_count": 1
        },
        {
            "sequence_id": "seq2",
            "length": 6,
            "gc_percent": 100.0,
            "at_percent": 0.0,
            "n_count": 2,
            "a_count": 0,
            "c_count": 2,
            "g_count": 2,
            "t_count": 0
        },
    ]


def test_invalid_fasta(tmp_path):
    fasta_file = tmp_path / "invalid.fasta"
    fasta_file.write_text("ACGT\n")

    with pytest.raises(ValueError):
        list(fasta_stats(fasta_file))


def test_write_fasta_stats(tmp_path):
    input_file = tmp_path / "input.fasta"
    output_file = tmp_path / "results" / "stats.csv"

    input_file.write_text(
        ">seq1\nACGT\n"
        ">seq2\nGGCCNN\n"
    )

    write_fasta_stats(input_file, output_file)

    assert output_file.exists()

    with output_file.open(newline="") as file:
        rows = list(csv.DictReader(file))

    assert len(rows) == 2

    assert rows[0] == {
        "sequence_id": "seq1",
        "length": "4",
        "gc_percent": "50.0",
        "at_percent": "50.0",
        "n_count": "0",
        "a_count": "1",
        "c_count": "1",
        "g_count": "1",
        "t_count": "1"
    }

    assert rows[1] == {
        "sequence_id": "seq2",
        "length": "6",
        "gc_percent": "100.0",
        "at_percent": "0.0",
        "n_count": "2",
        "a_count": "0",
        "c_count": "2",
        "g_count": "2",
        "t_count": "0"
    }

def test_run_workflow_with_custom_paths(tmp_path):
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text("rule all:\n    input: []\n")

    input_file = tmp_path / "sample.fasta"
    input_file.write_text(">seq1\nACGT\n")

    output_file = tmp_path / "stats.csv"

    with patch(
        "genomics_toolkit.workflow_runner.subprocess.run"
    ) as mock_run:
        mock_run.return_value.returncode = 0

        exit_code = run_workflow(
            cores=4,
            snakefile=snakefile,
            input_path=input_file,
            output_path=output_file,
        )

    assert exit_code == 0

    command = mock_run.call_args.args[0]

    assert "--cores" in command
    assert "4" in command
    assert "--config" in command
    assert f"input={input_file.resolve()}" in command
    assert f"output={output_file.resolve()}" in command

def test_lowercase_sequence():
    """Lowercase DNA should produce correct statistics."""
    result = calculate_stats("acgtacgtnn")

    assert result["length"] == 10
    assert result["gc_percent"] == 50.0
    assert result["at_percent"] == 50.0
    assert result["n_count"] == 2


def test_multiline_fasta(tmp_path):
    """Sequence lines should be concatenated correctly."""
    fasta_file = tmp_path / "multiline.fasta"

    fasta_file.write_text(
        ">seq1\n"
        "ACGT\n"
        "GGCC\n"
        ">seq2\n"
        "ATAT\n"
        "NNNN\n"
    )

    results = list(fasta_stats(fasta_file))

    assert len(results) == 2

    assert results[0]["sequence_id"] == "seq1"
    assert results[0]["length"] == 8
    assert results[0]["gc_percent"] == 75.0
    assert results[0]["n_count"] == 0

    assert results[1]["sequence_id"] == "seq2"
    assert results[1]["length"] == 8
    assert results[1]["gc_percent"] == 0.0
    assert results[1]["n_count"] == 4


def test_only_ambiguous_bases():
    """An all-N sequence should not cause division by zero."""
    result = calculate_stats("NNNNNN")

    assert result["length"] == 6
    assert result["gc_percent"] == 0.0
    assert result["at_percent"] == 0.0
    assert result["n_count"] == 6


def test_sequence_without_header(tmp_path):
    """Sequence data before a FASTA header should be rejected."""
    fasta_file = tmp_path / "invalid.fasta"

    fasta_file.write_text("ACGTACGT\n")

    with pytest.raises(ValueError, match="Invalid FASTA"):
        list(read_fasta(fasta_file))


def test_empty_fasta_sequence(tmp_path):
    """A FASTA record with no bases should be handled safely."""
    fasta_file = tmp_path / "empty_sequence.fasta"

    fasta_file.write_text(
        ">empty\n"
        ">valid\n"
        "ACGT\n"
    )

    results = list(fasta_stats(fasta_file))

    assert len(results) == 2

    assert results[0]["sequence_id"] == "empty"
    assert results[0]["length"] == 0
    assert results[0]["gc_percent"] == 0.0
    assert results[0]["n_count"] == 0

    assert results[1]["sequence_id"] == "valid"
    assert results[1]["length"] == 4
    assert results[1]["gc_percent"] == 50.0
    assert results[1]["n_count"] == 0

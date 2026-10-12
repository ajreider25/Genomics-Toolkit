
import csv

from genomics_toolkit.cli import main


def test_fasta_stats_cli(tmp_path):
    input_file = tmp_path / "sample.fasta"
    output_file = tmp_path / "stats.csv"

    input_file.write_text(">seq1\nACGTNN\n")

    exit_code = main([
        "fasta-stats",
        "--input", str(input_file),
        "--output", str(output_file),
    ])

    assert exit_code == 0
    assert output_file.exists()

    with output_file.open(newline="") as file:
        rows = list(csv.DictReader(file))

    assert rows == [{
        "sequence_id": "seq1",
        "length": "6",
        "gc_percent": "50.0",
        "n_count": "2",
    }]


def test_fasta_stats_cli_missing_input(tmp_path):
    output_file = tmp_path / "stats.csv"

    exit_code = main([
        "fasta-stats",
        "--input", str(tmp_path / "missing.fasta"),
        "--output", str(output_file),
    ])

    assert exit_code == 1
    assert not output_file.exists()

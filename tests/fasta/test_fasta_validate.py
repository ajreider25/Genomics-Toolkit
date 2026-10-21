
import gzip
import json
import pytest

from genomics_toolkit.operations.fasta.fasta_validate import (
    validate_fasta,
    write_validation_report
)

from genomics_toolkit.cli import main


def test_valid_fasta(tmp_path):
    fasta = tmp_path / "valid.fasta"
    fasta.write_text(">seq1\nACGTRYSWKMBDHVN\n")

    report = validate_fasta(fasta)

    assert report["valid"] is True
    assert report["sequence_count"] == 1
    assert report["error_count"] == 0


def test_duplicate_ids(tmp_path):
    fasta = tmp_path / "duplicate.fasta"
    fasta.write_text(">seq1\nACGT\n>seq1\nGGCC\n")

    report = validate_fasta(fasta)

    assert report["valid"] is False
    assert report["error_count"] == 1
    assert "Duplicate sequence ID" in report["errors"][0]


def test_invalid_nucleotide(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text(">seq1\nACGTZ\n")

    report = validate_fasta(fasta)

    assert report["valid"] is False
    assert "Z" in report["errors"][0]


def test_empty_sequence_permissive(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text(">seq1\n")

    report = validate_fasta(fasta)

    assert report["valid"] is True
    assert report["warning_count"] == 1


def test_empty_sequence_strict(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text(">seq1\n")

    report = validate_fasta(fasta, strict=True)

    assert report["valid"] is False
    assert report["error_count"] == 1


def test_empty_file(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("")

    report = validate_fasta(fasta)

    assert report["valid"] is True
    assert report["warning_count"] == 1


def test_missing_header(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text("ACGT\n")

    report = validate_fasta(fasta)

    assert report["valid"] is False
    assert "FASTA parsing error" in report["errors"][0]


def test_gzipped_fasta(tmp_path):
    fasta = tmp_path / "sample.fasta.gz"

    with gzip.open(fasta, "wt", encoding="utf-8") as handle:
        handle.write(">seq1\nACGT\n")

    report = validate_fasta(fasta)

    assert report["valid"] is True
    assert report["sequence_count"] == 1


def test_cli_valid_fasta(tmp_path):
    fasta = tmp_path / "valid.fasta"
    fasta.write_text(">seq1\nACGT\n")

    exit_code = main([
        "fasta-validate",
        "--input", str(fasta),
    ])

    assert exit_code == 0


def test_cli_invalid_fasta(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text(">seq1\nACGTZ\n")

    exit_code = main([
        "fasta-validate",
        "--input", str(fasta),
        "--strict",
    ])

    assert exit_code == 1

def test_malformed_fasta_returns_report(tmp_path):
    fasta = tmp_path / "malformed.fasta"
    fasta.write_text("ACGT\n")

    report = validate_fasta(fasta)

    assert report["valid"] is False
    assert report["error_count"] == 1
    assert "FASTA parsing error" in report["errors"][0]


def test_corrupted_gzip_returns_report(tmp_path):
    fasta = tmp_path / "corrupted.fasta.gz"
    fasta.write_bytes(b"not a gzip file")

    report = validate_fasta(fasta)

    assert report["valid"] is False
    assert report["error_count"] >= 1


def test_write_validation_report(tmp_path):
    fasta = tmp_path / "sample.fasta"
    output = tmp_path / "reports" / "validation.json"

    fasta.write_text(">seq1\nACGT\n")

    report = validate_fasta(fasta)
    write_validation_report(report, output)

    saved = json.loads(output.read_text())

    assert saved == report
    assert saved["valid"] is True


def test_cli_json_report(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    output = tmp_path / "validation.json"

    fasta.write_text(">seq1\nACGTZ\n")

    exit_code = main([
        "fasta-validate",
        "--input", str(fasta),
        "--report", str(output),
    ])

    assert exit_code == 1
    assert output.exists()

    saved = json.loads(output.read_text())
    assert saved["valid"] is False
    assert saved["error_count"] == 1

def test_cli_rejects_report_overwriting_input(tmp_path):
    fasta = tmp_path / "sample.fasta"
    original = ">seq1\nACGT\n"
    fasta.write_text(original)

    exit_code = main([
        "fasta-validate",
        "--input", str(fasta),
        "--report", str(fasta),
    ])

    assert exit_code == 1
    assert fasta.read_text() == original
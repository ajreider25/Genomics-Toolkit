
import gzip
import json

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fastq.fastq_validate import (
    validate_fastq,
    write_validation_report,
)


def test_valid_fastq(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "@read1\nACGT\n+\nIIII\n"
        "@read2\nGGNN\n+\n####\n"
    )

    report = validate_fastq(fastq)

    assert report["valid"] is True
    assert report["record_count"] == 2
    assert report["total_bases"] == 8
    assert report["error_count"] == 0


def test_invalid_header(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "read1\nACGT\n+\nIIII\n"
    )

    report = validate_fastq(fastq)

    assert report["valid"] is False
    assert report["error_count"] == 1


def test_quality_length_mismatch(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "@read1\nACGT\n+\nIII\n"
    )

    report = validate_fastq(fastq)

    assert report["valid"] is False
    assert any(
        "lengths do not match" in error
        for error in report["errors"]
    )


def test_invalid_nucleotide(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "@read1\nACGX\n+\nIIII\n"
    )

    report = validate_fastq(fastq)

    assert report["valid"] is False
    assert any(
        "invalid bases" in error
        for error in report["errors"]
    )


def test_invalid_separator(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "@read1\nACGT\n-\nIIII\n"
    )

    report = validate_fastq(fastq)

    assert report["valid"] is False


def test_incomplete_record(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "@read1\nACGT\n+\n"
    )

    report = validate_fastq(fastq)

    assert report["valid"] is False
    assert any(
        "incomplete" in error
        for error in report["errors"]
    )


def test_empty_file_warning(tmp_path):
    fastq = tmp_path / "empty.fastq"
    fastq.write_text("")

    report = validate_fastq(fastq)

    assert report["valid"] is True
    assert report["warning_count"] == 1


def test_empty_file_strict(tmp_path):
    fastq = tmp_path / "empty.fastq"
    fastq.write_text("")

    report = validate_fastq(fastq, strict=True)

    assert report["valid"] is False
    assert report["error_count"] == 1


def test_gzipped_fastq(tmp_path):
    fastq = tmp_path / "reads.fastq.gz"

    with gzip.open(
        fastq, "wt", encoding="utf-8"
    ) as handle:
        handle.write("@read1\nACGT\n+\nIIII\n")

    report = validate_fastq(fastq)

    assert report["valid"] is True
    assert report["record_count"] == 1


def test_corrupted_gzip(tmp_path):
    fastq = tmp_path / "reads.fastq.gz"
    fastq.write_bytes(b"not a gzip file")

    report = validate_fastq(fastq)

    assert report["valid"] is False
    assert report["error_count"] >= 1


def test_json_report(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n")

    report = validate_fastq(fastq)
    output = tmp_path / "report.json"

    write_validation_report(
        report,
        output,
        input_path=fastq,
    )

    saved = json.loads(output.read_text())

    assert saved["valid"] is True
    assert saved["record_count"] == 1


def test_report_cannot_overwrite_input(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n")

    report = validate_fastq(fastq)

    with pytest.raises(ValueError):
        write_validation_report(
            report,
            fastq,
            input_path=fastq,
        )


def test_cli_valid_fastq(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n")

    exit_code = main([
        "fastq-validate",
        "--input", str(fastq),
    ])

    assert exit_code == 0


def test_cli_invalid_fastq(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@read1\nACGT\n+\nIII\n")

    exit_code = main([
        "fastq-validate",
        "--input", str(fastq),
    ])

    assert exit_code == 1

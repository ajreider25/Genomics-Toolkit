
import gzip

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.fastq.fastq_pair_validate import (
    normalize_read_id,
    read_fastq_pairs,
    validate_fastq_pairs,
)


@pytest.fixture
def paired_files(tmp_path):
    r1 = tmp_path / "reads_R1.fastq"
    r2 = tmp_path / "reads_R2.fastq"

    r1.write_text(
        "@read1/1\nACGT\n+\nIIII\n"
        "@read2/1\nGGCC\n+\n5555\n",
        encoding="utf-8",
    )

    r2.write_text(
        "@read1/2\nTGCA\n+\nIIII\n"
        "@read2/2\nCCGG\n+\n5555\n",
        encoding="utf-8",
    )

    return r1, r2


def test_normalize_slash_ids():
    assert normalize_read_id("read1/1") == ("read1", 1)
    assert normalize_read_id("read1/2") == ("read1", 2)


def test_normalize_illumina_ids():
    assert normalize_read_id("read1 1:N:0:ATCG") == ("read1", 1)
    assert normalize_read_id("read1 2:N:0:ATCG") == ("read1", 2)


def test_normalize_unannotated_id():
    assert normalize_read_id("read1") == ("read1", None)


def test_valid_pairs(paired_files):
    r1, r2 = paired_files

    result = validate_fastq_pairs(r1, r2)

    assert result["valid"] is True
    assert result["pair_count"] == 2
    assert result["r1_bases"] == 8
    assert result["r2_bases"] == 8
    assert result["total_bases"] == 16


def test_pair_reader(paired_files):
    r1, r2 = paired_files

    pairs = list(read_fastq_pairs(r1, r2))

    assert len(pairs) == 2
    assert pairs[0][0][0] == "read1/1"
    assert pairs[0][1][0] == "read1/2"


def test_mismatched_ids(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"

    r1.write_text("@read1/1\nACGT\n+\nIIII\n")
    r2.write_text("@read2/2\nTGCA\n+\nIIII\n")

    result = validate_fastq_pairs(r1, r2)

    assert result["valid"] is False
    assert result["error_count"] == 1
    assert "mismatched read IDs" in result["errors"][0]


def test_different_read_counts(paired_files):
    r1, r2 = paired_files

    r2.write_text("@read1/2\nTGCA\n+\nIIII\n")

    result = validate_fastq_pairs(r1, r2)

    assert result["valid"] is False
    assert result["pair_count"] == 1
    assert "different read counts" in result["errors"][0]


def test_wrong_mate_label(paired_files):
    r1, r2 = paired_files

    r2.write_text("@read1/1\nTGCA\n+\nIIII\n")

    result = validate_fastq_pairs(r1, r2)

    assert result["valid"] is False
    assert "incorrect mate label" in result["errors"][0]


def test_gzip_pairs(paired_files, tmp_path):
    r1, r2 = paired_files
    gz_r1 = tmp_path / "reads_R1.fastq.gz"
    gz_r2 = tmp_path / "reads_R2.fastq.gz"

    for source, destination in [(r1, gz_r1), (r2, gz_r2)]:
        with gzip.open(destination, "wt", encoding="utf-8") as handle:
            handle.write(source.read_text())

    result = validate_fastq_pairs(gz_r1, gz_r2)

    assert result["valid"] is True
    assert result["pair_count"] == 2


def test_same_file_rejected(paired_files):
    r1, _ = paired_files

    result = validate_fastq_pairs(r1, r1)

    assert result["valid"] is False


def test_missing_file(paired_files, tmp_path):
    r1, _ = paired_files

    result = validate_fastq_pairs(r1, tmp_path / "missing.fastq")

    assert result["valid"] is False


def test_malformed_fastq(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"

    r1.write_text("@read1/1\nACGT\n+\nIII\n")
    r2.write_text("@read1/2\nTGCA\n+\nIIII\n")

    result = validate_fastq_pairs(r1, r2)

    assert result["valid"] is False


def test_empty_paired_files(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"

    r1.write_text("")
    r2.write_text("")

    result = validate_fastq_pairs(r1, r2)

    assert result["valid"] is True
    assert result["pair_count"] == 0


def test_cli_valid_pairs(paired_files):
    r1, r2 = paired_files

    exit_code = main([
        "fastq-pair-validate",
        "--r1", str(r1),
        "--r2", str(r2),
    ])

    assert exit_code == 0


def test_cli_invalid_pairs(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"

    r1.write_text("@read1/1\nACGT\n+\nIIII\n")
    r2.write_text("@other/2\nTGCA\n+\nIIII\n")

    exit_code = main([
        "fastq-pair-validate",
        "--r1", str(r1),
        "--r2", str(r2),
    ])

    assert exit_code == 1

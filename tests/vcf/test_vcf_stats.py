
import json

import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.vcf.vcf_stats import (
    calculate_vcf_stats,
    classify_variant,
    write_vcf_stats,
)


@pytest.fixture
def sample_vcf(tmp_path):
    path = tmp_path / "sample.vcf"

    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tG\t60\tPASS\t.\n"
        "chr1\t200\t.\tC\tT\t40\tPASS\t.\n"
        "chr1\t300\t.\tA\tC\t50\tPASS\t.\n"
        "chr2\t400\t.\tA\tAT\t30\tPASS\t.\n"
        "chr2\t500\t.\tAT\tA\t20\tq10\t.\n",
        encoding="utf-8",
    )

    return path


@pytest.mark.parametrize(
    "ref,alt,expected",
    [
        ("A", "G", "snv"),
        ("AC", "GT", "mnv"),
        ("A", "AT", "insertion"),
        ("AT", "A", "deletion"),
        ("AC", "GTG", "complex"),
        ("A", "<DEL>", "symbolic"),
        ("A", "*", "symbolic"),
        ("A", ".", "non_variant"),
        ("A", "A", "non_variant"),
    ],
)
def test_classify_variant(ref, alt, expected):
    assert classify_variant(ref, alt) == expected


def test_basic_stats(sample_vcf):
    stats = calculate_vcf_stats(sample_vcf)

    assert stats["total_records"] == 5
    assert stats["total_alt_alleles"] == 5
    assert stats["snvs"] == 3
    assert stats["insertions"] == 1
    assert stats["deletions"] == 1
    assert stats["transitions"] == 2
    assert stats["transversions"] == 1
    assert stats["titv_ratio"] == 2.0
    assert stats["pass_records"] == 4
    assert stats["mean_quality"] == 40.0


def test_multiallelic_record(tmp_path):
    path = tmp_path / "multi.vcf"

    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tG,C\t60\tPASS\t.\n"
    )

    stats = calculate_vcf_stats(path)

    assert stats["total_records"] == 1
    assert stats["multiallelic_records"] == 1
    assert stats["total_alt_alleles"] == 2
    assert stats["snvs"] == 2


def test_empty_vcf(tmp_path):
    path = tmp_path / "empty.vcf"

    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    stats = calculate_vcf_stats(path)

    assert stats["total_records"] == 0
    assert stats["mean_quality"] is None
    assert stats["titv_ratio"] is None


def test_invalid_vcf(tmp_path):
    path = tmp_path / "invalid.vcf"

    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t0\t.\tA\tG\t60\tPASS\t.\n"
    )

    with pytest.raises(ValueError):
        calculate_vcf_stats(path)


def test_json_output(sample_vcf, tmp_path):
    output = tmp_path / "stats.json"

    stats = calculate_vcf_stats(sample_vcf)
    write_vcf_stats(stats, output)

    saved = json.loads(output.read_text())

    assert saved["total_records"] == 5
    assert saved["titv_ratio"] == 2.0


def test_cli_stats(sample_vcf, tmp_path):
    output = tmp_path / "stats.json"

    exit_code = main([
        "vcf-stats",
        "--input", str(sample_vcf),
        "--output", str(output),
    ])

    assert exit_code == 0
    assert output.exists()

    saved = json.loads(output.read_text())
    assert saved["total_records"] == 5

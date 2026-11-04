from genomics_toolkit.operations.vcf.vcf_validate import (
    validate_vcf,
)
from genomics_toolkit.cli import main


def test_valid_vcf(tmp_path):
    path = tmp_path / "sample.vcf"

    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tG\t60\tPASS\t.\n"
        "chr2\t250\t.\tAT\tA\t50\tPASS\t.\n"
    )

    result = validate_vcf(path)

    assert result["valid"] is True
    assert result["record_count"] == 2
    assert result["error_count"] == 0


def test_invalid_position(tmp_path):
    path = tmp_path / "invalid.vcf"

    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t0\t.\tA\tG\t60\tPASS\t.\n"
    )

    result = validate_vcf(path)

    assert result["valid"] is False
    assert any("POS" in error for error in result["errors"])


def test_missing_header(tmp_path):
    path = tmp_path / "missing_header.vcf"

    path.write_text("##fileformat=VCFv4.2\n")

    result = validate_vcf(path)

    assert result["valid"] is False
    assert any(
        "Missing #CHROM" in error
        for error in result["errors"]
    )


def test_incorrect_column_count(tmp_path):
    path = tmp_path / "bad_columns.vcf"

    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tG\n"
    )

    result = validate_vcf(path)

    assert result["valid"] is False
    assert result["error_count"] == 1

def test_cli_valid_vcf(tmp_path):
    path = tmp_path / "valid.vcf"

    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tG\t60\tPASS\t.\n",
        encoding="utf-8",
    )

    exit_code = main([
        "vcf-validate",
        "--input",
        str(path),
    ])

    assert exit_code == 0


def test_cli_invalid_vcf(tmp_path):
    path = tmp_path / "invalid.vcf"

    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t0\t.\tA\tG\t60\tPASS\t.\n",
        encoding="utf-8",
    )

    exit_code = main([
        "vcf-validate",
        "--input",
        str(path),
    ])

    assert exit_code == 1


def test_cli_missing_file(tmp_path):
    exit_code = main([
        "vcf-validate",
        "--input",
        str(tmp_path / "missing.vcf"),
    ])

    assert exit_code == 1
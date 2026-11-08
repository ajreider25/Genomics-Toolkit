
import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.vcf.vcf_extract import (
    extract_vcf,
    load_variant_ids,
)


@pytest.fixture
def sample_vcf(tmp_path):
    path = tmp_path / "sample.vcf"

    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\trs1\tA\tG\t60\tPASS\t.\n"
        "chr1\t200\trs2\tC\tT\t50\tPASS\t.\n"
        "chr1\t300\trs3;custom3\tG\tA\t40\tPASS\t.\n"
        "chr2\t150\trs4\tA\tAT\t35\tPASS\t.\n"
        "chr2\t400\t.\tAT\tA\t30\tPASS\t.\n",
        encoding="utf-8",
    )

    return path


def positions(path):
    return [
        int(line.split("\t")[1])
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]


def test_extract_chromosome(sample_vcf, tmp_path):
    output = tmp_path / "chr1.vcf"

    result = extract_vcf(
        sample_vcf, output, chrom="chr1"
    )

    assert result["extracted_records"] == 3
    assert positions(output) == [100, 200, 300]


def test_extract_interval(sample_vcf, tmp_path):
    output = tmp_path / "region.vcf"

    extract_vcf(
        sample_vcf,
        output,
        chrom="chr1",
        start=150,
        end=300,
    )

    assert positions(output) == [200, 300]


def test_interval_boundaries(sample_vcf, tmp_path):
    output = tmp_path / "boundary.vcf"

    extract_vcf(
        sample_vcf,
        output,
        chrom="chr1",
        start=200,
        end=200,
    )

    assert positions(output) == [200]


def test_extract_id(sample_vcf, tmp_path):
    output = tmp_path / "rs2.vcf"

    extract_vcf(
        sample_vcf, output, variant_id="rs2"
    )

    assert positions(output) == [200]


def test_semicolon_separated_id(sample_vcf, tmp_path):
    output = tmp_path / "custom.vcf"

    extract_vcf(
        sample_vcf, output, variant_id="custom3"
    )

    assert positions(output) == [300]


def test_extract_ids_file(sample_vcf, tmp_path):
    ids = tmp_path / "ids.txt"
    ids.write_text("rs1\nrs4\n")

    output = tmp_path / "selected.vcf"

    extract_vcf(
        sample_vcf, output, ids_file=ids
    )

    assert positions(output) == [100, 150]


def test_combined_criteria(sample_vcf, tmp_path):
    output = tmp_path / "combined.vcf"

    extract_vcf(
        sample_vcf,
        output,
        chrom="chr1",
        start=150,
        variant_id="rs2",
    )

    assert positions(output) == [200]


def test_no_matches_preserves_header(sample_vcf, tmp_path):
    output = tmp_path / "empty.vcf"

    result = extract_vcf(
        sample_vcf, output, chrom="chr99"
    )

    assert result["extracted_records"] == 0

    lines = output.read_text().splitlines()

    assert len(lines) == 2
    assert lines[0].startswith("##fileformat")
    assert lines[1].startswith("#CHROM")


def test_invalid_interval(sample_vcf, tmp_path):
    with pytest.raises(ValueError):
        extract_vcf(
            sample_vcf,
            tmp_path / "bad.vcf",
            chrom="chr1",
            start=300,
            end=100,
        )


def test_interval_requires_chrom(sample_vcf, tmp_path):
    with pytest.raises(ValueError):
        extract_vcf(
            sample_vcf,
            tmp_path / "bad.vcf",
            start=100,
        )


def test_requires_selection(sample_vcf, tmp_path):
    with pytest.raises(ValueError):
        extract_vcf(
            sample_vcf,
            tmp_path / "bad.vcf",
        )


def test_same_input_output(sample_vcf):
    with pytest.raises(ValueError):
        extract_vcf(
            sample_vcf,
            sample_vcf,
            chrom="chr1",
        )


def test_load_variant_ids(tmp_path):
    path = tmp_path / "ids.txt"
    path.write_text(
        "# candidate variants\n"
        "rs1\n"
        "\n"
        "rs2\n"
    )

    assert load_variant_ids(path) == {"rs1", "rs2"}


def test_cli_extract(sample_vcf, tmp_path):
    output = tmp_path / "cli.vcf"

    exit_code = main([
        "vcf-extract",
        "--input", str(sample_vcf),
        "--output", str(output),
        "--chrom", "chr1",
        "--start", "150",
        "--end", "300",
    ])

    assert exit_code == 0
    assert positions(output) == [200, 300]


def test_cli_invalid_interval(sample_vcf, tmp_path):
    exit_code = main([
        "vcf-extract",
        "--input", str(sample_vcf),
        "--output", str(tmp_path / "bad.vcf"),
        "--chrom", "chr1",
        "--start", "300",
        "--end", "100",
    ])

    assert exit_code == 1

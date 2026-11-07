
import pytest

from genomics_toolkit.cli import main
from genomics_toolkit.operations.vcf.vcf_filter import (
    filter_vcf,
    get_info_depth,
)


@pytest.fixture
def sample_vcf(tmp_path):
    path = tmp_path / "sample.vcf"

    path.write_text(
        '##fileformat=VCFv4.2\n'
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        'chr1\t100\t.\tA\tG\t60\tPASS\tDP=35\n'
        'chr1\t200\t.\tC\tT\t15\tq10\tDP=12\n'
        'chr2\t300\t.\tA\tAT\t45\tPASS\tDP=40\n'
        'chr2\t400\t.\tAT\tA\t.\tPASS\tDP=8\n'
        'chr3\t500\t.\tG\tC\t80\tPASS\t.\n',
        encoding="utf-8",
    )

    return path


def positions(path):
    """Return positions from non-header VCF records."""
    return [
        int(line.split("\t")[1])
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]


def test_info_depth():
    assert get_info_depth("DP=35;AF=0.1") == 35
    assert get_info_depth("AF=0.1;DP=12") == 12
    assert get_info_depth(".") is None


def test_invalid_info_depth():
    with pytest.raises(ValueError):
        get_info_depth("DP=invalid")


def test_no_filters(sample_vcf, tmp_path):
    output = tmp_path / "output.vcf"

    result = filter_vcf(sample_vcf, output)

    assert result["total_records"] == 5
    assert result["kept_records"] == 5
    assert positions(output) == [100, 200, 300, 400, 500]


def test_min_quality(sample_vcf, tmp_path):
    output = tmp_path / "output.vcf"

    filter_vcf(sample_vcf, output, min_qual=40)

    assert positions(output) == [100, 300, 500]


def test_pass_only(sample_vcf, tmp_path):
    output = tmp_path / "output.vcf"

    filter_vcf(sample_vcf, output, pass_only=True)

    assert positions(output) == [100, 300, 400, 500]


def test_chromosome(sample_vcf, tmp_path):
    output = tmp_path / "output.vcf"

    filter_vcf(sample_vcf, output, chrom="chr2")

    assert positions(output) == [300, 400]


def test_min_depth(sample_vcf, tmp_path):
    output = tmp_path / "output.vcf"

    filter_vcf(sample_vcf, output, min_depth=20)

    assert positions(output) == [100, 300]


def test_snv_filter(sample_vcf, tmp_path):
    output = tmp_path / "output.vcf"

    filter_vcf(sample_vcf, output, variant_type="snv")

    assert positions(output) == [100, 200, 500]


def test_insertion_filter(sample_vcf, tmp_path):
    output = tmp_path / "output.vcf"

    filter_vcf(sample_vcf, output, variant_type="insertion")

    assert positions(output) == [300]


def test_combined_filters(sample_vcf, tmp_path):
    output = tmp_path / "output.vcf"

    result = filter_vcf(
        sample_vcf,
        output,
        min_qual=30,
        pass_only=True,
        min_depth=20,
        variant_type="snv",
    )

    assert positions(output) == [100]
    assert result["kept_records"] == 1
    assert result["filtered_records"] == 4


def test_headers_preserved(sample_vcf, tmp_path):
    output = tmp_path / "output.vcf"

    filter_vcf(sample_vcf, output, min_qual=100)

    lines = output.read_text().splitlines()

    assert lines[0] == "##fileformat=VCFv4.2"
    assert lines[2].startswith("#CHROM")
    assert len(lines) == 3


def test_same_input_output(sample_vcf):
    with pytest.raises(ValueError):
        filter_vcf(sample_vcf, sample_vcf)


def test_invalid_min_depth(sample_vcf, tmp_path):
    with pytest.raises(ValueError):
        filter_vcf(
            sample_vcf,
            tmp_path / "output.vcf",
            min_depth=-1,
        )


def test_invalid_vcf(tmp_path):
    source = tmp_path / "invalid.vcf"
    output = tmp_path / "output.vcf"

    source.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t0\t.\tA\tG\t60\tPASS\t.\n"
    )

    with pytest.raises(ValueError):
        filter_vcf(source, output)

    assert not output.exists()


def test_cli_filter(sample_vcf, tmp_path):
    output = tmp_path / "filtered.vcf"

    exit_code = main([
        "vcf-filter",
        "--input", str(sample_vcf),
        "--output", str(output),
        "--min-qual", "30",
        "--pass-only",
        "--variant-type", "snv",
    ])

    assert exit_code == 0
    assert positions(output) == [100, 500]


import pytest

from genomics_toolkit.operations.vcf.vcf_reader import (
    get_sample_genotypes,
    get_variant_type,
    iter_variants,
    open_variant_file,
    summarize_variant,
)


@pytest.fixture
def sample_vcf(tmp_path):
    path = tmp_path / "sample.vcf"

    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100000>\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,'
        'Description="Site depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,'
        'Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,'
        'Description="Sample depth">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,'
        'Description="Genotype quality">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
        "FILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "chr1\t100\trs1\tA\tG\t60\tPASS\t"
        "DP=35\tGT:DP:GQ\t0/1:18:99\n"
        "chr1\t200\trs2\tAT\tA\t50\tPASS\t"
        "DP=40\tGT:DP:GQ\t1/1:22:80\n"
        "chr1\t300\trs3\tC\tT,G\t45\tPASS\t"
        "DP=25\tGT:DP:GQ\t1/2:12:60\n",
        encoding="utf-8",
    )

    return path


def test_open_variant_file(sample_vcf):
    with open_variant_file(sample_vcf) as handle:
        assert list(handle.header.samples) == ["SAMPLE1"]


def test_iter_variants(sample_vcf):
    records = list(iter_variants(sample_vcf))

    assert len(records) == 3
    assert records[0].chrom == "chr1"
    assert records[0].pos == 100
    assert records[0].ref == "A"
    assert records[0].alts == ("G",)


def test_variant_types(sample_vcf):
    with open_variant_file(sample_vcf) as handle:
        records = list(handle)

        assert get_variant_type(records[0]) == ["snv"]
        assert get_variant_type(records[1]) == ["deletion"]
        assert get_variant_type(records[2]) == ["snv", "snv"]


def test_sample_genotypes(sample_vcf):
    with open_variant_file(sample_vcf) as handle:
        record = next(iter(handle))

        samples = get_sample_genotypes(record)

        assert samples["SAMPLE1"]["GT"] == (0, 1)
        assert samples["SAMPLE1"]["DP"] == 18
        assert samples["SAMPLE1"]["GQ"] == 99
        assert samples["SAMPLE1"]["phased"] is False


def test_multiallelic_genotype(sample_vcf):
    with open_variant_file(sample_vcf) as handle:
        records = list(handle)

        samples = get_sample_genotypes(records[2])

        assert samples["SAMPLE1"]["GT"] == (1, 2)


def test_summarize_variant(sample_vcf):
    with open_variant_file(sample_vcf) as handle:
        record = next(iter(handle))

        summary = summarize_variant(record)

        assert summary["chrom"] == "chr1"
        assert summary["pos"] == 100
        assert summary["id"] == "rs1"
        assert summary["ref"] == "A"
        assert summary["alts"] == ["G"]
        assert summary["qual"] == 60
        assert summary["types"] == ["snv"]


def test_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        with open_variant_file(
            tmp_path / "missing.vcf"
        ):
            pass


def test_compressed_vcf(sample_vcf, tmp_path):
    import pysam

    compressed = tmp_path / "sample.vcf.gz"

    pysam.tabix_compress(
        str(sample_vcf),
        str(compressed),
        force=True,
    )

    records = list(iter_variants(compressed))

    assert len(records) == 3
    assert records[0].pos == 100


def test_bcf(sample_vcf, tmp_path):
    import pysam

    bcf_path = tmp_path / "sample.bcf"

    with pysam.VariantFile(str(sample_vcf)) as source:
        with pysam.VariantFile(
            str(bcf_path),
            mode="wb",
            header=source.header,
        ) as destination:
            for record in source:
                destination.write(record)

    records = list(iter_variants(bcf_path))

    assert len(records) == 3
    assert records[0].pos == 100
    assert records[1].ref == "AT"
    assert records[2].alts == ("T", "G")

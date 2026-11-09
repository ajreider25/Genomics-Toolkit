
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import pysam

from .vcf_stats import classify_variant


@contextmanager
def open_variant_file(path: str | Path):
    """Open a VCF, VCF.gz, or BCF file for reading."""
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(path)

    with pysam.VariantFile(str(path), mode="r") as handle:
        yield handle


def iter_variants(
    path: str | Path,
) -> Iterator[pysam.VariantRecord]:
    """Stream variants while keeping the file open."""
    with open_variant_file(path) as handle:
        for record in handle:
            yield record


def get_variant_type(
    record: pysam.VariantRecord,
) -> list[str]:
    """Classify each alternate allele independently."""
    if not record.alts:
        return ["non_variant"]

    return [
        classify_variant(record.ref, alt)
        for alt in record.alts
    ]


def get_sample_genotypes(
    record: pysam.VariantRecord,
) -> dict:
    """Extract sample-level genotype information."""
    result = {}

    for sample_name, sample in record.samples.items():
        genotype = sample.get("GT")

        result[sample_name] = {
            "GT": (
                tuple(genotype)
                if genotype is not None
                else None
            ),
            "phased": sample.phased,
            "DP": sample.get("DP"),
            "GQ": sample.get("GQ"),
        }

    return result


def summarize_variant(
    record: pysam.VariantRecord,
) -> dict:
    """Return a simple, serializable variant summary."""
    quality = record.qual

    return {
        "chrom": record.chrom,
        "pos": record.pos,
        "id": record.id,
        "ref": record.ref,
        "alts": list(record.alts or []),
        "qual": quality,
        "filters": list(record.filter.keys()),
        "types": get_variant_type(record),
        "samples": get_sample_genotypes(record),
    }

import json
import math
from pathlib import Path

from .vcf_validate import validate_vcf

TRANSITIONS = {
    ("A", "G"),
    ("G", "A"),
    ("C", "T"),
    ("T", "C")
}


def classify_variant(ref: str, alt: str) -> str:
    """Classify a single REF/ALT allele pair."""

    ref = ref.upper()
    alt = alt.upper()

    if alt == ".":
        return "non_variant"

    if (
        alt == "*"
        or alt.startswith("<")
        or "[" in alt
        or "]" in alt
    ):
        return "symbolic"

    valid_bases = set("ACGTN")

    if (
        not ref
        or not alt
        or not set(ref) <= valid_bases
        or not set(alt) <= valid_bases
    ):
        return "other"

    if ref == alt:
        return "non_variant"

    if len(ref) == 1 and len(alt) == 1:
        return "snv"

    if len(ref) == len(alt):
        return "mnv"

    if len(ref) < len(alt) and alt.startswith(ref):
        return "insertion"

    if len(ref) > len(alt) and ref.startswith(alt):
        return "deletion"

    return "complex"


def calculate_vcf_stats(input_path: str | Path) -> dict:
    """Calculate summary statistics for a text VCF file."""

    input_path = Path(input_path)

    validation = validate_vcf(input_path)

    if not validation["valid"]:
        raise ValueError(
            "Cannot calculate stats for an invalid VCF: "
            + "; ".join(validation["errors"][:5])
        )

    stats = {
        "total_records": 0,
        "total_alt_alleles": 0,
        "multiallelic_records": 0,
        "snvs": 0,
        "mnvs": 0,
        "insertions": 0,
        "deletions": 0,
        "complex_variants": 0,
        "symbolic_alleles": 0,
        "non_variant_alleles": 0,
        "other_alleles": 0,
        "transitions": 0,
        "transversions": 0,
        "titv_ratio": None,
        "pass_records": 0,
        "mean_quality": None
    }

    quality_sum = 0.0
    quality_count = 0

    category_keys = {
        "snv": "snvs",
        "mnv": "mnvs",
        "insertion": "insertions",
        "deletion": "deletions",
        "complex": "complex_variants",
        "symbolic": "symbolic_alleles",
        "non_variant": "non_variant_alleles",
        "other": "other_alleles"
    }

    with input_path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\r\n").split("\t")

            ref = fields[3]
            alt_field = fields[4]
            qual = fields[5]
            filter_value = fields[6]

            stats["total_records"] += 1
            
            if filter_value == "PASS":
                stats["pass_records"] += 1

            if qual != ".":
                quality = float(qual)

                if math.isfinite(quality):
                    quality_sum += quality
                    quality_count += 1

            alternates = alt_field.split(",")

            if len(alternates) > 1:
                stats["multiallelic_records"] += 1

            for alt in alternates:
                category = classify_variant(ref, alt)

                stats[category_keys[category]] += 1

                if category != "non_variant":
                    stats["total_alt_alleles"] += 1

                if (
                    category == "snv"
                    and ref.upper() in "ACGT"
                    and alt.upper() in "ACGT"
                ):
                    pair = (ref.upper(), alt.upper())

                    if pair in TRANSITIONS:
                        stats["transitions"] += 1
                    else:
                        stats["transversions"] += 1

    if stats["transversions"] > 0:
        stats["titv_ratio"] = round(
            stats["transitions"] / stats["transversions"],
            4
        )

    if quality_count > 0:
        stats["mean_quality"] = round(
            quality_sum / quality_count,
            2
        )

    return stats


def write_vcf_stats(
        stats: dict,
        output_path: str | Path,
) -> None:
    """Write VCF stats as JSON."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("wt", encoding="utf-8") as handle:
        json.dump(stats, handle, indent=2)
        handle.write("\n")
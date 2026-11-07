
import math
import os
import tempfile
from pathlib import Path

from .vcf_stats import classify_variant
from .vcf_validate import validate_vcf


VARIANT_TYPES = (
    "snv",
    "mnv",
    "insertion",
    "deletion",
    "complex",
)


def get_info_depth(info_field: str) -> int | None:
    """Extract site-level DP from the VCF INFO field."""

    for entry in info_field.split(";"):
        if entry.startswith("DP="):
            value = entry[3:]

            try:
                depth = int(value)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid INFO/DP value: {value!r}"
                ) from exc

            if depth < 0:
                raise ValueError(
                    "INFO/DP cannot be negative."
                )

            return depth

    return None


def filter_vcf(
    input_path: str | Path,
    output_path: str | Path,
    min_qual: float | None = None,
    pass_only: bool = False,
    chrom: str | None = None,
    min_depth: int | None = None,
    variant_type: str | None = None,
) -> dict:
    """Filter VCF records while preserving headers."""

    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    if input_path.resolve() == output_path.resolve():
        raise ValueError(
            "Input and output must be different files."
        )

    if min_qual is not None and not math.isfinite(min_qual):
        raise ValueError("Minimum QUAL must be finite.")

    if min_depth is not None and min_depth < 0:
        raise ValueError("Minimum depth cannot be negative.")

    if variant_type is not None and variant_type not in VARIANT_TYPES:
        raise ValueError(
            f"Unsupported variant type: {variant_type}"
        )

    if chrom is not None and not chrom:
        raise ValueError("Chromosome cannot be empty.")

    validation = validate_vcf(input_path)

    if not validation["valid"]:
        raise ValueError(
            "Invalid VCF: " + "; ".join(validation["errors"][:5])
        )

    total_records = 0
    kept_records = 0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            dir=output_path.parent,
            prefix=".vcf_filter_",
            delete=False,
        ) as temp_file:
            temp_path = Path(temp_file.name)

            with input_path.open(
                "rt", encoding="utf-8"
            ) as source:
                for line_number, line in enumerate(source, start=1):
                    if line.startswith("#"):
                        temp_file.write(line)
                        continue

                    fields = line.rstrip("\r\n").split("\t")
                    total_records += 1

                    record_chrom = fields[0]
                    ref = fields[3]
                    alts = fields[4].split(",")
                    qual = fields[5]
                    filter_status = fields[6]
                    info = fields[7]

                    if chrom is not None and record_chrom != chrom:
                        continue

                    if pass_only and filter_status != "PASS":
                        continue

                    if min_qual is not None:
                        if qual == ".":
                            continue

                        try:
                            quality = float(qual)
                        except ValueError as exc:
                            raise ValueError(
                                f"Line {line_number}: invalid QUAL."
                            ) from exc

                        if not math.isfinite(quality):
                            raise ValueError(
                                f"Line {line_number}: non-finite QUAL."
                            )

                        if quality < min_qual:
                            continue

                    if min_depth is not None:
                        try:
                            depth = get_info_depth(info)
                        except ValueError as exc:
                            raise ValueError(
                                f"Line {line_number}: {exc}"
                            ) from exc

                        if depth is None or depth < min_depth:
                            continue

                    if variant_type is not None:
                        if not any(
                            classify_variant(ref, alt) == variant_type
                            for alt in alts
                        ):
                            continue

                    temp_file.write(line)
                    kept_records += 1

        os.replace(temp_path, output_path)

    except Exception:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
        raise

    return {
        "total_records": total_records,
        "kept_records": kept_records,
        "filtered_records": total_records - kept_records,
        "output_path": str(output_path),
    }
            
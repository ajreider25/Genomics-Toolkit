
import os
import tempfile
from pathlib import Path

from .vcf_validate import validate_vcf


def load_variant_ids(ids_file: str | Path) -> set[str]:
    """Load variant IDs, one per line."""

    ids = set()

    with Path(ids_file).open("rt", encoding="utf-8") as handle:
        for line in handle:
            variant_id = line.strip()

            if variant_id and not variant_id.startswith("#"):
                ids.add(variant_id)

    return ids


def extract_vcf(
    input_path: str | Path,
    output_path: str | Path,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    variant_id: str | None = None,
    ids_file: str | Path | None = None,
) -> dict:
    """Extract VCF records by chromosome, position, or ID."""

    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    if input_path.resolve() == output_path.resolve():
        raise ValueError(
            "Input and output must be different files."
        )

    if start is not None and start < 1:
        raise ValueError("Start position must be >= 1.")

    if end is not None and end < 1:
        raise ValueError("End position must be >= 1.")

    if (
        start is not None
        and end is not None
        and start > end
    ):
        raise ValueError(
            "Start position cannot exceed end position."
        )

    if (start is not None or end is not None) and not chrom:
        raise ValueError(
            "Chromosome is required for interval extraction."
        )

    if chrom == "":
        raise ValueError("Chromosome cannot be empty.")

    if variant_id == "":
        raise ValueError("Variant ID cannot be empty.")

    if not any([
        chrom is not None,
        variant_id is not None,
        ids_file is not None,
    ]):
        raise ValueError(
            "Provide --chrom, --id, or --ids-file."
        )

    validation = validate_vcf(input_path)

    if not validation["valid"]:
        raise ValueError(
            "Invalid VCF: "
            + "; ".join(validation["errors"][:5])
        )

    selected_ids = set()

    if variant_id is not None:
        selected_ids.add(variant_id)

    if ids_file is not None:
        selected_ids.update(load_variant_ids(ids_file))

    use_id_filter = variant_id is not None or ids_file is not None

    total_records = 0
    extracted_records = 0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            dir=output_path.parent,
            prefix=".vcf_extract_",
            delete=False,
        ) as temp_file:
            temp_path = Path(temp_file.name)

            with input_path.open(
                "rt", encoding="utf-8"
            ) as source:
                for line in source:
                    if line.startswith("#"):
                        temp_file.write(line)
                        continue

                    fields = line.rstrip("\r\n").split("\t")
                    total_records += 1

                    record_chrom = fields[0]
                    record_pos = int(fields[1])
                    record_ids = set(fields[2].split(";"))

                    if chrom is not None and record_chrom != chrom:
                        continue

                    if start is not None and record_pos < start:
                        continue

                    if end is not None and record_pos > end:
                        continue

                    if use_id_filter and not (
                        record_ids & selected_ids
                    ):
                        continue

                    temp_file.write(line)
                    extracted_records += 1

        os.replace(temp_path, output_path)

    except Exception:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
        raise

    return {
        "total_records": total_records,
        "extracted_records": extracted_records,
        "excluded_records": total_records - extracted_records,
        "output_path": str(output_path),
    }
from pathlib import Path

def validate_vcf(input_path: str | Path) -> dict:
    """Perform basic structural validation of a text VCF file."""

    input_path = Path(input_path)

    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    errors = []
    record_count = 0
    header_found = False
    fileformat_found = False
    expected_columns = None

    with input_path.open("rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")

            if line.startswith("##fileformat="):
                if record_count or header_found:
                    errors.append(
                        f"Line {line_number}: misplaced fileformat."
                    )
                fileformat_found = True
                continue

            if line.startswith('##'):
                if header_found:
                    errors.append(
                        f"Line {line_number}: metadata after header."
                    )
                continue

            if line.startswith("#CHROM"):
                if header_found:
                    errors.append(
                        f"Line {line_number}: duplicate column header."
                    )

                columns = line.split("\t")

                if columns[:8] != [
                    "#CHROM", "POS", "ID", "REF",
                    "ALT", "QUAL", "FILTER", "INFO"
                ]:
                    errors.append(
                        f"Line {line_number}: invalid VCF columns."
                    )

                if len(columns) > 8:
                    if columns[8] != "FORMAT":
                        errors.append(
                            f"Line {line_number}: expected FORMAT column."
                        )
                    if len(columns) == 9:
                        errors.append(
                            f"Line {line_number}: FORMAT without samples."
                        )

                expected_columns = len(columns)
                header_found = True
                continue

            if line.startswith("#"):
                errors.append(
                    f"Line {line_number}: unexpected header line."
                )
                continue

            if not line:
                errors.append(
                    f"Line {line_number}: empty record line."
                )

            if not header_found:
                errors.append(
                    f"Line {line_number}: record before column header."
                )
                continue

            fields = line.split("\t")

            if len(fields) != expected_columns:
                errors.append(
                    f"Line {line_number}: incorrect column count."
                )
                continue

            record_count += 1

            try:
                position = int(fields[1])
                if position < 1:
                    raise ValueError
            except ValueError:
                errors.append(
                    f"Line {line_number}: invalid POS value."
                )

            if not fields[0] or fields[0] == ".":
                errors.append(
                    f"Line {line_number}: invalid CHROM value."
                )

            if not fields[3] or fields[3] == ".":
                errors.append(
                    f"Line {line_number}: invalid REF allele."
                )

            if not fields[4]:
                errors.append(
                    f"Line {line_number}: empty ALT field."
                )

            if fields[5] != ".":
                try:
                    float(fields[5])
                except ValueError:
                    errors.append(
                        f"Line {line_number}: invalid QUAL value."
                    )

    if not fileformat_found:
        errors.append("Missing VCF fileformat declaration.")

    if not header_found:
        errors.append("Missing #CHROM column header.")

    return {
        "valid": len(errors) == 0,
        "record_count": record_count,
        "error_count": len(errors),
        "errors": errors
    }
    
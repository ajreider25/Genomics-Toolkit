
import gzip
import json
from pathlib import Path


DNA_BASES = set("ACGTN")


def validate_fastq(
    input_path: str | Path,
    strict: bool = False,
) -> dict:
    """Validate a four-line FASTQ file using Phred+33 encoding."""

    input_path = Path(input_path)

    if not input_path.is_file():
        raise FileNotFoundError(
            f"FASTQ file not found: {input_path}"
        )

    errors = []
    warnings = []
    record_count = 0
    total_bases = 0

    opener = (
        gzip.open
        if input_path.suffix.lower() == ".gz"
        else open
    )

    try:
        with opener(
            input_path, "rt", encoding="utf-8"
        ) as handle:

            while True:
                header = handle.readline()

                if not header:
                    break

                sequence = handle.readline()
                separator = handle.readline()
                quality = handle.readline()

                record_count += 1

                if not all((sequence, separator, quality)):
                    errors.append(
                        f"Record {record_count}: incomplete FASTQ record."
                    )
                    break

                header = header.rstrip("\r\n")
                sequence = sequence.rstrip("\r\n")
                separator = separator.rstrip("\r\n")
                quality = quality.rstrip("\r\n")

                if not header.startswith("@"):
                    errors.append(
                        f"Record {record_count}: header must start with '@'."
                    )

                if header == "@":
                    errors.append(
                        f"Record {record_count}: empty sequence identifier."
                    )

                if not separator.startswith("+"):
                    errors.append(
                        f"Record {record_count}: separator must start with '+'."
                    )

                if not sequence:
                    message = (
                        f"Record {record_count}: empty sequence."
                    )
                    (errors if strict else warnings).append(message)

                invalid_bases = sorted(
                    set(sequence.upper()) - DNA_BASES
                )

                if invalid_bases:
                    errors.append(
                        f"Record {record_count}: invalid bases: "
                        + ", ".join(invalid_bases)
                    )

                if len(sequence) != len(quality):
                    errors.append(
                        f"Record {record_count}: sequence and quality "
                        "lengths do not match."
                    )

                # Printable Phred+33 quality characters:
                # ASCII 33 through 126.
                if any(
                    ord(char) < 33 or ord(char) > 126
                    for char in quality
                ):
                    errors.append(
                        f"Record {record_count}: invalid quality characters."
                    )

                # If the separator repeats the identifier,
                # verify it matches the header.
                if separator != "+":
                    if separator[1:] != header[1:]:
                        errors.append(
                            f"Record {record_count}: separator "
                            "identifier does not match header."
                        )

                total_bases += len(sequence)

    except (OSError, EOFError, UnicodeError) as exc:
        errors.append(f"Unable to read FASTQ file: {exc}")

    if record_count == 0:
        message = "FASTQ file contains no records."
        (errors if strict else warnings).append(message)

    return {
        "valid": len(errors) == 0,
        "record_count": record_count,
        "total_bases": total_bases,
        "error_count": len(errors),
        "warning_count": len(warnings),
        "errors": errors,
        "warnings": warnings,
    }


def write_validation_report(
    report: dict,
    output_path: str | Path,
    input_path: str | Path | None = None,
) -> None:
    """Write a FASTQ validation report as JSON."""

    output_path = Path(output_path)

    if (
        input_path is not None
        and output_path.resolve() == Path(input_path).resolve()
    ):
        raise ValueError(
            "Report output cannot overwrite the FASTQ input."
        )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    
import json
from pathlib import Path

from .fasta_stats import read_fasta

VALID_DNA_BASES = set("ACGTRYSWKMBDHVN")

def validate_fasta(
        input_path: str | Path,
        strict: bool = False
) -> dict:
    """Validate a DNA FASTA file and return a structured report."""
    input_path = Path(input_path)

    if not input_path.is_file():
        raise FileExistsError(
            f"FASTA input file npt found: {input_path}"
        )

    errors = []
    warnings = []
    seen_ids = set()
    sequence_count = 0

    try:
        for sequence_id, sequence in read_fasta(input_path):
            sequence_count += 1

            if sequence_id in seen_ids:
                errors.append(f"Duplicate sequence ID: {sequence_id}")
            seen_ids.add(sequence_id)

            if not sequence:
                message = f"Empty sequence: {sequence_id}"
                if strict:
                    errors.append(message)
                else:
                    warnings.append(message)

            invalid_bases = sorted(
                set(sequence.upper()) - VALID_DNA_BASES
            )

            if invalid_bases:
                errors.append(
                    f"Invalid nucleotide(s) in {sequence_id}: "
                    f"{' ,'.join(invalid_bases)}"
                )

    except (ValueError, EOFError, OSError) as exc:
        errors.append(f"FASTA parsing error: {exc}")

    if sequence_count == 0 and not errors:
        message = "FASTA file contains no sequences"

        if strict:
            errors.append(message)
        else:
            warnings.append(message)

    return {
        "valid": len(errors) == 0,
        "sequence_count": sequence_count,
        "error_count": len(errors),
        "warning_count": len(warnings),
        "errors": errors,
        "warnings": warnings
    }

def write_validation_report(
    report: dict,
    output_path: str | Path,
    input_path: str | Path | None = None,
) -> None:
    """Write a FASTA validation report as JSON."""
    output_path = Path(output_path)

    if input_path is not None:
        if Path(input_path).resolve() == output_path.resolve():
            raise ValueError(
                "Input FASTA and JSON report must be different files."
            )
        
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
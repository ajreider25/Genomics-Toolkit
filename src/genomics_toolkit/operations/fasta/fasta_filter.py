import gzip
import tempfile
import os
from pathlib import Path

from .fasta_stats import read_fasta

def filter_fasta(
    input_path: str | Path,
    output_path: str | Path,
    min_length: int | None = None,
    max_length: int | None = None,
    ids: set[str] | None = None
) -> dict:
    """Filter FASTA records and write matching sequences."""
    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(
            f"FASTA input file not found: {input_path}"
        )

    if input_path.resolve() == output_path.resolve():
        raise ValueError(
            "Input and output cannot be the same file."
        )

    if min_length is not None and min_length < 0:
        raise ValueError("Maximum length cannot be negative.")

    if (
        min_length is not None
        and max_length is not None
        and min_length > max_length
    ):
        raise ValueError(
            "Minimum length cannot exceed maximum length."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    total_records = 0
    retained_records = 0
    temporary_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False
        ) as temporary_file:
            temporary_path = Path(temporary_file.name)

        opener = (
            gzip.open
            if output_path.suffix == ".gz"
            else open
        )

        with opener(
            temporary_path,
            "wt",
            encoding="utf-8"
        ) as handle:
            for sequence_id, sequence in read_fasta(input_path):
                total_records += 1

                if ids is not None and sequence_id not in ids:
                    continue

                if (
                    min_length is not None
                    and len(sequence) < min_length
                ):
                    continue

                if (
                    max_length is not None 
                    and len(sequence) > max_length
                ):
                    continue

                handle.write(f">{sequence_id}\n")

                for start in range(0, len(sequence), 80):
                    handle.write(sequence[start:start + 80] + "\n")

                retained_records += 1

        os.replace(temporary_path, output_path)
        temporary_path = None

    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    return {
        "total_records": total_records,
        "retained_records": retained_records,
        "removed_records": total_records - retained_records
    }
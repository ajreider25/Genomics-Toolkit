import gzip
import os
import tempfile
from pathlib import Path

from .fasta_stats import read_fasta

def read_ids_file(path: str| Path) -> set[str]:
    """Read sequence identifiers from a text file."""
    path = Path(path)
    selected_ids = set()

    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            sequence_id = line.strip()

            if not sequence_id or sequence_id.startswith("#"):
                continue

            if sequence_id in selected_ids:
                raise ValueError(
                    f"Duplicate ID in selection file "
                    f"at line {line_number}: {sequence_id}"
                )

            selected_ids.add(sequence_id)

    if not selected_ids:
        raise ValueError(
            "The ID selection file contains no sequence IDs."
        )

    return selected_ids

def extract_fasta(
        input_path: str | Path,
        ids_file: str | Path,
        output_path: str | Path
) -> dict:
    """Extract selected FASTA records and write output atomically."""
    input_path = Path(input_path)
    ids_file = Path(ids_file)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(
            f"FASTA input file not found: {input_path}"
        )

    if not ids_file.is_file():
        raise FileNotFoundError(
            f"ID selection file not found: {ids_file}"
        )

    if output_path.resolve() in {
        input_path.resolve(),
        ids_file.resolve()
    }:
        raise ValueError(
            "Output cannot be the same file as an input."
        )

    selected_ids = read_ids_file(ids_file)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    total_records = 0
    extracted_records = 0
    found_ids = set()
    temporary_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary_file:
            temporary_path = Path(temporary_file.name)

        opener = gzip.open if output_path.suffix == ".gz" else open

        with opener(
            temporary_path,
            "wt",
            encoding="utf-8",
        ) as handle:
            for sequence_id, sequence in read_fasta(input_path):
                total_records += 1

                if sequence_id not in selected_ids:
                    continue

                handle.write(f">{sequence_id}\n")

                for start in range(0, len(sequence), 80):
                    handle.write(sequence[start:start + 80] + "\n")

                found_ids.add(sequence_id)
                extracted_records += 1

        os.replace(temporary_path, output_path)
        temporary_path = None

    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    missing_ids = sorted(selected_ids - found_ids)

    return {
        "total_records": total_records,
        "extracted_records": extracted_records,
        "requested_ids": len(selected_ids),
        "found_ids": len(found_ids),
        "missing_ids": missing_ids
    }
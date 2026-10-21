
import gzip
import os
import tempfile
from pathlib import Path
from typing import Iterable

from .fasta_stats import read_fasta


def merge_fasta(
    input_paths: Iterable[str | Path],
    output_path: str | Path,
    allow_duplicates: bool = False,
) -> dict:
    """Merge FASTA files while preserving input order."""

    input_paths = [Path(path) for path in input_paths]
    output_path = Path(output_path)

    if not input_paths:
        raise ValueError("At least one input FASTA is required.")

    for path in input_paths:
        if not path.is_file():
            raise FileNotFoundError(
                f"FASTA input file not found: {path}"
            )

    if output_path.resolve() in {
        path.resolve() for path in input_paths
    }:
        raise ValueError(
            "Output cannot be the same file as an input."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    total_records = 0
    total_bases = 0
    duplicate_count = 0
    seen_ids = set()

    temporary_path = None

    try:
        # Create temporary output in the destination directory.
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary_file:
            temporary_path = Path(temporary_file.name)

        opener = (
            gzip.open
            if output_path.suffix.lower() == ".gz"
            else open
        )

        with opener(
            temporary_path,
            "wt",
            encoding="utf-8",
        ) as handle:
            for input_path in input_paths:
                for sequence_id, sequence in read_fasta(
                    input_path
                ):
                    if sequence_id in seen_ids:
                        duplicate_count += 1

                        if not allow_duplicates:
                            raise ValueError(
                                "Duplicate FASTA ID: "
                                f"{sequence_id}"
                            )

                    seen_ids.add(sequence_id)

                    handle.write(f">{sequence_id}\n")

                    for start in range(0, len(sequence), 80):
                        handle.write(
                            sequence[start:start + 80] + "\n"
                        )

                    total_records += 1
                    total_bases += len(sequence)

        # Replace the destination only after successful processing.
        os.replace(temporary_path, output_path)
        temporary_path = None

    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    return {
        "input_files": len(input_paths),
        "total_records": total_records,
        "total_bases": total_bases,
        "duplicate_ids": duplicate_count,
        "output_path": str(output_path),
    }

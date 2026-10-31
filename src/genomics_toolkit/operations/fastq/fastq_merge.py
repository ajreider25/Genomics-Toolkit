
import gzip
import os
import tempfile
from pathlib import Path

from .fastq_stats import read_fastq


def merge_fastq(
    input_paths: list[str | Path],
    output_path: str | Path,
) -> dict:
    """Merge FASTQ files in input order, preserving read data."""

    if not input_paths:
        raise ValueError("At least one input FASTQ file is required.")

    inputs = [Path(path) for path in input_paths]
    output_path = Path(output_path)

    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(f"FASTQ file not found: {path}")

        if path.resolve() == output_path.resolve():
            raise ValueError(
                "Output path cannot be one of the input files."
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    total_reads = 0
    total_bases = 0
    temp_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output_path.parent,
            prefix=".fastq_merge_",
            delete=False,
        ) as temp_file:
            temp_path = Path(temp_file.name)

        opener = (
            gzip.open
            if output_path.suffix.lower() == ".gz"
            else open
        )

        with opener(
            temp_path,
            "wt",
            encoding="utf-8",
            newline="\n",
        ) as handle:
            for input_path in inputs:
                for read_id, sequence, quality in read_fastq(input_path):
                    handle.write(f"@{read_id}\n")
                    handle.write(f"{sequence}\n")
                    handle.write("+\n")
                    handle.write(f"{quality}\n")

                    total_reads += 1
                    total_bases += len(sequence)

        os.replace(temp_path, output_path)

    except Exception:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
        raise

    return {
        "input_files": len(inputs),
        "total_reads": total_reads,
        "total_bases": total_bases,
        "output_path": str(output_path),
    }

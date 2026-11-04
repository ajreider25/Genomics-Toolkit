import gzip
import os
import tempfile
from pathlib import Path

from .fastq_stats import read_fastq


def filter_fastq(
    input_path: str | Path,
    output_path: str | Path,
    min_length: int | None = None,
    max_length: int | None = None,
    min_quality: float | None = None,
    max_n_percent: float | None = None,
) -> dict:
    """Filter FASTQ reads by length, quality, and N content."""

    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(f"FASTQ file not found: {input_path}")

    if input_path.resolve() == output_path.resolve():
        raise ValueError("Input and output paths must be different.")

    if min_length is not None and min_length < 0:
        raise ValueError("Minimum length cannot be negative.")

    if max_length is not None and max_length < 0:
        raise ValueError("Maximum length cannot be negative.")

    if (
        min_length is not None
        and max_length is not None
        and min_length > max_length
    ):
        raise ValueError("Minimum legnth cannot exceed maximum length.")

    if min_quality is not None and not 0 <= min_quality <= 93:
        raise ValueError("Minimum Phred quality must be between 0 and 93.")

    if max_n_percent is not None and not 0 <= max_n_percent <= 100:
        raise ValueError("Maximum N percentage must be between 0 and 100.")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    total_reads = 0
    kept_reads = 0
    filtered_reads = 0

    # Write to a temporary file so errors do not leave partial output
    temp_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output_path.parent,
            prefix=".fasta_filter_",
            delete=False
        ) as temp_file:
            temp_path = Path(temp_file.name)

        opener = gzip.open if output_path.suffix.lower() == ".gz" else open

        with opener(temp_path, "wt", encoding="utf-8", newline="\n") as handle:
            for read_id, sequence, quality in read_fastq(input_path):
                total_reads += 1
                length = len(sequence)

                mean_quality = (
                    sum(ord(char) - 33 for char in quality) / length
                    if length else 0.0
                )

                n_percent = (
                    100 * sequence.upper().count("N") / length
                    if length else 0.0
                )

                keep = True

                if min_length is not None and length < min_length:
                    keep = False

                if max_length is not None and length > max_length:
                    keep = False

                if min_quality is not None and mean_quality < min_quality:
                    keep = False

                if max_n_percent is not None and n_percent > max_n_percent:
                    keep = False

                if keep:
                    handle.write(f"@{read_id}\n")
                    handle.write(f"{sequence}\n")
                    handle.write("+\n")
                    handle.write(f"{quality}\n")
                    kept_reads += 1
                else:
                    filtered_reads += 1

        os.replace(temp_path, output_path)

    except Exception:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
        raise

    return {
        "total_reads": total_reads,
        "kept_reads": kept_reads,
        "filtered_reads": filtered_reads,
        "output_path": str(output_path),
    }

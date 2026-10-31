import gzip
import os
import tempfile
from pathlib import Path

from .fastq_stats import read_fastq


def trim_read(
    sequence: str,
    quality: str,
    min_quality: int,
) -> tuple[str, str]:
    """Remove low-quality bases from the 3' end of a read."""

    if not 0 <= min_quality <= 93:
        raise ValueError("Minimum Phred quality must be between 0 and 93.")

    if len(sequence) != len(quality):
        raise ValueError("Sequence and quality lengths do not match.")

    end = len(sequence)

    while end > 0 and ord(quality[end - 1]) - 33 < min_quality:
        end -= 1

    return sequence[:end], quality[:end]


def trim_fastq(
    input_path: str | Path,
    output_path: str | Path,
    min_quality: int = 20,
    min_length: int = 0,
) -> dict:
    """Trim low-quality 3' ends and discard reads below min_length."""

    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(f"FASTQ file not found: {input_path}")

    if input_path.resolve() == output_path.resolve():
        raise ValueError("Input and output paths must be different.")

    if not 0 <= min_quality <= 93:
        raise ValueError("Minimum Phred quality must be between 0 and 93.")

    if min_length < 0:
        raise ValueError("Minimum read length cannot be negative.")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    total_reads = 0
    trimmed_reads = 0
    discarded_reads = 0
    kept_reads = 0
    bases_removed = 0
    temp_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output_path.parent,
            prefix=".fastq_trim_",
            delete=False,
        ) as temp_file:
            temp_path = Path(temp_file.name)

        opener = gzip.open if output_path.suffix.lower() == ".gz" else open

        with opener(temp_path, "wt", encoding="utf-8", newline="\n") as handle:
            for read_id, sequence, quality in read_fastq(input_path):
                total_reads += 1

                trimmed_sequence, trimmed_quality = trim_read(
                    sequence,
                    quality,
                    min_quality,
                )

                removed = len(sequence) - len(trimmed_sequence)
                bases_removed += removed

                if removed > 0:
                    trimmed_reads += 1

                if len(trimmed_sequence) < min_length:
                    discarded_reads += 1
                    continue

                handle.write(f"@{read_id}\n")
                handle.write(f"{trimmed_sequence}\n")
                handle.write("+\n")
                handle.write(f"{trimmed_quality}\n")

                kept_reads += 1

        os.replace(temp_path, output_path)

    except Exception:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
        raise

    return {
        "total_reads": total_reads,
        "kept_reads": kept_reads,
        "trimmed_reads": trimmed_reads,
        "discarded_reads": discarded_reads,
        "bases_removed": bases_removed,
        "output_path": str(output_path),
    }
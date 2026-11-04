
import gzip
import json
from pathlib import Path


VALID_BASES = set("ACGTN")


def read_fastq(path: str | Path):
    """Yield (read_id, sequence, quality) from a four-line FASTQ.

    Quality strings use Phred+33 encoding.
    """
    path = Path(path)

    opener = gzip.open if path.suffix.lower() == ".gz" else open

    with opener(path, "rt", encoding="utf-8") as handle:
        record_number = 0

        while True:
            header = handle.readline()

            if not header:
                break

            record_number += 1

            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()

            if not all((sequence, separator, quality)):
                raise ValueError(
                    f"Record {record_number}: incomplete FASTQ record."
                )

            header = header.rstrip("\r\n")
            sequence = sequence.rstrip("\r\n")
            separator = separator.rstrip("\r\n")
            quality = quality.rstrip("\r\n")

            if not header.startswith("@") or len(header) == 1:
                raise ValueError(
                    f"Record {record_number}: invalid FASTQ header."
                )

            if not separator.startswith("+"):
                raise ValueError(
                    f"Record {record_number}: invalid separator."
                )

            if separator != "+" and separator[1:] != header[1:]:
                raise ValueError(
                    f"Record {record_number}: separator ID mismatch."
                )

            if len(sequence) != len(quality):
                raise ValueError(
                    f"Record {record_number}: sequence and quality "
                    "lengths do not match."
                )

            if set(sequence.upper()) - VALID_BASES:
                raise ValueError(
                    f"Record {record_number}: invalid nucleotide."
                )

            if any(
                ord(char) < 33 or ord(char) > 126
                for char in quality
            ):
                raise ValueError(
                    f"Record {record_number}: invalid quality character."
                )

            yield header[1:], sequence, quality


def calculate_fastq_stats(input_path: str | Path) -> dict:
    """Calculate read-level and base-quality statistics."""

    input_path = Path(input_path)

    if not input_path.is_file():
        raise FileNotFoundError(
            f"FASTQ file not found: {input_path}"
        )

    read_count = 0
    total_bases = 0
    gc_bases = 0
    called_bases = 0

    min_length = None
    max_length = 0

    quality_sum = 0
    q20_bases = 0
    q30_bases = 0

    position_quality_sums = []
    position_counts = []

    for _, sequence, quality in read_fastq(input_path):
        read_count += 1
        length = len(sequence)
        total_bases += length

        min_length = (
            length if min_length is None
            else min(min_length, length)
        )
        max_length = max(max_length, length)

        upper_sequence = sequence.upper()

        gc_bases += sum(
            base in "GC" for base in upper_sequence
        )
        called_bases += sum(
            base in "ACGT" for base in upper_sequence
        )

        for position, char in enumerate(quality):
            phred = ord(char) - 33

            quality_sum += phred

            if phred >= 20:
                q20_bases += 1

            if phred >= 30:
                q30_bases += 1

            if position >= len(position_quality_sums):
                position_quality_sums.append(0)
                position_counts.append(0)

            position_quality_sums[position] += phred
            position_counts[position] += 1

    mean_read_length = (
        total_bases / read_count
        if read_count else 0.0
    )

    gc_percent = (
        100 * gc_bases / called_bases
        if called_bases else 0.0
    )

    mean_quality = (
        quality_sum / total_bases
        if total_bases else 0.0
    )

    q20_percent = (
        100 * q20_bases / total_bases
        if total_bases else 0.0
    )

    q30_percent = (
        100 * q30_bases / total_bases
        if total_bases else 0.0
    )

    per_base_quality = [
        {
            "position": index + 1,
            "mean_quality": quality_total / count,
            "read_count": count,
        }
        for index, (quality_total, count) in enumerate(
            zip(position_quality_sums, position_counts)
        )
    ]

    return {
        "read_count": read_count,
        "total_bases": total_bases,
        "min_read_length": min_length if min_length is not None else 0,
        "max_read_length": max_length,
        "mean_read_length": mean_read_length,
        "gc_percent": gc_percent,
        "mean_quality": mean_quality,
        "q20_percent": q20_percent,
        "q30_percent": q30_percent,
        "per_base_quality": per_base_quality,
    }


def write_fastq_stats(
    stats: dict,
    output_path: str | Path,
    input_path: str | Path | None = None,
) -> None:
    """Write FASTQ statistics to a JSON file."""

    output_path = Path(output_path)

    if (
        input_path is not None
        and output_path.resolve() == Path(input_path).resolve()
    ):
        raise ValueError(
            "Statistics output cannot overwrite FASTQ input."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(stats, handle, indent=2)
        handle.write("\n")

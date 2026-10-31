
import gzip
import os
import random
import tempfile
from pathlib import Path

from .fastq_stats import read_fastq


def sample_fastq(
    input_path: str | Path,
    output_path: str | Path,
    count: int | None = None,
    fraction: float | None = None,
    seed: int = 42,
) -> dict:
    """Randomly sample FASTQ reads by exact count or probability."""

    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(f"FASTQ file not found: {input_path}")

    if input_path.resolve() == output_path.resolve():
        raise ValueError("Input and output paths must be different.")

    if (count is None) == (fraction is None):
        raise ValueError(
            "Specify exactly one of count or fraction."
        )

    if count is not None and count < 0:
        raise ValueError("Sample count cannot be negative.")

    if fraction is not None and not 0 <= fraction <= 1:
        raise ValueError("Sample fraction must be between 0 and 1.")

    rng = random.Random(seed)

    total_reads = 0
    sampled_reads = 0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output_path.parent,
            prefix=".fastq_sample_",
            delete=False,
        ) as temp_file:
            temp_path = Path(temp_file.name)

        opener = gzip.open if output_path.suffix.lower() == ".gz" else open

        if count is not None:
            # Reservoir sampling: uniform sample of up to count reads.
            # Store original indices to restore input order afterward.
            reservoir = []

            for read_id, sequence, quality in read_fastq(input_path):
                total_reads += 1

                if count == 0:
                    continue

                record = (
                    total_reads - 1,
                    read_id,
                    sequence,
                    quality,
                )

                if len(reservoir) < count:
                    reservoir.append(record)
                else:
                    index = rng.randrange(total_reads)

                    if index < count:
                        reservoir[index] = record

            reservoir.sort(key=lambda record: record[0])

            with opener(temp_path, "wt", encoding="utf-8", newline="\n") as handle:
                for _, read_id, sequence, quality in reservoir:
                    handle.write(f"@{read_id}\n")
                    handle.write(f"{sequence}\n")
                    handle.write("+\n")
                    handle.write(f"{quality}\n")

                    sampled_reads += 1

        else:
            # Bernoulli sampling: independently select each read.
            with opener(temp_path, "wt", encoding="utf-8", newline="\n") as handle:
                for read_id, sequence, quality in read_fastq(input_path):
                    total_reads += 1

                    if rng.random() < fraction:
                        handle.write(f"@{read_id}\n")
                        handle.write(f"{sequence}\n")
                        handle.write("+\n")
                        handle.write(f"{quality}\n")

                        sampled_reads += 1

        os.replace(temp_path, output_path)

    except Exception:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
        raise

    return {
        "total_reads": total_reads,
        "sampled_reads": sampled_reads,
        "discarded_reads": total_reads - sampled_reads,
        "seed": seed,
        "method": "count" if count is not None else "fraction",
        "output_path": str(output_path),
    }

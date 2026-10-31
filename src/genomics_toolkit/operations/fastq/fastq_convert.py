
import gzip
import os
import tempfile
from pathlib import Path

from .fastq_stats import read_fastq


def fastq_to_fasta(
    input_path: str | Path,
    output_path: str | Path,
    line_width: int = 80,
) -> dict:
    """Convert FASTQ records to FASTA, discarding quality scores."""

    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(
            f"FASTQ file not found: {input_path}"
        )

    if input_path.resolve() == output_path.resolve():
        raise ValueError("Input and output paths must be different.")

    if line_width < 1:
        raise ValueError("Line width must be at least 1.")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    record_count = 0
    total_bases = 0
    temp_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output_path.parent,
            prefix=".fastq_convert_",
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
            for read_id, sequence, quality in read_fastq(input_path):
                handle.write(f">{read_id}\n")

                for start in range(0, len(sequence), line_width):
                    handle.write(
                        sequence[start:start + line_width] + "\n"
                    )

                # FASTA representation of an empty sequence.
                if not sequence:
                    handle.write("\n")

                record_count += 1
                total_bases += len(sequence)

        os.replace(temp_path, output_path)

    except Exception:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
        raise

    return {
        "record_count": record_count,
        "total_bases": total_bases,
        "output_path": str(output_path),
    }

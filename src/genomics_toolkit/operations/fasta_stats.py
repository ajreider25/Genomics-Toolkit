import csv
from collections import Counter
from pathlib import Path
from typing import Iterator

def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield (sequence_id, sequence) pairs from a FASTA file."""

    sequence_id = None
    sequence_parts = []

    with Path(path).open("r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if sequence_id is not None:
                    yield sequence_id, "".join(sequence_parts)

                sequence_id = line[1:].split()[0]
                sequence_parts = []
            else:
                if sequence_id is None:
                    raise ValueError(
                        "Invalid FASTA: sequence before header"
                    )
                sequence_parts.append(line.upper())

    if sequence_id is not None:
        yield sequence_id, "".join(sequence_parts)

def calculate_stats(sequence: str) -> dict:
    """Calculate basic stats for a DNA sequence"""

    counts = Counter(sequence.upper())

    length = len(sequence)
    gc_count = counts["G"] + counts["C"]

    canonical_bases = sum(
        counts[base] for base in "ACGT"
    )

    gc_percent = (
        100 * gc_count / canonical_bases
        if canonical_bases
        else 0.0
    )

    return {
        "length": length,
        "gc_percent": round(gc_percent, 2),
        "n_count": counts["N"]
    }

def fasta_stats(path: str | Path) -> Iterator[dict]:
    """Generate stats for every FASTA record."""

    for sequence_id, sequence in read_fasta(path):
        yield {
            "sequence_id": sequence_id,
            **calculate_stats(sequence)
        }

def write_fasta_stats(
        input_path: str | Path,
        output_path: str | Path
) -> None:
    """Write FASTA sequence stats to a CSV file"""

    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(
            f"FASTA input file not found: {input_path}"
        )
    
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    fieldnames = [
        "sequence_id",
        "length",
        "gc_percent",
        "n_count"
    ]

    with output_path.open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames
        )

        writer.writeheader()

        for record in fasta_stats(input_path):
            writer.writerow(record)